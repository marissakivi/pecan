
#/fs/data2/output//PEcAn_1000009519 Lorimer
#/fs/data2/output/PEcAn_10000010299 HF
#/fs/data2/output/PEcAn_1000010320 NRP
#/fs/data2/output/PEcAn_1000010323 RH
#/fs/data2/output//PEcAn_1000010510 HF new spp.

setwd('/fs/data2/output/PEcAn_1000010510')

library(plyr)
library(PEcAn.settings)


input <- readRDS('/data/dbfiles/NPP_STAT_MODEL_HF.RDS')

#levels(input$taxon)[3] <- 'BEAL2' # Harvard Forest
#levels(input$taxon)[8] <- 'FRAM2'#7 Harvard Forest

levels(input$taxon)[4] <- 'PIGL' #HACK for Rooster Hill

levels(input$type)[1] <- 'ab'
levels(input$type)[2] <- 'abi'

settings = read.settings("pecan.SDA.xml")

var.names = 'AGB.pft'

get_treering_agb_pft <- function(settings, input, var.names){
  
  d <- settings$database$bety[c("dbname", "password", "host", "user")]
  bety <- src_postgres(host = d$host, user = d$user, password = d$password, dbname = d$dbname)
  
  #load_data_paleon_sda(settings)
  obvs<-list()
  obvs[[1]] <- input
  
  start_date <- settings$state.data.assimilation$start.date
  end_date   <- settings$state.data.assimilation$end.date
  
  obs.mean <- obs.mean.tmp <- list()
  obs.cov <- obs.cov.tmp <- list()
  
  obs.times <- seq(as.Date(start_date), as.Date(end_date), by = settings$state.data.assimilation$forecast.time.step)
  obs.times <- formatC(lubridate::year(obs.times), width = 4, format = "d", flag = "0")
  
  time.type <- 'year'
  
  i=1
  biomass2carbon <- 0.48
  
  obvs[[i]] <- obvs[[i]][obvs[[i]]$model=='Model RW + Census',] # + Census for HF
  #browser()
  if(!is.null(obvs[[i]]$value)) obvs[[i]]$value <- udunits2::ud.convert(obvs[[i]]$value,'Mg/ha','kg/m^2') * biomass2carbon #* kgm2Mgha 
  #if(!is.null(obvs[[i]]$abi)) obvs[[i]]$GWBI <- obvs[[i]]$GWBI * biomass2carbon  #* kgms2Mghayr 
  arguments <- list(.(year, iter, site_id), .(variable)) #, site_id
  arguments2 <- list(.(year), .(variable))
  arguments3 <- list(.(iter), .(variable), .(year))
  
  dataset <- obvs[[i]]
  
  ### Map species to model specific PFTs
  
  spp_id <- match_species_id(unique(dataset$taxon),format_name = 'usda',bety)
  pft_mat <- match_pft(spp_id$bety_species_id, settings$pfts,
                       con = bety$con, allow_missing = TRUE)
  
  x <- paste0('AGB.pft.', pft_mat$pft)
  names(x) <- spp_id$input_code
  
  PEcAn.logger::logger.info('Now, mapping data species to model PFTs')
  dataset$pft.cat <- x[as.character(dataset$taxon)]
  dataset <- dataset[dataset$pft.cat!='AGB.pft.NA',]
  
  variable <- c('value') #'ab'
  arguments <- list(.(site_id,year, iter,  pft.cat, site_id), .(variable))#, site_id
  arguments2 <- list(.(year, pft.cat), .(variable))
  arguments3 <- list(.(iter), .(pft.cat, variable), .(year))
  
  
  PEcAn.logger::logger.info('Now, aggregating data and creating SDA input lists')
  melt_id <- colnames(dataset)[-which(colnames(dataset) %in% variable)]
  melt.test <- reshape2::melt(dataset, id = melt_id, na.rm = TRUE)
  cast.test <- reshape2::dcast(melt.test, arguments, sum, margins = variable)
  
  melt_id_next <- colnames(cast.test)[-which(colnames(cast.test) %in% variable)]
  melt.next <- reshape2::melt(cast.test, id = melt_id_next)
  mean_mat <<- reshape2::dcast(melt.next, arguments2, mean)
  
  #browser()
  
  iter_mat <<- reshape2::acast(melt.next, arguments3, mean)
  cov.test <- apply(iter_mat,3,function(x){cov(x)})
  
  for(t in seq_along(obs.times)){
    obs.mean.tmp[[t]] <- mean_mat[mean_mat[,time.type]==obs.times[t], -c(1)] #THIS WONT WORK IF TIMESTEP ISNT ANNUAL
    
    if(any(var.names == 'AGB.pft')){
      obs.mean.tmp[[t]] <- rep(NA, length(unique(dataset$pft.cat)))
      names(obs.mean.tmp[[t]]) <- sort(unique(dataset$pft.cat))
      for(r in seq_along(unique(dataset$pft.cat))){
        k <- mean_mat[mean_mat$year==obs.times[t] & mean_mat$pft.cat==names(obs.mean.tmp[[t]][r]), variable]
        if(any(k)){
          obs.mean.tmp[[t]][r] <- k
        }
      }
    }
    
    obs.cov.tmp[[t]] <- matrix(cov.test[,which(colnames(cov.test) %in% obs.times[t])],
                               ncol = sqrt(dim(cov.test)[1]),
                               nrow = sqrt(dim(cov.test)[1]))
    
    colnames(obs.cov.tmp[[t]]) <- names(iter_mat[1,,t]) 
    
  }
  
  obs.mean <- obs.mean.tmp
  obs.cov <- obs.cov.tmp
  
  names(obs.mean) <- paste0(obs.times,'/12/31')
  names(obs.cov) <- paste0(obs.times,'/12/31')
  
  obs.list <- list(obs.mean = obs.mean, obs.cov = obs.cov)
  save(obs.list,file=file.path(settings$outdir,'sda.obs.Rdata'))
  
}

get_treering_agb_pft(settings = settings, input = input,
                     var.names = var.names)

load(file=file.path(settings$outdir,'sda.obs.Rdata'))


#check

obs.list$obs.mean

udunits2::ud.convert(sum(input[input$year==2010 & input$taxon=='QURU' & input$model == 'Model RW + Census' & input$type == 'ab','value'],na.rm=T)/250,'Mg/ha','kg/m^2')*.48 #600 iterations

matplot(do.call(rbind,obs.list$obs.mean),main='Means')
matplot(do.call(rbind,lapply(obs.list$obs.cov,diag)),main='Variances')

hist(input[input$year==2010 & input$taxon=='QURU' & input$model == 'Model RW + Census' & input$type == 'ab','value'])
hist(input[input$year==1965 & input$taxon=='QURU' & input$model == 'Model RW + Census' & input$type == 'ab','value'])


test <- input[input$year==2010 & input$taxon=='QURU' & input$model == 'Model RW + Census' & input$type == 'ab',]

site1 <- test[test$site_id==1,]
site2 <- test[test$site_id==2,]
site3 <- test[test$site_id==3,]

sum_site1 <- sum_site2 <- sum_site3 <- numeric(250)

for(i in 1:250){
  sum_site1[i] <- sum(site1[site1$iter==i,'value'],na.rm=T)
  sum_site2[i] <- sum(site1[site2$iter==i,'value'],na.rm=T)
  sum_site3[i] <- sum(site1[site3$iter==i,'value'],na.rm=T)
}

mean(c(sum_site1,sum_site2,sum_site3))

quants <- matrix(NA,2010,3)

for(i in 1960:2010){
  quants[i,] <- quantile(input[input$year==i & input$taxon=='QURU' & input$model == 'Model RW + Census' & input$type == 'ab','value'],c(.025,.5,.975),na.rm=T)
}
matplot(quants[1960:2010,])

