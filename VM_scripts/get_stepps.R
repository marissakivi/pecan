library(rgdal) # need to put in assim.sequential
library(ncdf4) # need to put in assim.sequential
library(PEcAn.all)
library(lubridate)

stepps2 <- readRDS('/fs/data2/output/paleon_data_products/STEPPS_predictions_median_sd.RDS')

get_stepps2 <- function(albers_coors, stepps){
  ten.get <- albers_coors
  dist <- fields::rdist(ten.get/10000000,stepps[,c('x','y')])
  stepps_coors_get <- stepps[which.min(dist),c('x','y')]
  stepps_use_x <- stepps[which(stepps[,c('x')]==as.numeric(stepps_coors_get[1])),]
  stepps_use_y <- stepps_use_x[which(stepps_use_x[,c('y')]==as.numeric(stepps_coors_get[2])),]
  
  stepps_use <- stepps_use_y[-grep('OTHER',stepps_use_y$taxon),] #getting rid of 'other' for linkages #would want to put back in for other models?
  
  stepps_obs_mean <- stepps_obs_cov <- list()
  for(t in 1:9){
    stepps_obs_mean[[t]] <- stepps_use[stepps_use$time==seq(200,2100,100)[t],'median']
    stepps_obs_cov[[t]] <- diag(stepps_use[stepps_use$time==seq(200,2100,100)[t],'sd'])
    
    names(stepps_obs_mean[[t]]) <- names(stepps_obs_cov[[t]]) <- stepps_use[stepps_use$time==seq(200,2100,100)[t],'taxon']
  }
  
  names(stepps_obs_mean) <- names(stepps_obs_cov) <- paste0(rev(seq(950,1750,100)),'/12/31')
  
  return(list(stepps_obs_mean=rev(stepps_obs_mean),
              stepps_obs_cov=rev(stepps_obs_cov)))
  
}

ten.get <-c(101608.9, 1247775.3) ### TENSION
obs.list <- get_stepps2(ten.get, stepps2)

#STEPPS1
#/fs/data2/output/paleon_data_products/1000000650/Pollen/fcomp/stepps1_v1.0.nc

get_stepps1 <- function(data.path,settings,time_step){
  d <- settings$database$bety[c("dbname", "password", "host", "user")]
  bety <- src_postgres(host = d$host, user = d$user, password = d$password, dbname = d$dbname)
  site <- PEcAn.DB::query.site(settings$run$site$id, bety$con)
  start_date <- settings$state.data.assimilation$start.date
  end_date   <- settings$state.data.assimilation$end.date
  
  obs.times <- seq(as.Date(start_date), as.Date(end_date), by = 'year')
  obs.times <- obs.times[seq(1,length(obs.times),time_step)]
  obs.times <- formatC(lubridate::year(obs.times), width = 4, format = "d", flag = "0")
  
  
  ### Pollen Data Product (STEPPS)
    ncin <- ncdf4::nc_open(data.path)
    
    coords <- data.frame(x=site$lon,y=site$lat)
    sp::coordinates(coords) <- ~ x + y
    sp::proj4string(coords) <- sp::CRS('+proj=longlat +ellps=WGS84')
    
    ### site utm coordinates
    utm <- sp::spTransform(coords, CRS("+proj=utm +zone=18N ellps=WGS84"))
    utm <- as.matrix(data.frame(utm))
    
    ### find grid cell
    site.x <- which(min(abs(ncvar_get(ncin, 'x') - utm[1])) == abs(ncvar_get(ncin, 'x') - utm[1]))
    site.y <- which(min(abs(ncvar_get(ncin, 'y') - utm[2])) == abs(ncvar_get(ncin, 'y') - utm[2]))
    years <- formatC(ncvar_get(ncin, 'year'), width = 4, format = "d", flag = "0")
    
    taxa <- names(ncin$var)
    if('other'%in%taxa) taxa <- taxa[-c(grep('other',taxa))]
    
    sims.keep <- array(NA,dim=c(length(taxa),length(ncin$dim$year$vals),length(ncin$dim$sample$vals)))
    for(n in seq_along(taxa)){
      taxa.start <- ncvar_get(ncin, taxa[n])
      
      # input is a matrix 'sims', with rows as time and columns as MCMC samples
      sims.keep[n,,] <- taxa.start[site.x,site.y,,]
    }
    
    mean.keep <- list()
    pecan.pfts <- as.character(lapply(settings$pfts, function(x) x[["name"]]))
    
    for(n in taxa){
      sims.start <- ncvar_get(ncin,n)
      # input is a matrix 'sims', with rows as time and columns as MCMC samples
      sims <- sims.start[site.x,site.y,,]
      mean.keep[[n]] <- rowMeans(sims)
      
    }
    
    var.inf <- (time_step/100)
    
    mean.mat <- as.data.frame(mean.keep)
    
    for(n in seq_len(ncol(mean.mat))){
      new.name <- pecan.pfts[grep(taxa[n],pecan.pfts,ignore.case = T)]
      if(any(nchar(new.name))){
        colnames(mean.mat)[n] <- paste0('Fcomp.',new.name)
      }
    }
    
    ##### Ordering to match pecan pft order because last one gets taken in ALR transform
    put_order <- order(colnames(mean.mat))
    mean.mat <- mean.mat[,put_order]
    sims.keep <- sims.keep[put_order,,]
    
    #####
    ##### Calculating Mean and Covariance
    #####
    
    obs.mean <- obs.mean.tmp <- list()
    for(n in 1:nrow(mean.mat)){
      obs.mean[[n]]<- log(mean.mat[n,1:8]/mean.mat[n,9])
    }
    
    sims.alr <- array(NA,dim=c(8,length(ncin$dim$year$vals),250))
    for(ii in 1:250){
      for(tt in 1:length(ncin$dim$year$vals)){
        sims.prop <- sims.keep[1:9,tt,ii]/sum(sims.keep[1:9,tt,ii])
        sims.alr[1:8,tt,ii] <- log(sims.prop[1:8]/sims.prop[9])
      }
    }
    
    alr.look <- sims.alr
    
    names(obs.mean) <- paste0(years,'/12/31')
    
    rownames(sims.alr) <- colnames(mean.mat)[1:8]
    obs.cov <- obs.cov.tmp <- list()
    for(n in 1:length(ncin$dim$year$vals)){
      
      sims.use <- sims.alr[,n,]
      obs.cov[[n]] <- cov(t(sims.use)) * var.inf
    }
    
    names(obs.cov) <- paste0(years,'/12/31')
    
    #### Interpolate over all years
    which.keep <- list()
    
    for(n in obs.times){
      min.vec <- na.omit(as.numeric(n) - year(as.Date(names(obs.mean))))
      which.keep[[n]] <- which(min(abs(min.vec))==abs(min.vec))
      obs.mean.tmp[[n]] <- obs.mean[[which.keep[[n]][1]]]
      obs.cov.tmp[[n]] <- obs.cov[[which.keep[[n]][1]]]
    }
    
    names(obs.mean.tmp)<-paste0(obs.times,'/12/31')
    names(obs.cov.tmp)<-paste0(obs.times,'/12/31')
    
    return(list(obs.mean = obs.mean.tmp, obs.cov = obs.cov.tmp))
    
}

settings <- read.settings('/fs/data2/output/PEcAn_1000008588/pecan.SDA_refactored_pollenonly.xml')

obs.list <- get_stepps1(data.path = c('/fs/data2/output/paleon_data_products/1000000650/Pollen/fcomp/stepps1_v1.0.nc'),
            settings = settings,
            time_step = 50)

corrplot::corrplot(cov2cor(obs.list$obs.cov[[1]]))

inv.alr <<-  nimble::nimbleFunction(
  run = function(alr = double(1)) {
    returnType(double(1))
    
    y = exp(c(alr, 0)) / sum(exp(c(alr, 0)))
    
    return(y)
  })

#check that it's what you think
rbind(sapply(settings$pfts,'[[','name'),
      inv.alr(as.numeric(obs.list$obs.mean[[1]])))

obs.mean <- obs.list$obs.mean


#remember:: different time step
plot(unlist(lapply(obs.mean,function(x)inv.alr(as.numeric(x))[4])))
plot(rev(mean.mat[,3]))


save(obs.list, file='obs.list.50.pollen.Rdata')
