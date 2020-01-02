
library(plyr)
library(PEcAn.settings)
library(PEcAn.DB)
library(PEcAn.data.land)
library(dplyr)

# ON DESKTOP:
# Step 1: Load and prepare data
input <- readRDS('~/Downloads/NPP_STAT_MODEL_RH.RDS')

# change for desired species only to make smaller data frame
# make sure that the codes are also recognized by LINKAGES
# for example, LINKAGES recognizes species code of red spruce as PIRU not PCRU
unique(input$taxon)
species <- c('ACRU','QURU','PIRU','PIST') 

# ON VM: 
# adjust species list to be for species at site
setwd('~/data_files')
species <- c('ACRU','QURU','PIRU','PIST')

# make bety connection
bety <- list(user = 'bety',
             password = 'bety',
             host = 'postgres',
             dbname = 'bety',
             driver = "PostgreSQL",
             write = TRUE)
dbcon <- db.open(bety)

# Step 2: Obtain bety pft and species information
spp_id <- match_species_id(species,format_name = 'usda',dbcon)
pft_mat = matrix(NA, length(spp_id$genus), 4)
pft_mat = as.data.frame(pft_mat)
names(pft_mat) = c('bety_pft_id', 'pft','bety_species_id','latin')
i = 0

# the following matches the species to a pft in BETY 
# upon completion, check to make sure the pfts are correct 
for(sppID in spp_id$bety_species_id){
  i = i + 1
  genus = spp_id$genus[i]
  spec = spp_id$species[i]
  pft.now <- PEcAn.DB::db.query(query = paste("SELECT pfts.id, pfts.name, species.id FROM pfts",
                                              "JOIN pfts_species ON pfts_species.pft_id = pfts.id",
                                              "JOIN species ON species.id = pfts_species.specie_id",
                                              "WHERE species.id =",
                                              sppID), 
                                con = dbcon)
  genus_check = sapply(pft.now$name, function(x){grepl(tolower(genus),tolower(x))})
  species_check = sapply(pft.now$name, function(x){grepl(tolower(spec),tolower(x))})
  if (!any(genus_check) & !any(species_check)){
    print('No pft available that reasonably matches!')
  }else{
    pft.now = pft.now[(genus_check & species_check),]
  }
  
  pft_mat[i,1] = pft.now$id
  pft_mat[i,2] = pft.now$name
  pft_mat[i,3] = sppID
  pft_mat[i,4] = paste(genus,spec)
}
x <- paste0('AGB.pft.', pft_mat$pft)
names(x) <- spp_id$input_code
save(x, file = 'rh_obs_VM.Rdata')

# BACK TO DESKTOP: 
load('~/Downloads/rh_obs_VM.Rdata')

# Step 3: Correct formatting of input 

# check here to make sure inputs from RDS file are correct! (e.g. check species codes, capitalization, etc.)
#levels(input$taxon)[3] <- 'BEAL2' #HF
#levels(input$taxon)[8] <- 'FRAM2' #HF
#levels(input$taxon)[4] <- 'PIGL'    #RH
levels(input$taxon)[which(levels(input$taxon) == 'PCRU')] = 'PIRU' #RH

# adjust "type" variable for consistency across sites
levels(input$type)
levels(input$type)[1] <- 'ab'
levels(input$type)[2] <- 'abi'

# condense input to only include species we are interested in
input = input %>% filter(taxon %in% species, type == 'ab', !is.na(value))
unique(input$taxon)input$pft.cat <- x[as.character(input$taxon)]
input$pft.cat <- x[as.character(input$taxon)]

# set up variables for aggregation step
start_date = min(input$year)
end_date = max(input$year)
obs.mean <- obs.mean.tmp <- list()
obs.cov <- obs.cov.tmp <- list()
obs.times <- c(start_date:end_date)
time.type <- 'year'
i=1
biomass2carbon <- 0.48
var.names = 'AGB.pft'

# convert from Mg/ha to kgC/m2
if(!is.null(input$value)) input$value <- udunits2::ud.convert(input$value,'Mg/ha','kg/m^2') * biomass2carbon #* kgm2Mgha 

# check to make sure no NAs persist
input <- input[input$pft.cat!='AGB.pft.NA',]

variable <- c('value') #'ab'
arguments <- list(.(year, iter,  pft.cat), .(variable))
arguments2 <- list(.(year, pft.cat), .(variable))
arguments3 <- list(.(iter), .(pft.cat, variable), .(year))

# melting and aggregating
melt_id <- colnames(input)[-which(colnames(input) %in% variable)]
melt.test <- reshape2::melt(input, id = melt_id, na.rm = TRUE)
cast.test <- reshape2::dcast(melt.test, arguments, sum, margins = variable)

melt_id_next <- colnames(cast.test)[-which(colnames(cast.test) %in% variable)]
melt.next <- reshape2::melt(cast.test, id = melt_id_next)
mean_mat <<- reshape2::dcast(melt.next, arguments2, mean)

iter_mat <<- reshape2::acast(melt.next, arguments3, mean)
cov.test <- apply(iter_mat,3,function(x){cov(x)})

for(t in seq_along(obs.times)){
  obs.mean.tmp[[t]] <- mean_mat[mean_mat[,time.type]==obs.times[t], -c(1)] #THIS WONT WORK IF TIMESTEP ISNT ANNUAL
  
  if(any(var.names == 'AGB.pft')){
    obs.mean.tmp[[t]] <- rep(NA, length(unique(input$pft.cat)))
    names(obs.mean.tmp[[t]]) <- sort(unique(input$pft.cat))
    for(r in seq_along(unique(input$pft.cat))){
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
save(obs.list,file=file.path(settings$outdir,'sda.obs.RH.Rdata'))

######################################################

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

