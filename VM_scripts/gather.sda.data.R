# SDA Site Data Compilation and Summarization
# 01 November 2019
# Marissa Kivi 
# Function to gather and organize data for SDA runs with LINKAGES

#### Inputs
# id :: PEcAn workflow ID
# obs.list :: observational data list containing obs.mean and obs.cov used in SDA workflow
# init :: site initials for tracking site in data frames

#### Outputs
# param.melt :: dataframe containing species name, parameter name, and parameter value for ever year of each ensemble
# growth.melt :: dataframe containing total and PFT-level predicted relative annual growth for each year of each ensemble
# life.melt :: dataframe containing birth, growth, and death rates for each species of each year of each ensemble 
# bias.melt :: dataframe containing total and PFT-level relative model bias for predicting growth
# weight.melt :: dataframe containing the annual ensemble weights for all years and ensembles 
# run.list :: list of dataframes that track the environmental variables, PFT-level growth, and bias for all yearly 
#     forecasts across all ensembles
# spp :: list of species in correct order for model 


gather_sda_data = function(id, obs.list, disturbance.yrs, init){
  
  require(reshape)
  require(mvtnorm)
  require(PEcAn.settings)
  
  # load workflow settings
  settings = read.settings(paste0('/data/workflows/PEcAn_',toString(id),'/pecan.SDA.xml'))
  
  # get SDA-specific variables from settings file
  nens = as.numeric(settings$ensemble$size)
  nspec = length(settings$pfts)
  site = settings$run$site$name
  sda.years = c(lubridate::year(settings$state.data.assimilation$start.date):lubridate::year(settings$state.data.assimilation$end.date))
  sda.years = sda.years[-1] # remove first year because the bias is generally so large
  nyr = length(sda.years)
  pfts = as.vector(sapply(settings$pfts, function(pft){pft$name}))
  spp = sapply(pfts, function(x){stringr::word(x,2,sep="_")})
  names(spp) = NULL
  lat = as.numeric(settings$run$site$lat)
  lon = as.numeric(settings$run$site$lon)
  
  # organize observed data 
  obs.mean = obs.list$obs.mean
  obs.cov = obs.list$obs.cov
  obs.mean <- obs.mean[sapply(lubridate::year(names(obs.mean)), function(obs.year) obs.year %in% (sda.years))]
  obs.cov <- obs.list$obs.cov[sapply(lubridate::year(names(obs.cov)), function(obs.year) obs.year %in% (sda.years))]
  
  # param.melt 
  param.names = c('DMAX','DMIN','B3','B2','AGEMX','G','SPRTND','SPRTMN','SPRTMX','MPLANT','D3','FROST',
                  'CM1','CM2','CM3','CM4','CM5','FWT','SLTA','SLTB')
  nparam = length(param.names)
  no_param = c(1,2,7,16,25,26)
  param.array = array(NA, dim= c(nspec, nparam, nens))
  
  # life.melt
  birth.array = array(NA, dim=c(nens, nyr, nspec))
  growth.array = array(NA, dim =c(nens, nyr, nspec))
  death.array = array(NA, dim=c(nens, nyr, nspec))

  # bias.melt 
  tot.bias.mat = matrix(NA,nens,nyr)
  err.agb.pft.array = array(NA, dim = c(nens,nyr,nspec))
  
  # growth.melt
  tot.growth.mat = matrix(NA,nens,nyr)

  # run.list
  obs = nens * nyr
  run.mat = matrix(NA, obs, 20)
  colnames(run.mat) = c('ens','year','lat','lon','pred','bias','weight','summer.temp','winter.temp','summer.precip','winter.precip',
                        'g.season','basal.area','algf25', 'algf50', 'algf75','sngf','stand.age','dominant','disturb')
  run.mat = as.data.frame(run.mat)
  run.mat$lat = rep(lat,obs)
  run.mat$lon = rep(lon,obs)
  run.list = list()

  # weight.melt
  # function for calculating weight of each ensemble forecast
  weight = function(X,mu.a,Pa){

    # calculate the likelihood of the ensemble members given mu.a and Pa for one given year
    nens <- nrow(X)
    wt.mat <- matrix(NA,nrow=nens)
    for (i in 1:nens){
      wt.mat[i] = dmvnorm(X[i,], mean=mu.a, sigma=Pa)
    }

    # if all predictions are insane, then give equal weight
    if (sum(wt.mat)<=0){
      return(rep(1/nens,nens))
    }

    # put into weights table
    wt.props <- t(prop.table(wt.mat))
    return(wt.props)
  }

 weights = matrix(0,nens,nyr)
 load(paste0(settings$outdir,'/SDA/sda.output.Rdata')) # met ensemble information
 for (j in 1:nyr){
   weights[,j] = weight(X = as.matrix(FORECAST[[j+1]]),
                        mu.a = as.matrix(obs.mean[[j]], nrow = 1),
                        Pa = as.matrix(obs.cov[[j]]))
 }
 
 load(paste0(settings$outdir,'/SDA/outconfig.Rdata')) # met ensemble information
 ids = list.dirs(settings$rundir, full.names = FALSE)[-1]

 # loop through ensembles
 for (i in 1:nens){

    cid = ids[(i+nens)]

    # load input
    load(file.path(settings$rundir,cid,'linkages.input.Rdata'))
    param.array[,,i] = as.matrix(spp.params[,-c(no_param)])

    # obtain met ensemble data
    load(outconfig$samples$met$samples[[i]])
    temp.array = as.matrix(temp.mat[rownames(temp.mat) %in% sda.years,], ncol=12)
    precip.array = as.matrix(precip.mat[rownames(precip.mat) %in% sda.years,], ncol=12)

    # loop through years (last year of SDA requires special exceptions)
    for (j in 1:nyr){

      # find index of observation for run.list
      ind = ((i-1)*nyr) + j
      yr = sda.years[j]
      disturb.update = yr - disturbance.yrs

      # collect error information
      tot.bias.mat[i,j] = (sum(FORECAST[[j+1]][i,]) - sum(ANALYSIS[[j+1]][i,])) / sum(ANALYSIS[[j+1]][i,])
      err.agb.pft.array[i,j,] = as.matrix((FORECAST[[j+1]][i,] - ANALYSIS[[j+1]][i,]) / ANALYSIS[[j+1]][i,])

      # collect ensemble-level run.list data 
      run.mat$ens[ind] = paste0(init,i)
      run.mat$year[ind] = yr
      run.mat$weight[ind] = weights[i,j]
      run.mat$dominant[ind] = spp[which.max(as.matrix((ANALYSIS[[j]][i,])))] # starting dominant species in stand
      run.mat$disturb[ind] = min(disturb.update[disturb.update >= 0]) # years since most recent disturbance 

      # load same year restart file, last year doesn't get renamed with year information
      if (j == nyr){
        load(file.path(settings$rundir,cid,'linkages.restart.Rdata'))
      }else{
        load(file.path(settings$rundir,cid,paste0(toString(yr),'-12-31 23:59:59linkages.restart.Rdata')))
      }

      run.mat$basal.area[ind] = sum(pi * dbh^2)
      run.mat$stand.age[ind] = quantile(iage[iage != 0], 0.5)

      # starting number of trees for calculating birth
      ntrees.res = ntrees

      # load output file
      if (j == nyr){
        load(file.path(settings$outdir,'out',cid,'linkages.out.Rdata'))
      }else{
        load(file.path(settings$outdir,'out',cid,paste0(toString(yr),'-12-31 23:59:59linkages.out.Rdata')))
      }

      # total growth
      tot.growth.mat[i,j] = (sum(agb.pft[,1,1]) - sum(ANALYSIS[[j]][i,])) / sum(ANALYSIS[[j]][i,])

      # species-level life process information  
      birth.array[i,j,] = ntrees.birth[,1,1] - ntrees.res
      for (entry in seq_along(growth.array[i,j,])){
        if (ANALYSIS[[j]][i,entry] == 0 & agb.pft[entry,1,1] != 0) growth.array[i,j,entry] <- NA
        if (ANALYSIS[[j]][i,entry] == 0 & agb.pft[entry,1,1] == 0) growth.array[i,j,entry] <- 0
        if (ANALYSIS[[j]][i,entry] != 0) growth.array[i,j,entry] = (agb.pft[entry,1,1] - ANALYSIS[[j]][i,entry]) / ANALYSIS[[j]][i,entry]
      }
      death.array[i,j,] = (ntrees.kill[,1,1] - ntrees.birth[,1,1]) / ntrees.birth[,1,1]

      # some climate variables to characterize model run
      run.mat$g.season[ind] = length(which(temp.array[j,]>10))*30
      run.mat$summer.temp[ind] = mean(temp.array[j,(6:8)])
      run.mat$winter.temp[ind] = mean(temp.array[j,c(1,11,12)])
      run.mat$summer.precip[ind] = mean(precip.array[j,(6:8)])
      run.mat$winter.precip[ind] = mean(precip.array[j,c(1,11,12)])

      # now, collect specific for each species growth and add to run.list
      for (k in 1:nspec){

        # add information collected until this point
        if (ind == 1){ # need to add a whole new part of the list the first iteration
          run.list[[k]] = run.mat
        }else{ # otherwise just add the newest row of information
          run.list[[k]][ind,] = run.mat[ind,]
        }

        # add growth factor information to run.mat
        run.list[[k]]$sngf[ind] = gf.vec.save[k,3,1,1]
        run.list[[k]]$algf25[ind] = quantile(algf.save.keep[,k,1,1], 0.25, na.rm = T)
        run.list[[k]]$algf50[ind] = quantile(algf.save.keep[,k,1,1], 0.50, na.rm = T)
        run.list[[k]]$algf75[ind] = quantile(algf.save.keep[,k,1,1], 0.75, na.rm = T)

        # add biomass information (predicted biomass relative increase and prediction error)
        run.list[[k]]$pred[ind] = growth.array[i,j,k]
        run.list[[k]]$bias[ind] = err.agb.pft.array[i,j,k]
      }
    }
  }

  # melt arrays into dataframes
 
  # arranging param.melt
  param.melt = melt(param.array)
  colnames(param.melt) = c('species','p.name', 'ensemble','p.value')
  param.melt$species = plyr::mapvalues(param.melt$species, from = c(1:nspec), to = spp)
  param.melt$p.name = plyr::mapvalues(param.melt$p.name, from = c(1:nparam), to = param.names)
  param.melt$ensemble = paste0(init,param.melt$ensemble)
  param.melt$site = rep(site, length(param.melt$species))

  # arranging life.melt 
  pft.growth.melt = melt(growth.array)
  colnames(pft.growth.melt) = c('ensemble','year','species','growth')
  birth.melt = melt(birth.array)
  colnames(birth.melt) = c('ensemble','year','species','birth')
  death.melt = melt(death.array)
  colnames(death.melt) = c('ensemble','year','species','death')
  
  life.melt = full_join(birth.melt, pft.growth.melt, id = c('ensemble','year','species')) %>%
    full_join(death.melt, id = c('ensemble','year','species'))
  life.melt$species = plyr::mapvalues(life.melt$species, from = c(1:nspec), to = spp)
  life.melt$ensemble = paste0(init,life.melt$ensemble)
  life.melt$site = rep(site, length(life.melt$year))
  
  # arranging growth.melt 
  growth.melt = melt(tot.growth.mat)
  colnames(growth.melt) = c('ensemble','year','total')
  pft.growth.melt$species = plyr::mapvalues(pft.growth.melt$species, from = c(1:nspec), to = spp)
  growth.melt = left_join(growth.melt, 
                          dcast(pft.growth.melt, ensemble + year ~ species, value.var = 'growth'), 
                          by = c('ensemble','year'))
  growth.melt$ensemble = paste0(init,growth.melt$ensemble)
  growth.melt$site = rep(site, length(growth.melt$year))
  
  # arranging bias.melt 
  bias.melt = melt(tot.bias.mat)
  colnames(bias.melt) = c('ensemble','year','total')
  error.melt = melt(err.agb.pft.array)
  colnames(error.melt) = c('ensemble','year','species','error')
  error.melt$species = plyr::mapvalues(error.melt$species, from = c(1:nspec), to = spp)
  bias.melt = left_join(bias.melt, 
                        dcast(error.melt, ensemble + year ~ species, value.var = 'error'),
                        by = c('ensemble','year'))
  bias.melt$ensemble = paste0(init,bias.melt$ensemble)
  bias.melt$site = rep(site,length(bias.melt$year))
  
  # arranging weight.melt 
  weight.melt = melt(weights)
  colnames(weight.melt) = c('ensemble','year','weight')
  weight.melt$ensemble = paste0(init,weight.melt$ensemble)
  weight.melt$site = rep(site,length(weight.melt$year))
  
  return(list(param.melt = param.melt, life.melt = life.melt, growth.melt = growth.melt,
         weight.melt = weight.melt, bias.melt = bias.melt, run.list = run.list, spp = spp))

}

