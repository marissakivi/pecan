# SDA Site Data Compilation and Summarization
# 01 November 2019
# Marissa Kivi 
# Function to gather and organize data for SDA runs with LINKAGES

#### Inputs
# settings :: PEcAn SDA settings file for site run
# lat :: site latitude
# lon :: site longitude
# obs.list :: observational data list containing obs.mean and obs.cov from converted npp stat model data 
# init :: site initials for tracking site in data frames

#### Outputs
# param.melt :: df that contains the parameter name and value for all species for each ensemble 
# bias.melt :: df that contains total agb prediction error for each year of each ensemble
# life.melt :: df that contains the growth, death, and birth rate for all species for each year of each ensemble
# run.list :: large list of data frames (one for each species) that tracks the environmental conditions for each yearly spp-level forecast

gather_sda_data = function(settings, lat, lon, obs.list, init){
  
  require(reshape)
  require(mvtnorm)
  
  # set SDA-specific variables from settings file
  nens = as.numeric(settings$ensemble$size)
  nspec = length(settings$pfts)
  site = settings$run$site$name
  sda.years = c(lubridate::year(settings$state.data.assimilation$start.date):lubridate::year(settings$state.data.assimilation$end.date))
  sda.years = sda.years[-1]
  nyr = length(sda.years)
  nparam = 20 # maybe want to make determined 
  no_param = c(1,2,7,16,25,26)
  runs = list.dirs(file.path(settings$outdir,'out'), full.names = FALSE)[-1]
  spp = as.vector(sapply(settings$pfts, function(pft){pft$name}))
  
  # storage for first step
  # param.melt
  param.array = array(NA, dim= c(nspec, nparam, nens))
  param.names = c('DMAX','DMIN','B3','B2','AGEMX','G','SPRTND','SPRTMN','SPRTMX','MPLANT','D3','FROST',
                  'CM1','CM2','CM3','CM4','CM5','FWT','SLTA','SLTB')
  
  # life.melt 
  birth.array = array(NA, dim=c(nens, nyr, nspec))
  growth.array = array(NA, dim =c(nens, nyr, nspec))
  death.array = array(NA, dim=c(nens, nyr, nspec))
  err.agb.pft.array = array(NA, dim = c(nens,nyr,nspec))
  
  # bias.melt 
  tot.bias.mat = matrix(NA,nens,nyr)
  
  # run.list 
  obs = nens * nyr
  run.mat = matrix(NA, obs, 19)
  colnames(run.mat) = c('ens','year','lat','lon','pred','bias','weight','summer.temp','winter.temp','summer.precip','winter.precip',
                        'g.season','basal.area','algf25', 'algf50', 'algf75','sngf','stand.age','dominant')
  run.mat = as.data.frame(run.mat)
  run.mat$lat = rep(lat,obs)
  run.mat$lon = rep(lon,obs)
  run.list = list()
  
  # other storage
  precip.array = array(NA, c(nyr,12,nens))
  temp.array = array(NA, c(nyr,12,nens))
  
  # get weights
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
  
  # contains forecast and analysis matrix
  weights = matrix(0,nens,nyr)
  load(paste0(settings$outdir,'/SDA/sda.output.Rdata')) # met ensemble information
  for (j in 1:nyr){
    weights[,j] = weight(X = as.matrix(FORECAST[[j+1]]), 
                         mu.a = as.matrix(obs.list$obs.mean[[j+1]], nrow = 1),
                         Pa = as.matrix(obs.list$obs.cov[[j+1]]))
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
    temp.array[,,i] = as.matrix(temp.mat[rownames(temp.mat) %in% sda.years,], ncol=12)
    precip.array[,,i] = as.matrix(precip.mat[rownames(precip.mat) %in% sda.years,], ncol=12)
    
    # loop through years (last year of SDA requires special exceptions)
    for (j in 1:nyr){
      
      # find index of observation
      ind = ((i-1)*nyr) + j
      yr = sda.years[j]
      
      run.mat$ens[ind] = paste0(init,i) 
      run.mat$year[ind] = yr
      run.mat$weight[ind] = weights[i,j]
      
      # error 
      tot.bias.mat[i,j] = (sum(FORECAST[[j+1]][i,]) - sum(ANALYSIS[[j+1]][i,])) / sum(ANALYSIS[[j+1]][i,])
      err.agb.pft.array[i,j,] = as.matrix((FORECAST[[j+1]][i,] - ANALYSIS[[j+1]][i,]) / ANALYSIS[[j+1]][i,])
      
      # starting dominant species in stand 
      run.mat$dominant[ind] = spp[which.max(as.matrix((ANALYSIS[[j]][i,])))]
      
      # load same year restart file
      if (j == nyr){
        load(file.path(settings$rundir,cid,'linkages.restart.Rdata'))
      }else{
        load(file.path(settings$rundir,cid,paste0(toString(yr),'-12-31 23:59:59linkages.restart.Rdata')))
      }
      
      run.mat$basal.area[ind] = sum(pi * dbh^2)
      run.mat$stand.age[ind] = mean(iage[iage != 0])
      
      # starting number of trees for calculating birth
      ntrees.res = ntrees
      
      # load output file
      if (j == nyr){
        load(file.path(settings$outdir,'out',cid,'linkages.out.Rdata'))
      }else{
        load(file.path(settings$outdir,'out',cid,paste0(toString(yr),'-12-31 23:59:59linkages.out.Rdata')))
      }
      
      # life process rates for averaging
      birth.array[i,j,] = (ntrees.birth[,1,1] - ntrees.res) / ntrees.res
      growth.array[i,j,] = (agb.pft[,1,1] - as.matrix(ANALYSIS[[j]][i,])) / as.matrix(ANALYSIS[[j]][i,])
      death.array[i,j,] = (ntrees.kill[,1,1] - ntrees.birth[,1,1]) / ntrees.birth[,1,1]
      
      # some climate variables to characterize model run 
      run.mat$g.season[ind] = length(which(temp.array[j,,i]>10))*30
      run.mat$summer.temp[ind] = mean(temp.array[j,(6:8),i])
      run.mat$winter.temp[ind] = mean(temp.array[j,c(1,11,12),i])
      run.mat$summer.precip[ind] = mean(precip.array[j,(6:8),i])
      run.mat$winter.precip[ind] = mean(precip.array[j,c(1,11,12),i])
      
      # now, collect specific for each species growth and add to run.list
      for (k in 1:nspec){
        
        # add information collected until this point
        if (i == 1 & j == 1){ # need to add a whole new part of the list the first iteration 
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
  param.melt = melt(param.array)
  colnames(param.melt) = c('species','p.name', 'ensemble','p.value')
  for (k in 1:nspec){
    param.melt$species[which(param.melt$species == k)] = spp[k]
  }
  param.melt$ensemble = paste0(init,param.melt$ensemble)
  param.melt$site = rep(site,length(param.melt$species))
  
  birth.melt = melt(birth.array)
  colnames(birth.melt) = c('ensemble','year','species','birth')
  growth.melt = melt(growth.array)
  colnames(growth.melt) = c('ensemble','year','species','growth')
  death.melt = melt(death.array)
  colnames(death.melt) = c('ensemble','year','species','death')
  error.melt = melt(err.agb.pft.array)
  colnames(error.melt) = c('ensemble','year','species','error')
  life.melt = full_join(error.melt, birth.melt, id = c('ensemble','year','species')) %>%
    full_join(growth.melt, id = c('ensemble','year','species')) %>%
    full_join(death.melt, id = c('ensemble','year','species'))
  for (k in 1:nspec){
    life.melt$species[which(life.melt$species == k)] = spp[k]
  }
  life.melt$ensemble = paste0(init,life.melt$ensemble)
  life.melt$site = rep(site, length(life.melt$year))
  
  bias.melt = melt(tot.bias.mat)
  colnames(bias.melt) = c('ensemble','year','bias')
  bias.melt$ensemble = paste0(init,bias.melt$ensemble)
  bias.melt$site = rep(site,length(bias.melt$year))
  
  return(list(param.melt=param.melt, bias.melt=bias.melt, life.melt=life.melt, run.list=run.list))
}

#save(param.melt=param.melt,
#          bias.melt=bias.melt,
#          life.melt=life.melt,
#          run.list=run.list,
#     file = '~/VM_scripts/test_data.Rdata')

