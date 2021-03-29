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


gather_sda_data = function(id, obs.list, init){
  
  require(reshape)
  require(mvtnorm)
  require(PEcAn.settings)
  
  # load workflow settings
  settings = read.settings(paste0('/save/workflows/PEcAn_',toString(id),'/pecan.SDA.xml'))
  outdir = paste0('/save/workflows/PEcAn_',toString(id),'/')
  
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
  
  # pft-level life rates
  birth.array = array(NA, dim=c(nens, nyr, nspec))
  death.array = array(NA, dim=c(nens, nyr, nspec))
  
  # pft-level growth rates
  growth.array = array(NA, dim =c(nens, nyr, nspec))
  adj.growth.array = array(NA, dim=c(nens,nyr,nspec))
  err.agb.pft.array = array(NA, dim=c(nens,nyr,nspec))
  
  # total growth rates
  tot.growth.mat = matrix(NA,nens,nyr)
  adj.growth.mat = matrix(NA,nens,nyr)
  tot.bias.mat = matrix(NA, nens, nyr)
  
  # lowest growth factor tracking
  gf.mat = matrix(NA, 1, 6)
  colnames(gf.mat) = c('ens','year','tree','species','value','ind')
  gf.mat = as.data.frame(gf.mat)
  
  # stand structure 
  dbh.skew.array = array(NA, dim = c(nens,nyr,nspec))
  dbh.kurt.array = array(NA, dim = c(nens,nyr,nspec))
  dbh.med.array = array(NA, dim = c(nens,nyr,nspec))
  ht.skew.array = array(NA, dim = c(nens,nyr,nspec))
  ht.kurt.array = array(NA, dim = c(nens,nyr,nspec))
  ht.med.array = array(NA, dim = c(nens,nyr,nspec))

  # run.list
  obs = nens * nyr
  run.mat = matrix(NA, obs, 16)
  colnames(run.mat) = c('ens','year','pred','adj','bias','weight','summer.temp','winter.temp','summer.precip','winter.precip',
                        'med.dbh','skew.dbh','kurt.dbh','contrib','avln','gap')
  run.mat = as.data.frame(run.mat)
  run.list = list()
  
  # weight tracking 
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
  load(paste0(outdir,'/SDA/sda.output.Rdata')) # met ensemble information
  for (j in 1:nyr){
    weights[,j] = weight(X = as.matrix(FORECAST[[j+1]]),
                         mu.a = as.matrix(obs.mean[[j]], nrow = 1),
                         Pa = as.matrix(obs.cov[[j]]))
  }
  
  load(paste0(outdir,'/SDA/outconfig.Rdata')) # met ensemble information
  ids = list.dirs(file.path(outdir,'run'), full.names = FALSE)[-1]
  
  # loop through ensembles
  for (i in 1:nens){
    
    print(i)
  
    # get current ID 
    cid = ids[(i+nens)]
    
    # load input
    load(file.path(outdir,'run',cid,'linkages.input.Rdata'))
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
      
      # collect ensemble-level run.list data 
      run.mat$ens[ind] = paste0(init,i)
      run.mat$year[ind] = yr
      run.mat$weight[ind] = weights[i,j]
  
      # error 
      tot.bias.mat[i,j] = (sum(FORECAST[[j+1]][i,]) - sum(ANALYSIS[[j+1]][i,])) / sum(ANALYSIS[[j+1]][i,])
      err.agb.pft.array[i,j,] = as.matrix((FORECAST[[j+1]][i,] - ANALYSIS[[j+1]][i,]) / ANALYSIS[[j+1]][i,])
      
      # load same year restart file, last year doesn't get renamed with year information
      if (j == nyr){
        load(file.path(outdir,'run',cid,'linkages.restart.Rdata'))
      }else{
        load(file.path(outdir,'run',cid,paste0(toString(yr),'-12-31 23:59:59linkages.restart.Rdata')))
      }
      
      # median diameter in stand 
      run.mat$med.dbh[ind] = quantile(dbh[dbh > 0], 0.5)
      run.mat$skew.dbh[ind] = e1071::skewness(dbh[dbh>0])
      run.mat$kurt.dbh[ind] = e1071::kurtosis(dbh[dbh>0])
      
      # starting number of trees for calculating birth
      ntrees.res = ntrees
      
      # load output file
      if (j == nyr){
        load(file.path(outdir,'out',cid,'linkages.out.Rdata'))
      }else{
        load(file.path(outdir,'out',cid,paste0(toString(yr),'-12-31 23:59:59linkages.out.Rdata')))
      }
      
      # total growth
      tot.growth.mat[i,j] = (sum(agb.pft[,1,1]) - sum(ANALYSIS[[j]][i,])) / sum(ANALYSIS[[j]][i,])
      adj.growth.mat[i,j] = (sum(ANALYSIS[[j+1]][i,]) - sum(ANALYSIS[[j]][i,])) / sum(ANALYSIS[[j]][i,])
      
      # species-level life process information  
      birth.array[i,j,] = ntrees.birth[,1,1] - ntrees.res
      for (entry in seq_along(growth.array[i,j,])){
        if (ANALYSIS[[j]][i,entry] == 0 & agb.pft[entry,1,1] != 0) growth.array[i,j,entry] <- NA
        if (ANALYSIS[[j]][i,entry] == 0 & agb.pft[entry,1,1] == 0) growth.array[i,j,entry] <- 0
        if (ANALYSIS[[j]][i,entry] != 0) growth.array[i,j,entry] = (agb.pft[entry,1,1] - ANALYSIS[[j]][i,entry]) / ANALYSIS[[j]][i,entry]
      }
      # ... and adjusted
      for (entry in seq_along(adj.growth.array[i,j,])){
        if (ANALYSIS[[j]][i,entry] == 0 & ANALYSIS[[j+1]][i,entry] != 0) adj.growth.array[i,j,entry] <- NA
        if (ANALYSIS[[j]][i,entry] == 0 & ANALYSIS[[j+1]][i,entry] == 0) adj.growth.array[i,j,entry] <- 0
        if (ANALYSIS[[j]][i,entry] != 0) adj.growth.array[i,j,entry] = (ANALYSIS[[j+1]][i,entry] - ANALYSIS[[j]][i,entry]) / ANALYSIS[[j]][i,entry]
      }
      death.array[i,j,] = (ntrees.kill[,1,1] - ntrees.birth[,1,1]) / ntrees.birth[,1,1]
      
      # some climate variables to characterize model run
      #run.mat$g.season[ind] = length(which(temp.array[j,]>10))*30
      run.mat$summer.temp[ind] = mean(temp.array[j,(6:8)])
      run.mat$winter.temp[ind] = mean(temp.array[j,c(1:4)])
      run.mat$summer.precip[ind] = mean(precip.array[j,(6:8)])
      run.mat$winter.precip[ind] = mean(precip.array[j,c(1:4)])
      run.mat$avln[ind] = avln[,1]
      
      # set up lgf tracking 
      lgf.mat = matrix(NA, sum(ntrees.kill), 6)
      colnames(lgf.mat) = c('ens','year','tree','species','value','ind')
      lgf.mat = as.data.frame(lgf.mat)
      spp.ind = c()
      # gf.vec.save:: algf, smgf, sngf, degdgf
      for (k in 1:nspec){
        spp.ind = c(spp.ind, rep(k,ntrees.kill[k]))
      }
      lgf.mat$species = spp.ind
      lgf.mat$tree = c(1:sum(ntrees.kill))
      lgf.mat$ens = rep(paste0(init,i), sum(ntrees.kill))
      lgf.mat$year = rep(yr, sum(ntrees.kill))
      
      for (tree in seq_along(lgf.mat$ens)){
        this.spp = lgf.mat$species[tree]
        lgf.mat$ind[tree] = which.min(c(algf.save.keep[tree,this.spp,1,1],gf.vec.save[this.spp,2:4,1,1]))
        lgf.mat$value[tree] = min(c(algf.save.keep[tree,this.spp,1,1],gf.vec.save[this.spp,2:4,1,1]))
      }
      gf.mat = rbind(gf.mat, lgf.mat)
      
      # calculation of light available to gap-height trees 
      dbh.melt = as.data.frame(dbh[dbh > 0])
      names(dbh.melt) = c('dbh')
      dbh.melt$age = iage[iage > 0]
      spp.ind = c()
      for (k in 1:nspec){
        spp.ind = c(spp.ind, rep(k,ntrees[k]))
      }
      dbh.melt$species = spp.ind
      
      # collect parameters for calculating canopy LAI 
      param.temp = data.frame(species = c(1:nspec), 
                              b2 = spp.params$B2,
                              b3 = spp.params$B3,
                              sltb = spp.params$SLTB,
                              slta = spp.params$SLTA,
                              fwt = spp.params$FWT,
                              frt = spp.params$FRT)
      canopy.melt = left_join(dbh.melt, param.temp, by = c('species'))
      
      # perform calculation of leaf biomass and, then, % of full sunlight for each tree 
      # height in this equation is calculated to the nearest 0.1 meters 
      canopy.melt$ht = ((canopy.melt$b2*canopy.melt$dbh-canopy.melt$b3*canopy.melt$dbh^2)/10)
      canopy.melt$frt[canopy.melt$age < canopy.melt$frt] = canopy.melt$age[canopy.melt$age < canopy.melt$frt]
      canopy.melt$leaf = ((((canopy.melt$slta + canopy.melt$sltb * canopy.melt$dbh) / 2) ^ 2) * 3.14 * canopy.melt$fwt * canopy.melt$frt)
      all.leaf = sum(canopy.melt$leaf)
      
      slar = vector()
      for (tree in seq_along(canopy.melt$dbh)){
        this.ht = canopy.melt$ht[tree]
        slar[tree] = ifelse(length(which(canopy.melt$ht > this.ht)) > 0, sum(canopy.melt %>% filter(ht > this.ht) %>% dplyr::select(leaf)), 0)
      }
      canopy.melt$slar = slar
      canopy.melt = canopy.melt %>% mutate(al = exp(-slar/93750))
      
      # calculate average available light to all trees of height less than or equal to 2.5 meters
      run.mat$gap[ind] = mean((canopy.melt %>% filter(ht <= 25) %>% dplyr::select(al))$al)
      if (is.na(run.mat$gap[ind])) run.mat$gap[ind] = exp(-all.leaf/93750)
      
      # convert leaf biomass to kg/tree, calculate biomass contributions
      canopy.melt = canopy.melt %>% mutate(leaf = leaf * 0.001) %>%
        mutate(biomass = (0.1193 * dbh^2.393) + leaf) %>% 
        arrange(desc(biomass))
      run.mat$contrib[ind] = round((min(which((cumsum(canopy.melt$biomass)/sum(canopy.melt$biomass)) > 0.5)) / length(canopy.melt$biomass)) * 100, 2)
      
      # now, collect specific for each species growth and add to run.list
      for (k in 1:nspec){
        
        # collect stand structure information
        dbh.skew.array[i,j,k] = skewness((canopy.melt %>% filter(species == k))$dbh)
        dbh.kurt.array[i,j,k] = kurtosis((canopy.melt %>% filter(species == k))$dbh)
        dbh.med.array[i,j,k] = median((canopy.melt %>% filter(species == k))$dbh)
        ht.skew.array[i,j,k] = skewness((canopy.melt %>% filter(species == k))$ht)
        ht.kurt.array[i,j,k] = kurtosis((canopy.melt %>% filter(species == k))$ht)
        ht.med.array[i,j,k] = median((canopy.melt %>% filter(species == k))$ht)
        
        # add information collected until this point
        if (ind == 1){ # need to add a whole new part of the list the first iteration
          run.list[[k]] = run.mat
        }else{ # otherwise just add the newest row of information
          run.list[[k]][ind,] = run.mat[ind,]
        }
        
        # add biomass information (predicted biomass relative increase and prediction error)
        run.list[[k]]$pred[ind] = growth.array[i,j,k]
        run.list[[k]]$adj[ind] = adj.growth.array[i,j,k]
        run.list[[k]]$bias[ind] = as.matrix((FORECAST[[j+1]][i,k] - ANALYSIS[[j+1]][i,k]) / ANALYSIS[[j+1]][i,k])
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
  
  # arranging life.melt 
  pft.growth.melt = melt(growth.array)
  colnames(pft.growth.melt) = c('ensemble','year','species','growth')
  birth.melt = melt(birth.array)
  colnames(birth.melt) = c('ensemble','year','species','birth')
  death.melt = melt(death.array)
  colnames(death.melt) = c('ensemble','year','species','death')
  adj.pft.growth.melt = melt(adj.growth.array)
  colnames(adj.pft.growth.melt) = c('ensemble','year','species','adj.growth')
  life.melt = full_join(birth.melt, pft.growth.melt, id = c('ensemble','year','species')) %>%
    full_join(death.melt, id = c('ensemble','year','species')) %>% 
    full_join(adj.pft.growth.melt, id = c('ensemble','year','species'))
  life.melt$species = plyr::mapvalues(life.melt$species, from = c(1:nspec), to = spp)
  life.melt$ensemble = paste0(init,life.melt$ensemble)
  
  # arranging growth.melt 
  growth.melt = melt(tot.growth.mat)
  colnames(growth.melt) = c('ensemble','year','total')
  adj.growth.melt = melt(adj.growth.mat)
  colnames(adj.growth.melt) = c('ensemble','year', 'adj.total')
  
  #pft.growth.melt$species = plyr::mapvalues(pft.growth.melt$species, from = c(1:nspec), to = spp)
  growth.melt = full_join(growth.melt, 
                          adj.growth.melt, 
                          by = c('ensemble','year'))
  growth.melt$ensemble = paste0(init,growth.melt$ensemble)
  
  # arranging weight.melt 
  weight.melt = melt(weights)
  colnames(weight.melt) = c('ensemble','year','weight')
  weight.melt$ensemble = paste0(init,weight.melt$ensemble)

  # arranging bias melt
  bias.melt = melt(tot.bias.mat)
  colnames(bias.melt) = c('ensemble','year','bias')
  bias.melt$ensemble = paste0(init,bias.melt$ensemble)
  bias.melt$site = rep(site,length(bias.melt$year))
  
  # arranging dbh and ht distribution information
  dist.melt = melt(dbh.kurt.array)
  colnames(dist.melt) = c('ensemble','year','species','dbh.kurtosis')
  temp1.melt = melt(dbh.med.array)
  colnames(temp1.melt) = c('ensemble','year','species','dbh.median')
  temp2.melt = melt(dbh.skew.array)
  colnames(temp2.melt) = c('ensemble','year','species','dbh.skew')
  temp3.melt = melt(ht.kurt.array)
  colnames(temp3.melt) = c('ensemble','year','species','ht.kurtosis')
  temp4.melt = melt(ht.med.array)
  colnames(temp4.melt) = c('ensemble','year','species','ht.median')
  temp5.melt = melt(ht.skew.array)
  colnames(temp5.melt) = c('ensemble','year','species','ht.skew')
  dist.melt = left_join(dist.melt, temp1.melt, by = c('ensemble','year','species')) %>% 
    left_join(temp2.melt, by = c('ensemble','year','species')) %>%
    left_join(temp3.melt, by = c('ensemble','year','species')) %>%
    left_join(temp4.melt, by = c('ensemble','year','species')) %>% 
    left_join(temp5.melt, by = c('ensemble','year','species'))
  
  # arranging pft-level error melt 
  error.melt = melt(err.agb.pft.array)
  colnames(error.melt) = c('ensemble','year','species','error')
  
  # remove first entry in gf.mat 
  gf.mat = gf.mat[-1,]
  
  return(list(param.melt = param.melt, 
              life.melt = life.melt, 
              growth.melt = growth.melt,
              weight.melt = weight.melt, 
              bias.melt = bias.melt, 
              error.melt = error.melt, 
              run.list = run.list, 
              spp = spp, 
              gf.mat = gf.mat,
              dist.melt = dist.melt))
}

