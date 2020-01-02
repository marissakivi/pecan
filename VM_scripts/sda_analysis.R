## SDA Analysis Script

# What years do analysis and forecast correspond to ???

# go through and find which variables we need to keep after collection step... we don't need everything

##########################################
## Step 1: Prepare environment and data ##
##########################################

rm(list=ls())
library(mvtnorm)
library(PEcAn.workflow)
library(PEcAn.settings)
library(ggplot2)
library(gridExtra)
library(boot)
library(reshape)
library(dplyr)
workflowID = '14000000018'
setwd(paste0('/data/workflows/PEcAn_',workflowID))

# read in settings file
settings = read.settings('pecan.SDA.xml')

# set SDA-specific variables
nens = as.numeric(settings$ensemble$size)
nspec = length(settings$pfts)
site = settings$run$site$name
sda.years = c(lubridate::year(settings$state.data.assimilation$start.date):lubridate::year(settings$state.data.assimilation$end.date))
sda.years = sda.years[-1]
nyr = length(sda.years)
nparam = 20 # maybe want to make determined 
no_param = c(1,2,7,16,25,26)
runs = list.dirs('./out', full.names = FALSE)[-1]

# storage for first step
# param.melt
param.array = array(NA, dim= c(nspec, nparam, nens))
param.names = c('DMAX','DMIN','B3','B2','AGEMX','G','SPRTND','SPRTMN','SPRTMX','MPLANT','D3','FROST',
                'CM1','CM2','CM3','CM4','CM5','FWT','SLTA','SLTB')

# coordinates for sites 
coords = matrix(NA, 4, 2)
coords[1,] = c() # HF
coords[2,] = c() # RH
coords[3,] = c() # NRP
coords[4,] = c() # GE

# life.melt 
birth.array = array(NA, dim=c(nens, nyr, nspec))
growth.array = array(NA, dim =c(nens, nyr, nspec))
death.array = array(NA, dim=c(nens, nyr, nspec))
err.agb.pft.array = array(NA, dim = c(nens,nyr,nspec))

# storage for second step
# bias.melt 
tot.bias.mat = matrix(NA,nens,nyr)

# storage for third step
# run.list 
obs = nens * nyr
run.mat = matrix(NA, obs, 17)
colnames(run.mat) = c('ens','year','lat','long','pred','bias','weight','summer.temp','winter.temp','summer.precip','winter.precip',
                      'g.season','basal.area','algf25', 'algf50', 'algf75','sngf','stand.age',
                      'dominant')
run.mat = as.data.frame(run.mat)
run.list = list()

# other storage
precip.array = array(NA, c(nyr,12,nens))
temp.array = array(NA, c(nyr,12,nens))

# get weights
load(paste0('./SDA/sda.output.Rdata')) # contains forecast and analysis matrix
load(paste0('/data/dbfiles/sda.obs.Rdata')) # weight calculation
weights = matrix(0,nens,nyr)
for (j in 1:nyr){
  weights[,j] = weight(X = as.matrix(FORECAST[[j+1]]), 
                       mu.a = as.matrix(obs.list$obs.mean[[j+1]], nrow = 1),
                       Pa = as.matrix(obs.list$obs.cov[[j+1]]))
}
rm(obs.list)
load(paste0('./SDA/outconfig.Rdata')) # met ensemble information
ids = list.dirs('./run/', full.names = FALSE)[-1]

# loop through ensembles
for (i in 1:nens){

  cid = ids[(i+nens)]

  # load input
  load(paste0('./run/',cid,'/linkages.input.Rdata'))
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
    
    run.mat$ens[ind] = i 
    run.mat$year[ind] = yr
    run.mat$weight[ind] = weights[i,j]
  
    # error 
    tot.bias.mat[i,j] = (sum(FORECAST[[j+1]][i,]) - sum(ANALYSIS[[j+1]][i,])) / sum(ANALYSIS[[j+1]][i,])
    err.agb.pft.array[i,j,] = as.matrix((FORECAST[[j+1]][i,] - ANALYSIS[[j+1]][i,]) / ANALYSIS[[j+1]][i,])
    
    # starting dominant species in stand 
    run.mat$dominant[ind] = which.max(as.matrix((ANALYSIS[[j]][i,])))
    
    # load same year restart file
    if (j == nyr){
      load(paste0('./run/',cid,'/linkages.restart.Rdata'))
    }else{
      load(paste0('./run/',cid,'/',toString(yr),'-12-31 23:59:59linkages.restart.Rdata'))
    }

    run.mat$basal.area[ind] = sum(pi * dbh^2)
    run.mat$stand.age[ind] = mean(iage[iage != 0])
    
    # starting number of trees for calculating birth
    ntrees.res = ntrees
    
    # load output file
    if (j == nyr){
      load(paste0('./out/',cid,'/linkages.out.Rdata'))
    }else{
      load(paste0('./out/',cid,'/',toString(yr),'-12-31 23:59:59linkages.out.Rdata'))
    }

    # life process rates for averaging
    birth.array[i,j,] = ntrees.birth[,1,1] - ntrees.res
    growth.array[i,j,] = agb.pft[,1,1] - as.matrix(ANALYSIS[[j]][i,]) / as.matrix(ANALYSIS[[j]][i,])
    death.array[i,j,] = ntrees.kill[,1,1] - ntrees.birth[,1,1]

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

bias.melt = melt(tot.bias.mat)
colnames(bias.melt) = c('ensemble','year','bias')

# after this point, we need: 
# nens, nspec, site, sda.years, nyr, nparam, runs

rm(birth.melt, growth.melt, death.melt, error.melt)

############################################
## Step 2: Determine sensitive parameters ##
############################################




############################################
## Step 2: Determine important parameters ##
############################################

# determine top <np> critical parameters for each species' growth
# only considers ensembles which have error rates less than or equal to 25% 
np = 3
cor.params = array(0, dim = c(np,4,nspec)) # s, p, real cor, abs cor 
for (k in 1:nspec){
  data = life.melt %>% filter(species == k) %>%
  #data = life.melt %>% filter(species == k, abs(error) <= 0.25) %>%
    select(ensemble,year,growth,error) %>% left_join(param.melt, id='ensemble') # change here
  for (s in 1:nspec){
    for (p in 1:nparam){
      
      x = data %>% filter(p.name == p, species == s) %>% select(p.value)
      y = data %>% filter(p.name == p, species == s) %>% select(error) # change here
      cor.c = cor(x, y)[1,1] # coefficient of correlation
      
      # check to make sure not NA or infinite
      if (is.na(cor.c) | abs(cor.c) == Inf){
        print(paste(s, p, 'is constant parameter'))
        next
      }
      
      # if greater than one listed, reorganize and crop 
      if (any(abs(cor.c) >= cor.params[,4,k])){
        
        if (any(abs(cor.c) == cor.params[,3,k])){
          print(paste('equal correlations for',k,s,p,'-',cor.c))
        }
        
        temp = rbind(cor.params[,,k], c(s,p,cor.c, abs(cor.c)))
        cor.params[,,k] = temp[order(temp[,4], decreasing=T)[1:np],]
      }
    }
  }
}

# get rid of absolute value column
cor.params = cor.params[,1:3,]
par.names = colnames(param.list[[1]])

# loop through and print ggplots for each sig parameter
plot.list = list()
now = 0
for (k in 1:nspec){
  for (i in 1:np){
    now = now + 1
    s = cor.params[i,1,k]
    p = cor.params[i,2,k]
    #df = life.melt %>% filter(species == k, abs(error) <= 0.25) %>%
    df = life.melt %>% filter(species == k) %>%
      select(ensemble,year,growth,error) %>% left_join(param.melt, id='ensemble') %>% # change here
      filter(p.name == p, species == s)
    pl = ggplot(df, aes(x = p.value, y = error, col = error)) + # change here
      geom_point()  + 
      geom_smooth(method='lm') + 
      labs(x=paste0(settings$pfts[s]$pft$name,'-', param.names[p]), 
           y = paste(settings$pfts[k]$pft$name,'-','rel. model error'))
    print(pl)
    plot.list[[now]] = ggplotGrob(pl)
  }
}
#grid.arrange(grobs = plot.list, nrow = 4)


#################################################
## Step 3: Look at temporal gaps in prediction ##
#################################################

ggplot(bias.melt, aes(x = year + 1960, y = bias, col = ensemble)) + 
  ylim(-1,1) +
  geom_point() + 
  labs(title = 'model bias over Time', 
       x = 'year', y = 'relative model bias')


###############################################################
## Step 4: Multidimensional analysis of species productivity ##
###############################################################

## THIS STEP SHOULD INCLUDE OBSERVATIONS FROM ALL SITES
## ONLY CONSIDERING ONE SPECIES AT A TIME RIGHT NOW

k = 1

########################################
# A. Predictors of prediction accuracy #
########################################

# Where does LINKAGES get predictions wrong about each species? # 
# Basic idea here is to use logistic regression to identify environmental variables which are significant to 
# correct predictions (binary variable: correct or incorrect) 

# get data from run.list for species
data = run.list[[k]]

# define good and bad predictions with set accuracy threshold
lim = 0.25 # 1 is good, 0 is bad
quality = sapply(data$bias, function(x){ifelse(abs(x) > lim, 0, 1)})
data$bias = quality

# look at all variables in relationship to classifying bias 
vars = names(data)[-(1:5)]
plot.list = list()
i = 0
for (c in vars){
  i = i + 1
  pl = ggplot(data, aes_string(x=c)) + geom_point(aes(y = bias, col = bias)) +
    labs(x = c, y = 'prediction accuracy')
  plot.list[[i]] = ggplotGrob(pl)
}
grid.arrange(grobs = plot.list, nrow = 4)

# run full logistic regression to determine significant variables 
# use summary to determine which variables to put in reduced model 
full.logreg = glm(bias~summer.temp+winter.temp+summer.precip+winter.precip+g.season+basal.area+
                  algf25+algf50+algf75+sngf+stand.age, data = data, family = 'binomial')
summary(full.logreg)

# fit reduced model to check - if variables are not significant here, remove them and run reduced model again
red.logreg = glm(bias~algf50+sngf, 
                 family = 'binomial', data = data)
summary(red.logreg)

# plot the significant features to look for patterns - where are predictions wrong?
ggplot(data) + geom_point(aes(x = sngf, y = pred, col = as.factor(bias)))
ggplot(data) + geom_point(aes(x = algf50, y = pred, col = as.factor(bias)))
ggplot(data) + geom_point(aes(x = basal.area, y = pred, col = as.factor(bias)))

# matrix of histograms to look for patterns in terms of bias that may have been missed
plot.list2 = list()
i = 1
for (c in vars){
  pl = ggplot(data = data, aes_string(x = c)) + 
    geom_histogram(aes(fill = as.factor(bias))) +
    labs(x = c) + 
    theme(legend.position = 'none')
  plot.list2[[i]] = ggplotGrob(pl)
  i = i + 1
}
grid.arrange(grobs = plot.list2, nrow = 4)

###################################
# A. Predictors of species growth #
###################################




plot.list = list()
for (c in 1:length(vars)){
  pl = ggplot(data) + geom_point(aes(x = data[,(5+c)], y = pred, col = bias)) +
    labs(x = vars[c], y = 'predicted growth')
  plot.list[[c]] = ggplotGrob(pl)
}
grid.arrange(grobs = plot.list, nrow = 4)

full.model = lm(pred~summer.temp+winter.temp+summer.precip+winter.precip+g.season+basal.area+
                  algf25+algf50+algf75+sngf+stand.age, data = data)
summary(full.model)

red.model = lm(pred~summer.temp+summer.precip+basal.area+algf25+stand.age, data = data)
summary(red.model)

# matrix of scatterplots
cols = rep(NA,length(data$bias))
cols[data$bias == 0] <- 'red'
cols[data$bias == 1] <- 'blue'
pairs(~summer.temp+summer.precip+basal.area+algf25+stand.age,
      lower.panel=panel.smooth, pch=20, cex = 0.5, col = cols, 
      data = data, main="Iris Scatterplot Matrix")




## FUNCTIONS ## 

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


