###-------------------------------------------------------------------###
### ensemble adjustment plots                                         ###
###-------------------------------------------------------------------###

library(mvtnorm)

# Calculate weights of ensemble members based on aboveground biomass by species
## Need Pa, mu.a, Pf, mu.f, X

## Inputs :
nspp = 4
site = 'HF'
num.runs = 200
first.id = 1002350792

# load observed data
directory = paste0('~/Desktop/Capstone_Spring/',site,'/')
load(paste0(directory,'Data/',site,'.gwbi.agb.Rdata'))

# calculate number of years
#data.start <- as.integer(substr(names(obs.list$obs.mean)[1],1,4))
#data.end <- as.integer(substr(names(obs.list$obs.mean)[length(obs.list$obs.mean)],1,4))
data.start = 1962
data.end = 2009
num.years <- data.end-data.start+1

# create storage matrices
ens.weights = matrix(NA,num.years,num.runs)

# iterate through years
for (i in 1:num.years){
  print(i)
  data.yr = data.start + i - 1

  # pinpoint analysis mean and covariance matrix (only first four entries)
  Pa = obs.list$obs.cov[[i]][1:4,1:4]
  mu.a = obs.list$obs.mean[[i]][1:4]

  X = matrix(0,num.runs,nspp)

  # load forecast data for all ensemble members
  for (j in 1:num.runs){
    id = first.id + j - 1
    load(paste0('~/Documents/Notre\ Dame/Senior/Capstone/LINKAGES/forecast.hypotheses/out.forecast/',id,'/',data.yr,'-12-31\ 23:59:59linkages.out.Rdata'))
    X[j,] = agb.pft[,1,1]
  }

  Pf = cov(X)
  mu.f = rep(NA,nspp)
  for (k in 1:nspp){
    mu.f[k] = mean(X[,k])
  }

  ens.weights[i,] = weight(X=X,mu.a=mu.a,Pa=Pa)
}

weight = function(X,mu.a,Pa){

  # calculate the likelihood of the ensemble members given mu.a and Pa for one given year
  nens <- nrow(X)
  wt.mat <- matrix(NA,nrow=nens)
  for (i in 1:nens){
    wt.mat[i] = dmvnorm(X[i,], mean=mu.a, sigma=Pa)
  }

  # put into weights table
  wt.props <- t(prop.table(wt.mat))
  return(wt.props)
}

