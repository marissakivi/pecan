## This script provides an alternative to the spinup method to start SDA with LINKAGES
## Instead of running spin-up, you will run one year of LINKAGES and then re-write the state variables
## that are propagated through time (by using the restart function)

## The script is broken up between aboveground and belowground processes because these are estimated differently

## Author: AM Willson
## Date modified: 01 July 2020

############
## Set-up ##
############

## Before running this script you must run LINKAGES for two years using the PEcAn web interface
## This is required because this process creates the infrastructure needed to run SDA in the PEcAn framework
## For example, we get input parameter values and site-specific inputs, as well as xml files in the correct format

## Choose the site you are interested in running and the species you want to include in your analysis
## The start and end dates should be chosen as follows:
  ## The end date should be 12/31 of the first year you plan to run SDA
  ## The start date should be 1/1 of the previous year
  ## These are chosen because they are the minimum time required to get an output RData object for each ensemble member
## Use the advanced setup option to specify the number of ensemble members you would like in your analysis

## Once this is done, we can overwrite the outputs for each ensemble member using independent estimates of each state

# Prepare workspace
rm(list = ls())
workflowID = 219
setwd(paste0('/data/workflows/PEcAn_14000000',workflowID,'/'))
nens = 20
nspec = 4

# Define moment matching function for distribution
beta_func = function(mu, sd){
  alpha = (mu^2 - mu^3 - mu * sd^2) / sd^2
  beta = (mu - 2 * mu^2 + mu^3 - sd^2 + mu * sd^2) / sd^2
  return(c(alpha, beta))
}

negbin_func = function(mu, sd){
  var = sd^2
  mu = mu
  size = (mu + mu^2) / var
  return(c(mu, size))
}

####################
## Re-write files ##
####################
## This part of the workflow is identical to that found and explained in workflow.spninup.R
## In short, the PEcAn interface does not use the correct met inputs, so we have to manually re-run everything
## See workflow.spinup.R for details

options(warn=1)
options(error=quote({
  PEcAn.utils::status.end("ERROR")
  PEcAn.remote::kill.tunnel(settings)
  if (!interactive()) {
    q(status = 1)
  }
}))

ensemble_location = '/data/dbfiles/met_data/SYLVANIA/weights/ensemble-weights-SYLVANIA-prism.csv'
metdir <- '/data/dbfiles/met_data/SYLVANIA/linkages/'
n = 20

ens_wts <- read.csv(ensemble_location)
clim_mods <- ens_wts$climate_model
avg_wt <- ens_wts$wts

clim_use <- sample(x = clim_mods, size = n, prob = avg_wt, replace = T)

wts_use <- vector()
for (i in 1:n) wts_use[i] = avg_wt[which(clim_mods == clim_use[i])]

wts_use_1 = wts_use / sum(wts_use)
save(ens_wts = wts_use_1, file = './weights.Rdata')

name <- numeric(length(clim_use))
for(i in 1:length(clim_use)){
  name[i] <- paste0('\t\t <path',i,'>',metdir,clim_use[i],'.Rdata</path',i,'>')
}
writeLines(name)

## Insert these lines into pecan.TRAIT.xml following the instructions in workflow.spinup.R
## Change method from sampling to looping
## Delete all lines below META in STATUS
## Delete run & out folders

settings <- PEcAn.settings::read.settings("pecan.TRAIT.xml")

if(!is.null(settings$meta.analysis)) {
  if (PEcAn.utils::status.check("META") == 0){
    PEcAn.utils::status.start("META")
    PEcAn.MA::runModule.run.meta.analysis(settings)
    PEcAn.utils::status.end()
  }
}

if (PEcAn.utils::status.check("CONFIG") == 0){
  PEcAn.utils::status.start("CONFIG")
  settings <- PEcAn.workflow::runModule.run.write.configs(settings)
  PEcAn.settings::write.settings(settings, outputfile='pecan.CONFIGS.xml')
  PEcAn.utils::status.end()
} else if (file.exists(file.path(settings$outdir, 'pecan.CONFIGS.xml'))) {
  settings <- PEcAn.settings::read.settings(file.path(settings$outdir, 'pecan.CONFIGS.xml'))
}

if ((length(which(commandArgs() == "--advanced")) != 0) && (PEcAn.utils::status.check("ADVANCED") == 0)) {
  PEcAn.utils::status.start("ADVANCED")
  q();
}

if (PEcAn.utils::status.check("MODEL") == 0) {
  PEcAn.utils::status.start("MODEL")
  PEcAn.remote::runModule.start.model.runs(settings, stop.on.error = FALSE)
  PEcAn.utils::status.end()
}

if (PEcAn.utils::status.check("OUTPUT") == 0) {
  PEcAn.utils::status.start("OUTPUT")
  runModule.get.results(settings)
  PEcAn.utils::status.end()
}

if ('ensemble' %in% names(settings) & PEcAn.utils::status.check("ENSEMBLE") == 0) {
  PEcAn.utils::status.start("ENSEMBLE")
  PEcAn.uncertainty::runModule.run.ensemble.analysis(settings, TRUE)    
  PEcAn.utils::status.end()
}

if ('sensitivity.analysis' %in% names(settings) & PEcAn.utils::status.check("SENSITIVITY") == 0) {
  PEcAn.utils::status.start("SENSITIVITY")
  runModule.run.sensitivity.analysis(settings)
  PEcAn.utils::status.end()
}

if ('assim.batch' %in% names(settings)) {
  if (PEcAn.utils::status.check("PDA") == 0) {
    PEcAn.utils::status.start("PDA")
    settings <- PEcAn.assim.batch::runModule.assim.batch(settings)
    PEcAn.utils::status.end()
  }
}

if("benchmarking" %in% names(settings) & "benchmark" %in% names(settings$benchmarking)){
  PEcAn.utils::status.start("BENCHMARKING")
  results <- papply(settings, function(x) calc_benchmark(x, bety))
  PEcAn.utils::status.end()
}

if (PEcAn.utils::status.check("FINISHED") == 0) {
  PEcAn.utils::status.start("FINISHED")
  PEcAn.remote::kill.tunnel(settings)
  db.query(paste("UPDATE workflows SET finished_at=NOW() WHERE id=", 
                 settings$workflow$id, "AND finished_at IS NULL"), 
           params = settings$database$bety)
  PEcAn.utils::status.end()
}

db.print.connections()
print("---------- PEcAn Workflow Complete ----------")

## This is the end of the fake "spinup" that has to be done to get files in the correct format
## The rest of the script is specific to the alternative method and is not copied from workflow.spinup.R

#####################
## Load in IC data ##
#####################

## Load in an RData file with five components:
## 1) ICs: contains parameter names and distribution parameters for all ICs except those in other components
## 2) tyl.store: contains TYL data for the final year of a nens-ensemble member, 1000-year (for old growth site) dry run at the site
## 3) C.mat.store: contains C.mat data for the final year of a nens-ensemble member, 1000-year (for old growth site) dry run at the site
## 4) ICs_DBH_cat: contains proportions of each species in each DBH class
## 5) DBH_cat: contains the minimum and maximum DBH for each DBH class

## Importantly, for this to work, the ICs must be formatted correctly. AW can provide an example.

# Load data
load('/home/acer/amw/ICs.RData')

####################################
## 1) Number of trees per species ##
####################################

# For each ensemble member, draw the number of trees for each species using ICs
# A normal distribution was chosen for convenience, but the data are rounded because this is a discrete quantity

# Storage
ntree = matrix(, nrow = nspec, ncol = nens)
ntree_params = matrix(, nrow = nspec, ncol = 2)

# Extract information from ICs
for(i in 1:nspec){
  ntree_params[i,1] = as.numeric(ICs$parama[which(ICs$Variable == paste0('ntrees_spec',i))])
  ntree_params[i,2] = as.numeric(ICs$paramb[which(ICs$Variable == paste0('ntrees_spec',i))])
}

# Draw number of trees from negative binomial distribution
# 1 is added to each draw to make sure that at least one of each tree occurs in each ensemble member

for(i in 1:nspec){
  ntree[i,] = rnegbin(nens, mu = negbin_func(ntree_params[i,1], ntree_params[i,2])[1], theta = negbin_func(ntree_params[i,1], ntree_params[i,2])[2]) + 1
}

#######################
## 2) Record species ##
#######################

## Assign each tree a species
## We also create a data frame to track the characteristics of each tree here

# Start by creating each variable for one ensemble member

# Individual tree ID number
id = 1:sum(ntree[,1])
ens = rep(1, times = length(id))
for(i in 1:nspec){
  if(i == 1){
    spec = rep(i, times = ntree[i,1])
  }else{
    spec = c(spec, rep(i, times = ntree[i,1]))
  }
}

# Then fill in for the other ensemble members
for(i in 2:nens){
  id = c(id, 1:sum(ntree[,i]))
  ens = c(ens, rep(i, times = sum(ntree[,i])))
  for(j in 1:nspec){
    spec = c(spec, rep(j, times = ntree[j,i]))
  }
}

track_trees = cbind(id, ens, spec)
track_trees = as.data.frame(track_trees)
colnames(track_trees) = c('id', 'ens', 'pft')

#####################
## 3) Assign NOGRO ##
#####################

## Assign trees NOGRO flags (0 or 1) based on species-specific probabilities
## Assumes that no trees will have 2 flags, so no trees will immediately die

# Storage
prob_nogro = matrix(, nrow = nspec, ncol = nens)

# Get species-specific NOGRO probability per ensemble member
for(i in 1:nspec){
  mu = as.numeric(ICs$parama[which(ICs$Variable == paste0('nogro_spec',i))])
  sd = as.numeric(ICs$paramb[which(ICs$Variable == paste0('nogro_spec',i))])
  prob_nogro[i,] = rbeta(nens, beta_func(mu, sd)[1], beta_func(mu, sd)[2])
}

# Storage
nogro = matrix(, nrow = 200, ncol = nens)

# Assign each tree a NOGRO number according to species and ensemble member
for(i in 1:nens){
  for(j in 1:sum(ntree[,i])){
    sp = track_trees$pft[which(track_trees$id == j & track_trees$ens == i)]
    nogro[j,i] = rbinom(1, 1, prob_nogro[sp,i])
  }
}

# Track nogro to add to track trees
ng = nogro[,1]
ng = ng[which(!is.na(ng))]

for(i in 2:nens){
  n = nogro[,i]
  n = n[which(!is.na(n))]
  ng = c(ng, n)
}

track_trees$nogro = ng

###################
## 4) Assign DBH ##
###################
## This step assigns a DBH based on species-specific probabilities of falling in DBH bins
## This information is then added to track_trees
## This part does not use ICs, but uses ICs_DBH_cat and DBH_cat
## Importantly, this assumes that are 12 DBH bins -- must change loop if this is not the case!!

# Storage
ICs_DBH_cumsum = matrix(, nrow = nspec, ncol = 12)

# Get ICs_DBH_cat cumulative sums
for(i in 1:nrow(ICs_DBH_cat)){
  ICs_DBH_cumsum[i,] = cumsum(ICs_DBH_cat[i,])
}

# Storage
dbh = matrix(, nrow = 200, ncol = nens)

# Draw DBH for each tree based on species-specific probabilities of falling into given DBH bins
for(i in 1:nens){
  for(j in 1:sum(ntree[,i])){
    rand = runif(1, 0, 1)
    sp = track_trees$pft[which(track_trees$id == j & track_trees$ens == i)]
    if(rand < ICs_DBH_cumsum[sp,1] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[1], DBH_cat$max[1])
    }
    if(rand >= ICs_DBH_cumsum[sp,1] & rand < ICs_DBH_cumsum[sp,2] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[2], DBH_cat$max[2])
    }
    if(rand >= ICs_DBH_cumsum[sp,2] & rand < ICs_DBH_cumsum[sp,3] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[3], DBH_cat$max[3])
    }
    if(rand >= ICs_DBH_cumsum[sp,3] & rand < ICs_DBH_cumsum[sp,4] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[4], DBH_cat$max[4])
    }
    if(rand >= ICs_DBH_cumsum[sp,4] & rand < ICs_DBH_cumsum[sp,5] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[5], DBH_cat$max[5])
    }
    if(rand >= ICs_DBH_cumsum[sp,5] & rand < ICs_DBH_cumsum[sp,6] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[6], DBH_cat$max[6])
    }
    if(rand >= ICs_DBH_cumsum[sp,6] & rand < ICs_DBH_cumsum[sp,7] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[7], DBH_cat$max[7])
    }
    if(rand >= ICs_DBH_cumsum[sp,7] & rand < ICs_DBH_cumsum[sp,8] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[8], DBH_cat$max[8])
    }
    if(rand >= ICs_DBH_cumsum[sp,8] & rand < ICs_DBH_cumsum[sp,9] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[9], DBH_cat$max[9])
    }
    if(rand >= ICs_DBH_cumsum[sp,9] & rand < ICs_DBH_cumsum[sp,10] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[10], DBH_cat$max[10])
    }
    if(rand >= ICs_DBH_cumsum[sp,10] & rand < ICs_DBH_cumsum[sp,11] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[11], DBH_cat$max[11])
    }
    if(rand >= ICs_DBH_cumsum[sp,11] & is.na(dbh[j,i])){
      dbh[j,i] = runif(1, DBH_cat$min[12], DBH_cat$max[12])
    }
  }
}

# Track DBH and add to track trees
d = dbh[,1]
d = d[which(!is.na(d))]

for(i in 2:nens){
  dd = dbh[,i]
  dd = dd[which(!is.na(dd))]
  d = c(d, dd)
}

track_trees$dbh = d

###################
## 5) Assign age ##
###################
## Age can be assigned in two ways
## 1) Age of each tree is deterministically based on species-specific linear relationships between DBH and age
## 2) Age is drawn from species-specific averages
## Both methods are processed here, but you must choose one to include in the output replacement below

## Method 1

# Storage
coeff = matrix(, nrow = nspec, ncol = 2)

# Extract coefficients of linear model
# We only use means because the variances result in impossible ages
for(i in 1:nspec){
  coeff[i,1] = as.numeric(ICs$parama[which(ICs$Variable == paste0('age_b0_spec',i))])
  coeff[i,2] = as.numeric(ICs$parama[which(ICs$Variable == paste0('age_b1_spec',i))])
}

# Storage
age1 = matrix(, nrow = 200, ncol = nens)

# Compute age from DBH
for(i in 1:nens){
  for(j in 1:sum(ntree[,i])){
    sp = track_trees$pft[which(track_trees$id == j & track_trees$ens == i)]
    age1[j,i] = round(coeff[sp,1] + coeff[sp,2] * dbh[j,i])
  }
}

# Check for negative values
which(age1 < 0)

## Method 2

# Storage
age_est = matrix(, nrow = nspec, ncol = 2)

# Extract species-specific mean and variance
for(i in 1:nspec){
  age_est[i,1] = as.numeric(ICs$parama[which(ICs$Variable == paste0('age_spec',i))])
  age_est[i,2] = as.numeric(ICs$paramb[which(ICs$Variable == paste0('age_spec',i))])
}

# Storage
age2 = matrix(, nrow = 200, ncol = nens)

# Make random draws for age
for(i in 1:nens){
  for(j in 1:sum(ntree[,i])){
    sp = track_trees$pft[which(track_trees$id == j & track_trees$ens == i)]
    age2[j,i] = rnegbin(1, mu = negbin_func(age_est[sp,1], age_est[sp,2])[1], theta = negbin_func(age_est[sp,1], age_est[sp,2])[2])
  }
}

# Check for zeros
which(age2 == 0)

# Track age (both methods) and add to track trees
ag = age1[,1]
ag = ag[which(!is.na(ag))]

for(i in 2:nens){
  a = age1[,i]
  a = a[which(!is.na(a))]
  ag = c(ag, a)
}

track_trees$age1 = ag

ag = age2[,1]
ag = ag[which(!is.na(ag))]

for(i in 2:nens){
  a = age2[,i]
  a = a[which(!is.na(a))]
  ag = c(ag, a)
}

track_trees$age2 = ag

# Those are all the aboveground quantities we need to estimate! 
# We will come back to this when we save these in the output files

###########################
## Belowground variables ##
###########################

## The belowground variables are the number of litter cohorts (ncohrt) and belowground C matrix (C.mat)
## Here, we don't need to create individuals, so each component can be +/- independent

## There is little information on soil properties and the information that does exist has large variance
## Additionally, the concept of litter cohorts in this regard is fairly LINKAGES-specific and would be difficult to quantify in a forest
## Therefore, the soil ICs are estimated from a long spinup run of LINKAGES and collated here

###############
## 6) NCOHRT ##
###############

## Number of litter cohorts taken from the number of existing cohorts in TYL

# Storage
ncohort = NA

for(i in 1:nens){
  ncohort[i] = length(which(tyl.store[,i] > 0))
}

####################
## 7) C.mat & TYL ##
####################

## C.mat & TYL were just taken from the output of a spinup run, so they do not need to be processed
newC = C.mat.store
newtyl = tyl.store

#################################
## Replace in ensemble outputs ##
#################################
# This part loads the output data for each ensemble member, replaces the quantities of interest, and re-saves the outputs in the same file

# Choose age option
age = age2

# List ensemble members
en = list.files(path = 'out/')

for(i in 1:nens){
  # Load output
  load(paste0('out/',en[i],'/linkages.out.Rdata'))
  
  # Replace
  ntrees.kill[,ncol(ntrees.kill),] = ntree[,i]
  ncohrt = ncohort[i]
  C.mat = newC[,,i]
  C.mat[which(is.na(C.mat))] = 0
  tyl = tyl.store[,i]
  tyl.save[,ncol(tyl.save),] = tyl.store[,i]
  nogro.save[,ncol(nogro.save),] = nogro[,i]
  nogro.save[which(is.na(nogro.save))] = 0
  iage.save[, ncol(iage.save),] = age[,i]
  iage.save[which(is.na(iage.save))] = 0
  dbh.save[,ncol(dbh.save),] = dbh[,i]
  dbh.save[which(is.na(dbh.save))] = 0
  
  # Resave
  save(abvgrnwood, abvgroundwood.biomass, aet.save, ag.biomass, ag.npp, area, avln, C.mat, cn, et, f.comp, ff, fl, hetero.resp, ksprt, leaf.litter, nee, sco2c, som, tab, tnap, total.soil.carbon, totl, tstem, tyl, water, agb.pft, algf.save.keep, awp.save, bar, dbh.save, gf.vec.save, iage.save, ncohrt, nogro.save, npp.spp.save, ntrees.birth, ntrees.kill, tyl.save, year,
       file = paste0('out/',en[i],'/linkages.out.Rdata'))
}

