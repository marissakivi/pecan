## Quick script to extract soil conditions from a long spinup to set ICs

## Author: AM Willson
## Date modified: 25 June 2020

rm(list = ls())
setwd('/data/workflows/PEcAn_14000000216/')
library(reshape2)
library(dplyr)
library(ggplot2)

# List run IDs
runs = list.files('out/')

# Number of ensembles
nens = length(runs)
# Number of years of interest
nyear = 10

# Storage
C.mat.store = array(, dim = c(100, 15, nens))
tyl.store = array(, dim = c(20, nyear, nens))

# Store C.mat for each ensemble member
for(i in 1:nens){
  load(paste0('out/',runs[i],'/linkages.out.Rdata'))
  C = C.mat
  for(j in 1:nrow(C)){
    if(C[j,1] == 0 & C[j,2] == 0 & C[j,3] == 0 & C[j,4] == 0 & C[j,5] == 0 & C[j,6] == 0 & C[j,7] == 0 & C[j,8] == 0 & C[j,9] == 0 & C[j,10] == 0 & C[j,11] == 0 & C[j,12] == 0){
      C[j,] = rep(NA, times = 15)
    }
  }
  C.mat.store[,,i] = C
  print(C[,5])
  rm(C.mat, C)
}

# Store last 10 years of tyl
for(i in 1:nens){
  load(paste0('out/',runs[i],'/linkages.out.Rdata'))
  for(j in 992:max(year)){
    tyl.store[,j-991,i] = tyl.save[,j,]
  }
  print(i)
}

###################
## Process C.mat ##
###################

# Melt to dataframe
C.mat_melt = melt(C.mat.store)

colnames(C.mat_melt) = c('Cohort', 'Column', 'Ensemble', 'Value')

# Subset out non-existent cohorts
C.mat_melt = subset(C.mat_melt, !is.na(C.mat_melt$Value))

# Maximum & minimum number of cohorts for this year
ncohort = NA
for(i in 1:nens){
  ncohort[i] = length(unique(C.mat_melt$Cohort[which(C.mat_melt$Ensemble == i)]))
}
min(ncohort)
max(ncohort)

# Separate by column of the original matrix
C.mat_1 = subset(C.mat_melt, C.mat_melt$Column == 1, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_2 = subset(C.mat_melt, C.mat_melt$Column == 2, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_3 = subset(C.mat_melt, C.mat_melt$Column == 3, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_4 = subset(C.mat_melt, C.mat_melt$Column == 4, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_5 = subset(C.mat_melt, C.mat_melt$Column == 5, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_6 = subset(C.mat_melt, C.mat_melt$Column == 6, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_7 = subset(C.mat_melt, C.mat_melt$Column == 7, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_8 = subset(C.mat_melt, C.mat_melt$Column == 8, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_9 = subset(C.mat_melt, C.mat_melt$Column == 9, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_10 = subset(C.mat_melt, C.mat_melt$Column == 10, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_11 = subset(C.mat_melt, C.mat_melt$Column == 11, select = c('Cohort', 'Ensemble', 'Value'))
C.mat_12 = subset(C.mat_melt, C.mat_melt$Column == 12, select = c('Cohort', 'Ensemble', 'Value'))

# Make column names refer to the quantity from each original column
colnames(C.mat_1) = c('Cohort', 'Ensemble', 'Weight')
colnames(C.mat_2) = c('Cohort', 'Ensemble', 'N.content')
colnames(C.mat_3) = c('Cohort', 'Ensemble', 'Nchange1')
colnames(C.mat_4) = c('Cohort', 'Ensemble', 'Nchange2')
colnames(C.mat_5) = c('Cohort', 'Ensemble', 'Litter.type')
colnames(C.mat_6) = c('Cohort', 'Ensemble', 'Destination')
colnames(C.mat_7) = c('Cohort', 'Ensemble', 'perc.lignin')
colnames(C.mat_8) = c('Cohort', 'Ensemble', 'Ligchange1')
colnames(C.mat_9) = c('Cohort', 'Ensemble', 'Ligchange2')
colnames(C.mat_10) = c('Cohort', 'Ensemble', 'Original.weight')
colnames(C.mat_11) = c('Cohort', 'Ensemble', 'perc.N')
colnames(C.mat_12) = c('Cohort', 'Ensemble', 'Frac.to.humus')

## We don't actually care about all of these columns
## N change parameters, destination, lignin change parameters, fraction to humus all don't change

rm(C.mat_3, C.mat_4, C.mat_6, C.mat_8, C.mat_9, C.mat_12)

# Combine each column with litter type (column 5)
C.mat_1 = C.mat_1 %>%
  full_join(C.mat_5, by = c('Cohort', 'Ensemble'))
C.mat_2 = C.mat_2 %>%
  full_join(C.mat_5, by = c('Cohort', 'Ensemble'))
C.mat_7 = C.mat_7 %>%
  full_join(C.mat_5, by = c('Cohort', 'Ensemble'))
C.mat_10 = C.mat_10 %>%
  full_join(C.mat_5, by = c('Cohort', 'Ensemble'))
C.mat_11 = C.mat_11 %>%
  full_join(C.mat_5, by = c('Cohort', 'Ensemble'))

###########################
## Look at C.mat columns ##
###########################

## 1) What cohorts exist in each ensemble
unique(C.mat_1$Litter.type)

## 2) Weight of each cohort for each ensemble
C.mat_1$Litter.type = as.character(C.mat_1$Litter.type)
# Remove small wood because there are few observations
C.mat_1 = subset(C.mat_1, C.mat_1$Litter.type != 14)

labs = c('1' = 'Dogwood leaves', '2' = 'Maple leaves', '3' = 'Cherry leaves', '4' = 'Birch leaves', '5' = 'White oak leaves', '6' = 'Hemlock & cedar leaves', '7' = 'Aspen leaves', '8' = 'Beech leaves', '9' = 'Red oak leaves', '10' = 'Fir leaves', '11' = 'Spruce leaves', '12' = 'Pine leaves', '13' = 'Roots', '14' = 'Small fresh wood', '15' = 'Large fresh wood', '16' = 'Twigs', '17' = 'Well-decayed wood', '18' = 'Humus')

ggplot(C.mat_1, aes(x = Weight, group = Litter.type)) +
  geom_histogram() +
  facet_wrap(~Litter.type, scales = 'free', labeller = as_labeller(labs)) +
  xlab('Cohort weight (t/ha)')

## 3) Current N content of each cohort  
C.mat_2$Litter.type = as.character(C.mat_2$Litter.type)
C.mat_2 = subset(C.mat_2, C.mat_2$Litter.type != 14)

ggplot(C.mat_2, aes(x = N.content, group = Litter.type)) +
  geom_histogram() +
  facet_wrap(~Litter.type, scales = 'free', labeller = as_labeller(labs)) +
  xlab('N content (t/ha)')

## 4) Current % lignin of each cohort
C.mat_7$Litter.type = as.character(C.mat_7$Litter.type)
C.mat_7 = subset(C.mat_7, C.mat_7$Litter.type != 14)

ggplot(C.mat_7, aes(x = perc.lignin, group = Litter.type)) +
  geom_histogram() +
  facet_wrap(~Litter.type, scales = 'free', labeller = as_labeller(labs)) +
  xlab('Lignin content (%)')

## 5) Original weight
C.mat_10$Litter.type = as.character(C.mat_10$Litter.type)
C.mat_10 = subset(C.mat_10, C.mat_10$Litter.type != 14)

ggplot(C.mat_10, aes(x = Original.weight, group = Litter.type)) +
  geom_histogram() +
  facet_wrap(~Litter.type, scales = 'free', labeller = as_labeller(labs)) +
  xlab('Original weight (t/ha)')

## 6) Current % lignin
C.mat_11$Litter.type = as.character(C.mat_11$Litter.type)
C.mat_11 = subset(C.mat_11, C.mat_11$Litter.type != 14)

ggplot(C.mat_11, aes(x = perc.N, group = Litter.type)) +
  geom_histogram() +
  facet_wrap(~Litter.type, scales = 'free', labeller = as_labeller(labs)) +
  xlab('Nitrogen content (%)')

#################
## Process tyl ##
#################

## Two things I want to look at: TYL in the last ~10 years for each ensemble member and 
## summary stats for TYL through time

# Take last 10 years
tyl_last10 = tyl.store

# Melt to data frame
tyl_last10_melt = melt(tyl_last10)
colnames(tyl_last10_melt) = c('Type', 'Year', 'Ensemble', 'Value')

# For some reason there are 20 rows, but only 18 of them can get filled. Remove last two
tyl_last10_melt = subset(tyl_last10_melt, tyl_last10_melt$Type %in% c(1:18))

# Make litter types descriptors
for(i in 1:nrow(tyl_last10_melt)){
  if(tyl_last10_melt$Type[i] == 1){
    tyl_last10_melt$Type[i] = 'Dogwood leaves'
  }
  if(tyl_last10_melt$Type[i] == 2){
    tyl_last10_melt$Type[i] = 'Maple/ash/basswood leaves'
  }
  if(tyl_last10_melt$Type[i] == 3){
    tyl_last10_melt$Type[i] = 'Cherry leaves'
  }
  if(tyl_last10_melt$Type[i] == 4){
    tyl_last10_melt$Type[i] = 'Birch leaves'
  }
  if(tyl_last10_melt$Type[i] == 5){
    tyl_last10_melt$Type[i] = 'White oak leaves'
  }
  if(tyl_last10_melt$Type[i] == 6){
    tyl_last10_melt$Type[i] = 'Hemlock & cedar leaves'
  }
  if(tyl_last10_melt$Type[i] == 7){
    tyl_last10_melt$Type[i] = 'Aspen leaves'
  }
  if(tyl_last10_melt$Type[i] == 8){
    tyl_last10_melt$Type[i] = 'Beech leaves'
  }
  if(tyl_last10_melt$Type[i] == 9){
    tyl_last10_melt$Type[i] = 'Red oak leaves'
  }
  if(tyl_last10_melt$Type[i] == 10){
    tyl_last10_melt$Type[i] = 'Fir leaves'
  }
  if(tyl_last10_melt$Type[i] == 11){
    tyl_last10_melt$Type[i] = 'Spruce leaves'
  }
  if(tyl_last10_melt$Type[i] == 12){
    tyl_last10_melt$Type[i] = 'Pine leaves'
  }
  if(tyl_last10_melt$Type[i] == 13){
    tyl_last10_melt$Type[i] = 'Roots'
  }
  if(tyl_last10_melt$Type[i] == 14){
    tyl_last10_melt$Type[i] = 'Small fresh wood'
  }
  if(tyl_last10_melt$Type[i] == 15){
    tyl_last10_melt$Type[i] = 'Large fresh wood'
  }
  if(tyl_last10_melt$Type[i] == 16){
    tyl_last10_melt$Type[i] = 'Twigs'
  }
  if(tyl_last10_melt$Type[i] == 17){
    tyl_last10_melt$Type[i] = 'Total leaf litter'
  }
  if(tyl_last10_melt$Type[i] == 18){
    tyl_last10_melt$Type[i] = 'Total litter'
  }
}

# Copy tyl.store for full time series analysis
tyl_full = tyl.store

# Storage
tyl_full_quants = array(, dim = c(20, nyear, 3))

# Get quants for each litter type and each year
for(i in 1:nyear){
  for(j in 1:20){
    tyl_full_quants[j,i,] = quantile(tyl_full[j,i,], probs = c(0.025, 0.5, 0.975))
  }
}

# Melt data frame
tyl_quants_melt = melt(tyl_full_quants)
colnames(tyl_quants_melt) = c('Type', 'Year', 'Quantile', 'Value')

# Subset out weird random columns
tyl_quants_melt = subset(tyl_quants_melt, tyl_quants_melt$Type %in% c(1:18))

# Make columns more readable
for(i in 1:nrow(tyl_quants_melt)){
  if(tyl_quants_melt$Type[i] == 1){
    tyl_quants_melt$Type[i] = 'Dogwood leaves'
  }
  if(tyl_quants_melt$Type[i] == 2){
    tyl_quants_melt$Type[i] = 'Maple/ash/basswood leaves'
  }
  if(tyl_quants_melt$Type[i] == 3){
    tyl_quants_melt$Type[i] = 'Cherry leaves'
  }
  if(tyl_quants_melt$Type[i] == 4){
    tyl_quants_melt$Type[i] = 'Birch leaves'
  }
  if(tyl_quants_melt$Type[i] == 5){
    tyl_quants_melt$Type[i] = 'White oak leaves'
  }
  if(tyl_quants_melt$Type[i] == 6){
    tyl_quants_melt$Type[i] = 'Hemlock & cedar leaves'
  }
  if(tyl_quants_melt$Type[i] == 7){
    tyl_quants_melt$Type[i] = 'Aspen leaves'
  }
  if(tyl_quants_melt$Type[i] == 8){
    tyl_quants_melt$Type[i] = 'Beech leaves'
  }
  if(tyl_quants_melt$Type[i] == 9){
    tyl_quants_melt$Type[i] = 'Red oak leaves'
  }
  if(tyl_quants_melt$Type[i] == 10){
    tyl_quants_melt$Type[i] = 'Fir leaves'
  }
  if(tyl_quants_melt$Type[i] == 11){
    tyl_quants_melt$Type[i] = 'Spruce leaves'
  }
  if(tyl_quants_melt$Type[i] == 12){
    tyl_quants_melt$Type[i] = 'Pine leaves'
  }
  if(tyl_quants_melt$Type[i] == 13){
    tyl_quants_melt$Type[i] = 'Roots'
  }
  if(tyl_quants_melt$Type[i] == 14){
    tyl_quants_melt$Type[i] = 'Small fresh wood'
  }
  if(tyl_quants_melt$Type[i] == 15){
    tyl_quants_melt$Type[i] = 'Large fresh wood'
  }
  if(tyl_quants_melt$Type[i] == 16){
    tyl_quants_melt$Type[i] = 'Twigs'
  }
  if(tyl_quants_melt$Type[i] == 17){
    tyl_quants_melt$Type[i] = 'Total leaf litter'
  }
  if(tyl_quants_melt$Type[i] == 18){
    tyl_quants_melt$Type[i] = 'Total litter'
  }
  if(tyl_quants_melt$Quantile[i] == 1){
    tyl_quants_melt$Quantile[i] = '2.5%'
  }
  if(tyl_quants_melt$Quantile[i] == 2){
    tyl_quants_melt$Quantile[i] = '50%'
  }
  if(tyl_quants_melt$Quantile[i] == 3){
    tyl_quants_melt$Quantile[i] = '97.5%'
  }
  print(i)
}

#######################
## Last 10 years TYL ##
#######################

ggplot(tyl_last10_melt, aes(x = Value, group = Type)) +
  geom_histogram() +
  facet_wrap(~Type, scales = 'free')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Birch leaves')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Birch leaves')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Hemlock & cedar leaves')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Hemlock & cedar leaves')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Maple/ash/basswood leaves')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Maple leaves')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Roots')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Roots')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Twigs')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Twigs')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Large fresh wood')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Large fresh wood')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Total leaf litter')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Total leaf litter')

ggplot() +
  geom_histogram(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Total litter')])) +
  xlab('Weight (t/ha)') +
  ggtitle('Total litter')

ggplot() +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Birch leaves')], color = 'Birch'), alpha = 0.5) +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Hemlock & cedar leaves')], color = 'Hemlock & cedar'), alpha = 0.5) +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Maple/ash/basswood leaves')], color = 'Maple'), alpha = 0.5) +
  xlab('Weight (t/ha)') +
  ggtitle('Leaf litter')

ggplot() +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Roots')], color = 'Roots'), alpha = 0.5) +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Twigs')], color = 'Twigs'), alpha = 0.5) +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Large fresh wood')], color = 'Wood'), alpha = 0.5) +
  geom_density(aes(x = tyl_last10_melt$Value[which(tyl_last10_melt$Type == 'Total leaf litter')], color = 'Leaves'), alpha = 0.5) +
  xlab('Weight (t/ha)') +
  ggtitle('Litter by organ')

#####################
## Time series TYL ##
#####################

quant50 = subset(tyl_quants_melt, tyl_quants_melt$Quantile == '50%', select = c('Type', 'Year', 'Value'))
quant25 = subset(tyl_quants_melt, tyl_quants_melt$Quantile == '2.5%', select = c('Type', 'Year', 'Value'))
quant97 = subset(tyl_quants_melt, tyl_quants_melt$Quantile == '97.5%', select = c('Type', 'Year', 'Value'))

ggplot(quant50, aes(x = Year, y = Value, group = Type)) +
  geom_line() +
  facet_wrap(~Type, scales = 'free')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Birch leaves')], y = quant50$Value[which(quant50$Type == 'Birch leaves')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Birch leaves')], ymin = quant25$Value[which(quant25$Type == 'Birch leaves')], ymax = quant97$Value[which(quant97$Type == 'Birch leaves')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Birch leaves')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Hemlock & cedar leaves')], y = quant50$Value[which(quant50$Type == 'Hemlock & cedar leaves')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Hemlock & cedar leaves')], ymin = quant25$Value[which(quant25$Type == 'Hemlock & cedar leaves')], ymax = quant97$Value[which(quant97$Type == 'Hemlock & cedar leaves')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Hemlock & cedar leaves')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Maple/ash/basswood leaves')], y = quant50$Value[which(quant50$Type == 'Maple/ash/basswood leaves')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Maple/ash/basswood leaves')], ymin = quant25$Value[which(quant25$Type == 'Maple/ash/basswood leaves')], ymax = quant97$Value[which(quant97$Type == 'Maple/ash/basswood leaves')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Maple leaves')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Roots')], y = quant50$Value[which(quant50$Type == 'Roots')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Roots')], ymin = quant25$Value[which(quant25$Type == 'Roots')], ymax = quant97$Value[which(quant97$Type == 'Roots')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Roots')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Twigs')], y = quant50$Value[which(quant50$Type == 'Twigs')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Twigs')], ymin = quant25$Value[which(quant25$Type == 'Twigs')], ymax = quant97$Value[which(quant97$Type == 'Twigs')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Twigs')

ggplot() +
  geom_line(aes(x = quant50$Year[which(quant50$Type == 'Large fresh wood')], y = quant50$Value[which(quant50$Type == 'Large fresh wood')])) +
  geom_ribbon(aes(x = quant50$Year[which(quant50$Type == 'Large fresh wood')], ymin = quant25$Value[which(quant25$Type == 'Large fresh wood')], ymax = quant97$Value[which(quant97$Type == 'Large fresh wood')]), alpha = 0.5) +
  xlab('Year') + ylab('Weight (t/ha)') +
  ggtitle('Large fresh wood')
  
##################
## Save for ICs ##
##################

## After looking at this analysis, it does not seem feasible to accurately estimate the weights of litter or
## number of litter cohorts
## Therefore, I will save the existing Cmat and last year of tyl to use as ICs with SDA
## This needs to be done for each ensemble member

tyl.store = tyl.store[,10,]

save(C.mat.store, tyl.store, file = 'IC_dat/soil_est.RData')

## This file needs to be downloaded and aggregated with the rest of the IC data for processing in another script