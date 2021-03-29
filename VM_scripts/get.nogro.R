## Quick script to extract NOGRO from a long spinup for Sylvania to set ICs

## Author: AM Willson
## Date modified: 19 June 2020

rm(list = ls())
setwd('/data/workflows/PEcAn_14000000216/')
library(reshape2)

# List run IDs
runs = list.files('out/')

# Storage
nogro.array = array(0,dim=c(4,100,20))

# Extract NOGRO and number of trees of each species and find fraction of trees of each species flagged for nogro each year
for (i in 1:20){
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  for (j in 902:1001){
    nu = 0
    for (k in 1:4){
      nl = nu + ntrees.kill[k,j,1] - 1
      nogro = sum(nogro.save[nu:nl,j,1] < 0)
      nogro.array[k,j-901,i] = nogro/(ntrees.kill[k,j,1])
      nu = nl + 1
    }
  }
}

# Make data frame
nogros_melt = melt(nogro.array)
colnames(nogros_melt) = c('species', 'year', 'ens', 'count')

# Subset by species
nogro_1 = subset(nogros_melt, nogros_melt$species == 1)
nogro_2 = subset(nogros_melt, nogros_melt$species == 2)
nogro_3 = subset(nogros_melt, nogros_melt$species == 3)
nogro_4 = subset(nogros_melt, nogros_melt$species == 4)

# Storage
ng_mean_1 = c()
ng_mean_2 = c()
ng_mean_3 = c()
ng_mean_4 = c()

# Average over years
for(i in 1:20){
  ng_mean_1[i] = mean(nogro_1$count[which(nogro_1$ens == i)])
  ng_mean_2[i] = mean(nogro_2$count[which(nogro_2$ens == i)])
  ng_mean_3[i] = mean(nogro_3$count[which(nogro_3$ens == i)])
  ng_mean_4[i] = mean(nogro_4$count[which(nogro_4$ens == i)])
}

# Get mean and sd for each species
nogro_mean_1 = mean(ng_mean_1,na.rm = T)
nogro_mean_2 = mean(ng_mean_2, na.rm = T)
nogro_mean_3 = mean(ng_mean_3, na.rm = T)
nogro_mean_4 = mean(ng_mean_4, na.rm = T)
nogro_sd_1 = sd(ng_mean_1, na.rm = T)
nogro_sd_2 = sd(ng_mean_2, na.rm = T)
nogro_sd_3 = sd(ng_mean_3, na.rm = T)
nogro_sd_4 = sd(ng_mean_4, na.rm = T)

# Save to export to local machine
nogro_est = data.frame(c(1:4))
colnames(nogro_est) = c('Species')
nogro_est$parama = c(nogro_mean_1, nogro_mean_2, nogro_mean_3, nogro_mean_4)
nogro_est$paramb = c(nogro_sd_1, nogro_sd_2, nogro_sd_3, nogro_sd_4)

dir.create('IC_dat')
save(nogro_est, file = 'IC_dat/nogro_est.RData')
