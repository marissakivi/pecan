
## Tree_ring_input.R 
# Author: Ann Raiho, Marissa Kivi 
# Last modified: January 2020

# This script has two purposes: 
# (1) Determine the 98% species for the site based on the tree ring data. 
# (2) Process the estimated biomass estimates into species-level mean and covariance
# matrices. 

################################################################################
################################################################################
################################################################################

## Part I: Data formatting

# set up working environment
library(plyr)
library(dplyr)
#library(PEcAn.settings)
#library(PEcAn.DB)
#library(PEcAn.data.land)
library(ggplot2)

rm(list=ls())

### ADJUST VARIABLES HERE
site = 'HARVARD'

# input RDS data file location and name
input = readRDS(paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/NPP_STAT_MODEL_',site,'.RDS'))

# 98% plot save location and name
plot.loc = paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/sda_priority_spp_plot12.jpg')

# final data product save location
output = paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/sda.obs.',site,'.plots34.Rdata')
###

# Before anything else, let's remove all NA values in order to reduce variable size.
input <- input %>% filter(!is.na(value))

########### 1. Are there only 250 iterations included? ###########
# The NPP stat model produces thousands of iterations of AGB estimates for each site. However, for
# our purposes, there was an agreement that we really only need to keep the last 250 iterations for
# each site. Here we make sure this stays consistent for all sites.
max.iter = max(input$iter)
if (max.iter > 250){
  first = max.iter - 250
  input <- input %>% filter(iter > first) %>% mutate(iter = iter - first)
}

########### 2. Does the file include the right data types? ###########
# Here we adjust "type" variable for consistency across sites.
levels(input$type) # Do we see two data types? Some variation of AB (1) and ABI (2) ?
levels(input$type)[1] <- 'ab'
levels(input$type)[2] <- 'abi'
# Then we reduce our data to only include "ab" data.
input = input %>% filter(type == 'ab') %>% select(-type)

########### 3. Which models are available for the site? ###########
# First, look at the number of models used. If more than one, which model do you want to use?
levels(input$model)
model.pick = 'Model RW + Census'
# Then, reduce the input data to be just that model.
input = input %>% filter(model == model.pick) %>% select(-model)

########### 4. Which plots do we want to include? ###########
# On rare occasions, we want to only include data for some of the plots
# with available data at a site. For example, North Round Pond has plots 
# with two distinct stand structures and disturbance histories, so we split 
# them up. If you want to do all the plots, you can skip this step.
# First, we need to consider the name of the column that gives the plot number for the estimate.
names(input)
names(input)[3] = 'plot'
# Choose which plots you want to include. 
unique(input$plot)
plts = c(1,2,3,4)
input <- input %>% filter(plot %in% plts) %>% select(-plot)

################################################################################
################################################################################
################################################################################

## Part II: Determining 98% species

# set up variables for aggregation step
start_date = min(input$year)
end_date = max(input$year)
obs.mean <- obs.mean.tmp <- list()
obs.cov <- obs.cov.tmp <- list()
obs.times <- c(start_date:end_date)
time.type <- 'year'
biomass2carbon <- 0.48
var.names = 'AGB.pft'

# convert all data points from Mg/ha of biomass to kg/m2 of carbon
if(!is.null(input$value)) input$value <- udunits2::ud.convert(input$value,'Mg/ha','kg/m^2') * biomass2carbon #* kgm2Mgha 

# check to make sure no NA-taxon entries persist
input <- input[!is.na(input$taxon),]

# start aggregation step
variable <- c('value') #'ab'
arguments <- list(.(year, iter,  taxon), .(variable))
arguments2 <- list(.(year, taxon), .(variable))

# melt
melt_id <- colnames(input)[-which(colnames(input) %in% variable)]
melt.test <- reshape2::melt(input, id = melt_id, na.rm = TRUE)

# find annual sum of carbon biomass by species for each iteration
cast.test <- reshape2::dcast(melt.test, arguments, sum, margins = variable)

# find mean annual sum of carbon biomass by species across all iterations
melt_id_next <- colnames(cast.test)[-which(colnames(cast.test) %in% variable)]
melt.next <- reshape2::melt(cast.test, id = melt_id_next)
mean_mat <- reshape2::dcast(melt.next, arguments2, mean)

# determine overall cumulative biomass contribution of each species across all species 
prior_mat <- mean_mat %>% group_by(taxon) %>% 
  summarize(contr = sum(value)) %>%
  arrange(desc(contr))
prior_mat$perc = prior_mat$contr/sum(prior_mat$contr,na.rm=T)
prior_mat$cumsum = cumsum(prior_mat$perc)
prior_mat$taxon = factor(prior_mat$taxon, levels = prior_mat$taxon)
pl = ggplot(prior_mat) +
  geom_point(aes(x=taxon, y=cumsum)) + 
  geom_hline(yintercept = 0.98, col = 'red') + 
  labs(title = 'Overall Cumulative Proportion of Biomass by Species', 
       y = 'Cumulative Proportion of Biomass')
pl

######################## STOP ########################

# looking at the plot, record in the following vector the codes of the species whose points fall
# under the red line 
species = c('QURU','ACRU','FAGR','BEAL','TSCA')

# save plot
ggsave(pl, filename = plot.loc)

################################################################################
################################################################################
################################################################################

## Part III: Reduce and reformat results for running SDA in PEcAn

# First, reduce data frame to only include top 98% taxa. 
melt.next = melt.next %>% filter(taxon %in% species)
unique(melt.next$taxon)

# Now, we need to label the data for each species in a way that PEcAn will recognize so we 
# can successfully match the data to the model results at the species-level.
# If you're unsure what species the code refers to, you can checkout the ReadMe file for the RDS.

# We will label each species with a PFT name that will remain consistent for each species
# throughout your analysis. Traditionally, for LINKAGES, they follow the pattern: 
# Genus.Species_Common.Name. However, as you will learn later, we need to slightly vary this 
# pattern so as to keep your PFTs/priors distinct from the work of others. A few possibilities 
# include Genus.Species_Common.Name.YOURINITIALS or Genus.Species_Common.Name.VERSIONNUMBER. 
# Please note, I already used Genus.Species_Common.Name.2 for my analysis. You will see this later.

# Once you choose your pattern, you will need to fill in the following vector with the chosen PFT
# name for each of your species IN THE SAME ORDER AS THE SPECIES VECTOR ABOVE. 

species # for order 
x = c('Quercus.Rubra_Northern.Red.Oak.2',
      'Acer.Rubrum_Red.Maple.2',
      'Fagus.Grandifolia_American.Beech.2',
      'Betula.Alleghaniensis_Yellow.Birch.2',
      'Tsuga.Canadensis_Hemlock.2')
names(x) = species
x

################################################################################
################################################################################
################################################################################

## Part IV. Finish up processing of species-level biomass mean and covariance matrices

# adjust names of taxa to be with PFT name from BETY 
melt.next$pft.cat <- x[as.character(melt.next$taxon)]

# find annual covariance matrices aross iterations for species 
arguments3 <- list(.(iter), .(pft.cat, variable), .(year))
iter_mat <- reshape2::acast(melt.next, arguments3, mean)
cov.test <- apply(iter_mat,3,function(x){cov(x)})

arguments2 <- list(.(year, pft.cat), .(variable))
mean_mat <- reshape2::dcast(melt.next, arguments2, mean)

# Organize results into SDA usable format. Note: you will get some warnings after running
# this section. Don't worry that's normal. 
for(t in seq_along(obs.times)){
  
  # first mean matrices
  if(any(var.names == 'AGB.pft')){
    obs.mean.tmp[[t]] <- rep(NA, length(unique(melt.next$pft.cat)))
    names(obs.mean.tmp[[t]]) <- sort(unique(melt.next$pft.cat))
    for(r in seq_along(unique(melt.next$pft.cat))){
      k <- mean_mat[mean_mat$year==obs.times[t] & mean_mat$pft.cat==names(obs.mean.tmp[[t]][r]), variable]
      if(any(k)){
        obs.mean.tmp[[t]][r] <- k
      }
    }
  }
  
  # then, covariance matrices
  obs.cov.tmp[[t]] <- matrix(cov.test[,which(colnames(cov.test) %in% obs.times[t])],
                             ncol = sqrt(dim(cov.test)[1]),
                             nrow = sqrt(dim(cov.test)[1]))
  colnames(obs.cov.tmp[[t]]) <- names(iter_mat[1,,t]) 
  
}

obs.mean <- obs.mean.tmp
obs.cov <- obs.cov.tmp

# This data format is necessary for the output to work with the SDA workflow
names(obs.mean) <- paste0(obs.times,'-12-31')
names(obs.cov) <- paste0(obs.times,'-12-31')

# Save the final obs.list 
obs.list <- list(obs.mean = obs.mean, obs.cov = obs.cov)
save(obs.list,file=output)

################################################################################
################################################################################
################################################################################

## Visualize mean AGB estimates for all species at site over time
ggplot(mean_mat, aes(x = year, y = value, col = pft.cat)) +
  geom_point() + 
  geom_smooth(aes(group=pft.cat)) + 
  labs(title = 'Biomass by PFT', 
       x = 'Year', 
       y = 'Aboveground Biomass (kgC/m2)', 
       col = 'Species') 

