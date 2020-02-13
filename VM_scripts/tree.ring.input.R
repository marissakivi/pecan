
## Tree_ring_input.R 
# Author: Ann Raiho, Marissa Kivi 
# Last modified: January 2020

# This script is composed of two different parts. The first part takes data for all species and summarizes the data, so as to determine
# the top 98% species for the site in question. The second part reduces the results to only include results for the top 98% species and formats
# results into an SDA-usable format. 

################################################################################
################################################################################
################################################################################

## Part I: Determine 98% species for site

# set up working environment
library(plyr)
library(PEcAn.settings)
library(PEcAn.DB)
library(PEcAn.data.land)
library(dplyr)
library(ggplot2)

rm(list=ls())
site = 'NORTHROUND'

# load and prepare data for all species
input <- readRDS(paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/NPP_STAT_MODEL_',site,'.RDS'))

# keep only last 250 iterations - this step isn't always necessary
sort(unique(input$iter))
if (site == 'GOOSE') input <- input %>% filter(iter >= 1250)

# adjust "type" variable for consistency across sites
levels(input$type)
levels(input$type)[1] <- 'ab'
levels(input$type)[2] <- 'abi'

# condense input to only include species we are interested in and simplify dataframe
# non-NA values, one model, total aboveground biomass as variable
if (site == 'HARVARD'){
  input = input %>% filter(type == 'ab', model == 'Model RW + Census', !is.na(value)) %>% select(-model, -site_id, -type)
}else{
  input = input %>% filter(type == 'ab', !is.na(value)) %>% select(-model, -type)
}

# for North Round Pond site, we need to separate the data for Plots 1 & 2 and Plots 3 & 4 
if (site == 'NORTHROUND') input <- input %>% filter(plot %in% c(1,2)) %>% dplyr::select(-plot)

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
       ylab = 'Cumulative Proportion of Biomass')
ggsave(pl, filename = paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/sda_priority_spp_plot12.jpg'))
pl

# determine which species to keep for modeling (top 98% species are under red line)
species = c('TSCA','QURU','PIST','ACRU')

################################################################################
################################################################################
################################################################################

## Part II: Reduce and reformat results for running SDA at site 

# first, reduce data frame to only include top 98% taxa 
melt.next = melt.next %>% filter(taxon %in% species)
unique(melt.next$taxon)

# are all of these species in the BETY database for LINKAGES? 

# make sure taxa codes are same as those found in LINKAGES spp_matrix.csv file 
# 1.) this only needs to be done for the species in top 98% and is necessary to ensure there
# are matches for the BETY database
# 2.) check the README file for RDS file for site to double check common species names
# 3.) THIS SECTION WILL NEED TO BE ADJUSTED FOR ALL SITES 
levels(melt.next$taxon)
#levels(melt.next$taxon)[3] = 'BEAL2' # HF, NRP
#levels(melt.next$taxon)[2] = 'ACSA3' # NRP
#levels(melt.next$taxon)[7] = 'FRAM2' # NRP
#levels(melt.next$taxon)[4] = 'PIRU' # RH
species = unique(melt.next$taxon)
species

################################################################################
################################################################################
################################################################################

## Part III. Match up species to PFT on BETY database

# You will need to look at the PFT names that you are going to be using in the BETY 
# database and match them to the species. These names should be in your priors CSV 
# file! Traditionally, for LINKAGES, they follow the pattern: Genus.Species_Common.Name

x = c('Quercus.Rubra_Northern.Red.Oak.2',
      'Pinus.Strobus_White.Pine.2',
      'Acer.Rubrum_Red.Maple.2',
      'Tsuga.Canadensis_Hemlock.2')
names(x) = species

# adjust names of taxa to be with PFT name from BETY 
melt.next$pft.cat <- x[as.character(melt.next$taxon)]

# find annual covariance matrices aross iterations for species 
arguments3 <- list(.(iter), .(pft.cat, variable), .(year))
iter_mat <- reshape2::acast(melt.next, arguments3, mean)
cov.test <- apply(iter_mat,3,function(x){cov(x)})

arguments2 <- list(.(year, pft.cat), .(variable))
mean_mat <- reshape2::dcast(melt.next, arguments2, mean)

# Step 3: Correct formatting of input 

# organize results into SDA usable format
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

names(obs.mean) <- paste0(obs.times,'-12-31')
names(obs.cov) <- paste0(obs.times,'-12-31')

obs.list <- list(obs.mean = obs.mean, obs.cov = obs.cov)
save(obs.list,file=paste0('/Users/marissakivi/Desktop/PalEON/SDA/sites/',site,'/sda.obs.',site,'.plots34.Rdata'))

################################################################################
################################################################################
################################################################################

## Visualize results
ggplot(mean_mat, aes(x = year, y = value, col = pft.cat)) +
  geom_point() + 
  geom_smooth(aes(group=pft.cat)) + 
  labs(title = 'Biomass by PFT', 
       x = 'Year', 
       y = 'Aboveground Biomass (kgC/m2)', 
       col = 'Species') 

