
## RUNNING SPIN-UP FOR ENSEMBLE-MET SDA MODEL RUNS

## Marissa Kivi 
## September 2019

# -------------------------------
# 0. Run workflow on interface
# -------------------------------

# Due to some VM quirks, we are unable to run the basic PEcAn workflow from the start; there are a few set-up 
# steps within the basic PEcAn workflow that we have to let the interface handle. If you're interested in what caused
# the issue, read the paragraph below. If you're not, then skip it and proceed to the steps that follow.

# The issue lies in how localhost is different between the Executor container and the RStudio Server. Thus, 
# if we have the host set up to be "docker", the workflow attempts to create a folder on the "remote docker server" 
# and fails doesn't create pecan.CHECKED.xml (see check.workflow.settings). If we set host to be "localhost", 
# it can't find the dbfiles for our input (see check.input) and pecan.CHECKED.xml remains incomplete. Therefore, 
# we can just let the interface do this step and that way, we do not have to mess with the main PEcAn functions. :) 

# First, we need to run the PEcAn web interface (paleon-pecan.virtual.crc.nd.edu:8000/pecan) with the exact same 
# characteristics that we want in our final spin-up model runs, except for the met. Choose the correct model, site, 
# species, number of ensembles, and dates. The dates in the interface should be the spin-up dates, and you can adjust 
# number of ensembles by checking the "Advanced Set-Up" box at the bottom of the page where you select species. 

# You should have inputted an input and file record for one of the site's met ensembles you will be using. Select that 
# met input for your interface run so that the dates are compatible and so that you do not have to change the dates
# in the XML files. We will add the actual ensembles later. 

# Run the interface once all of the desired settings are selected!

# --------------------------------
# 1. Set up working environment
# --------------------------------

library(PEcAn.settings)
library(PEcAn.uncertainty)
library(PEcAn.LINKAGES)
library(PEcAn.visualization)
library(PEcAn.assim.sequential)
library(PEcAn.remote)
library(PEcAn.visualization)
library(PEcAn.utils)
library(PEcAn.DB)
library(RCurl)
rm(list=ls())

# end whatever PEcAn processes might still be running
options(warn=1)
options(error=quote({
  PEcAn.utils::status.end("ERROR")
  PEcAn.remote::kill.tunnel(settings)
  if (!interactive()) {
    q(status = 1)
  }
}))

# set working directory to workflow folder from spinup 
workflowID = '14000000036'
setwd(paste0('/data/workflows/PEcAn_',workflowID))

# --------------------------------------------------
# 2. Enable writing and editing of workflow files
# --------------------------------------------------

# Next, we have to go to to the terminal tab in RStudio (below, next to "Console"). We are going to adjust the
# permissions of the files that the interface created so that we are able to edit and use these files within RStudio.
# This is another weird quirk about using PEcAn with Docker. The Executor container creates files as user "www-data" 
# and RStudio operates as user "carya", and they do not like to share apparently. 

### A. Go to the Terminal tab. 
### B. Paste the following commands in order, running each one separately. Pay attention to where you should be 
###    adjusting text for your own personal run. 
    # cd /data/workflows/PEcAn_<insert model ID> 
    # sudo chmod g+rw -R .
    # ---> type illinois for password 
### C. Return to the Console tab. 

# ----------------------------------
# 3. Edit pecan.TRAIT.xml file
# ----------------------------------

# Now, since we have enabled our ability to edit the files, it is time for us to edit the XML file so we can rerun the 
# spin-up with the appropriate met data. For this step, we will need the location of the ensemble weights of the met 
# ensembles we are going to be using, as well as the location of the met ensembles themselves. We should have placed 
# them in the /data folder in a previous step. We also need to know how many ensembles you are running. Adjust the 
# variables below.

ensemble_location = '/data/dbfiles/met_data/ROOSTER/weights/ensemble-weights-ROOSTER-prism.csv'
metdir <- '/data/dbfiles/met_data/ROOSTER/linkages/'
n = 25

### A. Get our sampled list of met ensembles. 

# load met ensemble weight files which contains model names, as well as weights and extract needed data
ens_wts <- read.csv(ensemble_location)
clim_mods <- ens_wts$climate_model
avg_wt <- ens_wts$wts

# randomly select n models from list of climate models based on their weights
# ANN: do we want to allow replacement here? do we just want to take the n most highly weighted ensembles?
clim_use <- sample(x = clim_mods, size = n, prob = avg_wt, replace = T)

# write the text to be added to the XML file; tabs have been added automatically so that the paste looks nicer
name <- numeric(length(clim_use))
for(i in 1:length(clim_use)){
  name[i] <- paste0('\t\t <path',i,'>',metdir,clim_use[i],'.Rdata</path',i,'>')
}
writeLines(name)

### B. Copy list of met paths into pecan.TRAIT.XML file near the bottom under <run>, <inputs>, <met>, and <path>. 
###    You don't have to remove the input ID number, but you should paste over the path that is currently in the file. 

### C. Adjust the ensemble method in the XML file to be "looping" and not "sampling." We already sampled the met and 
###    do not need to do it again. We just want the workflow to write the met paths to the ensembles in order. See below
###    for where you should edit. 

# <ensemble>
#   <size>25</size>
#   <variable>NPP</variable>
#   <samplingspace>
#   <parameters>
#   <method>uniform</method>
#   </parameters>
#   <met>
#   <method>sampling</method>  ### replace sampling with looping here 
#   </met>
#   </samplingspace>
#   <start.year>1860</start.year>
#   <end.year>1960</end.year>
# </ensemble>

### D. Now, we have to adjust the STATUS file so that the workflow will write over the current config files. In the 
###    workflow folder, open the STATUS file in RStudio and erase all of the rows except TRAIT and META. There is an 
###    example below. Once done, save the file and close it. 

# TRAIT	2019-09-04 19:28:32	2019-09-04 19:28:34	DONE	
# META	2019-09-04 19:28:34	2019-09-04 19:28:34	DONE	
# CONFIG	2019-09-04 19:28:34	2019-09-04 19:28:37	DONE	# Erase this row and all of the rows below it. 
# MODEL	2019-09-04 19:28:37	2019-09-04 19:28:48	DONE	
# OUTPUT	2019-09-04 19:28:48	2019-09-04 19:28:48	DONE	
# ENSEMBLE	2019-09-04 19:28:48	2019-09-04 19:28:48	DONE	
# FINISHED	2019-09-04 19:28:48	2019-09-04 19:28:48	DONE	

### E. Delete the "out" and "run" folders in the workflow ID so that they are correctly remade.

### F. Rerun the last parts of the PEcAn workflow 

# read in new adjusted settings
settings <- PEcAn.settings::read.settings("pecan.TRAIT.xml")

# write model specific configs with met ensembles
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

# start ecosystem model runs
if (PEcAn.utils::status.check("MODEL") == 0) {
  PEcAn.utils::status.start("MODEL")
  PEcAn.remote::runModule.start.model.runs(settings, stop.on.error = FALSE)
  PEcAn.utils::status.end()
}

# get results of model runs
if (PEcAn.utils::status.check("OUTPUT") == 0) {
  PEcAn.utils::status.start("OUTPUT")
  runModule.get.results(settings)
  PEcAn.utils::status.end()
}

# run ensemble analysis on model output. 
if ('ensemble' %in% names(settings) & PEcAn.utils::status.check("ENSEMBLE") == 0) {
  PEcAn.utils::status.start("ENSEMBLE")
  PEcAn.uncertainty::runModule.run.ensemble.analysis(settings, TRUE)    
  PEcAn.utils::status.end()
}

# the following two sections won't run unless you add the appropriate tags in the XML file
# run sensitivity analysis and variance decomposition on model output
if ('sensitivity.analysis' %in% names(settings) & PEcAn.utils::status.check("SENSITIVITY") == 0) {
  PEcAn.utils::status.start("SENSITIVITY")
  runModule.run.sensitivity.analysis(settings)
  PEcAn.utils::status.end()
}

# run parameter data assimilation
if ('assim.batch' %in% names(settings)) {
  if (PEcAn.utils::status.check("PDA") == 0) {
    PEcAn.utils::status.start("PDA")
    settings <- PEcAn.assim.batch::runModule.assim.batch(settings)
    PEcAn.utils::status.end()
  }
}

# run benchmarking
if("benchmarking" %in% names(settings) & "benchmark" %in% names(settings$benchmarking)){
  PEcAn.utils::status.start("BENCHMARKING")
  results <- papply(settings, function(x) calc_benchmark(x, bety))
  PEcAn.utils::status.end()
}

# workflow complete
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

### G. Visualize spin-up results. 

library(dplyr)

# read in and organize data 
n = as.numeric(settings$ensemble$size)
spp = length(settings$pfts)
years = lubridate::year(settings$run$start.date):lubridate::year(settings$run$end.date)

agb.total = matrix(0,n,length(years))
agb.pft.array = array(0,dim = c(n, length(years), spp))
temp.array = array(0, dim = c(n, length(years), 12))
precip.array = array(0, dim = c(n, length(years), 12))

for (i in 1:n){
  output_file = paste0(list.dirs(paste0(settings$outdir,'/out'), full.names = TRUE)[i+1],'/linkages.out.Rdata') 
  load(output_file)
  
  agb.total[i,] = ag.biomass[,1]
  agb.pft.array[i,,] = t(agb.pft[,,1])
  
  input_file = paste0(list.dirs(paste0(settings$rundir), full.names = TRUE)[i+1],'/linkages.input.Rdata')
  load(input_file)
  
  temp.array[i,,] = temp.mat
  precip.array[i,,] = precip.mat
}

# convert to data frames so we can use ggplot
agb.total.melt = reshape::melt(agb.total)
colnames(agb.total.melt) = c('ensemble','year','agb')
agb.pft.melt = reshape::melt(agb.pft.array)
colnames(agb.pft.melt) = c('ensemble','year','species','agb')
precip.melt = reshape::melt(precip.array)
colnames(precip.melt) = c('ensemble','year','month','precipitation')
precip.melt = precip.melt %>% mutate(date = year + 1859 + (1/month))
temp.melt = reshape::melt(temp.array)
colnames(temp.melt) = c('ensemble','year','month','temperature')
temp.melt = temp.melt %>% mutate(date = year + 1859 + (1/month))


# rename levels of species variable so easier to identify species
agb.pft.melt$species = as.factor(agb.pft.melt$species)
agb.pft.melt$species = plyr::mapvalues(agb.pft.melt$species, from = c(1,2,3,4), to = c('red maple', 'yellow birch', 'american beech', 'red oak'))

# plot 
library(ggplot2)
agb.pft.melt %>% ggplot(aes(x = year, y = agb, col = species)) + 
  geom_point() +
  facet_grid(species~., scales = 'free') + 
  labs(col = 'species', title = 'spin-up results: agb by pft')

agb.total.melt %>% ggplot(aes(x = year, y = agb, col = ensemble)) + 
  geom_point() +
  labs(title = 'spin-up results: total agb')

precip.melt %>% 
  ggplot(aes(x=date, y=precipitation, col = ensemble)) + 
  geom_line()

temp.melt %>%
  ggplot(aes(x=date, y=temperature, col = ensemble)) + 
  geom_line()
