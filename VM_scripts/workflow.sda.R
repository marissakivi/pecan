
## RUNNING SDA FOR ENSEMBLE-MET SDA MODEL RUNS

## Marissa Kivi 
## 12 September 2019

# ---------------------
# 0. Complete spin-up
# ---------------------

# Follow the R script for ensemble-based spin-up (workflow.spinup.R)

# If you completed the entire spin-up in the web interface, you will need to follow the second step in the spin-up 
# script so that you have write and read access to the files you will need to use in order to complete this workflow. 

rm(list=ls())

# -------------------------------
# 1. Set up working environment
# -------------------------------

# For this step, you will need the workflow ID from your spin-up. Adjust the variable below accordingly.  

ID = '14000000036'

# load necessary libraries
library(PEcAn.settings)
library(PEcAn.uncertainty)
library(PEcAn.SIPNET)
library(PEcAn.LINKAGES)
library(PEcAn.visualization)
library(PEcAn.ED2)
library(PEcAn.assim.sequential)
library(PEcAn.remote)
library(nimble)
library(lubridate)
library(PEcAn.visualization)
library(rgdal) # need to put in assim.sequential
library(ncdf4) # need to put in assim.sequential
library(dplyr)
library(dbplyr)

# set working directory to workflow info
setwd(paste0('/data/workflows/PEcAn_',ID))

# -------------------------------
# 2. Create pecan.SDA.xml file
# -------------------------------

### A. Open pecan.CONFIGS.xml from spin-up workflow folder.
### B. Immediately File > Save as > pecan.SDA.xml in the same folder. 
### C. Edit this newly-saved XML file for use in the SDA workflow.
###    Add <state.data.assimilation> tag with necessary settings information. You can just paste the following into your
###    xml file right after the very first <pecan> tag. You can uncomment it first by highlighting it and then by Code >
###    Comment/Uncomment Lines. Adjust the script as noted in the comments on the right side. Be sure to erase all of the 
###    comments!! 

  # <state.data.assimilation>
  #   <overwrite>TRUE</overwrite>
  #   <n.ensemble>30</n.ensemble>   # this needs to be changed to the correct number of ensembles members from spin-up
  #   <adjustment>TRUE</adjustment>
  #   <process.variance>FALSE</process.variance> # this should be TRUE!
  #   <sample.parameters>TRUE</sample.parameters>
  #   <inputs>
  #     <file>
  #       <input.id></input.id>  # we will load this input data later
  #       <path>
  #         <path></path>
  #       </path>
  #       <operator>direct</operator>
  #       <variable.id>1000000132</variable.id>  # if you are assimilating anything besides PFT aboveground biomass, change here
  #       <variable.name>
  #         <variable.name>AGB.pft</variable.name>
  #       </variable.name>
  #     </file>
  #   </inputs>
  #   <state.variables>
  #     <variable>
  #       <variable.name>AGB.pft</variable.name>  # ""
  #       <variable.id>1000000132</variable.id>  # ""
  #       <unit>MgC/ha/yr</unit>  # ""
  #       <min_value>0</min_value>   # the min and max values set the range for what is appropriate for state variable values
  #       <max_value>100000000</max_value>  # ""
  #     </variable>
  #   </state.variables>
  #   <spin.up>
  #     <start.date>1860/01/01</start.date>  # these dates should be the same as the start and end date of your spin-up (scroll down under <run>)
  #     <end.date>1960/12/31</end.date>
  #   </spin.up>
  #   <forecast.time.step>year</forecast.time.step>  # this setting controls how often the workflow is assimilating data
  #   <start.date>1960/01/01</start.date>  # first year of SDA - has to be the same year as last year of spin-up
  #   <end.date>2010/12/31</end.date> # last year of SDA
  # </state.data.assimilation>
        
### D. Save and close the file. 

# ---------------------------------------
# 3. Load settings and observation data
# ---------------------------------------

# Adjust the file location and name of the observation data in the following section to load the correct data
# to be assimilated into the model. 

# load transformed and reformatted observation data 
load('/data/dbfiles/sda.obs.RH.Rdata')
obs.mean <- obs.list$obs.mean
obs.cov <- obs.list$obs.cov

# load SDA xml as settings file
settings <- read.settings("pecan.SDA.xml")

# ----------------------
# 4. Run SDA workflow
# ----------------------

# steps of function
# 1. read settings 
# 2. splitting met for those models that need it 
# 3. tests for data assimilation 
# 4. set up for data assimilation
# 5. start loop through time


# set variables for workflow 
IC = NULL
Q = NULL
adjustment = TRUE
restart=F
control=list(trace=T,
             interactivePlot=T,
             TimeseriesPlot=T,
             BiasPlot=F,
             plot.title=NULL,
             debug=FALSE,
             pause = FALSE)

# sda.enkf(settings, obs.mean, obs.cov, Q = NULL, restart=F, 
#          control=list(trace=T,
#                       interactivePlot=T,
#                       TimeseriesPlot=T,
#                       BiasPlot=F,
#                       plot.title=NULL,
#                       debug=F))

if (control$debug) browser()

###-------------------------------------------------------------------###
### 1. read settings                                                  ###
###-------------------------------------------------------------------###
weight_list <- list()
adjustment <- settings$state.data.assimilation$adjustment
model      <- settings$model$type
write      <- settings$database$bety$write
defaults   <- settings$pfts
outdir     <- settings$modeloutdir # currently model runs locally, this will change if remote is enabled
rundir     <- settings$host$rundir
host       <- settings$host
forecast.time.step <- settings$state.data.assimilation$forecast.time.step  #idea for later generalizing
nens       <- as.numeric(settings$ensemble$size)
processvar <- as.logical(settings$state.data.assimilation$process.variance)
var.names <- sapply(settings$state.data.assimilation$state.variable, '[[', "variable.name")
names(var.names) <- NULL
input.vars <- sapply(settings$state.data.assimilation$inputs, '[[', "variable.name")
operators <- sapply(settings$state.data.assimilation$inputs, '[[', "operator")

# site details
# first col is the long second is the lat and row names are the site ids
site.ids <- settings$run$site$id
site.locs <- data.frame(Lon = as.numeric(settings$run$site$lon),
                        Lat = as.numeric(settings$run$site$lat))
colnames(site.locs) <- c("Lon","Lat")
rownames(site.locs) <- site.ids

# determine years for data assimilation
# start cut determines what is the best year to start spliting the met based on if we start with a restart or not.  
if (!is.null(restart)) {
  start.cut <-lubridate::ymd_hms(settings$state.data.assimilation$start.date, truncated = 3)-1
  Start.Year <-(lubridate::year(settings$state.data.assimilation$start.date)-1)
}else{
  start.cut <-lubridate::ymd_hms(settings$state.data.assimilation$start.date, truncated = 3)
  Start.Year <-(lubridate::year(settings$state.data.assimilation$start.date))
}
End.Year <-   lubridate::year(settings$state.data.assimilation$end.date) # years that assimilations will be done for - obs will be subsetted based on this

# filtering obs data based on years specifited in setting > state.data.assimilation
assimyears<-Start.Year:End.Year
obs.mean <- obs.mean[sapply(lubridate::year(names(obs.mean)), function(obs.year) obs.year %in% (assimyears))]
obs.cov <- obs.cov[sapply(lubridate::year(names(obs.cov)), function(obs.year) obs.year %in% (assimyears))]

# dir address based on the end date
if(!dir.exists("SDA")) dir.create("SDA",showWarnings = F)

# get model specific functions
do.call("library", list(paste0("PEcAn.", model)))
my.write_restart <- paste0("write_restart.", model)
my.read_restart <- paste0("read_restart.", model)
my.split_inputs  <- paste0("split_inputs.", model)

# double checking some of the inputs 
if (is.null(adjustment)) adjustment <- TRUE # unsure if this is necessary

###-------------------------------------------------------------------###
### 2. splitting/cutting the mets to the start and the end of SDA     ###
###-------------------------------------------------------------------### 

# models that don't need split_inputs, check register file for that
register.xml <- system.file(paste0("register.", model, ".xml"), package = paste0("PEcAn.", model))
register <- XML::xmlToList(XML::xmlParse(register.xml))
no_split <- !as.logical(register$exact.dates)
if (!exists(my.split_inputs)  &  !no_split) {
  PEcAn.logger::logger.warn(my.split_inputs, "does not exist")
  PEcAn.logger::logger.severe("please make sure that the PEcAn interface is loaded for", model)
  PEcAn.logger::logger.warn(
    my.split_inputs,
    "If your model does not need the split function you can specify that in register.Model.xml in model's inst folder by adding <exact.dates>FALSE</exact.dates> tag."
  )
}

if(!no_split){ 
  for(i in seq_along(settings$run$inputs$met$path)){
    
    ### model specific split inputs
    settings$run$inputs$met$path[[i]] <- do.call(my.split_inputs, 
                                                 args = list(settings = settings, 
                                                             start.time = start.cut, 
                                                             stop.time = lubridate::ymd_hms(settings$state.data.assimilation$end.date, truncated = 3, tz="UTC"),
                                                             inputs =  settings$run$inputs$met$path[[i]],
                                                             overwrite=T)) 
  }
}

###-------------------------------------------------------------------###
### 3. tests before data assimilation                                 ###
###-------------------------------------------------------------------###

# getting times where we have observations
# there will be an error if we need to fix them in next step
obs.times <- names(obs.mean)
obs.times.POSIX <- lubridate::ymd_hms(obs.times)

### TO DO: Need to find a way to deal with years before 1000 for paleon ### need a leading zero

for (i in seq_along(obs.times)) {
  if (is.na(obs.times.POSIX[i])) {
    if (is.na(lubridate::ymd(obs.times[i]))) {
      PEcAn.logger::logger.warn("Error: no dates associated with observations")
    } else {
      ### data does not have time associated with dates 
      ### adding 12:59:59PM assuming next time step starts one second later
      PEcAn.logger::logger.warn("Pumpkin Warning: adding one minute before midnight time assumption to dates associated with data")
      obs.times.POSIX[i] <- strptime(paste(obs.times[i], "23:59:59"),format="%Y-%m-%d %H:%M:%S",tz='UTC')#lubridate::ymd_hms(paste(obs.times[i], "23:59:59"))
    }
  }
}
obs.times <- obs.times.POSIX

###-------------------------------------------------------------------###
### 4. set up for data assimilation                                   ###
###-------------------------------------------------------------------###

# set up storage variables for DA 
nt          <- length(obs.times)
if (nt==0)     PEcAn.logger::logger.severe('There has to be at least one observation, before you can start the SDA code.')
FORECAST    <- ANALYSIS <- list()
enkf.params <- list()

# shape parameters estimated over time for process covariance
aqq         <- NULL
bqq         <- numeric(nt + 1)

# track range of state variables 
# interval remade everytime depending on data at time t
# state.interval stays constant and converts new.analysis to be within the correct bounds
interval    <- NULL 
state.interval <- cbind(as.numeric(lapply(settings$state.data.assimilation$state.variables, '[[', 'min_value')),
                        as.numeric(lapply(settings$state.data.assimilation$state.variables, '[[', 'max_value')))
rownames(state.interval) <- var.names

# parameters for ensembles should be created in main PEcAn workflow 
# set ensemble weights - either by Rdata file or sets all to constant 1 
if(!file.exists(file.path(settings$outdir, "ensemble_weights.Rdata"))){
  PEcAn.logger::logger.warn("ensemble_weights.Rdata cannot be found. Make sure you generate samples by running the get.ensemble.weights function before running SDA if you want the ensembles to be weighted.")
  #create null list
  for(tt in 1:length(obs.times)){
    weight_list[[tt]] <- rep(1,nens) #no weights
  }
} else{
  load(file.path(settings$outdir, "ensemble_weights.Rdata"))  ## loads ensemble.samples
}

# check for and load ensemble sample data
if(!file.exists(file.path(settings$outdir, "samples.Rdata"))) PEcAn.logger::logger.severe("samples.Rdata cannot be found. Make sure you generate samples by running the get.parameter.samples function before running SDA.")
load(file.path(settings$outdir, "samples.Rdata"))  ## loads ensemble.samples

# reformat parameters
new.params <- list()
for (i in seq_len(nens)) {
  new.params[[i]] <- lapply(ensemble.samples, function(x, n) {
    x[i, ] }, n = i)
} 

# define where time starts based on restart logical variable
# if this is a restart, pick up where we left last time    
if (restart){
  if(!file.exists(file.path(settings$outdir,"SDA", "sda.output.Rdata"))){
    PEcAn.logger::logger.warn("The SDA output from the older simulation doesn't exist.")
    t <- 1
  } else {
    load(file.path(settings$outdir,"SDA", "sda.output.Rdata"))
  }
  
  load(file.path(settings$outdir,"SDA", "outconfig.Rdata"))
  run.id <- outconfig$runs$id
  ensemble.id <- outconfig$ensemble.id
  
  # if you made it through the forecast and the analysis in t and failed on the analysis in t+1 so you didn't save t
  if(length(FORECAST) == length(ANALYSIS) && length(FORECAST) > 0) t = t + length(FORECAST) 
  
}else{
  t = 1
}

###-------------------------------------------------------------------###
### 5. start data assimilation loop                                   ###
###-------------------------------------------------------------------### 

for(t in t:nt){
  
  if (control$debug) browser()
  
  # check for observations in the given year and identify year
  obs <- which(!is.na(obs.mean[[t]]))
  obs.year <- lubridate::year(names(obs.mean)[t])
  
  ###----------------------------------------###
  ###  A. Checking if need to run forecast
  ###----------------------------------------###

  ## First: Do we have forecast output to compare to our data?
  
  # Check to see if SDA config file has been created and if there are run files for year 
  # If not, need to run forecast for year.  
  if(file.exists('run') & file.exists(file.path(settings$outdir,"SDA", "outconfig.Rdata"))){
    
    # need to load these in case during t==1 the analysis crashed so you have a forecast 
    # but didn't get to save the sda.output.Rdata
    load(file.path(settings$outdir,"SDA", "outconfig.Rdata")) 
    run.id <- outconfig$runs$id
    ensemble.id <- outconfig$ensemble.id
    if(t==1) inputs <- outconfig$samples$met 
    
    # looking for file for obs.year 
    sum_files <-
      sum(unlist(sapply(
        X = run.id,
        FUN = function(x){
          pattern = paste0(x, '/*.nc$')[1]
          grep(
            pattern = pattern,
            x = list.files(file.path(outdir,x), "*.nc$", recursive = F, full.names = T)
          )
        },
        simplify = T
      )))
    
  }else{
    sum_files <- 0 # if rundir or SDA outconfig hasn't been created yet 
  }
  
  ###----------------------------------###
  ###  B. Running forecast if needed
  ###----------------------------------###
  
  # if there are no files in the outdir for year, we need to run the forecast so we set up for that 
  if (sum_files == 0){
    
    ###--------------------------------------###
    ###  i. Setting up SV inputs if needed
    ###--------------------------------------###
    
    # WHY ARE INPUTS SPLIT IN THIS LOOP?
    
    # splitting the state variable input for models which need it 
    inputs.split <- list()
    if(!no_split & exists('outconfig')){
      for(i in seq_len(nens)){
        # model specific split inputs
        inputs.split$samples[i] <- do.call(
          my.split_inputs,
          args = list(
            settings = settings,
            start.time = (lubridate::ymd_hms(
              obs.times[t - 1], truncated = 3, tz = "UTC"
            )),
            stop.time = (lubridate::ymd_hms(
              obs.times[t], truncated = 3, tz = "UTC"
            )),
            inputs = inputs$samples[[i]]
          )
        )
      } 
      
    }else{
      # if t == 0 : we need to set SDA configs file
      if(t > 1) inputs.split <- inputs
    }
    
    ###-------------------------###
    ###  ii. Setting up restart
    ###-------------------------###
    
    ## Second: Has the analysis been run for the past year? 
    # If yes, then we restart from analysis. 
    # If not, we set it to null and start from the beginning. 

    if(exists('new.state')){
      restart.arg<-list(runid = run.id, 
                        start.time = lubridate::ymd_hms(obs.times[t - 1], truncated = 3),
                        stop.time = lubridate::ymd_hms(obs.times[t], truncated = 3), 
                        settings = settings,
                        new.state = new.state, 
                        new.params = new.params, 
                        inputs = inputs.split, 
                        RENAME = TRUE,
                        ensemble.id=ensemble.id)
    }else{ 
      restart.arg = NULL
    }
    
    ###-------------------------------------###
    ###  iii. Writing configs for model runs
    ###-------------------------------------###

    # Writing the config and submunning the model and reading the outputs for each ensemble
    outconfig <- write.ensemble.configs(defaults = settings$pfts, 
                                        ensemble.samples = ensemble.samples, 
                                        settings = settings,
                                        model = settings$model$type, 
                                        write.to.db = settings$database$bety$write,
                                        restart = restart.arg)
    
    # Save outconfig file for ensembles
    save(outconfig, file = file.path(settings$outdir,"SDA", "outconfig.Rdata"))
  
    # Extract config variables for analysis
    run.id <- outconfig$runs$id
    ensemble.id <- outconfig$ensemble.id
    # for any time after t==1, the met is the split met
    if(t==1) inputs <- outconfig$samples$met
    
    ###-------------------------------###
    ###  iv. Running model ensembles 
    ###-------------------------------###
    
    if(control$debug) browser()
    
    # submit model job submissions to RabbitMQ service
    # you will see a lot PEcAn output after running this command
        # - job submissions for each of the ensembles
        # - list of jobs that are to be completed with their folder names
        # - status bar tracking progress of those model runs until they are done
            # (the function continues to check statuses until they are all done)
    
    PEcAn.remote::start.model.runs(settings, settings$database$bety$write)
  } 
  
  ###--------------------------###
  ###  C. Reading the output
  ###--------------------------###
  
  # extract state variable forecasts from all ensembles for analysis (X)
  # extract ensemble parameters data (new.params)
  
  # the following chunk will print all of the state variable forecasts for each of the ensembles

  Sys.sleep(60)
  
  X_tmp <- vector("list", 2)
  X <- list()
  for (i in seq_len(nens)) {
    X_tmp[[i]] <- do.call(
      my.read_restart,
      args = list(
        outdir = outdir,
        runid = run.id[i],
        stop.time = obs.times[t],
        settings = settings,
        var.names = var.names,
        params = new.params[[i]]
      )
    )
    
    # we also want to carry some deterministic relationships to write_restart
    # these will be stored in params
    X[[i]]      <- X_tmp[[i]]$X
    if (!is.null(X_tmp[[i]]$params))
      new.params[[i]] <- X_tmp[[i]]$params
  }
  
  #changing the extension of nc files to a more specific date related name
  files <-  list.files(
    path = file.path(settings$outdir, "out"),
    "*.nc$",
    recursive = TRUE,
    full.names = TRUE)
  files <-  files[grep(pattern = "SDA*", files, invert = TRUE)]
  
  file.rename(files, 
              file.path(dirname(files), 
                        paste0("SDA_", basename(files), "_", gsub(" ", "", names(obs.mean)[t]), ".nc") ) )
  
  # set up variables for analysis
  X <- do.call(rbind, X)
  
  # check to make sure there are successful forecasts 
  if(sum(X,na.rm=T) == 0){
    logger.severe(paste('NO FORECAST for',obs.times[t],'Check outdir logfiles or read restart. Do you have the right variable names?'))
  }
  
  ###------------------------------###
  ###  D. Preparing observations
  ###------------------------------###
  
  if (any(obs)) {
    
    # finding obs data - which type of observation do we have at this time point?
    input.order <- sapply(input.vars, agrep, x=names(obs.mean[[t]]))
    names(input.order) <- operators 
    input.order.cov <- sapply(input.vars, agrep, x=colnames(obs.cov[[t]]))
    names(input.order.cov) <- operators
    choose <- unlist(sapply(colnames(X), agrep, x=names(obs.mean[[t]]), max=1, USE.NAMES = F))
    choose.cov <- unlist(sapply(colnames(X), agrep, x=colnames(obs.cov[[t]]), max=1, USE.NAMES = F))
    
    # I'm unsure what this section is for 
    if(!any(choose)){
      choose <- unlist(input.order)
      choose <- order(names(obs.mean[[t]]))
      choose.cov <- unlist(input.order.cov)
      choose.cov <- order(colnames(obs.cov[[t]]))
    }
    
    # dropping observations with NA mean values   
    na.obs.mean <- which(is.na(unlist(obs.mean[[t]][choose])))
    if (length(na.obs.mean) > 0) choose <- choose [-na.obs.mean]
    Y <- unlist(obs.mean[[t]][choose]) 
    
    # set up covariance matrix for desired species 
    R <- as.matrix(obs.cov[[t]][choose.cov,choose.cov])
    R[is.na(R)] <- 0.1
    
    if (control$debug) browser()
    
    # making the mapping matrix between observed data and their forecast state variables
    # TO DO: doesn't work unless it's one to one
    #if(length(operators)==0) H <- Construct_H(choose, Y, X)
    H <- Construct_H(choose, Y, X)
    
    ###------------------------###
    ###  E. Running analysis 
    ###------------------------###
    
    # define method of analysis 
    if(processvar == FALSE){
      an.method<-EnKF  
    }else{    
        an.method <- GEF
    }  
    
    # define extra args
    if (processvar && t > 1) {
      aqq <- enkf.params[[t-1]]$aqq
      bqq <- enkf.params[[t-1]]$bqq
      X.new<-enkf.params[[t-1]]$X.new
    }
    if(!exists('Cmcmc_tobit2space') | !exists('Cmcmc')) {
      recompile = TRUE
    }else{
      recompile = FALSE
    }
    
    # get weights 
    wts <- unlist(weight_list[[t]][outconfig$samples$met$ids])
    
    # run analysis function
    enkf.params[[t]] <- Analysis.sda(settings,
                                     FUN=an.method,
                                     Forecast=list(Q=Q, X=X),
                                     Observed=list(R=R, Y=Y),
                                     H=H,
                                     extraArg=list(aqq=aqq, bqq=bqq, t=t,
                                                   recompile=recompile,
                                                   wts = wts),
                                     nt=nt,
                                     obs.mean=obs.mean,
                                     obs.cov=obs.cov)
    
    # reading back analysis state variable variables for forecast . . 
    FORECAST[[t]] <- X
    mu.f <- enkf.params[[t]]$mu.f
    Pf <- enkf.params[[t]]$Pf
    # and analysis
    Pa <- enkf.params[[t]]$Pa
    mu.a <- enkf.params[[t]]$mu.a
    
    # hack for zero variance
    diag(Pf)[which(diag(Pf) == 0)] <- 0.1 
    
    # extracting extra outputs from analysis 
    if (processvar) {
      aqq <- enkf.params[[t]]$aqq
      bqq <- enkf.params[[t]]$bqq
      X.new <- enkf.params[[t]]$X.new
    }
    
    ###-------------------------------------------###
    ###  F. Writing trace for tracking progress
    ###-------------------------------------------###     

    if (control$trace) {
      PEcAn.logger::logger.info ("\n --------------------------- ",
                                 obs.year,
                                 " ---------------------------\n")
      PEcAn.logger::logger.info ("\n --------------Obs mean----------- \n")
      print(Y)
      PEcAn.logger::logger.info ("\n --------------Obs Cov ----------- \n")
      print(R)
      PEcAn.logger::logger.info ("\n --------------Forecast mean ----------- \n")
      print(enkf.params[[t]]$mu.f)
      PEcAn.logger::logger.info ("\n --------------Forecast Cov ----------- \n")
      print(enkf.params[[t]]$Pf)
      PEcAn.logger::logger.info ("\n --------------Analysis mean ----------- \n")
      print(t(enkf.params[[t]]$mu.a))
      PEcAn.logger::logger.info ("\n --------------Analysis Cov ----------- \n")
      print(enkf.params[[t]]$Pa)
      PEcAn.logger::logger.info ("\n ------------------------------------------------------\n")
    }
    if (control$debug) browser()
    if (control$pause) readline(prompt="Press [enter] to continue \n")
    
  } else {
    
    mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
    Pf <- cov(X)

    # if no process variance, forecast is the same as the analysis 
    # if yes, no data 
    if (!processvar) {
      mu.a <- mu.f
      Pa   <- Pf + Q
    } else {
      mu.a <- mu.f
      if(!exists('q.bar')){
        q.bar <- diag(ncol(X))
        PEcAn.logger::logger.info('Process variance not estimated. Analysis has been given uninformative process variance')
      } 
      Pa   <- Pf + solve(q.bar)
    }
    enkf.params[[t]] <- list(mu.f = mu.f, Pf = Pf, mu.a = mu.a, Pa = Pa)
  }
  
  ###-------------------------------------------###
  ###  G. Adjust/update state variable matrix
  ###-------------------------------------------### 
  
  # start adjustment of state variable matrix 
  # this step gives weights to the different ensemble members based on their likelihood given data
  # then it adjusts the analysis mean estimates based on these weights 
  if(adjustment){
    # if we have process var, then x is x.new
    if (processvar & exists('X.new')){
      X.adj.arg <- X.new 
    }else{ 
      X.adj.arg <- X
      print('Using X not X.new. GEF was skipped or process variance == FALSE')
      }
    analysis <- adj.ens(Pf, X.adj.arg, mu.f, mu.a, Pa)
  }else{
    analysis <- as.data.frame(rmvnorm(as.numeric(nrow(X)), mu.a, Pa, method = "svd"))
  }
  
  colnames(analysis) <- colnames(X)
  
  # map analysis vectors to be in bounds of state variables
  if(processvar){
    for(i in 1:ncol(analysis)){
      int.save <- state.interval[which(startsWith(colnames(analysis)[i],
                                                  var.names)),]
      analysis[analysis[,i] < int.save[1],i] <- int.save[1]
      analysis[analysis[,i] > int.save[2],i] <- int.save[2]
    }
  }
  
  new.state  <- as.data.frame(analysis)
  ANALYSIS[[t]] <- analysis
  
  ###---------------------###
  ###  H. Save outputs
  ###---------------------### 
  
  # keep obs data and settings for later visualization in Dashboard
  Viz.output <- list(settings, obs.mean, obs.cov) 
  save(site.locs, t, X, FORECAST, ANALYSIS, enkf.params, new.state, new.params, run.id,
       ensemble.id, ensemble.samples, inputs, Viz.output,  file = file.path(settings$outdir,"SDA", "sda.output.Rdata"))
  
  # interactive plotting 
  if (t > 1 & control$interactivePlot) { #
    print(interactive.plotting.sda(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS))
  }
} 
# end loop over time
 
if (control$debug) browser()

###-------------------------------------------------------------------###
### time series plots                                                 ###
###-------------------------------------------------------------------###----- 
if(control$TimeseriesPlot) post.analysis.ggplot(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS,plot.title=control$plot.title)
if(control$TimeseriesPlot) PEcAn.assim.sequential:::post.analysis.ggplot.violin(settings, t, obs.times, obs.mean, obs.cov, obs, X, FORECAST, ANALYSIS)
###-------------------------------------------------------------------###
### bias diagnostics                                                  ###
###-------------------------------------------------------------------###----
if(control$BiasPlot)   PEcAn.assim.sequential:::postana.bias.plotting.sda(settings,t,obs.times,obs.mean,obs.cov,obs,X,FORECAST,ANALYSIS)
###-------------------------------------------------------------------###
### process variance plots                                            ###
###-------------------------------------------------------------------###----- 
if (processvar & control$BiasPlot) postana.bias.plotting.sda.corr(t,obs.times,X,aqq,bqq)

# end of SDA analysis function


warnings()
