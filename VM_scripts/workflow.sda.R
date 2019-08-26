library(PEcAn.all)
library(PEcAn.SIPNET)
library(PEcAn.LINKAGES)
library(PEcAn.visualization)
library(PEcAn.ED2)
library(PEcAn.assim.sequential)
library(nimble)
library(lubridate)
library(PEcAn.visualization)
#PEcAn.assim.sequential::
library(rgdal) # need to put in assim.sequential
library(ncdf4) # need to put in assim.sequential

if(FALSE){
  ### sample parameters
  # Query the trait database for data and priors
  settings <- read.settings("pecan.CHECKED.xml")
  
  if (PEcAn.utils::status.check("TRAIT") == 0){
    PEcAn.utils::status.start("TRAIT")
    settings <- PEcAn.workflow::runModule.get.trait.data(settings)
    PEcAn.settings::write.settings(settings, outputfile='pecan.TRAIT.xml')
    PEcAn.utils::status.end()
  } else if (file.exists(file.path(settings$outdir, 'pecan.TRAIT.xml'))) {
    settings <- PEcAn.settings::read.settings(file.path(settings$outdir, 'pecan.TRAIT.xml'))
  }
  
  get.parameter.samples(settings, ens.sample.method = settings$ensemble$method) 
  
}

setwd('/home/carya/workflows/PEcAn_98000000048')
settings <- read.settings("pecan.SDA.xml")
load('/data/dbfiles/sda.obs.Rdata')
obs.mean <- obs.list$obs.mean
obs.cov <- obs.list$obs.cov

IC = NULL
Q = NULL
adjustment = TRUE
restart=F
control=list(trace=T,
             interactivePlot=T,
             TimeseriesPlot=T,
             BiasPlot=F,
             plot.title=NULL,
             debug=FALSE)

sda.enkf(settings, obs.mean, obs.cov, Q = NULL, restart=F, 
         control=list(trace=T,
                      interactivePlot=T,
                      TimeseriesPlot=T,
                      BiasPlot=F,
                      plot.title=NULL,
                      debug=F))

#xml changes
#need to change ens outside of sda
#need to change dates in run tags

## marissa's HF tree rings - /fs/data2/output/PEcAn_1000010299
## marissa's NRP tree rings - /fs/data2/output/PEcAn_1000010320
## marissa's RH tree rings - /fs/data2/output/PEcAn_1000010323

warnings()
