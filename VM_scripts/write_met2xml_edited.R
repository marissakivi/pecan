#write met to xml format
#example: <path2>/home/araiho/linkages_ens_hf_met//bcc-csm1-1_001.01/climate.Rdata</path2> 

# load met ensemble weight files which contains model names 
# and average weights of ensmebles across all iterations and years
ens_wts <- read.csv('/data/dbfiles/met_data/HARVARD/weights/ensemble-weights-HARVARD-prism.csv')

# extract needed data
clim_mods <- ens_wts$climate_model
avg_wt <- ens_wts$wts

# following cna be used to randomly select n models from list of climate models
# based on their probability of occurrence 
clim_use <- sample(x = clim_mods,size = 300,prob = avg_wt,replace = T)
rbind(sort(avg_wt),sort(table(clim_use)))

metdir <- '/data/dbfiles/met_data/HARVARD/linkages/'

if(FALSE){
  files_out <- list.files(metdir)
  name <- numeric(length(files_out))
}

# hack to just use all met ensmbles once 
clim_use = clim_mods
name <- numeric(length(clim_use))

for(i in 1:length(clim_use)){
  name[i] <- paste0('<path',i+1,'>',metdir,clim_use[i],'/climate.Rdata</path',i+1,'>')
}

writeLines(name)

#copy and paste from console to inputs tag in xml