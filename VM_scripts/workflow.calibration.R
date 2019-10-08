# SDA Process: Calibration
# Marissa Kivi 
# October 2019

# I. Use script to update priors in BETY database
# script => update_bety_priors.R

# II. Run dry run on PEcAn interface
library(PEcAn.all)
library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(RColorBrewer)

workflowID = '14000000022'
setwd(paste0('/data/workflows/PEcAn_',workflowID))

sppname = c('red maple', 'yellow birch','beech','red oak')
sppcol = c('maroon','gold','darkgreen','orange')

settings = read.settings('pecan.CONFIGS.xml')
n = as.numeric(settings$ensemble$size)
nspp = length(settings$pfts)
nyear = 200
runs = list.dirs('./out', full.names = FALSE)[-1]

# III. Visual check 

dbh.array = array(0,dim=c(200,nyear,n))
ntrees.array = array(0,dim=c(nspp,nyear,n))
spp.list = c()

# gather information
for (i in 1:n){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  dbh.array[,,i] = dbh.save[,,1]
  ntrees.array[,,i] = ntrees.kill[,,1]
  for (j in 1:nyear){
    for (s in 1:nspp){
      spp.list = c(spp.list, rep(s, ntrees.kill[s,j,1]))
    }
  }
}

# put into one dataframe
dbh.melt = reshape::melt(dbh.array) %>% filter(value != 0)
tree.melt = dbh.melt
colnames(tree.melt) = c('index','year','ensemble','dbh')
tree.melt$species = spp.list

# bin diameters
bins.melt = tree.melt %>% 
  mutate(dbh = as.factor(round(dbh))) %>% 
  group_by(year, ensemble, species, dbh) %>%
  summarize(num = n()) %>%
  ungroup() %>%
  mutate(dbh = as.numeric(dbh), 
         species = as.factor(species))

# 1. diameter plots
tree.melt %>% mutate(species = as.factor(species)) %>%
  ggplot(aes(x = year, y = dbh, col = species)) + 
  geom_point() + 
  facet_wrap(~ensemble) + 
  labs(title='Diameter vs. Year') +
  scale_color_manual(
    values = sppcol, 
    limits = c(1,2,3,4), 
    name = 'species',
    labels = sppname
  ) + 
  geom_hline(yintercept=10)

# 2. diameter bins plots
bins.melt %>% ggplot(aes(x = year, y = dbh, col = species, size = num)) +
  geom_point() + 
  facet_wrap(~ensemble) + 
  scale_size_continuous(
    range = c(0.01,1.5)
  ) +
  scale_color_manual(
    values = sppcol, 
    limits = c(1,2,3,4), 
    name = 'species',
    labels = sppname
  ) + 
  geom_hline(yintercept=10)

# 3. circle plots

j = 175 # you can change year of interest here
pl.diam = 3256.7 #cm
ids = c('origin', sapply(c(1:n),toString)) # initialize origin and ensemble nodes
parent = c(rep('origin',n)) 
child = c(sapply(c(1:n),toString))
spp.ids = c(NA, rep('plot', n))
dbh.ids = c(NA,rep(pl.diam, n))

# loop through each ensemble and add data 
for (i in 1:n){
  # get data for year from ensemble for all species
  data = tree.melt %>% filter(ensemble == i, year == j) %>% mutate(species = as.numeric(species)) 
  num = length(data$dbh)
  spp.ids = c(spp.ids, data$species)
  dbh.ids = c(dbh.ids, data$dbh)
  parent = c(parent, rep(i,num))
  child = c(child, paste0(i,'-',sapply(data$index, toString)))
  ids = c(ids, paste0(i,'-',sapply(data$index, toString)))
}
    
# organize data into hierarchical format
edges <- data.frame(parent = parent, child = child)
vertices <- data.frame(id = ids, species=spp.ids, size=dbh.ids)
hData <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
    
# create plot with hierarchical data
pl <- ggraph(hData, layout = 'circlepack', weight="size") +
  geom_node_circle(aes(fill=species)) +
  scale_fill_manual(
    values = sppcol, 
    limits = c(1,2,3,4), 
    name = 'species',
    labels = sppname
    ) +
  ggtitle(paste0('Stand Structure Plot in Year ',j)) + 
  theme_void()
pl

# IV. Investigate model output for weird behaviors 

# is the species unable to regenerate/seed? 

# check: 
birth.array = array(0,dim=c(nspp,nyear,n))
for (i in 1:n){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  birth.array[,1,i] = ntrees.birth[,1,1]
  for (j in 2:nyear){
    birth.array[,j,i] = ntrees.birth[,j,1] - ntrees.kill[,(j-1),1]
  }
}
birth.melt = reshape::melt(birth.array)
colnames(birth.melt) = c('species','year','ensemble','births')

birth.melt %>% 
  ggplot(aes(x=year,y=births,col=species, group=species)) +
  geom_line() + 
  facet_wrap(~ensemble)

# diagnose: 
gmult.array = array(0,dim=c(nspp,4,nyear,n))
jtemp.mat = matrix(0,nyear,n)
frost.mat = matrix(0,nspp,n)
mplant.mat = matrix(0,nspp,n)
al.mat = matrix(1,nyear,n)
nplant.array = array(0,dim=c(nspp,nyear,n))

# gather data
for (i in 1:n){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  load(paste0('./run/',runs[i],'/linkages.input.Rdata'))
  gmult.array[,,,i] = gf.vec.save[,,,1]
  jtemp.mat[,i] = temp.mat[,1]
  frost.mat[,i] = spp.params$FROST
  mplant.mat[,i] = spp.params$MPLANT
  
  slite = 1.5 * (1 - exp(-1.136*(1.0-.08)))
  yfl = runif(1,0,1)
  nplant = spp.params$MPLANT * slite * gmult.array[,2,,i] * gmult.array[,4,,i] * yfl

  # #calculate leaf weight in G/plot and leaf area index (from birth subroutine)
  # for (j in 2:nyear){
  #   
  #   #initialize foliage biomass (folw) and foliage area (fola)
  #   folw = 0
  #   fola = 0
  #   nl = 1
  #   
  #   for(s in 1:nspp){
  #     if(ntrees.kill[s,(j-1),1] == 0) next
  #     nu = nl + ntrees.kill[s,(j-1),1] - 1
  #     ret = spp.params$FRT[s]
  #     for(k in nl:nu){
  #       age = iage.save[k,(j-1),1]
  #       if(age<ret){
  #         ret = age
  #         folw = folw + ((spp.params$SLTA[s]+spp.params$SLTB[s]*dbh.save[k,(j-1),1])/2)^2*(3.14*spp.params$FWT[s]*ret)
  #         fola = fola + ((1.9283295) * 10^-4)*((dbh.save[k,(j-1),1])^2.129)
  #       }
  #     }
  #     nl = nl + ntrees.kill[s,(j-1),1]
  #   }
  #   
  #   fola[is.na(fola)] <- 0
  #   #calculate amount of light at forest floor
  #   al.mat[j,i] = 1 * exp(-folw/93750)
  # }

}

# prepare for plotting
gmult.melt = reshape::melt(gmult.array)
colnames(gmult.melt) = c('species','gmult','year','ensemble','value')
birth.limit.melt = birth.melt %>% 
  right_join(gmult.melt, by = c('species','year','ensemble')) %>% 
  filter(births < 1,
         gmult != 3)
birth.limit.melt$species = as.factor(birth.limit.melt$species)
birth.limit.melt$species = plyr::mapvalues(birth.limit.melt$species, from = c(1,2,3,4), to = sppname)
birth.limit.melt$gmult = as.factor(birth.limit.melt$gmult)
birth.limit.melt$gmult = plyr::mapvalues(birth.limit.melt$gmult, 
                                   from = c(1,2,4), 
                                   to = c('light', 'moisture','degree day'))

jtemp.melt = reshape::melt(jtemp.mat)
colnames(jtemp.melt) = c('year','ensemble','jan.temp')
frost.melt = reshape::melt(frost.mat)
colnames(frost.melt) = c('species','ensemble','frost')
limit.melt = left_join(frost.melt, jtemp.melt, by = c('ensemble')) %>% 
  mutate(limit = (frost > jan.temp))

# this plot will show the number of ensembles in which each species is unable to seed due to
# the frost condition each year 
limit.melt %>% 
  filter(limit) %>%
  ggplot(aes(x = year, fill = species)) +
  facet_wrap(~ensemble) + 
  geom_histogram(binwidth = 1)

# these plots will show each growth factor for each species across time for all the years 
# where birth was less than 1 birth/year 
birth.limit.melt %>% ggplot(aes(x = year, y = value, col = value)) + 
  geom_point() + 
  facet_grid(gmult~species) + 
  scale_color_continuous(low = 'red', high = 'green')

