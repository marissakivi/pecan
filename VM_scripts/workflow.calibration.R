# SDA Process: Calibration
# Marissa Kivi
# October 2019

# I. Use script to update priors in BETY database
# script => update_bety_priors.R

# II. Run dry run on PEcAn interface
rm(list=ls())

# III. Set run-specific variables
library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(PEcAn.all)

workflowID = '14000000085'
setwd(paste0('/data/workflows/PEcAn_',workflowID))

plot.dir = paste0('/data/workflows/PEcAn_',workflowID,'/calibration')
if (!dir.exists(plot.dir)) dir.create(plot.dir)

# HARVARD
#sppname = c('red maple', 'yellow birch','american beech','red oak','eastern hemlock') #HARVARD
#sppcol = c('purple','gold','blue','red','darkgreen') #HARVARD

# ROOSTER
#sppname = c('red maple','red spruce','white pine','red oak') #ROOSTER
#sppcol = c('purple','pink','green','red') #ROOSTER

# GOOSE
#sppname = c('red maple', 'american beech', 'white pine','white oak', 'red oak')
#sppcol = c('purple', 'blue', 'green', 'black', 'red')

# NORTHROUND 1 & 2
#sppname = c('sugar maple','yellow birch', 'american beech','white ash','eastern hemlock')
#sppcol = c('magenta','gold','blue','brown','darkgreen')

# NORTHROUND 3 & 4 
#sppname = c('red maple','white pine','red oak','eastern hemlock')
#sppcol = c('purple','green','red','darkgreen')

settings = read.settings('pecan.CONFIGS.xml')
nens = as.numeric(settings$ensemble$size)
n = nens
nspec = length(settings$pfts)
nyear = 100
runs = list.dirs('./out', full.names = FALSE)[-1]

# IV. Visual check

dbh.array = array(0,dim=c(200,nyear,nens))
ntrees.array = array(0,dim=c(nspec,nyear,nens))
spp.list = c()

# gather information
for (i in 1:nens){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  dbh.array[,,i] = dbh.save[,,1]
  ntrees.array[,,i] = ntrees.kill[,,1]
  for (j in 1:nyear){
    for (s in 1:nspec){
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
pl1 = tree.melt %>% mutate(species = as.factor(species)) %>%
  filter(species %in% c(4,3,2)) %>%
  ggplot(aes(x = year, y = dbh, col = species)) +
  geom_point(size = 0.1) +
  facet_wrap(~ensemble) +
  labs(title='Diameter vs. Year') +
  scale_color_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  ) +
  geom_hline(yintercept=10)

ggsave(plot = pl1, file.path(plot.dir,'diameter.plot.jpg'))

# 2. diameter bins plots
pl2 = bins.melt %>% ggplot(aes(x = year, y = dbh, col = species, size = num)) +
  geom_point() +
  facet_wrap(~ensemble) +
  scale_size_continuous(
    range = c(0.01,1.5)
  ) +
  scale_color_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  ) +
  geom_hline(yintercept=10)

ggsave(plot = pl2, file.path(plot.dir,'diameter.bins.plot.jpg'))

# 3. circle plots

j = 100 # you can change year of interest here
pl.diam = 3256.7 #cm
ids = c('origin', sapply(c(1:nens),toString)) # initialize origin and ensemble nodes
parent = c(rep('origin',nens))
child = c(sapply(c(1:nens),toString))
spp.ids = c(NA, rep('plot', nens))
dbh.ids = c(NA,rep(pl.diam, nens))

# loop through each ensemble and add data
for (i in 1:nens){
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
pl3 <- ggraph(hData, layout = 'circlepack', weight="size") +
  geom_node_circle(aes(fill=species)) +
  scale_fill_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
    ) +
  ggtitle(paste0('Stand Structure Plot in Year ',j)) +
  theme_void()
ggsave(plot = pl3, file.path(plot.dir,'circle.plot.jpg'))

# print maximum diameter reached by each species 
tree.melt %>% group_by(species) %>% summarize(maxdia = max(dbh)) 

# V. Investigate model output for weird behaviors

# is the species unable to regenerate/seed?

# check:
birth.array = array(0,dim=c(nspec,nyear,nens))
for (i in 1:nens){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  birth.array[,1,i] = ntrees.birth[,1,1]
  for (j in 2:nyear){
    birth.array[,j,i] = ntrees.birth[,j,1] - ntrees.kill[,(j-1),1]
  }
}
birth.melt = reshape::melt(birth.array)
colnames(birth.melt) = c('species','year','ensemble','births')
birth.melt$species = as.factor(birth.melt$species)

birth.melt %>%
  ggplot(aes(x=year,y=births,col=species, group=species)) +
  geom_line() +
  facet_wrap(~ensemble) +
  scale_color_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  )

# diagnose:
gmult.array = array(0,dim=c(nspec,4,nyear,n))
jtemp.mat = matrix(0,nyear,n)
frost.mat = matrix(0,nspec,n)
mplant.mat = matrix(0,nspec,n)
al.mat = matrix(1,nyear,n)
nplant.array = array(0,dim=c(nspec,nyear,n))

# gather data
for (i in 1:nens){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  load(paste0('./run/',runs[i],'/linkages.input.Rdata'))
  gmult.array[,,,i] = gf.vec.save[,,,1]
  jtemp.mat[,i] = temp.mat[,1]
  frost.mat[,i] = spp.params$FROST
  mplant.mat[,i] = spp.params$MPLANT

  slite = 1.5 * (1 - exp(-1.136*(1.0-.08)))
  yfl = runif(1,0,1)
  nplant = spp.params$MPLANT * slite * gmult.array[,2,,i] * gmult.array[,4,,i] * yfl
}

# prepare for plotting
gmult.melt = reshape::melt(gmult.array)
colnames(gmult.melt) = c('species','gmult','year','ensemble','value')
gmult.melt$species = as.integer(gmult.melt$species)
birth.melt$species = as.integer(birth.melt$species)
birth.limit.melt = birth.melt %>%
  right_join(gmult.melt, by = c('species','year','ensemble')) %>%
  filter(births < 1,
         gmult != 3)
birth.limit.melt$species = as.factor(birth.limit.melt$species)
birth.limit.melt$species = plyr::mapvalues(birth.limit.melt$species, from = c(1:nspec), to = sppname)
birth.limit.melt$gmult = as.factor(birth.limit.melt$gmult)
birth.limit.melt$gmult = plyr::mapvalues(birth.limit.melt$gmult,
                                   from = c(1,2,4),
                                   to = c('light', 'moisture','degree day'))

jtemp.melt = reshape::melt(jtemp.mat)
colnames(jtemp.melt) = c('year','ensemble','jan.temp')
frost.melt = reshape::melt(frost.mat)
colnames(frost.melt) = c('species','ensemble','frost')
frost.melt$species = as.factor(frost.melt$species)
limit.melt = left_join(frost.melt, jtemp.melt, by = c('ensemble')) %>%
  mutate(limit = (frost > jan.temp))

# this plot will show the number of ensembles in which each species is unable to seed due to
# the frost condition each year
limit.melt %>%
  filter(limit) %>%
  ggplot(aes(x = year, fill = species)) +
  facet_wrap(~ensemble) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  )

# these plots will show each growth factor for each species across time for all the years
# where birth was less than 1 birth/year
birth.limit.melt %>% ggplot(aes(x = year, y = value, col = value)) +
  geom_point() +
  facet_grid(gmult~species) +
  scale_color_continuous(low = 'red', high = 'green')

# is the species unable to grow?

# check:
nogro.array = array(0,dim=c(nspec,nyear,n))
for (i in 1:nens){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  for (j in 1:nyear){
    nu = 1
    for (k in 1:nspec){
      nl = nu + ntrees.kill[k,j,1] - 1
      nogro = sum(nogro.save[nu:nl,j,1] < 0)
      nogro.array[k,j,i] = nogro/ntrees.kill[k,j,1]
      nu = nl + 1
    }
  }
}

nogro.melt = reshape::melt(nogro.array)
colnames(nogro.melt) = c('species','year','ensemble','count')
nogro.melt$species = as.factor(nogro.melt$species)

nogro.melt %>%
  ggplot(aes(x=year,y=count,col=species, group=species)) +
  geom_line() +
  facet_wrap(~ensemble) +
  scale_color_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  )

# diagnose:

gf.array = array(NA,dim=c(200,4,nyear,nens))
tree.array = array(NA,dim=c(200,nyear,nens))

# growth dictated by growth factors
for (i in 1:nens){
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  for (j in 1:nyear){
    # get sppid and growth factors
    sppid = c(rep(1,ntrees.birth[1,j,1]))
    algf = c(algf.save.keep[which(!is.na(algf.save.keep[,1,j,1])),1,j,1])
    smgf = c(rep(gf.vec.save[1,2,j,1], ntrees.birth[1,j,1]))
    degdgf = c(rep(gf.vec.save[1,4,j,1], ntrees.birth[1,j,1]))
    sngf = c(rep(gf.vec.save[1,3,j,1], ntrees.birth[1,j,1]))

    for (k in 2:nspec){
      sppid = c(sppid, rep(k,ntrees.birth[k,j,1]))
      algf = c(algf, algf.save.keep[which(!is.na(algf.save.keep[,k,j,1])),k,j,1])
      smgf = c(smgf,rep(gf.vec.save[k,2,j,1], ntrees.birth[k,j,1]))
      degdgf = c(degdgf,rep(gf.vec.save[k,4,j,1], ntrees.birth[k,j,1]))
      sngf = c(sngf,rep(gf.vec.save[k,3,j,1], ntrees.birth[k,j,1]))
    }
    nogro.inds = which(nogro.save[,j,1] < 0)
    t = length(algf)
      algf[is.na(algf)] = smgf[is.na(smgf)] = sngf[is.na(sngf)] = degdgf[is.na(degdgf)] = 0
      gf.array[1:t,1,j,i] = algf#[nogro.inds] # algf
      gf.array[1:t,2,j,i] = smgf#[nogro.inds] # smgf
      gf.array[1:t,3,j,i] = sngf#[nogro.inds] # sngf
      gf.array[1:t,4,j,i] = degdgf#[nogro.inds] # degdgf
      tree.array[1:t,j,i] = sppid#[nogro.inds] # species id
    #}
  }
}

gf.melt = reshape::melt(gf.array)
colnames(gf.melt) = c('tree','gf','year','ensemble','value')
tree.melt = reshape::melt(tree.array)
colnames(tree.melt) = c('tree','year','ensemble','species')
lgf.melt = left_join(gf.melt,tree.melt) %>% 
  filter(!is.na(value)) %>%
  group_by(ensemble,year,tree,species) %>%
  summarize(lgf = ifelse(any(is.na(value)),NA,which.min(value)),
            lgf.val = ifelse(any(is.na(value)),NA,min(value)))

# check to make sure smgf=0 is not degdgf=0 
# in LINKAGES, if degdgf=0, than smgf and sngf= 0, too 
zeros.id = which(lgf.melt$lgf.val==0)

for (i in zeros.id){
  
  ens = lgf.melt$ensemble[i]
  yr = lgf.melt$year[i]
  tr = lgf.melt$tree[i]
  
  entry = gf.melt %>% filter(tree == tr, ensemble == ens, year == yr)
  
  # check if degdgf is zero, too
  if (entry %>% filter(gf==4) %>% dplyr::select(value) == 0){
    lgf.melt$lgf[i] = 4
  }
}

lgf.melt$lgf = as.factor(lgf.melt$lgf)
lgf.melt$species = as.factor(lgf.melt$species)

# this plot shows the lowest growth factor name and value for each species each year of the enesmblse 
lgf.melt %>% ggplot(aes(x=year, y=lgf.val, col=lgf)) +
  geom_point() +
  geom_jitter() +
  facet_wrap(~species) + 
  scale_color_manual(
    values = c('gold','blue','red','green'),
    limits = c(1,2,3,4),
    name = 'growth factor',
    labels = c("algf","smgf","sngf","degdgf")
  )

# is the species being killed too early/often?

# check:
death.array = array(0,dim=c(nspec,nyear,nens))
for (i in 1:nens){ # loop through ensembles
  load(paste0('./out/',runs[i],'/linkages.out.Rdata'))
  for (j in 1:nyear){
    death.array[,j,i] = (ntrees.birth[,j,1] - ntrees.kill[,j,1])/ntrees.birth[,j,1]
  }
}
death.melt = reshape::melt(death.array)
colnames(death.melt) = c('species','year','ensemble','deaths')
death.melt$species = as.factor(death.melt$species)

death.melt %>%
  ggplot(aes(x=year,y=deaths,col=species, group=species)) +
  geom_line() +
  facet_wrap(~ensemble) +
  scale_color_manual(
    values = sppcol,
    limits = c(1:nspec),
    name = 'species',
    labels = sppname
  )

# diagnose: 
 
# track number of percentage of nogro trees for each species and compare with death rates
nogro.death.melt = left_join(death.melt,nogro.melt, by = c('species','year','ensemble'))
nogro.death.melt %>% 
  ggplot(aes(x=year,y=deaths,size=count,col=ensemble,group=ensemble)) + 
  geom_point() +
  facet_wrap(~species)
