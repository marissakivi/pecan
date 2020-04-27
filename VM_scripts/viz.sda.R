# SDA Visualization Scripts
# Marissa Kivi 
# October 2019

sppname = c('sugar maple', 'yellow birch', 'american beech', 'white ash', 'eastern hemlock')

# read in settings
settings <- read.settings("pecan.SDA.xml")

# read in data 
load('/data/dbfiles/sda.obs.NORTHROUND.plots12.Rdata')
obs.mean <- obs.list$obs.mean
obs.cov <- obs.list$obs.cov

library(dplyr)
library(ggplot2)

#################
## collect data
#################

# read in and organize data 
n = as.numeric(settings$ensemble$size)
spp = length(settings$pfts)
years = lubridate::year(settings$state.data.assimilation$start.date):lubridate::year(settings$state.data.assimilation$end.date)

adj.agb.total = matrix(0,n,length(years))
adj.agb.pft.array = array(0,dim = c(n, length(years), spp))
pred.agb.total = matrix(0,n,length(years))
pred.agb.pft.array = array(0,dim = c(n, length(years), spp))
bias.total = matrix(0,n,length(years))
bias.pft.array = array(0,dim=c(n,length(years),spp))

# open sda.output.Rdata file
load('./SDA/sda.output.Rdata')
for (i in 1:length(years)){
  pred.agb.total[,i] = apply(FORECAST[[i]], 1, sum)
  pred.agb.pft.array[,i,] = FORECAST[[i]]
  adj.agb.total[,i] = apply(ANALYSIS[[i]],1,sum)
  adj.agb.pft.array[,i,] = as.matrix(ANALYSIS[[i]])
  bias.total[,i] = adj.agb.total[,i] - pred.agb.total[,i]
  bias.pft.array[,i,] = adj.agb.pft.array[,i,] - pred.agb.pft.array[,i,]
}

#################
## format data
#################

# convert to data frames so we can use ggplot
# melt, name columns, adjust, and add type column
# make adjusted values at a slightly later time point so not overlapping entirely if same
adj.agb.total.melt = reshape::melt(adj.agb.total)
colnames(adj.agb.total.melt) = c('ensemble','year','agb')
adj.agb.total.melt$year = adj.agb.total.melt$year + years[1] - 1
adj.agb.total.melt$type = rep('adjusted',length(adj.agb.total.melt$year))

adj.agb.pft.melt = reshape::melt(adj.agb.pft.array)
colnames(adj.agb.pft.melt) = c('ensemble','year','species','agb')
adj.agb.pft.melt$year = adj.agb.pft.melt$year + years[1] - 1
adj.agb.pft.melt$type = rep('adjusted',length(adj.agb.pft.melt$year))

pred.agb.total.melt = reshape::melt(pred.agb.total)
colnames(pred.agb.total.melt) = c('ensemble','year','agb')
pred.agb.total.melt$year = pred.agb.total.melt$year + years[1] - 1
pred.agb.total.melt$type = rep('predicted',length(pred.agb.total.melt$year))

pred.agb.pft.melt = reshape::melt(pred.agb.pft.array)
colnames(pred.agb.pft.melt) = c('ensemble','year','species','agb')
pred.agb.pft.melt$year = pred.agb.pft.melt$year + years[1] - 1
pred.agb.pft.melt$type = rep('predicted',length(pred.agb.pft.melt$year))

bias.total.melt = reshape::melt(bias.total)
colnames(bias.total.melt) = c('ensemble','year','bias')
bias.total.melt$year = bias.total.melt$year + years[1] - 1
bias.pft.melt = reshape::melt(bias.pft.array)
colnames(bias.pft.melt) = c('ensemble','year','species','bias')
bias.pft.melt$year = bias.pft.melt$year + years[1] - 1

# combine into two dataframes
agb.total.melt = rbind(adj.agb.total.melt, pred.agb.total.melt)
agb.pft.melt = rbind(adj.agb.pft.melt, pred.agb.pft.melt)

# rename levels of species variable so easier to identify species
agb.pft.melt$species = as.factor(agb.pft.melt$species)
agb.pft.melt$species = plyr::mapvalues(agb.pft.melt$species, from = c(1:spp), to = sppname)

# format observed data
inds = which(lubridate::year(names(obs.mean)) %in% years)
agb.pft.obs = matrix(0,length(inds),spp)
track = 1
for (yr in inds){
  agb.pft.obs[track,] = obs.mean[[yr]]
  track = track + 1
}
agb.pft.obs.melt = reshape2::melt(agb.pft.obs)
names(agb.pft.obs.melt) = c('year','species','agb.obs')
agb.pft.obs.melt$year = agb.pft.obs.melt$year + years[1] - 1
agb.pft.obs.melt$species = plyr::mapvalues(agb.pft.obs.melt$species, from = c(1:spp), to = sppname)
agb.pft.obs.melt$species = as.factor(agb.pft.obs.melt$species)
agb.obs.total = agb.pft.obs.melt %>% group_by(year) %>% summarize(total.obs = sum(agb.obs))

agb.pft.melt = left_join(agb.pft.melt,agb.pft.obs.melt, by = c('species','year'))

#################
## plot data
#################

# plot pft agb
pred.pft.plot = pred.agb.pft.melt %>% 
  group_by(year, species) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75), 
            mean.pred = mean(agb))
adj.pft.plot = adj.agb.pft.melt %>% 
  group_by(year, species) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75), 
            mean.adj = mean(agb))
adj.pft.plot$species = as.factor(adj.pft.plot$species)
adj.pft.plot$species = plyr::mapvalues(adj.pft.plot$species, from = c(1:spp), to = sppname)
pred.pft.plot$species = as.factor(pred.pft.plot$species)
pred.pft.plot$species = plyr::mapvalues(pred.pft.plot$species, from = c(1:spp), to = sppname)

# first plot looks at the range of yearly predictions by LINKAGES prior to adjustment
# the black line is the observed data 
# the red line is the adjusted mean
just.mean = adj.pft.plot %>% select(year,species,mean.adj)
pl1 = left_join(pred.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  left_join(just.mean, by = c('year','species')) %>%
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(x = year, y = agb.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.pred, col = 'prediction mean')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) + 
  facet_grid(species~., scales='free') + 
  labs(col = 'species', y = 'aboveground biomass (Kg C / m2)', 
       title = 'sda results: predicted biomass vs. data') + 
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1'),name ='predicted interval') + 
  scale_color_manual(values = c('data'='black','prediction mean'='blue', 'adjusted mean'='red'),name =NULL)
ggsave(plot = pl1, filename = 'predicted-by-pft.jpg')

# second plot looks at the range of yearly adjusted values for ensembles 
# the black line is the observed data 
pl2 = left_join(adj.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(x = year, y = agb.obs, col = 'data')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  facet_grid(species~., scales='free') + 
  labs(col = 'species', title = 'sda results: adjusted biomass vs. data', 
       y = 'aboveground biomass (Kg C/m2)') + 
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1'),name ='adjusted interval') + 
  scale_color_manual(values = c('data'='black','adjusted mean'='red'),name =NULL)
ggsave(plot = pl2, filename = 'adjusted-by-pft.jpg')

# plot total agb 
pred.intervals = pred.agb.total.melt %>% 
  group_by(year) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75), 
            mean.pred = mean(agb))
adj.intervals = adj.agb.total.melt %>% 
  group_by(year) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75), 
            mean.adj = mean(agb))
just.mean.tot = adj.intervals %>% select(year, mean.adj)

# the third plot looks at the range of total predicted biomass compared to the data and the adjusted values
# data is shown by the black line
pl3 = left_join(pred.intervals, agb.obs.total, by = c('year'))  %>%
  left_join(just.mean.tot, by = c('year')) %>%
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(x = year, y = total.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.pred, col = 'prediction mean')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) + 
  labs(y = 'aboveground biomass (Kg C / m2)', title = 'sda results: predicted total biomass vs. data') + 
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1'),name ='predicted interval') + 
  scale_color_manual(values = c('data'='black','prediction mean'='blue', 'adjusted mean'='red'),name =NULL)
ggsave(plot = pl3, filename = 'predicted-total.jpg')

  
# the fourth plot looks at the range of total adjusted biomass compared to the data
# data is shown by the black line
pl4 = left_join(adj.intervals, agb.obs.total, by = c('year'))  %>%
  ggplot() + 
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(x = year, y = total.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) + 
  labs(y = 'aboveground biomass (Kg C/m2)', title = 'sda results: adjusted total biomass vs. data') + 
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1'),name ='adjusted interval') + 
  scale_color_manual(values = c('data'='black','adjusted mean'='red'),name =NULL)
ggsave(plot = pl4, filename = 'adjusted-total.jpg')

# plot model bias in predicting total agb
bias.plot = bias.total.melt %>% group_by(year) %>%
  summarize(mean.bias = mean(bias), 
            b05 = quantile(bias,0.05), 
            b95 = quantile(bias,0.95), 
            b25 = quantile(bias,0.25),
            b75 = quantile(bias,0.75))
pl5 = bias.plot %>% ggplot(aes(x = year, y = mean.bias)) + 
  geom_ribbon(aes(ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(col = 'mean')) +
  scale_color_manual(values = c('mean'='black'),name =NULL) + 
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1', name = 'predicted interval')) + 
  labs(title = 'bias in total agb', y = 'bias', x = 'year') + 
  geom_hline(yintercept = 0)
ggsave(plot = pl5, filename = 'bias-total.jpg')

