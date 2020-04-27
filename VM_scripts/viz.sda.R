# SDA Visualization Scripts
# Marissa Kivi 
# October 2019

## 1. set working directory to your SDA run 
setwd('/data/workflows/PEcAn_14000000089')

## 2. write out species names for your run in the same order as you find them in the PFT folder of your working directory
## I didn't make this hard-coded because I wanted to keep things neat for nice plotting

sppname = c('red maple','yellow birch','american beech','red oak','eastern hemlock') # HF
#sppname = c('red maple','red spruce','white pine','red oak') # RH
#sppname = c('red maple','american beech','white pine','white oak','red oak') # GE 
#sppname = c('sugar maple', 'yellow birch', 'american beech', 'white ash', 'eastern hemlock') # NRP12

## 3. load your obs.list file for your run as well 

load('/data/dbfiles/sda.obs.HARVARD.RData')
#load('/data/dbfiles/sda.obs.ROOSTER.Rdata')
#load('/data/dbfiles/sda.obs.GOOSE.RData')
#load('/data/dbfiles/sda.obs.NORTHROUND.plots12.Rdata')

######################
## set up environment
######################

# read in settings
settings <- read.settings("pecan.SDA.xml")

# create plots folder
if (!dir.exists('./plots')) dir.create('./plots')

# read in data 
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
years = (lubridate::year(settings$state.data.assimilation$start.date)+1):lubridate::year(settings$state.data.assimilation$end.date)

adj.agb.total = matrix(0,n,length(years))
adj.agb.pft.array = array(0,dim = c(n, length(years), spp))
pred.agb.total = matrix(0,n,length(years))
pred.agb.pft.array = array(0,dim = c(n, length(years), spp))
bias.total = matrix(0,n,length(years))
# bias.pft.array = array(0,dim=c(n,length(years),spp))

# open sda.output.Rdata file
load('./SDA/sda.output.Rdata')
for (i in 1:length(years)){
  pred.agb.total[,i] = apply(FORECAST[[i+1]], 1, sum)
  pred.agb.pft.array[,i,] = FORECAST[[i+1]]
  adj.agb.total[,i] = apply(ANALYSIS[[i+1]],1,sum)
  adj.agb.pft.array[,i,] = as.matrix(ANALYSIS[[i+1]])
  bias.total[,i] = adj.agb.total[,i] - pred.agb.total[,i]
  # bias.pft.array[,i,] = adj.agb.pft.array[,i,] - pred.agb.pft.array[,i,]
}

######################
## melt to dataframes
######################

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

# bias.pft.melt = reshape::melt(bias.pft.array)
# colnames(bias.pft.melt) = c('ensemble','year','species','bias')
# bias.pft.melt$year = bias.pft.melt$year + years[1] - 1

# format observed data
inds = which(lubridate::year(names(obs.mean)) %in% years)
agb.pft.obs = matrix(0,length(inds),spp)
agb.pft.var = matrix(0,length(inds),spp)
agb.obs.tot.var = vector()
track = 1
for (yr in inds){
  agb.pft.obs[track,] = obs.mean[[yr]]
  agb.pft.var[track,] = diag(obs.cov[[yr]])
  agb.obs.tot.var[track] = sum(colSums(obs.cov[[yr]]))
  track = track + 1
}
agb.pft.obs.melt = reshape2::melt(agb.pft.obs)
names(agb.pft.obs.melt) = c('year','species','agb.obs')
agb.pft.obs.melt$year = agb.pft.obs.melt$year + years[1] - 1
agb.pft.obs.melt$species = plyr::mapvalues(agb.pft.obs.melt$species, from = c(1:spp), to = sppname)
agb.pft.obs.melt$species = as.factor(agb.pft.obs.melt$species)

agb.pft.var.melt = melt(agb.pft.var)
names(agb.pft.var.melt) = c('year','species','agb.var')
agb.pft.var.melt$species = plyr::mapvalues(agb.pft.var.melt$species, from = c(1:spp), to = sppname)
agb.pft.var.melt$year = agb.pft.var.melt$year + years[1] - 1
agb.pft.obs.melt = left_join(agb.pft.obs.melt, agb.pft.var.melt, by = c('year','species')) %>% 
  mutate(b05.data = qnorm(0.05, mean = agb.obs, sd = sqrt(agb.var)),
         b95.data = qnorm(0.95, mean = agb.obs,sd = sqrt(agb.var)), 
         b25.data = qnorm(0.25, mean = agb.obs,sd = sqrt(agb.var)),
         b75.data = qnorm(0.75, mean = agb.obs,sd = sqrt(agb.var)))

agb.obs.var = as.data.frame(agb.obs.tot.var)
colnames(agb.obs.var) = c('total.var')
agb.obs.var$year = c(1:length(inds)) + years[1] - 1

agb.obs.total = agb.pft.obs.melt %>% group_by(year) %>% summarize(total.obs = sum(agb.obs)) %>% 
  left_join(agb.obs.var, by = 'year') %>% mutate(b05.data = qnorm(0.05, mean = total.obs, sd = sqrt(total.var)),
                                                 b95.data = qnorm(0.95, mean = total.obs,sd = sqrt(total.var)), 
                                                 b25.data = qnorm(0.25, mean = total.obs,sd = sqrt(total.var)),
                                                 b75.data = qnorm(0.75, mean = total.obs,sd = sqrt(total.var))
                                                 )
 
##############################
## organize data for plotting
##############################

# plot-ready pft agb predicted
pred.pft.plot = pred.agb.pft.melt %>%
  group_by(year, species) %>%
  summarize(b05 = quantile(agb,0.05),
            b95 = quantile(agb,0.95),
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75),
            mean.pred = mean(agb))
pred.pft.plot$species = as.factor(pred.pft.plot$species)
pred.pft.plot$species = plyr::mapvalues(pred.pft.plot$species, from = c(1:spp), to = sppname)

# plot-ready pft agb adjusted
adj.pft.plot = adj.agb.pft.melt %>%
  group_by(year, species) %>%
  summarize(b05 = quantile(agb,0.05),
            b95 = quantile(agb,0.95),
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75),
            mean.adj = mean(agb))
adj.pft.plot$species = as.factor(adj.pft.plot$species)
adj.pft.plot$species = plyr::mapvalues(adj.pft.plot$species, from = c(1:spp), to = sppname)

# plot-ready total agb 
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

# plot-ready total model bias 
bias.plot = bias.total.melt %>% group_by(year) %>%
  summarize(mean.bias = mean(bias),
            b05 = quantile(bias,0.05),
            b95 = quantile(bias,0.95),
            b25 = quantile(bias,0.25),
            b75 = quantile(bias,0.75))

#########
## plot
#########

# first plot looks at the range of yearly predictions by LINKAGES prior to adjustment
# black line is the observed data
# red line is the adjusted mean
just.mean = adj.pft.plot %>% select(year,species,mean.adj)
pl1 = left_join(pred.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  left_join(just.mean, by = c('year','species')) %>%
  ggplot() +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90% predicted'), alpha = 0.6) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50% predicted'), alpha = 0.6) +
  geom_line(aes(x = year, y = agb.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.pred, col = 'prediction mean')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  facet_grid(species~., scales='free') +
  labs(col = 'species', y = 'aboveground biomass (Kg C / m2)',
       title = 'sda results: predicted biomass vs. data') +
  scale_fill_manual(values = c('90% predicted'='coral2','50% predicted'='rosybrown1'), name = 'intervals') +
  scale_color_manual(values = c('data'='black','prediction mean'='blue', 'adjusted mean'='red'),name =NULL)
pl1
ggsave(plot = pl1, filename = './plots/predicted-by-pft.jpg')

pl1b = left_join(pred.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  left_join(just.mean, by = c('year','species')) %>%
  ggplot() +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90% predicted'), alpha = 0.6) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50% predicted'), alpha = 0.6) +
  geom_line(aes(x = year, y = agb.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.pred, col = 'prediction mean')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  facet_grid(species~., scales='fixed') +
  labs(col = 'species', y = 'aboveground biomass (Kg C / m2)',
       title = 'sda results: predicted biomass vs. data') +
  scale_fill_manual(values = c('90% predicted'='coral2','50% predicted'='rosybrown1'), name = 'intervals') +
  scale_color_manual(values = c('data'='black','prediction mean'='blue', 'adjusted mean'='red'),name =NULL)
pl1b
ggsave(plot = pl1b, filename = './plots/predicted-by-pft.fixed.jpg')

# second plot looks at the range of yearly adjusted values for ensembles 
# the black line is the observed data 
pl2 = left_join(adj.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  ggplot() +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90% adjusted'), alpha = 0.6) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50% adjusted'), alpha = 0.6) +
  geom_line(aes(x = year, y = agb.obs, col = 'data')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  facet_grid(species~., scales='free') +
  labs(col = 'species', title = 'sda results: adjusted biomass vs. data',
       y = 'aboveground biomass (Kg C/m2)') +
  scale_fill_manual(values = c('90% adjusted'='coral2','50% adjusted'='rosybrown1'), name = 'intervals') +
  scale_color_manual(values = c('data'='black','adjusted mean'='red'),name =NULL)
pl2
ggsave(plot = pl2, filename = './plots/adjusted-by-pft.jpg')

pl2b = left_join(adj.pft.plot, agb.pft.obs.melt, by = c('year','species'))  %>%
  ggplot() +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90% adjusted'), alpha = 0.6) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50% adjusted'), alpha = 0.6) +
  geom_line(aes(x = year, y = agb.obs, col = 'data')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  facet_grid(species~., scales='fixed') +
  labs(col = 'species', title = 'sda results: adjusted biomass vs. data',
       y = 'aboveground biomass (Kg C/m2)') +
  scale_fill_manual(values = c('90% adjusted'='coral2','50% adjusted'='rosybrown1'), name = 'intervals') +
  scale_color_manual(values = c('data'='black','adjusted mean'='red'),name =NULL)
pl2b
ggsave(plot = pl2, filename = './plots/adjusted-by-pft.fixed.jpg')

# third plot looks at the range of total predicted biomass compared to the data and the adjusted values
# data is shown by the black line
pl3 = left_join(pred.intervals, agb.obs.total, by = c('year'))  %>%
  left_join(just.mean.tot, by = c('year')) %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90% predicted'), alpha = 0.6) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50% predicted'), alpha = 0.6) +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = total.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.pred, col = 'prediction mean')) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  labs(y = 'aboveground biomass (Kg C / m2)', title = 'sda results: predicted total biomass vs. data') +
  scale_fill_manual(values = c('90% predicted'='coral2','50% predicted'='rosybrown1'),name = 'intervals') +
  scale_color_manual(values = c('data'='black','prediction mean'='blue', 'adjusted mean'='red'),name =NULL)
pl3
ggsave(plot = pl3, filename = './plots/predicted-total.jpg')

# the fourth plot looks at the range of total adjusted biomass compared to the data
# data is shown by the black line
pl4 = left_join(adj.intervals, agb.obs.total, by = c('year'))  %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(x = year, ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(x = year, y = b05.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = b95.data), linetype = 'dashed') +
  geom_line(aes(x = year, y = total.obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean.adj, col = 'adjusted mean')) +
  labs(y = 'aboveground biomass (Kg C/m2)', title = 'sda results: adjusted total biomass vs. data') +
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1'),name ='adjusted interval') +
  scale_color_manual(values = c('data'='black','adjusted mean'='red'),name =NULL)
pl4
ggsave(plot = pl4, filename = './plots/adjusted-total.jpg')


pl5 = bias.plot %>% ggplot(aes(x = year, y = mean.bias)) +
  geom_ribbon(aes(ymin = b05, ymax = b95, fill = '90%')) +
  geom_ribbon(aes(ymin = b25, ymax = b75, fill = '50%')) +
  geom_line(aes(col = 'mean')) +
  scale_color_manual(values = c('mean'='black'),name =NULL) +
  scale_fill_manual(values = c('90%'='coral2','50%'='rosybrown1', name = 'predicted interval')) +
  labs(title = 'bias in total agb', y = 'bias', x = 'year') +
  geom_hline(yintercept = 0)
ggsave(plot = pl5, filename = './plots/bias-total.jpg')

