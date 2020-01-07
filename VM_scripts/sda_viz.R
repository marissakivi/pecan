# SDA Visualization Scripts
# Marissa Kivi 
# October 2019

# read in settings
settings <- read.settings("pecan.SDA.xml")

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
adj.agb.total.melt$type = rep('adjusted',length(adj.agb.total.melt$year))

adj.agb.pft.melt = reshape::melt(adj.agb.pft.array)
colnames(adj.agb.pft.melt) = c('ensemble','year','species','agb')
adj.agb.pft.melt$type = rep('adjusted',length(adj.agb.pft.melt$year))

pred.agb.total.melt = reshape::melt(pred.agb.total)
colnames(pred.agb.total.melt) = c('ensemble','year','agb')
pred.agb.total.melt$type = rep('predicted',length(pred.agb.total.melt$year))

pred.agb.pft.melt = reshape::melt(pred.agb.pft.array)
colnames(pred.agb.pft.melt) = c('ensemble','year','species','agb')
pred.agb.pft.melt$type = rep('predicted',length(pred.agb.pft.melt$year))

bias.total.melt = reshape::melt(bias.total)
colnames(bias.total.melt) = c('ensemble','year','bias')
bias.pft.melt = reshape::melt(bias.pft.array)
colnames(bias.pft.melt) = c('ensemble','year','species','bias')

# combine into two dataframes
agb.total.melt = rbind(adj.agb.total.melt, pred.agb.total.melt)
agb.pft.melt = rbind(adj.agb.pft.melt, pred.agb.pft.melt)

# rename levels of species variable so easier to identify species
agb.pft.melt$species = as.factor(agb.pft.melt$species)
agb.pft.melt$species = plyr::mapvalues(agb.pft.melt$species, from = c(1,2,3,4), to = c('red maple', 'yellow birch', 'american beech', 'red oak'))

#################
## plot data
#################

# plot pft agb
pred.pft.intervals = pred.agb.pft.melt %>% 
  group_by(year, species) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75))
agb.pft.plot = adj.agb.pft.melt %>% group_by(year, species) %>%
  summarize(mean.agb = mean(agb)) %>%
  ungroup() %>%
  right_join(pred.pft.intervals, c('year','species'))
agb.pft.plot$species = as.factor(agb.pft.plot$species)
agb.pft.plot$species = plyr::mapvalues(agb.pft.plot$species, from = c(1,2,3,4), to = c('red maple', 'yellow birch', 'american beech', 'red oak'))

agb.pft.plot %>%
  ggplot(aes(x = year, y = mean.agb)) + 
  geom_ribbon(aes(ymin = b05, ymax = b95), fill = 'lightcoral') +
  geom_ribbon(aes(ymin = b25, ymax = b75), fill = 'lightcyan', alpha = 0.50) +
  geom_line(col='black') +
  facet_grid(species~., scales='free') + 
  labs(col = 'species', title = 'sda results: agb by pft')

# plot total agb 
pred.intervals = pred.agb.total.melt %>% 
  group_by(year) %>% 
  summarize(b05 = quantile(agb,0.05), 
            b95 = quantile(agb,0.95), 
            b25 = quantile(agb,0.25),
            b75 = quantile(agb,0.75))
agb.plot = adj.agb.total.melt %>% group_by(year) %>%
  summarize(mean.agb = mean(agb)) %>%
  right_join(pred.intervals, 'year')
  
agb.plot %>% ggplot(aes(x = year, y = mean.agb)) + 
  geom_ribbon(aes(ymin = b05, ymax = b95), fill = 'lightcoral') +
  geom_ribbon(aes(ymin = b25, ymax = b75), fill = 'lightcyan', alpha = 0.50) +
  geom_line(col='black') +
  labs(title = 'sda results: total agb')

# plot model bias in predicting total agb
bias.plot = bias.total.melt %>% group_by(year) %>%
  summarize(mean.bias = mean(bias), 
            b05 = quantile(bias,0.05), 
            b95 = quantile(bias,0.95), 
            b25 = quantile(bias,0.25),
            b75 = quantile(bias,0.75))
bias.plot %>% ggplot(aes(x = year, y = mean.bias)) + 
  geom_ribbon(aes(ymin = b05, ymax = b95, fill = '90% Pred')) +
  geom_ribbon(aes(ymin = b25, ymax = b75, fill = '50% Pred'), alpha = 0.85) +
  geom_line(col = 'black') +
  labs(title = 'bias in total agb', ylab = 'bias', xlab = 'year of sda', fill = NULL,
       col = 'Mean Bias') + 
  geom_hline(yintercept = 0)

