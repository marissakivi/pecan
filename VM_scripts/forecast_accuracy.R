
## Comparing annual forecast accuracy for a single site analysis 
## Here we do not include the first year of assimilation because they're the same in terms of forecast 
## NOTE: all runs should be in the /save directory for this script

rm(list=ls())

sdaID = 'PEcAn_14000000233'
openID = 'PEcAn_14000000236'
nens = 150
nspp = 4
spp = c("red maple", 'red spruce', 'white pine', 'red oak')
years = c(1961:2010)

# organize observed data
load('/data/dbfiles/sda.obs.ROOSTER.Rdata')
obsAGB = matrix(0, nspp, length(years))
obsCOV = matrix(0, nspp, length(years))
for (i in seq_along(years)){
  date = paste0(years[i], '-12-31')
  ind = which(names(obs.list$obs.mean) == date)
  obsAGB[,i] = as.vector(obs.list$obs.mean[[ind]])
  obsCOV[,i] = as.vector(diag(obs.list$obs.cov[[ind]]))
}

# gather agb predictions from free run 
freeAGB = array(0, dim = c(nspp, nens, length(years)))
runs = list.dirs(paste0('/save/workflows/',openID,'/out/'))[-1]
for (i in seq_along(runs)){
  r = runs[i]
  load(paste0(r,'/linkages.out.Rdata'))
  freeAGB[,i,] = agb.pft[,((dim(agb.pft)[2]-length(years)+1):dim(agb.pft)[2]),1]
}

# gather agb predictions from SDA run 
sdaAGB = array(0, dim = c(nspp, nens, length(years)))
runs = list.dirs(paste0('/save/workflows/',sdaID,'/out/'))[(nens+2):(nens*2+1)]
for (i in seq_along(runs)){
  r = runs[i]
  for (j in seq_along(years)){
    if (j != length(years)){
      yr = years[j]
      load(paste0(r,'/',yr, '-12-31 23:59:59linkages.out.Rdata'))
    }else{
      load(paste0(r,'/linkages.out.Rdata'))
    }
    sdaAGB[,i,j] = agb.pft
  }
}

# melt all matrices/array into dataframes to combine 
free = reshape2::melt(freeAGB)
colnames(free) = c('species','ens','year', 'agb')
free$type = 'free'

sda = reshape2::melt(sdaAGB)
colnames(sda) = c('species','ens','year', 'agb')
sda$type = 'sda'

obs = reshape2::melt(obsAGB) 
colnames(obs) = c('species', 'year', 'obs')

cov = reshape2::melt(obsCOV)
colnames(cov) = c('species', 'year', 'obsVar')

obs = left_join(obs, cov, by = c('species', 'year'))


# join the dataframes into one 
data = rbind(free, sda) 
allData = left_join(data, obs, by = c('species','year'))
allData$year = allData$year + years[1] - 1 

# summarize info on forecasted and observed total AGB
totalObs = obs %>% 
  group_by(year) %>%
  summarize(obs = sum(obs)) %>% 
  ungroup()
totalData = data %>% 
  group_by(year, type, ens) %>% 
  summarize(total = sum(agb)) %>% 
  ungroup() %>%
  left_join(totalObs, by = c('year'))

# calculate PFT-level RMSE and normalized RMSE 
comparison = allData %>% 
  mutate(SE = (agb - obs)^2) %>%
  group_by(species, type, ens) %>% 
  summarize(RMSE = sqrt(mean(SE, na.rm=TRUE)), 
            nRMSE = (sqrt(mean(SE, na.rm=TRUE)) / mean(obs, na.rm=TRUE))*100) %>%
  ungroup() %>% 
  group_by(species, type) %>% 
  summarize(RMSE = mean(RMSE),
            nRMSE = mean(nRMSE))

# calculate total RMSE and normalized RMSE 
allComparison = totalData %>% 
  mutate(SE = (total-obs)^2) %>% 
  group_by(type, ens) %>% 
  summarize(RMSE = sqrt(mean(SE, na.rm = TRUE)), 
            nRMSE = (sqrt(mean(SE, na.rm = TRUE)) / mean(obs, na.rm = TRUE))*100) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  summarize(RMSE = mean(RMSE), 
            nRMSE = mean(nRMSE))


plotData = allData %>% 
  group_by(species, year, type) %>% 
  summarize(obs.lower = mean(obs) - 1.96 * sqrt(mean(obsVar)), 
            obs.upper = mean(obs) + 1.96 * sqrt(mean(obsVar)),
            obs = mean(obs), 
            mean = mean(agb, na.rm = TRUE), 
            lower = quantile(agb, 0.025, na.rm = TRUE), 
            upper = quantile(agb, 0.975, na.rm = TRUE)) %>% 
  ungroup()
plotData$species = plyr::mapvalues(plotData$species, 
                                   from = c(1:nspp), 
                                   to = spp)

ggplot(plotData) + 
  geom_line(aes(x = year, y = mean, group = type, color = type)) + 
  geom_line(aes(x = year, y = obs)) + 
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, group = type, fill = type), alpha = 0.25) + 
  facet_wrap(species~., scales = 'free', ncol = 1)

totalPlot = totalData %>% group_by(year, type) %>%
  summarize(mean = mean(total), 
            lower = quantile(total, 0.025), 
            upper = quantile(total, 0.975), 
            obs = mean(obs)) %>%
  ungroup()

ggplot(totalPlot) + 
  geom_line(aes(x = year, y = mean, group = type, color = type)) + 
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, group = type, fill = type), alpha = 0.25) + 
  geom_line(aes(x = year, y = obs)) + 
  labs(title = 'Total Aboveground Biomass', 
       y = 'aboveground biomass (Kg C / m2)')
  


freePlot = ggplot(plotData %>% filter(type == 'free')) +
  geom_line(aes(x = year, y = obs.lower), linetype = 'dashed') +
  geom_line(aes(x = year, y = obs.upper), linetype = 'dashed') +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper), fill = 'blue', alpha = 0.25) +
  geom_line(aes(x = year, y = obs, col = 'data'),size = 1.2) +
  geom_line(aes(x = year, y = mean, col = 'prediction mean')) +
  facet_grid(species~., scales='free') +
  labs(col = 'species', y = 'aboveground biomass (Kg C / m2)',
       title = 'free run results: predicted biomass vs. data') +
  scale_color_manual(values = c('data'='black','prediction mean'='blue'),name =NULL)
freePlot
