# Looking at adjusted stand structure for three years in an SDA run 

rm(list=ls())
settings <- read.settings("pecan.SDA.xml")
run.id = as.list(c(14000006687:14000006701))
time_to_rewind = 1961
sda_rewind(settings,run.id,time_to_rewind)

library(igraph)
library(ggraph)

runID = 14000000117

years = c(1960:1968)
spp = c('Red Maple','Yellow Birch','American Beech','Red Oak','Hemlock')
sppcol = c('blue'='Red Maple','gold','green','red','purple')
ens = c(14000006732:14000006746)

for (yr in years){
  
  pl.diam = 3256.7 #cm
  ids = c('origin', sapply(ens,toString)) # initialize origin and ensemble nodes
  parent = c(rep('origin',length(ens)))
  child = c(sapply(ens,toString))
  spp.inds = c(NA, rep('plot', length(ens)))
  dbh.inds = c(NA,rep(pl.diam, length(ens)))
  
  # loop through each ensemble and add data
  for (i in 1:length(ens)){
    
    # load data
    ens.now = ens[i]
    load(file.path('/data','workflows',paste0('PEcAn_',runID),'run',ens.now, paste0(toString(yr+1),'-12-31 23:59:59linkages.restart.Rdata')))
    num = sum(ntrees)

    spp.ids = c()
    for (s in seq_along(spp)){
      spp.ids = c(spp.ids, rep(s, ntrees[s]))
    }
    
    for (tt in 1:num){
      tree = paste0(i,'tree',tt)
      parent = c(parent, toString(ens.now))
      child = c(child,tree)
      sp = spp.ids[tt]
      ids = c(ids,tree)
      spp.inds = c(spp.inds, spp[sp])
      dbh.inds = c(dbh.inds, dbh[tt])
      
    }
  }
  
  # organize data into hierarchical format
  edges <- data.frame(parent = parent, child = child)
  vertices <- data.frame(id = ids, species=spp.inds, size=dbh.inds)
  hData <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
  
  # create plot with hierarchical data
  pl3 <- ggraph(hData, layout = 'circlepack', weight="size") +
    geom_node_circle(aes(fill=species)) +
    scale_fill_manual(
      values = c('Red Maple'='blue','plot'='gray',
                 'Yellow Birch'='gold','American Beech'='green','Red Oak'='red','Hemlock'='purple')) +
    ggtitle(paste0('50th percentile trees :: ',yr)) +
    theme_void()
  ggsave(plot = pl3, paste0('~/data_files/','SLT.',toString(yr),'.circle.plot.jpg'))
}
