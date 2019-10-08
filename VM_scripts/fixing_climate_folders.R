
outdir = '/data/dbfiles/met_data/HARVARD/linkages'
dirs = list.dirs('/data/dbfiles/met_data/HARVARD/linkages', full.names=FALSE)[-1]
for (dir in dirs){
  load(paste0(outdir,'/',dir,'/climate.Rdata'))
  filename = paste0(outdir,'/',dir,'.Rdata')
  save(temp.mat,precip.mat, file=filename)
}
