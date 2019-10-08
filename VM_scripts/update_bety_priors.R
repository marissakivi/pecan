## Mass updates on BETY species parameter priors 
## September 21 2019
## Marissa Kivi 

# inputs: settings file, updated prior matrix (use as reference the file created in workflow)

# might just be easier to overwrite all of them - assume they have been changed

# check to make sure all species are entered on database (should be in order to have on interface)
# determine which variables have been associated with species already 
# input updated priors for species variables 
# how to associate these with species? 

###############################
# set up working environment
###############################

library(dplyr)
library(PEcAn.all)

# load settings
settings <- read.settings("pecan.xml")
pfts <- settings$pfts
bety <- settings$database$bety
modeltype = settings$model$id

# make bety connection
dbcon <- db.open(bety)
on.exit(db.close(dbcon))
# OR 
#drv <- RPostgreSQL::PostgreSQL()
#con = DBI::dbConnect(drv,
#                     host = 'postgres',
#                     dbname = 'bety',
#                     user = 'bety',
#                     password = 'bety')

# variable names that are needed for LINKAGES species: 
pars = c('DMAX','DMIN','HTMAX','AGEMX','DBHMAX','Gmax','SPRTND','SPRTMN','SPRTMX','MPLANT','D3',
         'FROST','CM1','CM2','CM3','CM4','CM5','FWT','SLTA','SLTB')

# get variable ids for each desired parameter
query.text <- paste(
  "SELECT name, variables.id",
  "FROM variables") 
vars.id = db.query(query = query.text, con = dbcon) %>%
  filter(name %in% pars) 
vars.id$id = sapply(vars.id$id, toString)

#################################
# check what is available in db
#################################

# check pfts in database
if (!is.list(pfts)) {
  PEcAn.logger::logger.severe('pfts must be a list')
}
pft_names <- vapply(pfts, "[[", character(1), "name")
pft_ids <- sapply(query_pfts(dbcon, pft_names, strict = TRUE)[["id"]], toString)

for (pft in pft_names){
  
  # get appropriate pft 
  pftres <- query_pfts(dbcon, pft)
  pfttype <- pftres[["pft_type"]]
  pftid <- toString(pftres[["id"]])
  
  if (nrow(pftres) == 0) {
    PEcAn.logger::logger.severe("Could not find pft", pft[["name"]])
    break
  }
  
  # determine priors in database for pft
  query.text <- paste(
    "SELECT variables.name, distn, parama, paramb, pfts_priors.pft_id",
    "FROM priors",
    "JOIN variables ON priors.variable_id = variables.id",
    "JOIN pfts_priors ON pfts_priors.prior_id = priors.id",
    "JOIN pfts ON pfts.id = pfts_priors.pft_id") 
  priors <- db.query(query = query.text, con = dbcon) %>% 
    filter(name %in% pars,
           pft_id == pftid)
  
  # determine priors missing from database
  missing = pars[which(!(pars %in% priors$name))]
  
  # read in updated priors
  new.priors = read.csv('~/VM_scripts/updated_priors.csv', header = T) %>% 
    filter(pft == pftres$name) %>%
    dplyr::select(pft, varname, distn, parama, paramb)

  # adjust all necessary priors
  for (i in 1:length(new.priors$pft)){
    p = new.priors[i,]
    
    varname = p$varname # variable name not pft
    distn = p$distn
    parama = p$parama
    paramb = p$paramb
    varid = vars.id %>% filter(name == varname) %>% dplyr::select(id)
    pftid = pftid
    
    query.priors <- paste(
      "SELECT variables.name, pfts_priors.prior_id, parama, paramb, distn", 
      "FROM priors",
      "JOIN variables ON priors.variable_id = variables.id",
      "JOIN pfts_priors ON pfts_priors.prior_id = priors.id",
      "JOIN pfts ON pfts.id = pfts_priors.pft_id",
      "WHERE pfts.id =", pftid,
      'AND priors.variable_id =', varid) 
    priorid <- db.query(query = query.priors, con = dbcon)
    existing.prior <- priorid
    
    # prior in database already?
    
    now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    # we need to update the existing prior
    if (length(existing.prior$name) > 0){
      
      existing.prior$prior_id = toString(existing.prior$prior_id)
      
      # update database query
      db.query(paste0("UPDATE priors SET distn = '", distn,
                      "', parama = ", parama, ", paramb = ", paramb,
                      ", updated_at = '",now, "' WHERE id = ",existing.prior$prior_id), dbcon)
     
    }else{ # else we need to input a new prior for the species
      
      # input database query (first into priors table)
      db.query(paste0("INSERT INTO priors ",
                      "(variable_id, phylogeny, distn, parama, paramb, notes, created_at, updated_at)",
                      " VALUES (",varid, ", '', '", distn, "', ", parama, ", ", paramb, ", '", pftres$name, 
                      " linkages', '", now, "', '", now, "')"), dbcon)
      new.prior.id = db.query(paste0("SELECT id FROM priors where variable_id = ",varid, "AND notes = '", 
                                     pftres$name, " linkages'"), dbcon)
      
      # (then into pfts_priors)
      db.query(paste0("INSERT INTO pfts_priors ",
                      "(pft_id, prior_id, created_at, updated_at)", 
                      "VALUES (", pftid, ", ", new.prior.id, ", '", now, "', '", now, "')"), dbcon)
    }
  }
}




