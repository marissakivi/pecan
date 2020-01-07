## Mass updates on BETY species parameter priors 
## September 21 2019
## Marissa Kivi 

# this script takes as input a settings file (for pft information) and an updated prior matrix, where
# the matrix features are pft name, variable name, and distribution name (ALL GIVEN IN PECAN FORMAT), as 
# well as the parameters for the distribution - this function only works with two parameter distns right
# now but could be adapted to include more 

# script workflow: 
# - load settings and open DB connection
# FOR EACH PFT: 
# - read in updated priors for pft 
# FOR EACH UPDATED PRIOR: 
# check if existing prior is good:
#   if yes, go to next prior
#   else, check if prior exists, but isn't associated:
#       if yes, check if another prior is associated:
#             if yes, dissociate it
#             associate good prior with pft
#       else, input a new prior and associate it 

###############################
# set up working environment
###############################
rm(list=ls())
library(dplyr)
library(PEcAn.all)

all.priors <- read.csv('~/VM_scripts/reset_priors.csv', header = T)
pfts <- as.vector(unique(all.priors$pft))
bety <- list(user = 'bety',
             password = 'bety',
             host = 'postgres',
             dbname = 'bety',
             driver = "PostgreSQL",
             write = TRUE)
modeltype =  "1000000008" # linkages model 

# make bety connection
dbcon <- db.open(bety)
on.exit(db.close(dbcon))

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
  
  print(pft)
  
  # get appropriate pft 
  pftres <- query_pfts(dbcon, pft)
  pfttype <- pftres[["pft_type"]]
  pftid <- toString(pftres[["id"]])
  
  if (nrow(pftres) == 0) {
    PEcAn.logger::logger.severe("Could not find pft", pft[["name"]])
    break
  }

  # read in updated priors
  new.priors = all.priors %>% 
    filter(pft == pftres$name) %>%
    dplyr::select(pft, varname, distn, parama, paramb)

  # adjust all priors in table if they are different
  for (i in 1:length(new.priors$pft)){
    
    print(i)
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
    add = TRUE
    
    # check to see if existing priors is the same as what we need
    if (length(existing.prior$name) > 0){
      
      print('existing.prior is good')
      if (existing.prior$parama == parama & existing.prior$paramb == paramb | existing.prior$distn == distn){
        add = FALSE
      }
    }
    
    # if add is still TRUE, we need a new prior 
    if (add){
      
      print('need to add new prior')
      
      # let's see if there are any priors in the database that we can just associate with the pft
      query.priors <- paste0(
        "SELECT variables.name, priors.id, parama, paramb, distn ", 
        "FROM priors ",
        "JOIN variables ON priors.variable_id = variables.id ",
        "WHERE variables.name = '", varname,
        "' AND parama = ", parama,
        " AND paramb = ", paramb,
        " AND distn = '", distn,"'") 
      poss_prior <- db.query(query = query.priors, con = dbcon)
      
      # a match! we need to associate the prior with the pft
      if (length(poss_prior) > 0){
        
        print('there is a prior available')
        
        # remove all associations with pft for this variable if it has some
        if (length(existing.prior) > 0){
          # update database query
          db.query(paste0("DELETE from pfts_priors WHERE prior_id = ",existing.prior$prior_id), dbcon)
        }
          
        print('associating prior')
        # associate other prior with the pft 
        now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        db.query(paste0("INSERT into pfts_priors ",
                        "(pft_id, prior_id, created_at, updated_at)",
                        " VALUES (", pftid, ", ",  poss_prior$id, ", '", now, "', '", now, "')"), dbcon)
        
      }else{ # we need to input a new prior 
        
        print('no prior available - need to insert new')
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
}



