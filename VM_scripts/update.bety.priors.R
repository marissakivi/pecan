## Mass updates on BETY species parameter priors for LINKAGES
## September 21 2019
## Marissa Kivi 

# This script takes a specially-formatted prior matrix as input (example given in shared Google drive).
# The matrix features are pft name, variable name, and distribution name (ALL GIVEN IN PECAN FORMAT; see example), as 
# well as the parameters for the distribution. Right now, this function only works with two-parameter distributions. 

# The general workflow of the script is as follows: 
# (1) Set up the working environment and open BETY connection
# (2) Loop through the PFTs
#     (A) Loop through priors
#        (i) Check if the current prior for this species is OK. If yes, go to next prior. 
#        (ii) Else, check if there is another prior in the database that would works. If yes, associate it with the PFT. 
#        (iii) Else, create a new prior for our needs. 

###############################
# 1. Set up working environment
###############################
rm(list=ls())
library(dplyr)
library(PEcAn.DB)

# load prior information
all.priors <- read.csv('~/data_files/priors-04-2020.csv', header = T)
pfts <- as.vector(unique(all.priors$pft))
if (pfts[length(pfts)]=='') pfts = pfts[1:(length(pfts)-1)]

# make bety connection
bety <- list(user = 'bety',
             password = 'bety',
             host = 'postgres',
             dbname = 'bety',
             driver = "PostgreSQL",
             write = TRUE)
modeltype =  "1000000008" # linkages model 
dbcon <- db.open(bety)
on.exit(db.close(dbcon))

# variable names that are needed for LINKAGES species: 
pars = c('DMAX','DMIN','HTMAX','AGEMX','DBHMAX','Gmax','SPRTND','SPRTMN','SPRTMX','MPLANT','D3',
         'FROST','CM1','CM2','CM3','CM4','CM5','SLA','SLTA','SLTB')

# get variable id numbers for each LINKAGES parameter
query.text <- paste(
  "SELECT name, variables.id",
  "FROM variables") 
vars.id = db.query(query = query.text, con = dbcon) %>%
  filter(name %in% pars) 
vars.id$id = sapply(vars.id$id, toString)

##################################
# 2. Loop through the PFTs
##################################

# get pft id numbers for each pft in the matrix
pft_names = pfts
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

  # read in updated priors for this pft 
  new.priors = all.priors %>% 
    filter(pft == pftres$name) %>%
    dplyr::select(pft, varname, distn, parama, paramb)

########################################
# 3. Loop through the priors for the PFT
########################################
  
  # loop through the priors for this pft
  for (i in 1:length(new.priors$pft)){
    
    p = new.priors[i,]
    
    # pick out updated prior information
    varname = p$varname # variable name not pft
    print(toString(varname))
    distn = p$distn
    parama = p$parama
    paramb = p$paramb
    varid = vars.id %>% filter(name == varname) %>% dplyr::select(id)
    pftid = pftid
    
##########################################################
# 4. Is the current prior OK? (i.e. no changes necessary)
##########################################################
    
    # get current prior associated with this pft in BETY
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
    
    # if there is more than one prior, remove one
    # if one is actually correct, we will find it again later
    while (length(existing.prior$name) > 1){
      print('there is more than one prior for this PFT and variable combination, removing all of them')
      db.query(paste0("DELETE from pfts_priors WHERE prior_id = ",existing.prior$prior_id[1]), dbcon)
      
      # re-query
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
    }
    
    # check to see if existing priors is the same as what we need
    if (length(existing.prior$name) > 0){
      if (existing.prior$parama == parama & existing.prior$paramb == paramb & existing.prior$distn == distn){
        add = FALSE
        print('existing.prior is good')
      }
    }
    
##################################################
# 5. Current prior is not OK; need a different one
##################################################
    
    # if add is still TRUE, we need a new prior because the old one doesn't line up
    if (add){
      
      print('need to add new prior')
      
###################################################
# 6. Is there a good prior already in the database?
###################################################

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
        
######################################
# 7. Need to add an entirely new prior
######################################
        
      }else{ # we need to input a new prior 
        
        print('no prior available - need to insert new')
        
        # remove all associations with pft for this variable if it has some
        if (length(existing.prior) > 0){
          # update database query
          db.query(paste0("DELETE from pfts_priors WHERE prior_id = ",existing.prior$prior_id), dbcon)
        }
        
        # input database query (first into priors table)
        now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        db.query(paste0("INSERT INTO priors ",
                        "(variable_id, phylogeny, distn, parama, paramb, notes, created_at, updated_at)",
                        " VALUES (",varid, ", '', '", distn, "', ", parama, ", ", paramb, ", '", pftres$name, 
                        " linkages', '", now, "', '", now, "')"), dbcon)
        new.prior.id = db.query(paste0("SELECT id FROM priors where variable_id = ",varid, " AND distn = '", 
                                       distn,"' AND parama = ", parama, " AND paramb = ", paramb), dbcon)
        
        # (then into pfts_priors)
        db.query(paste0("INSERT INTO pfts_priors ",
                        "(pft_id, prior_id, created_at, updated_at)", 
                        "VALUES (", pftid, ", ", new.prior.id, ", '", now, "', '", now, "')"), dbcon)
      }
    }
  }
}



