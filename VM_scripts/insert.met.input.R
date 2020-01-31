########################################################################################
########################################################################################

rm(list=ls())
# restart R just to start with a clean slate 

## INSERTING MET DATA TO BETY
# adding met file record and input record to BETY database

#1. establish BETY connection
drv <- RPostgreSQL::PostgreSQL()
con = DBI::dbConnect(drv,
                     host = 'postgres',
                     dbname = 'bety',
                     user = 'bety',
                     password = 'bety'
                     )

#2. set variables for input file
in.path = '/data/dbfiles/met_data/NORTHROUND/linkages' ## path to file directory (not including file name)
in.prefix = 'bcc.csm1.1_001.01.Rdata'
siteid = 1000026710 ## site id number, directions on how to obtain in google doc
startdate = '0850-01-01 00:00:00' ## adjust date years as needed for available data
enddate = '2015-12-31 00:00:00' 

# do not change these if entering linkages input data
mimetype = 'text/plain'
formatname = 'LINKAGES met'

#3. insert input and file record in database (for met ensembles, only need to enter one of the ensembles for use in our workflow)
library(PEcAn.DB)
library(DBI)
file_input = PEcAn.DB::dbfile.input.insert(in.path = in.path,
                                 in.prefix = in.prefix,
                                 siteid = siteid,
                                 startdate = startdate,
                                 enddate = enddate, 
                                 mimetype = mimetype,
                                 formatname = formatname,
                                 con = con,
                                 hostname = 'docker')


