########################################################################################
########################################################################################

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
in.path = '/data/dbfiles/met_data/ROOSTER/linkages' ## path to file directory (not including file name)
in.prefix = 'bcc.csm1.1_001.01.Rdata'
siteid = 1000026714 ## site id number, directions on how to obtain in google doc
startdate = '1600-01-01 00:00:00' ## adjust date years as needed for available data
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

########################################################################################
########################################################################################

## finding soil 
library(XML)
library(dplyr)
library(PEcAn.data.land)
workdir = '/home/carya'
PEcAn.data.land::extract_soil_gssurgo(outdir = workdir,lat = 47.1950765,lon = -95.1648107)
ncin <- ncdf4::nc_open(file.path(workdir,'gSSURGO_soil_2.nc'))
print(paste('%clay =',ncdf4::ncvar_get(ncin,'fraction_of_clay_in_soil')[1]*100))
print(paste('%sand =',ncdf4::ncvar_get(ncin,'fraction_of_sand_in_soil')[1]*100))


