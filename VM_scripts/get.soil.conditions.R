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

# if you run into errors, just continue to run the next lines 
# if it still doesn't work, reduce the number of decimals in the lat and lon coordinates and try again 

