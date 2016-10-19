############ STACK RCP urban frac data     ##############
# Load urban frac netcdf data from Hurt et al (2011)   #
# And create stacks that can be used in further gridding#
#                                                       #
# Author: Niklas Boke-Olen                              #
# niklas.boke-olen@nateko.lu.se                         #
#                                                       #
#########################################################

##NETCDF RCP urban fraction data downloaded from https://daac.ornl.gov/VEGETATION/guides/Land_Use_Harmonization_V1.html
# Chini, L.P., G.C. Hurtt, and S. Frolking. 2014. Harmonized Global Land Use for Years 1500 – 2100, V1. Data set. Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA. http://dx.doi.org/10.3334/ORNLDAAC/1248

library(raster)
library(rgdal)
world <- readOGR(paste(data_dir,"WorldMap/TM_WORLD_BORDERS-0.3.shp",sep=""),layer="TM_WORLD_BORDERS-0.3") 
world.df <- data.frame(world)
africa <- world[world.df$REGION==2,] #take out africa

##1700-2005 urban fraction
past <- brick('/media/niklas/data/Satellite/grazing/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_gurbn.nc4')
past <- crop(past,africa)
#2000-2004
past <- past[[301:305]]


# MESSAGE (8.5 W/m2), AIM (6 W/m2), GCAM (4.5 W/m2), and IMAGE (2.6 W/m2).

# rcp8.5 2005 -2100
r <- brick('/media/niklas/data/Satellite/grazing/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2t1.v1_message.v1_gurbn.nc4')
r <- crop(r,africa)
combine <- stack(past,r)
names(combine) <- paste('y',c(2000:2100),sep='')
writeRaster(combine,'/media/niklas/data/Satellite/grazing/lu_out_africa/hurtt_urbanfrac_rcp_8_5.tif',format='GTiff')


# rcp 6 2005 -2100

r <- brick('/media/niklas/data/Satellite/grazing/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2.v1_aim.v1.1_gurbn.nc4')
r <- crop(r,africa)
combine <- stack(past,r)
names(combine) <- paste('y',c(2000:2100),sep='')
writeRaster(combine,'/media/niklas/data/Satellite/grazing/lu_out_africa/hurtt_urbanfrac_rcp_6.tif',format='GTiff')


# rcp 2.6 2005 -2100

r <- brick('/media/niklas/data/Satellite/grazing/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2.v1_image.v1.1_gurbn.nc4')
r <- crop(r,africa)
combine <- stack(past,r)
names(combine) <- paste('y',c(2000:2100),sep='')
writeRaster(combine,'/media/niklas/data/Satellite/grazing/lu_out_africa/hurtt_urbanfrac_rcp_2_6.tif',format='GTiff')


# rcp 4.5 2005 -2100

r <- brick('/media/niklas/data/Satellite/grazing/LAND_USE_HARMONIZATION_V1_1248/data/LUHa_u2.v1_minicam.v1_gurbn.nc4')
r <- crop(r,africa)
combine <- stack(past,r)
names(combine) <- paste('y',c(2000:2100),sep='')
writeRaster(combine,'/media/niklas/data/Satellite/grazing/lu_out_africa/hurtt_urbanfrac_rcp_4_5.tif',format='GTiff')



