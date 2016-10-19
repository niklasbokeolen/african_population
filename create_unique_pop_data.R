############ UNIQUE POPULATION 2000 AFRICA ##############
# Combines SSP population projections and               #
# RCP urban fractions to grid population for the future #
#                                                       #
# Author: Niklas Boke-Olen                              #
# niklas.boke-olen@nateko.lu.se                         #
#                                                       #
#########################################################


### INITIAL SETTINGS AND LIBRARY ####
rm(list=ls())
setwd('/media/niklas/data/Satellite/grazing/pop_scripts_publication/')
library(raster)
library(rgdal)
library(SDMTools)

##directory settings 
data_dir <- 'data/'   ## directory of data
out_dir <- 'output/' ##output directory


## LOAD DATA ###

## population count Africa 2000 from WorldPoP 
pop_2000_file <- paste(data_dir,'Worldpop/ap00v4_TOTAL_adj.tif',sep='')
pop_2000 <- raster(pop_2000_file)

## MATCH WATER MASK and apply to pop_2000
water.mask  <- raster(paste(data_dir,'gis_data/poly2raster_coast_arcmap.tif',sep='')) ## created with polygon2raster in arcmap. cell size set to match pop_2000.
water.mask <- crop(water.mask ,pop_2000)
water.mask  <- resample(water.mask ,pop_2000,method='ngb')
water.mask [!is.na(water.mask )] <- 1
pop_2000 <- pop_2000*water.mask

### rcp based urban frac from hurtt as sample for refernce use
urb_frac.stack <-stack(paste(data_dir,'lu_out_africa/hurtt_urbanfrac_rcp_8_5.tif',sep=''))
urb_frac <- crop(urb_frac.stack[[1]],extent(-25.36,63.5,-40.4,37.54)+0.5,snap='out') ## crop for africa extent +0.5 


## load distance to roads created in ArcMap to matcg pop_2000 pixels
dist_roads <- raster('/media/niklas/data/Satellite/grazing/pop/data/gis_data/dist_to_road_arcHighres.tif')

#extent pop_2000 so that it matches urban frac extent
pop <- extend(pop_2000,extent(urb_frac))

#Create empty grid to be filled with the unique population dataset.
pop_unique <- pop; pop_unique[,] <- NA

same_all <- NA; pp <- 1

## loop through each urb_frac pixel
for(row in 1:nrow(urb_frac)){ 
  for(col in 1:ncol(urb_frac)){ 
    
    ## locate gridcell
    x <-  xFromCol(urb_frac,col) - 0.25 
    y <- yFromRow(urb_frac,row) -0.25
    ext <- extent(x,x+0.5,y,y+0.5)
    
    #check if this should be processed. only if within area chosen
    a <- intersect(ext,extent(pop)) 
    b <- intersect(ext,extent(dist_roads))
    
    if(!is.null(a) & !is.null(b)){ 
      if(a == ext & b == ext){
        
        
        pop_cell <- crop(pop,ext)  ## crop population to gridcell.
        
        if(sum(pop_cell[!is.na(pop_cell)])>0){ ##perform only if we have population value inside the urban frac pixel. if all NA do not process, if all zero process.
          
          # distance to road
          ##rescaled to inverse distance to road.
          dist_roads_cell <- crop(dist_roads,ext)
          dist_roads_cell <- 0.00001 + (max(dist_roads_cell[,])- dist_roads_cell) / (max(dist_roads_cell[,])- min(dist_roads_cell[,])) / 10000 
          
          ## make sure grids match exactly
          pop_rsmp <- resample(pop_cell,dist_roads_cell,method='ngb')
          
          ##center of gravity
          ##rescaled to inverse distance to COG.
          cog <- as.numeric(COGravity(pop_rsmp, y = NULL, z = NULL, wt = NULL))
          cog.spdf <- SpatialPointsDataFrame(coords = data.frame(latitude=cog[1],longitude=cog[3],data=1), 
                                             data = data.frame(data=1), proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          cog_dist <- pop_rsmp;  cog_dist[,] <- NA
          cog_dist[cellFromXY(cog_dist, cog.spdf) ] <- 1
          cog_dist <- raster::distance(cog_dist)
          cog_dist <- 0.00001 + (max(cog_dist[,])- cog_dist) / (max(cog_dist[,])- min(cog_dist[,])) / 10000
          # 
          
          ## add together
          pop_w <- pop_rsmp + dist_roads_cell + cog_dist   
          pop_unique[ext] <-  pop_w[,] ##update resulting unique dataset

          ## PRINTOUT TO SCREEN
          max_same <- max(t(as.data.frame(table(pop_w[,]))[,2])) ## calculate maximum equal pixels for printing purposes
          uniq <- length(unique(pop_w[!is.na(pop_w)])) == length(pop_w[!is.na(pop_w )]) || max_same<3
          print(paste(row,col,Sys.time(),uniq,length(unique(pop_w[!is.na(pop_w)])),max_same,sum(is.na(pop_w[,])))) ## print progress
          
        }else if(!is.null(pop_cell[!is.na(pop_cell)])){ ##only zero population, add small values based on inverse distance to road
          
          # distance to road
          ##rescaled to inverse distance to road.
          dist_roads_cell <- crop(dist_roads,ext)
          dist_roads_cell <- 0.00001 + (max(dist_roads_cell[,])- dist_roads_cell) / (max(dist_roads_cell[,])- min(dist_roads_cell[,])) / 10000 
          
          ## make sure grids match exactly
          pop_rsmp <- resample(pop_cell,dist_roads_cell,method='ngb')
          
          ## add together
          pop_w <- pop_rsmp + dist_roads_cell  #

          pop_unique[ext] <-  pop_w[,] ##update rsulting unique dataset
          
          ## PRINTOUT TO SCREEN
          max_same <- max(t(as.data.frame(table(pop_w[,]))[,2])) ## calculate maximum equal pixels for printing purposes
          uniq <- length(unique(pop_w[!is.na(pop_w)])) == length(pop_w[!is.na(pop_w )]) || max_same<3
          print(paste(row,col,Sys.time(),uniq,length(unique(pop_w[!is.na(pop_w)])),max_same,sum(is.na(pop_w[,])))) ## print progress
          
        }
       
        
      }
    }
    
  }

  
}

## SAVE OUTPUT
writeRaster(pop_unique,paste(out_dir,'pop2000_unique.tif',sep=''),format='GTiff',overwrite=T,NAflag=-9999)

# 

