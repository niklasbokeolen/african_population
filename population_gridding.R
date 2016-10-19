############ POPULATION GRIDDING AFRICA #################
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
library(foreach)
library(doMC)
cores <- 4 #cores to use for parallell processing
registerDoMC(cores)
source('pop_functions.R')
##directory settings 
data_dir <- 'data/'   ## directory of data
out_dir <- 'output/' ##output directory


## LOAD DATA ###

## population density year 2000
pop2000_unique <- raster(paste(out_dir,'pop2000_unique.tif',sep='')) ## created in file create_unique_pop_data.R

##world map , take out africa and add index field
world <- readOGR(paste(data_dir,"WorldMap/TM_WORLD_BORDERS-0.3.shp",sep=""),layer="TM_WORLD_BORDERS-0.3") 
africa <- world[data.frame(world)$REGION==2,] #take out africa
africa.df <- data.frame(africa)
africa.df$ind <- seq(1,length(africa.df$FIPS),1)
africa@data <- africa.df

country_code <- read.csv(paste(data_dir,'popssps/countryCodes.csv',sep='')) ### SSP specific population data countrycodes

##settings
runyears <- 2000 #c(2000:2100) ## years to grid
ssp_list <- c(  1  , 1 , 2 , 3 , 4 ,  5  ,5  , 4   , 3   , 2   , 3   , 1   , 2   , 4   ,  5  ) ##ssp run list
rcp_list <- c('4_5','6','6','6','6','8_5','6','4_5','8_5','8_5','4_5','2_6','4_5','8_5','4_5') ##rcp run list
pixels <- 3600 #number of population pixels within each hurtt 0.5 degree cell



for(scen in 1:1){ #length(rcp_list)){
  rcp <- rcp_list[scen]  
  ssp <- ssp_list[scen]   
  
  ### rcp urban frac from hurtt
  urb_frac.stack <-stack(paste(data_dir,'lu_out_africa/hurtt_urbanfrac_rcp_',rcp,'.tif',sep=''))
  
  ## SSP population estimates.
  pop <- read.csv(paste(data_dir,'popssps/popSSP',ssp,'.csv',sep=''),header=TRUE,row.names=1)
  colnames(pop) <- c(1:length(pop[1,]))
  urban<- read.csv(paste(data_dir,'popssps/urbanShareSSP',ssp,'.csv',sep=''),header=TRUE,row.names=1)
  colnames(urban) <- c(1:length(urban[1,]))

  for(year in runyears){
    
    ## PRINT PROGRESS
    print(paste('Gridding','year ',year,'---',Sys.time(),'RCP=',rcp,'SSP=',ssp,'scen=',scen,'of',length(rcp_list)))
    
    ## define name of output files
    outfile_main <- paste(out_dir,'grid_pop_count',year,'_SSP',ssp,'_RCP',rcp,'.tif',sep='')
    
    ##Load Hurtt urban fraction 
    urb_frac <- crop(urb_frac.stack[[year-1999]],extent(africa)+0.5,snap='out')
    
    ## load population from last year (or for first year load unique data)
    if(!file.exists(outfile_main)){ #only run if file not exist allready, allows to stop script and restart
      if(year==2000){
        pop_grid <- crop(pop2000_unique,extent(urb_frac),snap='out') #country population density
      }else{
        pop_grid <-  raster(paste(out_dir,'grid_pop_count',year-1,'_SSP',ssp,'_RCP',rcp,'.tif',sep=''))
      }
      
      ## Create urban mask
      urb_mask <- to_urban_mask_parallell(urb_frac,pop_grid,pixels,cores)
      
      ## create country number raster from africa shape file
      cnum_file <- paste(out_dir,'africa_countrynum.tif',sep='')
      if(!file.exists(cnum_file)){
        africa_countrynum <- rasterize(africa,urb_mask,field='ind') ##slow portion. for speed-up it should be saved to tif file first time run
        writeRaster(africa_countrynum,cnum_file,format='GTiff',NAflag=-9999)
      }else if (!exists("africa_countrynum")){ ##load if not loaded.
        africa_countrynum <- raster(cnum_file )
      }
      

      ## distribute SSP population for each copuntry using the urban mask population for year before (unique for year 2000)
      c_pop_grid.list <- foreach(i=1:length(africa.df$ISO3)) %dopar% { #perform in parallell , returns a list with a grid per country
        
        ## extract country and check if it should be processed
        country <- africa[africa.df$ISO3==africa.df$ISO3[i],] 
        a <- intersect(extent(urb_frac),extent(country))
        b <-  intersect(extent(africa_countrynum),extent(country))
        c <- FALSE
        if(!is.null(b) && !is.null(a)) c <- floor(a) == floor(extent(country)) #make sure that the entire country is present in data
        
        nr <- as.numeric(country_code$nr[as.character(country_code$iso)==as.character(africa.df$ISO3[i])]) ##extract country number from country_code
        if(length(nr) > 0 && c){ #make sure country exist in SSP list
          
          c_pop <- pop[year-1999,paste(nr)] #country population total from ssp file
          c_pop_urb <- c_pop*urban[year-1999,paste(nr)]/100 #number of people in urban
          c_pop_rur <- c_pop - c_pop_urb #number of rural people
          
          
          
          cnum_rast <- crop(africa_countrynum,country) #crop countrynumber raster
          cnum_rast[!cnum_rast==africa.df$ind[i]] <- NA; cnum_rast[!is.na(cnum_rast)] <- 1 #create country mask
          c_pop_grid <- crop(pop_grid,country)  *  cnum_rast #get country population grid from year before
          c_urb_mask <- crop(urb_mask,country) * cnum_rast #country urban mask
          
          if(c_pop_urb==0){ c_urb_mask[c_urb_mask==1] <- 0 }# in the cases where urban population values are missing distirbute for all in same way
         
          
          
          c_pop_urb_grid <- c_pop_grid * c_urb_mask ## create a populaiton grid for urban pixels
          if(sum(c_pop_urb_grid[,],na.rm=T)>0){  #only if we have urban population
            
            c_pop_urb_grid <-   c_pop_urb_grid  *  c_pop_urb / sum(c_pop_urb_grid[,],na.rm=T)  #grid SSP urban population
            
          }else{ # no urban pixels
            if(c_pop_urb>0) c_pop_rur <- c_pop_rur + c_pop_urb #Special case if we do not have urban pixels but have urban population put that into rural.
            
            c_pop_urb_grid <- c_pop_grid; c_pop_urb_grid[!is.na(c_pop_urb_grid)] <- 0
            
          }
          
          ##distribute rural population
          c_pop_rur_grid <-  c_pop_grid * !c_urb_mask
          c_pop_rur_grid <- c_pop_rur_grid *  c_pop_rur / sum(c_pop_rur_grid[,],na.rm=T) #gridded rural pop
          c_pop_grid <- c_pop_urb_grid + c_pop_rur_grid #combine rural and urban grid.
          c_pop_grid #return population grid for country, foreach loop specific
          
          
        }
        
      }
      
      ##determine the order to combien the countries, sort by spatial location for speed up.
      p <- 1; loc <- data.frame(i=NA,x=NA,y=NA)
      for(i in 1:length(c_pop_grid.list)){
        if(!is.null(c_pop_grid.list[[i]])){
          ext <- extent(c_pop_grid.list[[i]])
          loc <- rbind(loc,data.frame(i=i,x = (ext[2]+ext[1])/2, y= (ext[4]+ext[3])/2))
          p <- p+1
        }
      }
      loc <- na.omit(loc)
      loc$group <- NA; loc$group[loc$y>=15] <- 1
      loc$group[loc$y>= 5 & loc$y< 15] <- 2
      loc$group[loc$y>= -10 & loc$y< 5] <- 3
      loc$group[loc$y>= -50 & loc$y< -10] <- 4
      loc.ordered <- loc[with(loc, order(group, x)), ]
      new_list <- list()
      for(i in 1:length(loc.ordered$i)){ new_list[[i]] <-  c_pop_grid.list[[loc.ordered$i[i]]] }##reorder list for speedup
      
      pop_result_grid <- do.call(merge,new_list) ##merge rasters in list
      
      ##write output to file
      writeRaster(pop_result_grid,outfile_main,format='GTiff',overwrite=T,NAflag=-9999)
      removeTmpFiles(h=0.01) ##remove temporary files to avoid out of memory
    }
  }
  

  removeTmpFiles(h=0.01) ##remove temporary files to avoid out of memory
} 


removeTmpFiles(h=0.001)  ##remove temporary files to avoid out of memory
