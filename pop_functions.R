############ POPULATION FUNCTIONS AFRICA   ##############
# Functions used for the population gridding of africa  #
#                                                       #
#                                                       #
# Author: Niklas Boke-Olen                              #
# niklas.boke-olen@nateko.lu.se                         #
# niklasbokeolen@gmail.com                              #
#                                                       #
#########################################################
to_urban_mask_parallell <- function(urb_frac,pop_grid,pixels,cores){
  urb_frac1 <- round(urb_frac*pixels)
  ind_rast <- urb_frac1; 
  ind_rast [,] <- 1:length(ind_rast[,]); 
  ind_mask <- urb_frac1; 
  ind_mask[ind_mask>0] <- 1; ind_mask[ind_mask==0] <- NA; 
  ind_rast <- ind_mask*ind_rast;
  ind.df <- na.omit(as.data.frame(ind_rast))
  
  
  registerDoMC(cores)  #change to CPU cores  
  
  
  cutval <- NA
  
  cutval1 <- foreach(i=1:length(ind.df$layer),.combine=rbind) %dopar% {
    
    urbP_cell <- urb_frac1[ind.df$layer[i]]
    xy <- xyFromCell(urb_frac1, ind.df$layer[i], spatial=FALSE) - 0.25
    ext <- extent(xy[1],xy[1]+0.5,xy[2],xy[2]+0.5)
    a <- intersect(ext,extent(pop_grid)) #check if  should be processed. only if within area chosen
    
    if(!is.null(a)){
      if(a == ext){
        sort_pix <- sort(extract(pop_grid,ext))
        if(length(sort_pix)>0){
          index_sort <- (length(sort_pix)- urbP_cell)
          if(index_sort<1){ #if all pixels should be included allow for this (could happend with alot of NA for example)
            cutval[i] <- max((sort_pix[1] - 1),0)
          }else{
            cutval[i] <- sort_pix[index_sort]
            
            if(sum(sort_pix>cutval[i])<urbP_cell){  #if we have equal values select all equal at the cutoff to be urban.(should only be pairs)
              diff <- sort_pix - sort_pix[index_sort]
              cutval[i] <- sort_pix[diff<0][length(sort_pix[diff<0])]
            }
          }
        }
      }
    }
    cutval[i]
  }
  
  beginCluster(cores)
  cutval_grid <- urb_frac1; cutval_grid[,] <- NA
  cutval_grid[ind.df$layer[1:length(ind.df$layer)]] <- cutval1
  cutval_rsmp <- resample(cutval_grid,pop_grid,method='ngb')
  urban_mask <- pop_grid > cutval_rsmp
  urban_mask <-  reclassify(urban_mask, cbind(NA, 0))
  endCluster()
  return(urban_mask)
}
