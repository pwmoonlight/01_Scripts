###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4", showWarnings = F)


### ------------
### Read in Mask

require(raster)
mask <- raster("000_GIS_LAYERS/brazil_mask_minus_centroids.tif")


for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4/", species[[x]], sep=""), showWarnings = F)
    
    distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
    dims <- dim(distribution_data)[1]

    if(length(distribution_data[,1])>0){
      distribution_data[,12] <- extract(mask, SpatialPoints(distribution_data[,2:1]))
      index <- which(is.na(distribution_data[,12]))
      if(length(index)>0){
        distribution_data <- distribution_data[-index,]
      }
    
      if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-length(distribution_data[,1]),"fewer points..."))}

      if(length(distribution_data[,1])>0){
        distribution_data <- distribution_data[,-12]
        write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4/", species[[x]], ".csv", sep=""))
      }
    }
  }
}
