###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3", showWarnings = F)



### ---------------------
### Read in Altitude Data

require(raster)
alt <- raster("000_GIS_LAYERS/Brazil_Masked_GIS_Layers/alt/1km.asc")


### -------------------------
### Cycle Through The Species

for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3/", species[[x]], sep=""), showWarnings = F)
    distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_2/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
    
    dims <- dim(distribution_data)[1]
    if(length(distribution_data[,1])>0){
      for(y in length(distribution_data[,1]):1){
        if(is.na(distribution_data[y,1])){distribution_data <- distribution_data[-y,]}
      }
    
      distribution_data[,8] <- gsub("[^0-9]","", distribution_data[,8])
      distribution_data[,8] <- as.numeric(distribution_data[,8])
      distribution_data[,9] <- gsub("[^0-9]","", distribution_data[,9])
      distribution_data[,9] <- as.numeric(distribution_data[,9])
      distribution_data[,12] <- extract(alt, SpatialPoints(distribution_data[,2:1]))
      distribution_data[,13] <- NA
      distribution_data[,14] <- NA
      
      for(y in 1:length(distribution_data[,1])){
        if(!is.na(distribution_data[y,12])){
          if(!is.na(distribution_data[y,8])){
            if(distribution_data[y,8] != 0){
              distribution_data[y,13] <- distribution_data[y,8] - distribution_data[y,12]
              if(distribution_data[y,13]<0){distribution_data[y,13] <- 0-distribution_data[y,13]}
            }
          }
          if(!is.na(distribution_data[y,9])){
            if(distribution_data[y,9] != 0){
              distribution_data[y,14] <- distribution_data[y,9] - distribution_data[y,12]
              if(distribution_data[y,14]<0){distribution_data[y,14] <- 0-distribution_data[y,14]}
            }
          }
          distribution_data[y,13] <- max(distribution_data[y,13], distribution_data[y,14], na.rm=T)
        }
      }
      
      distribution_data <- distribution_data[,-14]
    
      index <- which(distribution_data[,13]>500)
      if(length(index)>0){distribution_data <- distribution_data[-index,]}
      
      distribution_data <- distribution_data[,-13]
      distribution_data <- distribution_data[,-12]
    
      if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-length(distribution_data[,1]),"fewer points..."))}
      write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3/", species[[x]], ".csv", sep=""))
    }
  }
}

