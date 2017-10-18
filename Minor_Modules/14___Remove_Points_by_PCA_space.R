###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

# Load in required packages
require(dismo)
require(rgeos)
require(raster)
library(ade4)

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_6", showWarnings = F)

PCA <- c(1, 6, 12, 13, 14, 15, 21, 26, 27, 28, 30, 31, 33)
bg <- lapply(list.files(path="000_GIS_LAYERS/Brazil_Masked_GIS_Layers", pattern="*.tif$", full.names = T), raster)
bg <- bg[PCA]
bg <- stack(bg)







for(x in 1:length(species)){
  distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
  dims <- dim(distribution_data)[1]
  if(length(distribution_data[,1])>4){
    
    pca_data <- extract(bg, SpatialPoints(distribution_data[,2:1]))
    pca_results <- dudi.pca(pca_data,scannf = F, nf = 2)
    
    SD1 <- sd(pca_results$li[,1])
    SD1_times <- pca_results$li[,1]/(SD1*1.5)
    for(i in 1:length(SD1_times)){
      if(SD1_times[i]< 0){SD1_times[i] <- 0-SD1_times[i]}
    }
    
    SD1_index <- which(SD1_times > 1)
    
    if(length(pca_results$li[1,])>1){
      SD2 <- sd(pca_results$li[,2])
      SD2_times <- pca_results$li[,1]/(SD2*1.5)
      for(i in 1:length(SD2_times)){
        if(SD2_times[i]< 0){SD2_times[i] <- 0-SD2_times[i]}
      }
      
      SD2_index <- which(SD2_times > 1)
      
      SD1_index <- c(SD1_index, SD2_index)
      SD1_index <- unique(SD1_index)
    }
    
    distribution_data <- distribution_data[-SD1_index,]
  }
  
  if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-dim(distribution_data)[1],"fewer points..."))}
  if(length(distribution_data[,1])>0){
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_6/", species[[x]], ".csv", sep=""))
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]], ".csv", sep=""))
  }
}

rm(distribution_data, SD1_index, SD2_index, SD1, SD2, SD1_times, SD2_times, dims)








