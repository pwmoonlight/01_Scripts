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
library(rrcov)

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
    pca_results <- PcaClassic(pca_data)
    
    pca_index <- pca.distances(pca_results, pca_data, r=rankMM(pca_data))
    
    index <- which(as.data.frame(pca_index@flag)[,1] == FALSE)
    
    if(length(index)>0){distribution_data <- distribution_data[-index,]}
  }
  
  if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-dim(distribution_data)[1],"fewer points..."))}
  if(length(distribution_data[,1])>0){
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_6/", species[[x]], ".csv", sep=""))
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]], ".csv", sep=""))
  }
}

rm(distribution_data, dims, pca_data, pca_results, pca_index)



