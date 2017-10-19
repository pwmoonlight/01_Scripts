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
<<<<<<< HEAD
library(rrcov)
=======
<<<<<<< HEAD
library(rrcov)
=======
>>>>>>> ff303605a6e19477b6f08ce833057f75dc4dc1ab
>>>>>>> 1fe41a805103cf4e6d7b88e73c371acf9c27d52d

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
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 1fe41a805103cf4e6d7b88e73c371acf9c27d52d
    pca_results <- PcaClassic(pca_data)
    
    pca_index <- pca.distances(pca_results, pca_data, r=rankMM(pca_data))
    
    index <- which(as.data.frame(pca_index@flag)[,1] == FALSE)
    
    if(length(index)>0){distribution_data <- distribution_data[-index,]}
<<<<<<< HEAD
=======
=======
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
>>>>>>> ff303605a6e19477b6f08ce833057f75dc4dc1ab
>>>>>>> 1fe41a805103cf4e6d7b88e73c371acf9c27d52d
  }
  
  if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-dim(distribution_data)[1],"fewer points..."))}
  if(length(distribution_data[,1])>0){
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_6/", species[[x]], ".csv", sep=""))
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]], ".csv", sep=""))
  }
}

<<<<<<< HEAD
rm(distribution_data, dims, pca_data, pca_results, pca_index)
=======
<<<<<<< HEAD
rm(distribution_data, dims, pca_data, pca_results, pca_index)
=======
rm(distribution_data, SD1_index, SD2_index, SD1, SD2, SD1_times, SD2_times, dims)





>>>>>>> ff303605a6e19477b6f08ce833057f75dc4dc1ab
>>>>>>> 1fe41a805103cf4e6d7b88e73c371acf9c27d52d



