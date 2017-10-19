###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

# Load in the CBI function
source(paste(getwd(), "/01_Scripts/Minor_Modules/5___CBI_Function.R", sep=""))

require(dismo)
require(maptools)

kde_raster <- raster("000_GIS_LAYERS/Brazil_Masked_GIS_Layers/KDE_Raster/kde_raster.asc")

CBI_results <- as.data.frame(matrix(nrow=length(species), ncol=12))
rownames(CBI_results) <- species
colnames(CBI_results) <- c("CBI_1", "CBI_2","CBI_3","CBI_4","CBI_5", "CBI_mean", "AUC_1", "AUC_2","AUC_3","AUC_4","AUC_5", "AUC_mean")


#Load in the model AUCs
for(x in 1:length(species)){
  AUCs <- read.csv(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], "/maxentResults.csv", sep=""))
  CBI_results[x,7:12] <- AUCs[,9]
}





#Find the model CBIs
for(x in 1:length(species)){
  
  writeLines(paste("Finding CBI for", species[[x]]))
  
  background_data <- read.csv(paste("03_Modelling/10_Background_Data/Biased/", species[[x]], ".csv", sep=""))[,-1]
  species_data <- read.csv(paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
  
  for(y in 0:4){
    
    model <- raster(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], "_", y+1, ".tif", sep=""))
    
    
    #partitions <- read.csv(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], "/species_", y, "_samplePredictions.csv", sep=""))
    #partitions.test <- partitions[which(partitions[,3] == "test"),]
    
    #partitions.test <- species_data[partitions.test[,1],]
    
    #partitions.test <- extract(model, SpatialPoints(partitions.test[,2:1]))
    partitions.test <- extract(model, species_data[,2:1])
    partitions.bg <- extract(model, SpatialPoints(background_data[,3:4]))
    
    
    CBI <- contBoyce(predAtPres = partitions.test, predAtBg = partitions.bg, numClasses = 10)
    
    CBI_results[x,y+1] <- CBI
  }
  CBI_results[x,6] <- mean(as.numeric(CBI_results[x,1:5]))
  writeLines(paste("...", mean(as.numeric(CBI_results[x,1:5])), "\n"))
  write.csv(CBI_results, file="03_Modelling/CBI_results.csv")
}

rm(CBI, CBI_results, partitions, partitions.test, partitions.bg, x, y, background_data, species_data, AUCs)




