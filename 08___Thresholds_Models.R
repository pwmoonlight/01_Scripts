###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

require(raster)

dir.create("03_Modelling/12_Thresholded_Models", showWarnings=F)

## Define Threshold Function
rc <- function(x) {
  ifelse(x <=  threshold_value, 0,
         ifelse(x >  threshold_value, 1, NA)) }



CBIs <- read.csv("03_Modelling/CBI_results.csv", row.names = 1)
CBIs <- CBIs[which(CBIs[,6] >= 0.5),]
species <- rownames(CBIs)

<<<<<<< HEAD
example_model <- model <- raster(paste("11_models/Bias_Spatial_Filtering/", species[[1]], ".tif", sep=""))
=======
example_model <- model <- raster(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[1]], ".tif", sep=""))
>>>>>>> 05df30d7c27e26ba5ca13b66d316e101f9462003
threshold_value <- 1.1
example_model <- calc(example_model, fun=rc)

for(x in 1:length(species)){
<<<<<<< HEAD
  model <- raster(paste("11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
  threshold_value <- read.csv(paste("11_models/Bias_Spatial_Filtering/", species[[x]], "/maxentResults.csv", sep=""))
  threshold_value <- threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]
=======
  model <- raster(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
  threshold_value <- read.csv(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], "/maxentResults.csv", sep=""))
<<<<<<< HEAD
  if(length(threshold_value$X10.percentile.training.presence.logistic.threshold[6]) == 1){
    threshold_value <- threshold_value$X10.percentile.training.presence.logistic.threshold[6]
  }
  if(length(threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]) == 1){
    threshold_value <- threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]
  }
  #sp.circle <- readShapePoly(paste("03_Modelling/05_Species_Circles/", species[[x]], ".shp", sep=""))
=======
  threshold_value <- threshold_value$X10.percentile.training.presence.logistic.threshold[6]
>>>>>>> 05df30d7c27e26ba5ca13b66d316e101f9462003
  
  sp.circle <- readShapePoly(paste("03_Modelling/05_Species_Circles/", species[[x]], ".shp", sep=""))
>>>>>>> 0a325e238c5a9880c3629d2f114d65d86dde2fd2
  
  writeLines(paste("\nWorking on", species[[x]]))
  writeLines(paste("...The 10 percentile training presence logistic threshold  for", species[[x]], "is:", threshold_value))#
  
  ## Perform Thresholding
  thresholded_model <- calc(model, fun=rc)
  
  #thresholded_model <- mask(thresholded_model, sp.circle)
  
  thresholded_model <- merge(thresholded_model, example_model)
  
  plot(thresholded_model, main = species[[x]])
  
  writeRaster(thresholded_model,  paste("03_Modelling/12_Thresholded_Models/", species[[x]], ".tif", sep = ""))
}



rm(CBIs, model, thresholded_model, threshold_value, species, x)
