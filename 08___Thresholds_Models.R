###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("12_Thresholded_Models", showWarnings=F)

## Define Threshold Function
rc <- function(x) {
  ifelse(x <=  threshold_value, 0,
         ifelse(x >  threshold_value, 1, NA)) }



CBIs <- read.csv("CBI_results.csv", row.names = 1)
CBIs <- CBIs[which(CBIs[,6] >= 0.5),]
species <- rownames(CBIs)

example_model <- model <- raster(paste("11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
threshold_value <- 1.1
example_model <- calc(example_model, fun=rc)

for(x in 1:length(species)){
  model <- raster(paste("11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
  threshold_value <- read.csv(paste("11_models/Bias_Spatial_Filtering/", species[[x]], "/maxentResults.csv", sep=""))
  threshold_value <- threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]
  
  sp.circle <- readShapePoly(paste("05_Species_Circles/", species[[x]], ".shp", sep=""))
  
  writeLines(paste("\nWorking on", species[[x]]))
  writeLines(paste("...The 10 percentile training presence logistic threshold  for", species[[x]], "is:", threshold_value))#
  
  ## Perform Thresholding
  thresholded_model <- calc(model, fun=rc)
  
  thresholded_model <- mask(thresholded_model, sp.circle)
  
  thresholded_model <- merge(thresholded_model, example_model)
  
  plot(thresholded_model, main = species[[x]])
  
  writeRaster(thresholded_model,  paste("12_Thresholded_Models/", species[[x]], ".tif", sep = ""))
}



rm(CBIs, model, thresholded_model, threshold_value, species, x)
