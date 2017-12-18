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

example_model <- model <- raster(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[1]], ".tif", sep=""))
threshold_value <- 1.1
example_model <- calc(example_model, fun=rc)

species_2 <- substr(list.files("03_Modelling/12_Thresholded_Models", full.names=F, recursive=F),1, nchar(list.files("03_Modelling/12_Thresholded_Models", full.names=F, recursive=F))-4)
species <- species[!(species %in% species_2)]

for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/12_Thresholded_Models/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/12_Thresholded_Models/", species[[x]], sep=""))
    
    model <- raster(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
    threshold_value <- read.csv(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], "/maxentResults.csv", sep=""))
    #if(length(threshold_value$Fixed.cumulative.value.10.area[6]) == 1){
    #  threshold_value <- threshold_value$Fixed.cumulative.value.10.area[6]
    #}
    #if(length(threshold_value$X10.percentile.training.presence.logistic.threshold[6]) == 1){
    #  threshold_value <- threshold_value$X10.percentile.training.presence.logistic.threshold[6]
    #}
    if(length(threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]) == 1){
      threshold_value <- threshold_value$X10.percentile.training.presence.Cloglog.threshold[6]
    }
    #if(length(threshold_value$Equal.training.sensitivity.and.specificity.Cloglog.threshold [6]) == 1){
    #  threshold_value <- threshold_value$Equal.training.sensitivity.and.specificity.Cloglog.threshold [6]
    #}
    
    writeLines(paste("\nWorking on", species[[x]]))
    writeLines(paste("...The 10 percentile training presence logistic threshold  for", species[[x]], "is:", threshold_value))#
    
    ## Perform Thresholding
    thresholded_model <- calc(model, fun=rc)
    
    #thresholded_model <- mask(thresholded_model, sp.circle)
    
    thresholded_model <- merge(thresholded_model, example_model)
    
    plot(thresholded_model, main = species[[x]])
    
    writeRaster(thresholded_model,  paste("03_Modelling/12_Thresholded_Models/", species[[x]], ".tif", sep = ""), overwrite=T)
  }
}



rm(CBIs, model, thresholded_model, threshold_value, species, x, species_2)
