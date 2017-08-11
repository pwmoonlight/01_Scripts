dir.create("13_Models_Bioclim/Bias_Spatial_Filtering", showWarnings=F)

writeLines(paste("...Producing Model"))
model <- maxent(bg_bioclim,p=species_data[,2:1],a=background_data[,3:4], path=paste("13_models_Bioclim/Bias_Spatial_Filtering/", species[[x]], sep=""), args="replicates=5")

writeLines(paste("...Predicting Model Using Current Climate Variables"))
predict(model, bg_bioclim, filename=paste("13_Models_Bioclim/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))

rm(model)
