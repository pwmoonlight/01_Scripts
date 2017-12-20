dir.create("03_Modelling/11_Models/Bias_Spatial_Filtering", showWarnings=F)

writeLines(paste("...Producing Model"))
model <- maxent(bg,p=species_data[,2:1],a=background_data[,3:4], path=paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], sep=""), args="replicates=5")

writeLines(paste("...Predicting Model Using Current Climate Variables"))
predict(model, bg, filename=paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""), overwrite=T)

rm(model)