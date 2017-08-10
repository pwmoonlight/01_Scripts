dir.create("11_Models/Bias_Spatial_Filtering", showWarnings=F)

writeLines(paste("...Parititioning Presence and Absence Data"))

writeLines(paste("...Producing Model"))
model <- maxent(bg,p=species_data[,2:1],a=background_data[,3:4], path=paste("11_models/Bias_Spatial_Filtering/", species[[x]], sep=""), args="replicates=5")

writeLines(paste("...Predicting Model Using Current Climate Variables"))
predict(model, bg, filename=paste("11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))

rm(model)