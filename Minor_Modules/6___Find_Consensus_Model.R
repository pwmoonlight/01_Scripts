###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("03_Modelling/11_models/Consensus", showWarnings=F)
dir.create("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", showWarnings=F)

species <- sub(".csv", "", list.dirs("03_Modelling/11_models/Bias_Spatial_Filtering", full.names=F, recursive=F))

require(raster)

for(x in 1:length(species)){
  
  if(!dir.exists(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], sep=""), showWarnings = F)
    
    writeLines(paste("   ...Working With ", species[[x]], sep=""))
    
    model <- lapply(list.files(path="03_Modelling/11_models/Bias_Spatial_Filtering/", full.names = T, pattern=paste(species[[x]], "_", sep="")), raster)
    model <- stack(model)
  
    writeRaster(model, file=paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""), overwrite=T)
  }
}

rm(x, model)
