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
species_2 <- substr(list.files("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering", full.names=F, recursive=F),1, nchar(list.files("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering", full.names=F, recursive=F))-4)
species <- species[!(species %in% species_2)]

require(raster)

for(x in length(species):1){
  
  if(!file.exists(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))){
    
    writeLines(paste("   ...Working With ", species[[x]], sep=""))
    
    model <- lapply(list.files(path="03_Modelling/11_models/Bias_Spatial_Filtering/", full.names = T, pattern=paste(species[[x]], "_", sep="")), raster)
    model <- stack(model)
  
    writeRaster(model, file=paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""), overwrite=T)
  }
}

rm(x, model)
