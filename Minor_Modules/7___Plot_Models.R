###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("/000_Modeling_Working_Directory_000/11_models/Plots", showWarnings=F)

for(x in 1:length(species)){
  writeLines(paste("   ...Working With ", species[[x]], sep=""))
  
  model <- raster(paste("03_Modelling/11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
  sp.circle <- readShapePoly(paste("03_Modelling/05_Species_Circles/", species[[x]], ".shp", sep=""))
  masked.model <- mask(model, sp.circle)
  species.data <- read.csv(paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
  
  png(filename=paste("03_Modelling/11_models/Plots/", species[[x]], ".png"), width=8126, height=8126, units="px")
    plot(model, col=gray.colors(100, start = 0.7, end = 1, gamma = 2.2, alpha = NULL), legend=F)
    plot(masked.model, add=T, legend=F)
    plot(SpatialPoints(species.data[,2:1]), add=T, cex=5, pch=17, col="blue")
  dev.off()
}




rm(x, model, sp.circle, masked.model, species.data)
