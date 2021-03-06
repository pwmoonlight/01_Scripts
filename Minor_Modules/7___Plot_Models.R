###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

require(raster)
require(maptools)
require(rgdal)

dir.create("03_Modelling/11_models/Plots", showWarnings=F)
dir.create("03_Modelling/11_models/Plots/Bias_Spatial_Filtering", showWarnings=F)

BRA_ADM <- list.files("000_GIS_LAYERS/BRA_Adm_2/", pattern="[.]shp$", full.names=T)
BRA_ADM <- readOGR(BRA_ADM[1])

CBIs <- read.csv("03_Modelling/CBI_results.csv", row.names = 1)
CBIs <- CBIs[which(CBIs[,6] >= 0.5),]
species <- rownames(CBIs)


for(x in 1:length(species)){
  
  if(!dir.exists(paste("03_Modelling/11_models/Plots/Bias_Spatial_Filtering/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/11_models/Plots/Bias_Spatial_Filtering/", species[[x]], sep=""))
    png(filename=paste("03_Modelling/11_models/Plots/Bias_Spatial_Filtering/", species[[x]], ".png", sep=""), width=8126, height=8126, units="px")
    
    writeLines(paste("   ...Working With ", species[[x]], sep=""))
    
    model <- raster(paste("03_Modelling/11_models/Consensus/Bias_Spatial_Filtering/", species[[x]], ".tif", sep=""))
    sp.circle <- readShapePoly(paste("03_Modelling/05_Species_Circles/", species[[x]], ".shp", sep=""))
    masked.model <- mask(model, sp.circle)
    species.data <- read.csv(paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
    
    plot(model, col=gray.colors(100, start = 0.7, end = 1, gamma = 2.2, alpha = NULL), legend=F)
    plot(masked.model, add=T, legend=F)
    plot(BRA_ADM, add=T)
    plot(SpatialPoints(species.data[,2:1]), add=T, cex=5, pch=17, col="blue")
    dev.off()
  }
}




rm(x, model, sp.circle, masked.model, species.data, species, species_2, BRA_ADM)
