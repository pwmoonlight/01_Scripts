
###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
###############################################################################################################

require(maptools)
require(rgdal)
require(dismo)
require(rgeos)
require(raster)

dir.create("03_Modelling/05_Species_Circles", showWarnings=F)
dir.create("03_Modelling/10_Background_Data", showWarnings=F)
dir.create("03_Modelling/10_Background_Data/Biased", showWarnings=F)
dir.create("03_Modelling/10_Background_Data/Non_Biased", showWarnings=F)

nordeste <- readOGR("000_GIS_LAYERS/nordeste.shp", layer="nordeste")

for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/10_Background_Data/Biased/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/10_Background_Data/Biased/", species[[x]], sep=""), showWarnings = F)
    
  
    writeLines(paste("\nFinding Background data for ", species[[x]], ":", sep=""))
  
    circle <- read.csv(paste("03_Modelling/08_Species_To_Model_Non_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
    
    
    writeLines("...Producing Circle")
    circle <- circles(circle[,2:1], d=400000, lonlat=TRUE)
    #circle <- gUnaryUnion(circle@polygons)
    circle <- gUnion(circle@polygons, SpatialPolygons(nordeste@polygons))
    
    writePolyShape(as(circle, "SpatialPolygonsDataFrame"), paste("03_Modelling/05_Species_Circles/", species[[x]], ".shp", sep=""))
    
    writeLines("...Finding Cells in Circle")
    circle <- cellFromPolygon(kde_raster, circle)
  
    writeLines("...Finding Coordinates of Cells")
    circle <- xyFromCell(kde_raster, circle[[1]])
  
    writeLines("...Finding Values at Coordinates")
    background <- extract(kde_raster, circle, na.rm=T, cellnumbers=T)
  
    background <- cbind(background, circle)
  
    writeLines("...Removing NA Values")
    background <- background[which(!is.na(background[,2])),]
  
  
  
    writeLines("...Producing Background Samples")
    background_biased <- background[sample(1:length(background[,1]), prob=background[,2], replace=TRUE, size=1000),]
    background_not_biased <- background[sample(1:length(background[,1]), replace=TRUE, size=1000),]
  
    writeLines("...Saving Background Samples")
    write.csv(background_biased, file=paste("03_Modelling/10_Background_Data/Biased/", species[[x]], ".csv", sep=""))
    write.csv(background_not_biased, file=paste("03_Modelling/10_Background_Data/Non_Biased/", species[[x]], ".csv", sep=""))
  }
}


rm(bg, background, background_biased, background_not_biased, circle, x, nordeste)






#background <- cbind(background, extract(kde_raster, circle, na.rm=T))

#writeLines("...Removing NA Values")
#background <- background[which(!is.na(background[,2])),]

#colnames(background)[length(background[1,])] <- "kde_raster"




