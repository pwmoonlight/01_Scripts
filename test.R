setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))


require(dismo)
require(igraph)
require(rgeos)
require(rgdal)

nordeste <- readOGR("000_GIS_LAYERS/nordeste.shp", layer="nordeste")


species <- gsub(".tif$", "", list.files("03_Modelling/12_Thresholded_Models/", pattern="*.tif", full.names=F))
dir.create("03_Modelling/test", showWarnings = F)
for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/test/", species[[x]], sep=""))){
    writeLines(paste("Working on", species[[x]]))
    
    
    dir.create(paste("03_Modelling/test/", species[[x]], sep=""))
    
    model <- raster(paste("03_Modelling/12_Thresholded_Models/", species[[x]], ".tif", sep=""))
    model[which(model[1:length(model)] == 0)] <- NA

    data <- SpatialPoints(read.csv(paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,3:2])
    
    circle <- circles(data, d=400000, lonlat=TRUE)
    circle <- gUnaryUnion(circle@polygons)
    circle <- cellFromPolygon(model, circle)
    
    target <- which(model[1:length(model)] == 1)
    
    model[which(model[1:length(model)] == 1)[!which(model[1:length(model)] == 1) %in% circle[[1]]]] <- 2
    


#    adjacents <- as.data.frame(adjacent(model, index_twos, target=target, sorted=T))
    
#    carry_on <- 10000
#    while(carry_on != length(which(adjacents[,2] %in% index_ones))){
#      carry_on <- length(which(adjacents[,2] %in% index_ones))
#      model[adjacents[which(adjacents[,2] %in% index_ones),1]] <- 1
#      
#      index_ones <- c(index_ones, adjacents[which(adjacents[,2] %in% index_ones),1])
#    }

    
    carry_on <- 1
    
    index_twos <- which(model[1:length(model)] == 2)
    index_ones <- which(model[1:length(model)] == 1)
    
    while(carry_on > 0){
      circle <- xyFromCell(model, index_ones)
      temp <- model
      temp[which(1:length(temp) == 2)] <- NA
      circle <- circles(circle, d=100000, lonlat=TRUE)
      circle <- gUnaryUnion(circle@polygons)
      circle <- cellFromPolygon(model, circle)
      
      model[index_twos[which(index_twos %in% circle[[1]])]] <- 1
    }
    
        
    index_ones <- which(model[1:length(model)] == 1)
    circle <- xyFromCell(model, index_ones)
    circle <- circles(circle, d=125000, lonlat=TRUE)
    circle <- gUnaryUnion(circle@polygons)
    circle <- cellFromPolygon(model, circle)
    
    index_twos <- which(model[1:length(model)] == 2)
    model[index_twos[which(index_twos %in% circle[[1]])]] <- 1
    
    
    model <- crop(model, nordeste)
    model <- mask(model, nordeste)
    
    writeRaster(model, file=paste("03_Modelling/test/", species[[x]], ".tif", sep=""), overwrite=T)
  }
}

# Plot for a visual check
species <- gsub(".{4}$", "", list.files("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", pattern=".tif$", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/8___Plot_Thresholded_Models.R", sep=""))

