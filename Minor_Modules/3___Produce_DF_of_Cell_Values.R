###############################################################################################################
  ##########################################################################################lorena villanueva-almanza#####################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

require(dismo)
require(rgeos)
require(maptools)

writeLines(paste("\nProducing Circles Around ALL Data...", sep=""))
all_circles <- lapply(seq_along(1:length(species)), function(x){read.csv(paste("04_Species_To_Model_Distribution_Data/", species[[x]], ".csv", sep=""))})
all_circles <- do.call(rbind, all_circles)
all_circles <- circles(crop(SpatialPoints(all_circles[,3:2]), extent(kde_raster)), d=80000, lonlat=TRUE)
all_circles <- gUnaryUnion(all_circles@polygons)
dir.create("05_Species_Circles", showWarnings=F)
writePolyShape(as(all_circles, "SpatialPolygonsDataFrame"), paste("05_Species_Circles/all_circles.shp", sep=""))

writeLines(paste("Finding Cells Within Circles...", sep=""))
all_circles <- cellFromPolygon(kde_raster, all_circles)

bg_data <- rasterToPoints(kde_raster)
all_cells <- cellFromXY(kde_raster, bg_data[,1:2])


all_circles <- all_circles[[1]][all_circles[[1]] %in% all_cells]


rownames(bg_data) <- all_cells

bg_data <- bg_data[as.character(all_circles),]

writeLines(paste("Producing DF of Cell Values...", sep=""))
for(x in 1:length(bg@layers)){
  print(x)
  temp <- rasterToPoints(bg[[x]])
  print(x)
  rownames(temp) <- all_cells
  print(x)
  temp <- temp[as.character(all_circles),]
  print(x)
  bg_data <- cbind(bg_data, temp[,3])
}
colnames(bg_data) <- append(c("x", "y", "kde_raster"), names(bg))





dir.create("06_BG_data", showWarnings=F)

writeLines(paste("Saving DF of Cell Values...", sep=""))
write.csv(bg_data, file="06_BG_data/BG_data.csv")

rm(bg, temp, all_cells, all_circles)
