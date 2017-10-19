###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

# Load in Spatial Filtering Function

filterByProximity <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy,longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    discard <- as.matrix(discard)
    return(xy[-discard[,1],])
    #return <- discard
  }
  if (nrow(closePts) == 0) {
    return(xy)
    #return <- discard
  }
  return(return)
}

# Load in Data Frames
species_data_frames <- lapply(1:length(species), function(x){read.csv(paste("03_Modelling/08_Species_To_Model_Non_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]})
species_data_frames <- lapply(1:length(species), function(x){species_data_frames[[x]][!is.na(species_data_frames[[x]][,1]),]})


# Run Spatial Filtering
scale_corrected_species_data_frames <- lapply(species_data_frames, function(x){filterByProximity(as.matrix(x[,1:2]), scale_distance, mapUnits=F)})

# Count the number of cells left per species
scale_corrected_cell_number <- lapply(scale_corrected_species_data_frames, function(x){dim(as.data.frame(x))[1]})

writeLines(paste("\nThe following species have fewer than five points following Spatial Filtering:\n"))
dir.create("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", showWarnings = FALSE)
for(x in 1:length(species)){
  if(scale_corrected_cell_number[[x]]>(cutoff-1)){
    scale_corrected_species_data_frames[[x]] <- as.data.frame(scale_corrected_species_data_frames[[x]])
    scale_corrected_species_data_frames[[x]][,3] <- cellFromXY(kde_raster, scale_corrected_species_data_frames[[x]][,2:1])
    scale_corrected_species_data_frames[[x]] <- scale_corrected_species_data_frames[[x]][!duplicated(scale_corrected_species_data_frames[[x]][,3]),]
    
    #scale_corrected_species_data_frames[[x]] <- bg_data[as.character(scale_corrected_species_data_frames[[x]][,3]),]
    
    
    write.csv(scale_corrected_species_data_frames[[x]], paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))
  }
  if(scale_corrected_cell_number[[x]]<cutoff){
    writeLines(paste(species[[x]]))
  }
}

rm(scale_corrected_cell_number, scale_corrected_species_data_frames, species_data_frames, scale_distance, cutoff, x, species)
