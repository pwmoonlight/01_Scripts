###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

require(raster)

species_data <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]

writeLines("   ...Finding Cells At Points")
species_data[,4] <- cellFromXY(kde_raster, species_data[,2:1])

writeLines("   ...Finding Unique Cells")
species_data <- species_data[!duplicated(species_data[,4]),]
species_data <- species_data[!is.na(species_data[,4]),]

#writeLines("   ...Finding Values at Cells")
#species_data <- bg_data[as.character(species_data[,4]),]
#species_data <- species_data[!is.na(species_data[,4]),]
