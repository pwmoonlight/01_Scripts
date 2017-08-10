   #########################################################
  ###########################################################
 #############################################################
######## by Peter Moonlight, Tiina Sarkinen et al 2017 ########
 #############################################################
  ###########################################################
   #########################################################

### ---------------------
### Prepare the workspace
### ---------------------
  
# Set the working directory
setwd("C:/000_Modeling_Working_Directory_000")
   
# Clear the R workspace
rm(list = ls(all=T))

# Allocate larger than usual memory and 
memory.size(5000000)
options(java.parameters = "-Xmx12288m")



### -----------------------------
### Load in the distribution data
### -----------------------------

require("XLConnect")

# Distribution data should be in a folder called "02_Input Data" within the working directory
# The script currently works with excel 1997-2003 files as exported from PADME.
# This section will need altering for different input data formats.

distribution_data <- list.files('02_Input_Distribution_Data', pattern="[.]xls$", full.names=T)
distribution_data <- loadWorkbook(distribution_data)
distribution_data <- readWorksheet(distribution_data, sheet="sheet4", header=T)
distribution_data <- as.data.frame(distribution_data)

colnames(distribution_data) <- c("lat","long","CoordSource", "Specimen", "Raw_Species", "Species", "Altitude")

# Plot the input data on a crude world map for a visual check
plot(distribution_data$long, distribution_data$lat, col='orange', pch=20, cex=0.75)
points(distribution_data$long, distribution_data$lat, col='red', cex=0.75)
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, add=T)



### -----------------------------------------
### Separate the distribution data by species
### -----------------------------------------

# Find the species names
species <- as.vector(unique(distribution_data[,"Species"]))

# Create a data frame for each species
require("rgdal")
species.data.frames<-c()

#Separate the data by species
species_distribution_data <- lapply(1:length(species), function(x){distribution_data[distribution_data$Species == species[[x]],]})

#Save the data
dir.create("03_Input_Distribution_Data_by_Species")
lapply(species, function(x){dir.create(paste("03_Input_Distribution_Data_by_Species/", x, sep=""))})
lapply(1:length(species), function(x){write.csv(species_distribution_data[[x]], file=paste("03_Input_Distribution_Data_by_Species/", species[[x]], "/data.csv", sep=""))})



### -----------------------------------------------------------
### Find Which Species Have not Previously Been Modelled at all
### -----------------------------------------------------------

dir.create("99_Modelled_Species_Distribution_Data")
dir.create("04_Species_To_Model_Distribution_Data")

species <- list.dirs("03_Input_Distribution_Data_by_Species", full.names=F, recursive=F)

# List previously modelled species
previously_modelled_species <- list.dirs("99_Modelled_Species_Distribution_Data", full.names=F, recursive=F)
not_previously_modelled_species <- !(species %in% previously_modelled_species)

lapply(species[not_previously_modelled_species], function(x){dir.create(paste("04_Species_To_Model_Distribution_Data/", x, sep=""))})
lapply((1:length(species))[not_previously_modelled_species], function(x){write.csv(species_distribution_data[[x]], paste("04_Species_To_Model_Distribution_Data/", species[[x]], "/data.csv", sep=""))})



### -----------------------------------------------------------------------
### Find Which Species Have not Previously Been Modelled with the same data
### -----------------------------------------------------------------------

previously_modelled_species <- list.dirs("99_Modelled_Species_Distribution_Data", full.names=F, recursive=F)
previously_modelled_species_distribution_data <- lapply(1:length(previously_modelled_species), function(x){read.csv(paste("99_Modelled_Species_Distribution_Data/", previously_modelled_species[[x]], "/data.csv", sep=""))[,-1]})

species_to_crossref <- species[which(species %in% previously_modelled_species)]

for(x in 1:length(species_to_crossref)){
  old <- which(previously_modelled_species == species_to_crossref[[x]])
  old_data <- previously_modelled_species_distribution_data[[old]]

  new <- which(species == species_to_crossref[[x]])
  new_data <- read.csv(paste("03_Input_Distribution_Data_by_Species/", species[[new]], "/data.csv", sep=""))[,-1]

  if(!identical(old_data,new_data)){
    print(x)
    dir.create(paste("04_Species_To_Model_Distribution_Data/", species_to_crossref[[x]], sep=""))
    write.csv(species_distribution_data[[old]],paste("04_Species_To_Model_Distribution_Data/", species_to_crossref[[x]], "/data.csv", sep=""))
  }
}
