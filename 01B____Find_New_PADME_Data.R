###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################



### ----------------------------- ###
### Load in the distribution data ###
### ----------------------------- ###

require("XLConnect")

distribution_data <- list.files('02_Input_Distribution_Data/PADME', pattern="[.]xls$", full.names=T)
distribution_data <- loadWorkbook(distribution_data)
distribution_data <- readWorksheet(distribution_data, sheet="sheet4", header=T)
distribution_data <- as.data.frame(distribution_data)

colnames(distribution_data) <- c("lat","long","CoordSource", "Specimen", "Raw_Species", "Species", "Altitude")

distribution_data[,6] <- sub(" ", "_", distribution_data[,6])
for(x in 1:length(distribution_data[,6])){distribution_data[x,6] <- strsplit(distribution_data[x,6]," ")[[1]][1]}
distribution_data[,6] <- gsub("_", " ", distribution_data[,6])


### ----------------------------------------- ###
### Separate the distribution data by species ###
### ----------------------------------------- ###

# Find the species names
species <- as.vector(unique(distribution_data[,"Species"]))

# Create a data frame for each species
require("rgdal")

#Separate the data by species
species_distribution_data <- lapply(1:length(species), function(x){distribution_data[distribution_data$Species == species[[x]],]})

#Save the data
dir.create("03_Input_Distribution_Data_by_Species", showWarnings = F)
dir.create("03_Input_Distribution_Data_by_Species/PADME", showWarnings = F)
lapply(species, function(x){dir.create(paste("03_Input_Distribution_Data_by_Species/PADME/", x, sep=""), showWarnings = F)})
lapply(1:length(species), function(x){write.csv(species_distribution_data[[x]], file=paste("03_Input_Distribution_Data_by_Species/PADME/", species[[x]], "/data.csv", sep=""))})



### -----------------------------------------------------------
### Find Which Species Have not Previously Been Modelled at all
### -----------------------------------------------------------

dir.create("99_Modelled_Species_Distribution_Data", showWarnings = F)
dir.create("99_Modelled_Species_Distribution_Data/PADME", showWarnings = F)
dir.create("04_Species_To_Model_Distribution_Data", showWarnings = F)
dir.create("04_Species_To_Model_Distribution_Data/PADME", showWarnings = F)

# List previously modelled species
previously_modelled_species <- list.dirs("99_Modelled_Species_Distribution_Data/PADME", full.names=F, recursive=F)
not_previously_modelled_species <- !(species %in% previously_modelled_species)

lapply(species[not_previously_modelled_species], function(x){dir.create(paste("04_Species_To_Model_Distribution_Data/PADME/", x, sep=""), showWarnings = F)})
lapply((1:length(species))[not_previously_modelled_species], function(x){write.csv(species_distribution_data[[x]], paste("04_Species_To_Model_Distribution_Data/PADME/", species[[x]], "/data.csv", sep=""))})

writeLines("\n\n##############################################################")
writeLines("### The species never previously modelled with PADME data are:\n")
for(x in species[not_previously_modelled_species]){writeLines(x)}


### -----------------------------------------------------------------------
### Find Which Species Have not Previously Been Modelled with the same data
### -----------------------------------------------------------------------

previously_modelled_species <- list.dirs("99_Modelled_Species_Distribution_Data/PADME", full.names=F, recursive=F)
previously_modelled_species_distribution_data <- lapply(1:length(previously_modelled_species), function(x){read.csv(paste("99_Modelled_Species_Distribution_Data/PADME/", previously_modelled_species[[x]], "/data.csv", sep=""))[,-1]})

species_to_crossref <- species[which(species %in% previously_modelled_species)]
crossrefed_species <- c()
y <- 1

for(x in 1:length(species_to_crossref)){
  old <- which(previously_modelled_species == species_to_crossref[[x]])
  old_data <- previously_modelled_species_distribution_data[[old]]

  new <- which(species == species_to_crossref[[x]])
  new_data <- read.csv(paste("03_Input_Distribution_Data_by_Species/PADME/", species[[new]], "/data.csv", sep=""))[,-1]

  if(!identical(old_data,new_data)){
    dir.create(paste("04_Species_To_Model_Distribution_Data/PADME/", species_to_crossref[[x]], sep=""), showWarnings = F)
    write.csv(new_data,paste("04_Species_To_Model_Distribution_Data/PADME/", species_to_crossref[[x]], "/data.csv", sep=""))
    crossrefed_species[[y]] <- species_to_crossref[[x]]
    y <- y+1
  }
}

writeLines("\n\n##################################################################")
writeLines("### The species previously modelled with different PADME data are:\n")

for(x in crossrefed_species){writeLines(x)}

# Clear the R workspace
rm(distribution_data, new_data, old_data, crossrefed_species, new, not_previously_modelled_species, old, previously_modelled_species, previously_modelled_species_distribution_data, species_distribution_data, species_to_crossref, x, y)


