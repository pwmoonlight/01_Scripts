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

distribution_data <- list.files('03_Modelling/02_Input_Distribution_Data/REFLORA', pattern="[.]xlsx", full.names=T)
distribution_data <- loadWorkbook(distribution_data)
distribution_data <- readWorksheet(distribution_data, sheet="Relatório", header=T)
distribution_data <- as.data.frame(distribution_data)

distribution_data <- distribution_data[,c(30,32,17,19,10,11,27)]

# Remove specimens with no latitude or longitude
for(y in length(distribution_data[,6]):1){if(is.na(distribution_data[y,1])){distribution_data <- distribution_data[-y,]}}
for(y in length(distribution_data[,6]):1){if(distribution_data[y,1] ==""){distribution_data <- distribution_data[-y,]}}
for(y in length(distribution_data[,6]):1){if(is.na(distribution_data[y,2])){distribution_data <- distribution_data[-y,]}}
for(y in length(distribution_data[,6]):1){if(distribution_data[y,2] ==""){distribution_data <- distribution_data[-y,]}}


# Remove specimens not determined to species
for(y in length(distribution_data[,6]):1){if(is.na(distribution_data[y,6])){distribution_data <- distribution_data[-y,]}}
for(y in length(distribution_data[,6]):1){if(distribution_data[y,6] == ""){distribution_data <- distribution_data[-y,]}}
for(y in length(distribution_data[,6]):1){if(distribution_data[y,6] == "sp."){distribution_data <- distribution_data[-y,]}}

distribution_data[,8] <- paste(distribution_data[,5], distribution_data[,6])
colnames(distribution_data) <- c("lat","long","collector", "collectornumber", "genus", "species", "altitude", "Scientific_name")

# Convert DMS to DD
distribution_data[,1] <- sub("º", "", distribution_data[,1])
distribution_data[,2] <- sub("º", "", distribution_data[,2])
distribution_data[,1] <- sub("'", "", distribution_data[,1])
distribution_data[,2] <- sub("'", "", distribution_data[,2])
distribution_data[,1] <- sub("\"", "", distribution_data[,1])
distribution_data[,2] <- sub("\"", "", distribution_data[,2])

for(y in 1:length(distribution_data[,1])){
  split <- strsplit(distribution_data[y,1], " ")
  if(split[[1]][4] == "N"){
    distribution_data[y,1] <- as.numeric(as.numeric(split[[1]][1]) + (as.numeric(split[[1]][2]) + as.numeric(split[[1]][3])/60)/60)
  }
  if(split[[1]][4] == "S"){
    distribution_data[y,1] <- 0-as.numeric(as.numeric(split[[1]][1]) + (as.numeric(split[[1]][2]) + as.numeric(split[[1]][3])/60)/60)
  }
  
  split <- strsplit(distribution_data[y,2], " ")
  if(split[[1]][4] == "E"){
    distribution_data[y,2] <- as.numeric(as.numeric(split[[1]][1]) + (as.numeric(split[[1]][2]) + as.numeric(split[[1]][3])/60)/60)
  }
  if(split[[1]][4] == "W"){
    distribution_data[y,2] <- 0-as.numeric(as.numeric(split[[1]][1]) + (as.numeric(split[[1]][2]) + as.numeric(split[[1]][3])/60)/60)
  }
}

#write.csv(distribution_data, file="test.csv")

### ----------------------------------------- ###
### Separate the distribution data by species ###
### ----------------------------------------- ###

# Find the species names
species <- as.vector(unique(distribution_data[,"Scientific_name"]))

# Create a data frame for each species
require("rgdal")
species.data.frames<-c()

#Separate the data by species
species_distribution_data <- lapply(1:length(species), function(y){distribution_data[distribution_data$Scientific_name == species[[y]],]})

#Save the data
dir.create("03_Modelling/03_Input_Distribution_Data_by_Species", showWarnings = F)
dir.create("03_Modelling/03_Input_Distribution_Data_by_Species/REFLORA", showWarnings = F)
lapply(species, function(y){dir.create(paste("03_Modelling/03_Input_Distribution_Data_by_Species/REFLORA/", y, sep=""), showWarnings = F)})
lapply(1:length(species), function(y){write.csv(species_distribution_data[[y]], file=paste("03_Modelling/03_Input_Distribution_Data_by_Species/REFLORA/", species[[y]], "/data.csv", sep=""))})



### -----------------------------------------------------------
### Find Which Species Have not Previously Been Modelled at all
### -----------------------------------------------------------

dir.create("03_Modelling/99_Modelled_Species_Distribution_Data", showWarnings = F)
dir.create("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA", showWarnings = F)
dir.create("03_Modelling/04_Species_To_Model_Distribution_Data", showWarnings = F)
dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA", showWarnings = F)

#species <- list.dirs("03_Input_Distribution_Data_by_Species/REFLORA", full.names=F, recursive=F)

# List previously modelled species
previously_modelled_species <- list.dirs("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA", full.names=F, recursive=F)
not_previously_modelled_species <- !(species %in% previously_modelled_species)

lapply(species[not_previously_modelled_species], function(y){dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", y, sep=""), showWarnings = F)})
lapply((1:length(species)), function(y){
  if(not_previously_modelled_species[[y]]){
    write.csv(species_distribution_data[[y]], paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", species[y], "/data.csv", sep=""))
  }})

writeLines("\n\n##############################################################")
writeLines("### The species never previously modelled with REFLORA data are:\n")
for(y in previously_modelled_species){writeLines(y)}



### -----------------------------------------------------------------------
### Find Which Species Have not Previously Been Modelled with the same data
### -----------------------------------------------------------------------

previously_modelled_species <- list.dirs("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA", full.names=F, recursive=F)
previously_modelled_species_distribution_data <- lapply(1:length(previously_modelled_species), function(y){read.csv(paste("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA/", previously_modelled_species[[y]], "/data.csv", sep=""))[,-1]})

species_to_crossref <- species[which(species %in% previously_modelled_species)]
crossrefed_species <- c()
z <- 1

for(y in 1:length(species_to_crossref)){
  old <- which(previously_modelled_species == species_to_crossref[[y]])
  old_data <- previously_modelled_species_distribution_data[[old]]

  new <- which(species == species_to_crossref[[y]])
  new_data <- read.csv(paste("03_Modelling/03_Input_Distribution_Data_by_Species/REFLORA/", species[[new]], "/data.csv", sep=""))[,-1]

  if(!identical(old_data,new_data)){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", species_to_crossref[[y]], sep=""), showWarnings = F)
    write.csv(new_data,paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", species_to_crossref[[y]], "/data.csv", sep=""))
    crossrefed_species[[z]] <- species_to_crossref[[y]]
    z <- z+1
  }
}

writeLines("\n\n##################################################################")
writeLines("### The species previously modelled with different REFLORA data are:\n")

for(y in crossrefed_species){writeLines(y)}

# Clear the R workspace
rm(distribution_data, new_data, old_data, crossrefed_species, new, not_previously_modelled_species, old, previously_modelled_species, previously_modelled_species_distribution_data, species_distribution_data, species_to_crossref, z, y)

