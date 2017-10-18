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

#require("XLConnect")

distribution_data <- list.files('02_Input_Distribution_Data/CRIA', pattern="[.]txt$", full.names=T)
#distribution_data <- loadWorkbook(distribution_data, create=T)
#distribution_data <- readWorksheet(distribution_data, sheet="sheet1", header=T)
#distribution_data <- as.data.frame(distribution_data)

distribution_data <- read.table(distribution_data[2], header=T, sep="\t", fill=T, quote="")

distribution_data <- distribution_data[,c(25,26,17,18,1,2,3,23,24,19)]

# Remove specimens not determined to species
#for(x in length(distribution_data[,7]):1){
#  if(is.na(distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#  if(distribution_data[x,7] ==""){distribution_data <- distribution_data[-x,]}
#  if(distribution_data[x,7] =="sp"){distribution_data <- distribution_data[-x,]}
#  if(distribution_data[x,7] =="sp."){distribution_data <- distribution_data[-x,]}
#  
#  if(grepl("cf.", distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#  if(grepl("aff.", distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#  if(grepl("subsp.", distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#  if(grepl("var.", distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#  if(grepl(" ", distribution_data[x,7])){distribution_data <- distribution_data[-x,]}
#}

distribution_data[,11] <- paste(distribution_data[,6], distribution_data[,7])
colnames(distribution_data) <- c("lat","long","collector", "collectornumber","family" ,"genus", "species", "min_alt", "max_alt", "state", "Scientific_name")


### ----------------------------------------- ###
### Separate the distribution data by species ###
### ----------------------------------------- ###

# Find the species names
species <- as.vector(unique(distribution_data[,"Scientific_name"]))

# Find Nordeste Species
Reflora_Species <- read.csv("000_Species/Angiospermas_Flora do Brasil 2020_25 set 2017_Working.csv")
Reflora_Species <- Reflora_Species[which(Reflora_Species[,47]==1),9]
Reflora_Species <- as.vector(Reflora_Species)

species <- species[species %in% Reflora_Species]


# Create a data frame for each species
require("rgdal")

#Separate the data by species
species_distribution_data <- lapply(1:length(species), function(x){distribution_data[distribution_data$Scientific_name == species[[x]],]})

#Save the data
dir.create("03_Modelling/03_Input_Distribution_Data_by_Species", showWarnings = F)
dir.create("03_Modelling/03_Input_Distribution_Data_by_Species/CRIA", showWarnings = F)
for(y in 1:length(species)){dir.create(paste("03_Modelling/03_Input_Distribution_Data_by_Species/CRIA/", species[[y]], sep=""), showWarnings = F)}
for(y in 1:length(species)){write.csv(species_distribution_data[[y]], file=paste("03_Modelling/03_Input_Distribution_Data_by_Species/CRIA/", species[[y]], "/data.csv", sep=""))}

### -----------------------------------------------------------
### Find Which Species Have not Previously Been Modelled at all
### -----------------------------------------------------------

dir.create("03_Modelling/99_Modelled_Species_Distribution_Data", showWarnings = F)
dir.create("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA", showWarnings = F)
dir.create("03_Modelling/04_Species_To_Model_Distribution_Data", showWarnings = F)
dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned", showWarnings = F)
dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA", showWarnings = F)

#species <- list.dirs("03_Input_Distribution_Data_by_Species/CRIA", full.names=F, recursive=F)

# List previously modelled species
previously_modelled_species <- list.dirs("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA", full.names=F, recursive=F)
not_previously_modelled_species <- !(species %in% previously_modelled_species)

lapply(species[not_previously_modelled_species], function(y){dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", y, sep=""), showWarnings = F)})
for(y in 1:length(species)){
  if(not_previously_modelled_species[y])
  {
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", species[[y]], sep=""))
    write.csv(species_distribution_data[[y]], paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", species[[y]], "/data.csv", sep=""))}
  }

writeLines("\n\n##############################################################")
writeLines("### The species never previously modelled with CRIA data are:\n")
for(x in species[not_previously_modelled_species]){writeLines(x)}


### -----------------------------------------------------------------------
### Find Which Species Have not Previously Been Modelled with the same data
### -----------------------------------------------------------------------

previously_modelled_species <- list.dirs("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA", full.names=F, recursive=F)

if(length(previously_modelled_species)>0){
  previously_modelled_species_distribution_data <- lapply(1:length(previously_modelled_species), function(x){read.csv(paste("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA/", previously_modelled_species[[x]], "/data.csv", sep=""))[,-1]})
}

species_to_crossref <- species[which(species %in% previously_modelled_species)]
crossrefed_species <- c()
y <- 1

for(x in 1:length(species_to_crossref)){
  old <- which(previously_modelled_species == species_to_crossref[[x]])
  old_data <- previously_modelled_species_distribution_data[[old]]

  new <- which(species == species_to_crossref[[x]])
  new_data <- read.csv(paste("03_Modelling/03_Input_Distribution_Data_by_Species/CRIA/", species[[new]], "/data.csv", sep=""))[,-1]

  if(!identical(old_data,new_data)){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", species_to_crossref[[x]], sep=""), showWarnings = F)
    write.csv(new_data,paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", species_to_crossref[[x]], "/data.csv", sep=""))
    crossrefed_species[[y]] <- species_to_crossref[[x]]
    y <- y+1
  }
}

writeLines("\n\n##################################################################")
writeLines("### The species previously modelled with different CRIA data are:\n")

for(x in crossrefed_species){writeLines(x)}


# Clear the R workspace
rm(distribution_data, new_data, old_data, crossrefed_species, new, not_previously_modelled_species, old, previously_modelled_species, previously_modelled_species_distribution_data, species_distribution_data, species_to_crossref, x, y)
