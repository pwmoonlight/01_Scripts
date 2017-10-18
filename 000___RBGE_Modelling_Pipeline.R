
###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################


#####################################################################################################################
##                                                                                                                 ##
## This script calls various other scripts to run the whole RBGE modelling pipeline                                ##
##                                                                                                                 ##
## Scripts can be run separately or together. They can be restarted from any point, and run in (almost!) any order ##
##                                                                                                                 ##
## The modules are:                                                                                                ##
##    (1) When a new file is inputted:                                                                             ##
##        a) Compare input data from CRIA, PADME or REFLORA with previously outputted data                         ##
##        b) Select input data which differs to run fresh models                                                   ##
##        c) Combine CRIA, PADME and REFLORA data ready to run the models                                          ##
##                                                                                                                 ##
##    (2) Remove species with too few points                                                                       ##
##        note: This also removes NA presence points                                                               ##
##                                                                                                                 ##
##    (3) Perform Spatial Filtering                                                                                ##
##                                                                                                                 ##
##    (4) Create a background sample for each species                                                              ##
##        note: Both a biased and a non-biased sample are created                                                  ##
##                                                                                                                 ##
##    (4) Run a PCA to select environmental variables                                                              ##
##                                                                                                                 ##
##    (5) Run four models with CHIRPs/MODIS data:                                                                  ##
##        a) No spatial filtering, no bias file                                                                    ##
##        b) Spatial filtering, no bias file                                                                       ##
##        c) No Spatial filtering, bias file                                                                       ##
##        d) Spatial filtering, bias file                                                                          ##
##                                                                                                                 ##
##    (6) Evaluate the CHIRPs/MODIS models                                                                         ##
##                                                                                                                 ##
##    (7) Run four models with bioclim data and project into the future:                                           ##
##        a) No spatial filtering, no bias file                                                                    ##
##        b) Spatial filtering, no bias file                                                                       ##
##        c) No Spatial filtering, bias file                                                                       ##
##        d) Spatial filtering, bias file                                                                          ##
##                                                                                                                 ##
##    (8) Predict Models into Future Climate Scenarios                                                             ##
##                                                                                                                 ##
##    (9) Evaluate the bioclim models                                                                              ##
##                                                                                                                 ##
##    (10) Mask to the Northeast of Brazil:                                                                        ##
##        a) The CHIRPS/MODIS models                                                                               ##
##        b) The CHIRPS/MODIS models - thresholded                                                                 ##
##        c) The bioclim models                                                                                    ##
##        d) The bioclim models - thresholded                                                                      ##
##        e) The future projected models                                                                           ##
##        f) The future projected models - thresholded                                                             ##
##                                                                                                                 ##
##    (11) Produce a Predicted Niche Occupance (PNO) profile for each species                                      ##
##                                                                                                                 ##
#####################################################################################################################


########################################
####----------- MODULE 1 -----------####
########################################
###### - Finds species to model - ######
########################################


### Prepare the working space
### -------------------------

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))


### Find any new or altered CRIA data ###
### --------------------------------- ###

# CRIA distribution data should be in a folder called "02_Input Data/CRIA" within the working directory
# This module currently works with excel 1997-2003 files as downloaded from CRIA

source(paste(getwd(), "/01_Scripts/01A____Find_New_CRIA_Data.R", sep=""))


### Find any new or altered PADME data ###
### ---------------------------------- ###

# PADME distribution data should be in a folder called "02_Input Data/PADME" within the working directory
# This module currently works with excel 1997-2003 files as exported from PADME

source(paste(getwd(), "/01_Scripts/01B____Find_New_PADME_Data.R", sep=""))


### Find any new or altered REFLORA data ###
### ------------------------------------ ###

# REFLORA distribution data should be in a folder called "02_Input Data/REFLORA" within the working directory
# This module currently works with excel 1997-2003 files as downloaded from REFLORA

source(paste(getwd(), "/01_Scripts/01C____Find_New_REFLORA_Data.R", sep=""))


### Combine any new or altered CRIA, PADME, and REFLORA data ###
### -------------------------------------------------------- ###

source(paste(getwd(), "/01_Scripts/01D____Combine_Data.R", sep=""))

################################################################
################################################################



#######################################
####----------- MODULE 2 ----------####
#######################################
######### - Clean the Data - ##########
#######################################


### Remove Data Points in states the species is not supposed to live in ###
### ------------------------------------------------------------------- ###

species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/9___Remove_Points_In_Wrong_States.R", sep=""))



### Remove Data Points in states they are not stated to be in ###
### --------------------------------------------------------- ###

species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/10___Remove_Points_In_Incorrect_States.R", sep=""))



### Remove Data Points with incorrect altitudes ###
### ------------------------------------------- ###

species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_2", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/11___Remove_Points_With_Incorrect_Altitudes.R", sep=""))


### Remove Data Points with outwith Brazil ###
### -------------------------------------- ###

species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_3", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/12___Remove_Points_Outwith_Brazil.R", sep=""))


### Remove Data Points with outwith Species' Range ###
### ---------------------------------------------- ###
species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/13___Remove_Points_Outwith_Species_Range.R", sep=""))


### Remove Data Points with outwith Species' PCA Range ###
### -------------------------------------------------- ###
species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/13___Remove_Points_by_PCA_space.R", sep=""))



require(dismo)
require(rgeos)
require(raster)

### Remove Species with Fewer than X Distribution Points ###
### ---------------------------------------------------- ###

# Specify how many distribution points are necessary
cutoff <- 5
kde_raster <- raster("000_GIS_LAYERS/Brazil_Masked_GIS_Layers/KDE_Raster/kde_raster.asc")
#kde_raster <- aggregate(kde_raster, 2)

species <- sub(".csv", "", list.files("03_Modelling/04_Species_To_Model_Distribution_Data", pattern=".csv", full.names=F, recursive=F))
dir.create("03_Modelling/08_Species_To_Model_Non_Scale_Corrected_Distribution_Data", showWarnings = F)
for(x in 1:length(species)){
  writeLines(paste("\nChecking ", species[[x]], " ..."))
  
  # Read in that species' data and find unique points
  source(paste(getwd(), "/01_Scripts/Minor_Modules/2___Find_Unique_Species_Distribution_Data.R", sep=""))

  species_data <- species_data[rownames(species_data)[which(rownames(species_data)!="NA")],]
  
  # Remove species with fewer than 5 unique points
  source(paste(getwd(), "/01_Scripts/02A____Remove_Species_With_Too_Few_Points.R", sep=""))
  
  # Save a "model ready" csv file for species with enough data points
  if(length(species_data[,1])>cutoff){
    writeLines("   ...Sufficient Presence Points")
    unlink(paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]],".csv", sep=""))
    write.csv(species_data, file=paste("03_Modelling/08_Species_To_Model_Non_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))
  }
}

##########################################################################
##########################################################################



##########################################################################
####---------------------------- MODULE 3 ----------------------------####
##########################################################################
######### - Perform Spatial Filtering a.k.a. Scale Correction - ##########
##########################################################################



### Create a "Scale Corrected" dataset for each species ###
### --------------------------------------------------- ###
# This is also known as "Spatial Filtering"

# Set the distance over which to preform Spatial Filtering in km
scale_distance <- 10

species <- sub(".csv", "", list.files("03_Modelling/08_Species_To_Model_Non_Scale_Corrected_Distribution_Data", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/03____Spatial_Filtering.R", sep=""))



### Produce a biased and a non-biased background sample for each species ###
### -------------------------------------------------------------------- ###

require(raster)
kde_raster <- raster("000_GIS_LAYERS/Brazil_Masked_GIS_Layers/KDE_Raster/kde_raster.asc")

# Produce a background sample for each species
# This involves producing both a non-biased and a biased background sample
species <- sub(".csv", "", list.files("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data", pattern=".csv", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/04____Background_Samples.R", sep=""))
#source(paste(getwd(), "/01_Scripts/04____Background_Samples_backwards.R", sep=""))

########################################################################
########################################################################



########################################################################
####--------------------------- MODULE 4 ---------------------------####
########################################################################
######### - Perform a PCA to select environmental variables - ##########
########################################################################

bg <- lapply(list.files(path="000_GIS_LAYERS/Brazil_Masked_GIS_Layers", pattern="*.tif$", full.names = T), raster)
bg <- stack(bg)
bg <- aggregate(bg, 2)

source(paste(getwd(), "/01_Scripts/04____PCA.R", sep=""))

# This line will need manual editing depending upon which variables are selected
PCA <- c(1, 6, 12, 13, 14, 15, 21, 26, 27, 28, 30, 31, 33)



#####################################################
####----------------- MODULE 5 ------------------####
#####################################################
######### - Runs the CHIRPS/MODIS models - ##########
#####################################################

# These is the A* models, using MODIS/CHIRPs data
# The bioclim models, including projection into the future are produced below

# Load in required packages
require(dismo)
require(rgeos)
require(raster)

### Run 4 models ###
### ------------ ###

dir.create("03_Modelling/11_models", showWarnings = F)
require(dismo)

bg <- lapply(list.files(path="000_GIS_LAYERS/Brazil_Masked_GIS_Layers", pattern="*.tif$", full.names = T), raster)
# Usually you will have done a PCA by this point. Keep only those BG layers selected during the PCA  
bg <- bg[PCA]
bg <- stack(bg)
#bg <- aggregate(bg, 2)


species <- sub(".csv", "", list.files("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data", pattern=".csv", full.names=F, recursive=F))

for(x in 1:length(species)){
  writeLines(paste("\nWorking on ", species[[x]], " ..."))
  
  species_data <- read.csv(paste("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
  background_data <- read.csv(paste("03_Modelling/10_Background_Data/Biased/", species[[x]], ".csv", sep=""))[,-1]
  
  source(paste(getwd(), "/01_Scripts/05d____Run_Model_Bias_Spatial_Filtering.R", sep=""))
  
}


rm(background_data, species_data, bg, kde_raster, model, PCA, x)



#####################################################
####----------------- MODULE 6 ------------------####
#####################################################
####### - Evaluates the CHIRPS/MODIS models - #######
#####################################################



### Finds the model CBI ###
### ------------------- ### 

species <- sub(".csv", "", list.files("03_Modelling/09_Species_To_Model_Scale_Corrected_Distribution_Data", pattern=".csv", full.names=F, recursive=F))

source(paste(getwd(), "/01_Scripts/07___Find_Model_CBIs.R", sep=""))


### Finds the model mean and creates a consensus models ###
### --------------------------------------------------- ### 

species <- sub(".csv", "", list.dirs("03_Modelling/11_models/Bias_Spatial_Filtering", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/6___Find_Consensus_Model.R", sep=""))


### Plots the model for visual check ###
### -------------------------------- ### 


source(paste(getwd(), "/01_Scripts/Minor_Modules/7___Plot_Models.R", sep=""))


### Converts the Models to Binary (thresholding) ###
### -------------------------------------------- ###

# This is done with a 10% training threshold

source(paste(getwd(), "/01_Scripts/08___Thresholds_Models.R", sep=""))
species <- gsub(".{4}$", "", list.files("03_Modelling/12_Thresholded_Models", pattern=".tif", full.names=F, recursive=F))
source(paste(getwd(), "/01_Scripts/Minor_Modules/8___Plot_Thresholded_Models.R", sep=""))


################################################################################
####------------------------------ MODULE 5a -------------------------------####
################################################################################
######### - Perform a PCA to select bioclim environmental variables - ##########
################################################################################


source(paste(getwd(), "/01_Scripts/04____PCA_bioclim.R", sep=""))

# This line will need manual editing depending upon which variables are selected
PCA_bioclim <- c(4,9,10,11,13,17,18)



################################################
####--------------- MODULE 7 ---------------####
################################################
######### - Runs the bioclim models - ##########
################################################


dir.create("14_Models_Bioclim", showWarnings = F)
require(dismo)


### Read in the bioclim data ###
### ------------------------ ###

# The folder may need changing depending upon where you are running the model

#bg_bioclim <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM", pattern="*3_degrees.tif$", full.names = T), raster)
bg_bioclim <- lapply(list.files(path="000_GIS_LAYERS/Brazil Masked BIOCLIM", pattern="*3_degrees.tif$", full.names = T), raster)

# Usually you will have done a PCA by this point. Keep only those BG layers selected during the PCA  
bg_bioclim <- bg_bioclim[PCA_bioclim]
bg_bioclim <- stack(bg_bioclim)

# Prepare the future climate data
source(paste(getwd(), "/01_Scripts/Minor_Modules/4___Prepare_Future_Climate_Data.R", sep=""))

species <- sub(".csv", "", list.files("09_Species_To_Model_Scale_Corrected_Distribution_Data", pattern=".csv", full.names=F, recursive=F))

for(x in 1:length(species)){
  writeLines(paste("\nWorking on ", species[[x]], " ..."))
  
  species_data <- read.csv(paste("09_Species_To_Model_Scale_Corrected_Distribution_Data/", species[[x]], ".csv", sep=""))[,-1]
  background_data <- read.csv(paste("10_Background_Data/Biased/", species[[x]], ".csv", sep=""))[,-1]
  
  source(paste(getwd(), "/01_Scripts/11d____Run_Model_Bias_Spatial_Filtering_BIOCLIM.R", sep=""))
  
  ### Project Those Models onto Future Scenarios ###
  ### ------------------------------------------ ### 
  
  # This is a list of models chosen following via K-clustering. This is currently done in a separate script.
  chosen_models_45 <- c(1, 6, 8, 9, 12, 13)
  chosen_models_85 <- c(6, 8, 11, 12, 13)
  
  source(paste(getwd(), "/01_Scripts/06____Predict_Model_Into_Future_Scenarios.R", sep=""))
  
  rm(chosen_models_45, chosen_models_85)
}
rm(background_data, species_data, bg_bioclim, kde_raster, model, PCA, x, chosen_models_45, chosen_models_85)


##############################################
####-------------- MODULE 10 -------------####
##############################################
####### - Mask Models to the Nordeste - ######
##############################################



nordeste <- readOGR("000_GIS_layers/nordeste.shp", layer="nordeste")

### Mask the CHIRPS/MODIS models ###
### ---------------------------- ###

species <- gsub(".{4}$", "", list.files("12_Thresholded_models", pattern="*.tif", full.names=F))
dir.create("11a_Models_Masked_Nordeste", showWarnings = F)
for(x in 1:length(species)){
  model <- mask(raster(paste("11_models/Bias_Spatial_Filtering/", species[[x]], ".tif", sep="")), nordeste)
  writeRaster(model, file=paste("11a_Models_Masked_Nordeste/", species[[x]], ".tif", sep=""))
}

### Mask the CHIRPS/MODIS models - with thresholding ###
### ------------------------------------------------ ###

species <- gsub(".{4}$", "", list.files("12_Thresholded_models", pattern="*.tif", full.names=F))
dir.create("12a_Thresholded_Models_Masked_Nordeste", showWarnings = F)
for(x in 1:length(species)){
  model <- mask(raster(paste("12_Thresholded_Models/", species[[x]], ".tif", sep="")), nordeste)
  writeRaster(model, file=paste("12a_Thresholded_Models_Masked_Nordeste/", species[[x]], ".tif", sep=""))
}


### Mask the bioclim models ###
### ----------------------- ###

### Mask the bioclim models - with thresholding ###
### ------------------------------------------- ###

### Mask the future bioclim models ###
### ------------------------------ ###

### Mask the future bioclim models - with thresholding ###
### -------------------------------------------------- ###






########################################
####----------- MODULE 4 -----------####
########################################
####### - Produces PNO profiles - ######
########################################
