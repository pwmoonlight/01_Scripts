
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
##    (1) Correct names by taxonomy using the REFLORA R package                                                    ##
##                                                                                                                 ##
##    (2) Shift points that are just in the sea on to land                                                         ##
##                                                                                                                 ##
##    (3) Remove duplicate occurrance points                                                                       ##
##                                                                                                                 ##
##    (4) Remove points that are in states outside the species' known ranges as per REFLORA data                   ##
##                                                                                                                 ##
##    (5) Remove points whose altitude and altitude at coordinates differ by >250m                                 ##
##                                                                                                                 ##
##    (6) Remove points at known Brazil and State Centroids                                                        ##
##                                                                                                                 ##
##    (7) Plot species' distributions for visual check by taxonomic experts                                        ##
##                                                                                                                 ##
##    (8) Save cleaned results ready for modelling                                                                 ##
##                                                                                                                 ##
#####################################################################################################################



#########################################################################
####---------------------------- MODULE 1 ---------------------------####
#########################################################################
###### - Correct names by taxonomy using the REFLORA R package   - ######
#########################################################################


### Prepare the working space
### -------------------------

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))


### Find CRIA data ###
### -------------- ###

# CRIA distribution data should be in a folder called "02a_Input_Distribution_Data_DIRTY/CRIA" within the working directory
# This module currently works with excel 1997-2003 files as downloaded from CRIA

source(paste(getwd(), "/01_Scripts/01A____Find_New_CRIA_Data.R", sep=""))

# Get suggested accepted names for specimens
require(flora)
print(Sys.time())
distribution_data[,10] <-   as.character(lapply(distribution_data[,9], remove.authors))
distribution_data[,11:19] <- get.taxa(as.character(lapply(distribution_data[,10], suggest.names)), replace.synonyms = T, drop=c("id", "search.str", "threat.status", "notes", "name.status"))[,1:10]
print(Sys.time())
distribution_data[,19] <- paste(distribution_data[,5], distribution_data[,6])
distribution_data_test <- distribution_data[,-7]

which(is.na(distribution_data[,10]))

####################################################################
####------------------------- MODULE 2 -------------------------####
####################################################################
###### - Shift points that are just in the sea on to land   - ######
####################################################################



####################################################################
####------------------------- MODULE 3 -------------------------####
####################################################################
############## - Remove duplicate occurrance point - ###############
####################################################################



########################################################################################################
####------------------------------------------- MODULE 4 -------------------------------------------####
########################################################################################################
###### - Remove points that are in states outside the species' known ranges as per REFLORA data - ######
########################################################################################################



##########################################################################################
####------------------------------------ MODULE 5 ------------------------------------####
##########################################################################################
###### - Remove points whose altitude and altitude at coordinates differ by >250m - ######
##########################################################################################



#####################################################################
####-------------------------- MODULE 6 -------------------------####
#####################################################################
###### - Remove points at known Brazil and State Centroids   - ######
#####################################################################



###################################################################################
####--------------------------------- MODULE 7 --------------------------------####
###################################################################################
###### - Plot species' distributions for visual check by taxonomic experts - ######
###################################################################################



############################################################
####--------------------- MODULE 8 ---------------------####
############################################################
###### - Save cleaned results ready for modelling   - ######
############################################################
