###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################


dir.create("03_Modelling/07_Species_With_Too_Little_Data_To_Model", showWarnings = F)

CRIA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/CRIA", full.names=F, recursive=F)
PADME <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/PADME", full.names=F, recursive=F)
REFLORA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA", full.names=F, recursive=F)

if(!length(species_data[,1])>cutoff){

  writeLines("   ...Too Few Presence Points\n   ...Removing From Analysis")
  write.csv(species_data, paste("03_Modelling/07_Species_With_Too_Little_Data_To_Model/", species[[x]], ".csv", sep=""))

  c <- which(CRIA == species[[x]])
  p <- which(PADME == species[[x]])
  r <- which(REFLORA == species[[x]])
  
  if(length(c)>0){
    c <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/CRIA/", species[[x]], "/data.csv", sep=""))[,-1]
    dir.create(paste("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA/", species[[x]], sep=""), showWarnings = F)
    write.csv(c, paste("03_Modelling/99_Modelled_Species_Distribution_Data/CRIA/", species[[x]], "/data.csv", sep=""))
    unlink(paste("03_Modelling/04_Species_To_Model_Distribution_Data/CRIA/", species[[x]], sep=""), recursive=T)
    }
  
  if(length(p)>0){
    p <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/PADME/", species[[x]], "/data.csv", sep=""))[,-1]
    dir.create(paste("03_Modelling/99_Modelled_Species_Distribution_Data/PADME/", species[[x]], sep=""), showWarnings = F)
    write.csv(p, paste("03_Modelling/99_Modelled_Species_Distribution_Data/PADME/", species[[x]], "/data.csv", sep=""))
    unlink(paste("03_Modelling/04_Species_To_Model_Distribution_Data/PADME/", species[[x]], sep=""), recursive=T)
    }
  
  if(length(r)>0){
    r <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", species[[x]], "/data.csv", sep=""))[,-1]
    dir.create(paste("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA/", species[[x]], sep=""), showWarnings = F)
    write.csv(r, paste("03_Modelling/99_Modelled_Species_Distribution_Data/REFLORA/", species[[x]], "/data.csv", sep=""))
    unlink(paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", species[[x]], sep=""), recursive=T)
    }
  
  rm(c, p, r)
  
  unlink(paste("03_Modelling/04_Species_To_Model_Distribution_Data/", species[[x]],".csv", sep=""))
}


rm(CRIA, PADME, REFLORA)
