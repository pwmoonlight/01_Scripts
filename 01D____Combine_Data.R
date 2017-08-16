###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################



### ------------------------------------- ###
### Find a list of species to be modelled ###
### ------------------------------------- ###

CRIA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/CRIA", full.names=F, recursive=F)
PADME <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/PADME", full.names=F, recursive=F)
REFLORA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA", full.names=F, recursive=F)

species <- c(CRIA, PADME, REFLORA)
species <- as.vector(unique(species))



### -------------------------------------------------- ###
### Read in the distribution data, combine and save it ###
### -------------------------------------------------- ###

writeLines("\n\n#################################################")
writeLines("### The following species are now ready to model:\n")

for(x in species){
  c <- which(CRIA == x)
  p <- which(PADME == x)
  r <- which(REFLORA == x)
  
  cpr <- c()
  
  if(length(c)>0){
    c <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/CRIA/", x, "/data.csv", sep=""))[,-1]
    c <- c[,c(1,2,7)]
    colnames(c) <- c("lat", "long", "species")
    cpr[[1]] <- c}
  if(length(p)>0){
    p <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/PADME/", x, "/data.csv", sep=""))[,-1]
    p <- p[,c(1,2,6)]
    colnames(p) <- c("lat", "long", "species")
    cpr[[2]] <- p}
  if(length(r)>0){
    r <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/REFLORA/", x, "/data.csv", sep=""))[,-1]
    r <- r[,c(1,2,8)]
    colnames(r) <- c("lat", "long", "species")
    cpr[[3]] <- r}
  
  classes <- c(class(c), class(p), class(r))
  classes <- which(classes == "data.frame")
  
  write.csv(do.call(rbind, cpr[classes]), file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/", x, ".csv", sep=""))
  
  writeLines(x)
}

# Clear the R workspace
rm(CRIA, PADME, REFLORA, c, p, r, classes, species, x, distribution_data, cpr)
