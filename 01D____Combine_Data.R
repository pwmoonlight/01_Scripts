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

CRIA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA", full.names=F, recursive=F)
PADME <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/PADME", full.names=F, recursive=F)
REFLORA <- list.dirs("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/REFLORA", full.names=F, recursive=F)

species <- c(CRIA, REFLORA)
species <- as.vector(unique(species))


for(x in length(species):1){if(grepl("1", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("2", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("3", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("4", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("5", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("6", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("7", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("8", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("9", species[x])){species <- species[-x]}}
for(x in length(species):1){if(grepl("0", species[x])){species <- species[-x]}}

### -------------------------------------------------- ###
### Read in the distribution data, combine and save it ###
### -------------------------------------------------- ###

writeLines("\n\n#################################################")
writeLines("### The following species are now ready to model:\n")

for(x in species){
  c <- which(CRIA == x)
  #p <- which(PADME == x)
  r <- which(REFLORA == x)
  
  cpr <- c()
  
  if(length(c)>0){
    c <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/CRIA/", x, "/data.csv", sep=""))[,-1]
    #c <- c[,c(1,2,7)]
    #colnames(c) <- c("lat", "long", "species")
    cpr[[1]] <- c}
  #if(length(p)>0){
  #  p <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/PADME/", x, "/data.csv", sep=""))[,-1]
  #  p <- p[,c(1,2,6)]
  #  colnames(p) <- c("lat", "long", "species")
  #  cpr[[2]] <- p}
  if(length(r)>0){
    r <- read.csv(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/REFLORA/", x, "/data.csv", sep=""))[,-1]
    #r <- r[,c(1,2,8)]
    #colnames(r) <- c("lat", "long", "species")
    cpr[[3]] <- r}
  
  classes <- c(class(c), class(NULL), class(r))
  classes <- which(classes == "data.frame")
  
  write.csv(do.call(rbind, cpr[classes]), file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/", x, ".csv", sep=""))
  
  writeLines(x)
}

# Clear the R workspace
rm(CRIA, PADME, REFLORA, c, p, r, classes, species, x, distribution_data, cpr)
