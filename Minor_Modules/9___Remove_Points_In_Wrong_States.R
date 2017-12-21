###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

Reflora_Species <- read.csv("000_Species/Angiospermas_Flora do Brasil 2020_25 set 2017_Working.csv")
States <- read.table("000_Species/States.txt", header=F, sep="\t", fill=T, quote="")

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1", showWarnings = F)

for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1/", species[[x]], sep=""), showWarnings = F)
    
    distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Uncleaned/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
    dims <- dim(distribution_data)[1]
    for(y in length(distribution_data[,1]):1){
      if(is.na(distribution_data[y,10])){distribution_data[y,10] <- ""}
      if(nchar(distribution_data[y,10])>0){
        if(length(as.character(States[which(States[,1] == distribution_data[y,10]),2]))>0){
          distribution_data[y,10] <- as.character(States[which(States[,1]==distribution_data[y,10]),2])
          if(Reflora_Species[which(Reflora_Species$Species==species[[x]]),which(colnames(Reflora_Species)==distribution_data[y,10])]==0){
            distribution_data <- distribution_data[-y,]
          }
        }
        if(length(as.character(States[which(States[,1] == distribution_data[y,10]),2]))==0){
          distribution_data <- distribution_data[-y,]
        }
      
      }
    }
    if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has fewer points..."))}
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1/", species[[x]], ".csv", sep=""))
  }
}
