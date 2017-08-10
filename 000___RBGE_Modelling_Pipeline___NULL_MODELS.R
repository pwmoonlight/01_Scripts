
###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
###############################################################################################################
###############################################################################################################
###############################################################################################################


#####################################################################################################################
##                                                                                                                 ##
##                                                                                                                 ##
#####################################################################################################################



### Prepare the working space
### -------------------------

setwd("C:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))


dir.create("50_Null_Presence_data", showWarnings = F)
dir.create("51_Null_Absence_data", showWarnings = F)
dir.create("52_Null_CBI", showWarnings = F)



### Read in the Background data
### ---------------------------

bg_data <- read.csv("06_BG_data/BG_data.csv", header=T)
rownames(bg_data) <- bg_data[,1]
bg_data <- bg_data[,-1]
bg_data <- bg_data[,-3:(0-length(colnames(bg_data)))]



### Sample the BG Data to produce Null Species
### ------------------------------------------

sample_sizes <- append(append(c(seq(from=50, to=100, by=10)), c(seq(from=120, to=200, by=20))), append(c(seq(from=260, to=500, by=60)), c(seq(from=600, to=1000, by=100))))

for(x in sample_sizes){
  for(y in 1:1000){
    null_sample <- bg_data[sample(1:length(bg_data[,1]), size=x, replace=F),]
    write.csv(null_sample, file=paste("50_Null_Presence_data/", x, "_", y, ".csv", sep=""))
  }
}



### Produce BG data for each Null Species
### -------------------------------------

kde_raster <- raster("KDE_Raster_CRIA_Data_masked___PWM_10_5_2017.asc")
require(maptools)
require(dismo)
require(rgeos)

for(x in sample_sizes){
  for(y in 1:1000){
    writeLines(paste(x, ": ", y/10, "%", sep=""))
    
    circle <- read.csv(paste("50_Null_Presence_data/", x, "_", y, ".csv", sep=""))[,-1]
    
    circle <- circles(circle[,1:2], d=80000, lonlat=TRUE)
    circle <- gUnaryUnion(circle@polygons)
    
    circle <- cellFromPolygon(kde_raster, circle)
    circle <- xyFromCell(kde_raster, circle[[1]])

    background <- circle[sample(1:length(circle[,1]), replace=TRUE, size=1000),]
    
    write.csv(background, file=paste("51_Null_Absence_data/", x, "_", y, ".csv", sep=""))
  }
}
