
##-----########################################################-----#
###-----#########################################################-----#
####-----##########################################################-----#
#####-----### by Peter Moonlight, Tiina Sarkinen et al 2017 #########-----#
####-----##########################################################-----#
###-----#########################################################-----#
##-----########################################################-----#

#This script reads data as downloaded from CRIA species link
#It then compares the data a raster of altitude
#Specimens whose are >100m different to stated are saved separate file for checking



### ---------------------
### Prepare the workspace
### ---------------------

# Set the working directory
setwd("C:/000_Modeling_Working_Directory_000")

# Clear the R workspace
rm(list = ls(all=T))

# Allocate larger than usual memory and 
memory.size(5000000)
options(java.parameters = "-Xmx12288m")



### -----------------------------
### Load in the distribution data
### -----------------------------

require("XLConnect")

#Distribution data should be a xls file downloaded "as it comes" from CRIA Species Link
distribution_data <- list.files('02_Input_Distribution_Data', pattern="[.]xls$", full.names=T)
distribution_data <- loadWorkbook(distribution_data)
distribution_data <- readWorksheet(distribution_data, sheet="sheet1", header=T)
distribution_data <- as.data.frame(distribution_data)

distribution_data[,34] <- as.numeric(distribution_data[,34])
distribution_data[,35] <- as.numeric(distribution_data[,35])

# Create a Spatial Points version of the distribution data
distribution_data_SP <- SpatialPoints(distribution_data[,c(34,35)], proj4string=CRS("+proj=longlat +ellps=GRS80 +datum=WGS84"))

# Plot the input data on a crude world map for a visual check
plot(distribution_data$longitude, distribution_data$latitude, col='orange', pch=20, cex=0.75)
points(distribution_data$long, distribution_data$lat, col='red', cex=0.75)
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, add=T)

require(rgeos)



### -----------------------------
### Load in the altitude data
### -----------------------------

library(raster)

alt <- raster("Y:/World/elevation_1km_resolution/srtm_1km.asc")

#Plot for a visual check
plot(distribution_data_SP)
plot(alt, add=T)
plot(distribution_data_SP, add=T)

# Create a data frame into which to put the results
differences <- as.data.frame(matrix(nrow=length(distribution_data_SP), ncol=4))
colnames(differences) <- c("alt", "alt_min", "alt_max", "diff")



### --------------------#
### ----------------------#
### Run the analysis--------#
### ----------------------#
### --------------------#

# Choose a cutoff. This is the number of meters difference between the stated and "actual" altitude you wish to be highlighted
cutoff <- 100

# Extract the altitude at all points
differences[,1] <- extract(alt, distribution_data_SP)

# Extract the stated altitude
differences[,2] <- as.numeric(distribution_data[,40])
differences[,3] <- as.numeric(distribution_data[,41])
differences[,4] <- 0

# Find the difference
for(x in 1:length(differences[,1])){
  if(!is.na(differences[x,2])){
    if(differences[x,2]!=0){
      if(differences[x,2]!=9999){
        A <- max(differences[x,1], differences[x,2]) - min(differences[x,1], differences[x,2])
        B <- max(differences[x,1], differences[x,3]) - min(differences[x,1], differences[x,3])
        differences[x,4] <- min(c(A,B))
      }
    }
  }
}

# Save a csv of those records with a larger difference than the cutoff
above_cutoff <- distribution_data[which(differences[,4]>cutoff),]
above_cutoff[,52] <- differences[which(differences[,4]>cutoff),4]
colnames(above_cutoff)[52] <- "altitude_difference"
write.csv(above_cutoff, "distribution_data_not_matching_altitudes.csv")
