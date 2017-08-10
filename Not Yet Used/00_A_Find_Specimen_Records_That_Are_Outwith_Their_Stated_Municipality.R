
##-----########################################################-----#
###-----#########################################################-----#
####-----##########################################################-----#
#####-----### by Peter Moonlight, Tiina Sarkinen et al 2017 #########-----#
####-----##########################################################-----#
###-----#########################################################-----#
##-----########################################################-----#

#This script reads data as downloaded from CRIA species link
#It then compares the data to a shapefile of municipalities in the northeast of Brazil
#Specimens whose coordinates fall outside their stated municipality are saved in a separate file for checking

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



### -----------------------------------
### Load in Brazilian Municipality File
### -----------------------------------

require(rgdal)
municipalities <- readOGR("Y:/South America GIS/Brasil/divisao_municipal_sab/divisao_municipal_sab.SHP", layer="divisao_municipal_sab")
#If needed, this can be downloaded from http://www.insa.gov.br/sigsab/acervoDigital



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



### -----------------------------------------------------------------------
### Find municipalities in both the municipality file and distribution data
### -----------------------------------------------------------------------

distribution_data_municipalities <- toupper(iconv(distribution_data$county, to='ASCII//TRANSLIT'))
municipalities_municipalities <- toupper(iconv(municipalities@data$NM_MUNICIP, to='ASCII//TRANSLIT'))

match <- which(distribution_data_municipalities %in% municipalities_municipalities)



### -----------------------------------------------------------------------
### Load in a function for turning the municipality data into a raster file
### -----------------------------------------------------------------------

FunR<-function(r){
  ext<-raster(extent(r), nrow=(extent(r)[4]-extent(r)[3])/0.01, ncol=(extent(r)[2]-extent(r)[1])/0.01)
  crs(ext)<-crs(r)
  D<-rasterize(r,ext,field=1,update=T)
}



### ---------------------------------------------------------------------------------------
### Add a column into the distribution data to record any distance outside the municipality
### ---------------------------------------------------------------------------------------

distribution_data[,52] <- NA
colnames(distribution_data)[52] <- "dist_outside_munic"



### --------------------#
### ----------------------#
### Run the analysis--------#
### ----------------------#
### --------------------#

#Loop around every point where the distribution data municipality matches a municipality in the municipalities file
for(x in 1:length(distribution_data_municipalities[match])){
  #Determine which municipality in it matches
  y <- which(distribution_data_municipalities[match][x] == municipalities_municipalities)
  #Extract only that municipality 
  municipality <- municipalities[y,]
  #Turn that municipality into a raster file
  municipality <- FunR(municipality)
  #Extract a value. If the point is in the municipality this will be 1, if not it will be NA
  value <- extract(municipality, distribution_data_SP[match][x])
  
  #If the point is outside the municipality
  if(is.na(value)){
    #Find all points in the municipality
    cell.values <- extract(municipality, coordinates(municipality), cellnumbers=T)
    CoorNotNA <- xyFromCell(municipality, cell.values[!is.na(cell.values[,2]),1])
    #Add the distance to that municipality to the column in the distribution data
    distribution_data[,52][match][x] <- min(pointDistance(distribution_data_SP[match][x], CoorNotNA, lonlat=T))/1000
  }
}

#Write a csv with those data points outside the municipalities. The final line of that file is the distance to the municipality
write.csv(distribution_data[which(distribution_data[,52]>0),], "distribution_data_not_matching_municipalities.csv")
