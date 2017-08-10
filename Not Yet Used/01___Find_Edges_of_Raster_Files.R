
##-----########################################################-----#
###-----#########################################################-----#
####-----##########################################################-----#
#####-----### by Peter Moonlight, Tiina Sarkinen et al 2017 #########-----#
####-----##########################################################-----#
###-----#########################################################-----#
##-----########################################################-----#



#This script finds the edges of a raster


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



### ------------------
### Read in the Raster
### ------------------

require(raster)
raster <- raster("bio1.asc")

# Plot for a visual check
plot(raster)



### ----------------------------
### Find the edges of the raster
### ----------------------------

edges <- boundaries(raster, type="inner", classes=FALSE, directions=8, asNA=FALSE)

# Plot for a visual check
plot(edges)



### -------------------------------
### Remove the centre of the raster
### -------------------------------

edges_data_frame <- as.data.frame(edges)
centre <- which(edges_data_frame[,1] == 0)
edges[centre] <- NA

# Plot for a visual check
plot(edges)



### ------------------------
### Write the altered raster
### ------------------------

writeRaster(edges, file="edges.asc")
