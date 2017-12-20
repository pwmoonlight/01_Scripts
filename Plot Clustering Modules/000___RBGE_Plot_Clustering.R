
###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################




### Prepare the working space
### -------------------------

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))


### Read in the nordeste mask
### -------------------------

nordeste <- raster("000_GIS_LAYERS/nordeste.tif")
nordeste <- raster("000_GIS_LAYERS/nordeste.grd")


### Read in the Plot Localities
### ---------------------------

sites <- read.csv("07_Virtual_Plot_Checklists/sites.csv")[,c(2,4,10,11)]
#sites <- sites[which(!is.na(extract(nordeste, SpatialPoints(sites[,c(4:3)])))),]


### Read in the Models
### ------------------

species <- gsub(".tif$", "", list.files("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", pattern="*.tif$", full.names=F))
species.level.dir <- lapply(1:length(species), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)


### Expand the sites matrix to house the results
### --------------------------------------------
sites[,5:(length(species)+4)] <- NA
colnames(sites)[5:(length(species)+4)] <- species


### Populate the results
### --------------------

for(x in 1:nlayers(SDM.b)){
  print(x)  
  sites[,(x+4)] <- extract(SDM.b[[x]], SpatialPoints(sites[4:3]))
}

sites <- t(sites)
colnames(sites) <- sites[1,]
sites <- sites[-c(1:4),]
rownames(sites) <- species
write.csv(sites, "07_Virtual_Plot_Checklists/site_results.csv")




# a) Load all of the required spreadsheets into R
# Obs: I used the read.csv function in here, but there is a limit to the number of rows and columns that csv files can have. If data is being lost because of this, please load tables as text files

spp <- read.csv("07_Virtual_Plot_Checklists/site_results.csv", sep=",", head=TRUE, row.names=1)
colnames(spp) <- gsub("\\.", "", colnames(spp))
colnames(spp) <- gsub("X", "", colnames(spp))

sums <- colSums(spp, na.rm=T)
spp <- spp[,-which(sums==0)]


sites <- read.csv("07_Virtual_Plot_Checklists/sites.csv", fileEncoding = "latin1")[,-1]
sites <- sites[which(sites$AreaID %in% colnames(spp)),]

sppxsites <- as.data.frame(matrix(ncol=2))
colnames(sppxsites) <- c("AreaID", "SppID")
index <- 1
for(x in 1:dim(spp)[1]){
  for(y in 1:dim(spp)[2]){
    if(spp[x,y] == 1){
      print(index)
      sppxsites[index,] <- c(0,0)
      
      sppxsites[index,1] <- colnames(spp)[x]
      sppxsites[index,2] <- rownames(spp)[y]
        
      index <- index+1
    }
  }
}
