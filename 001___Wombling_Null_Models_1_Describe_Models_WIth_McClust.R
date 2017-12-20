
############################################################################################################################
############################################################################################################################
#
# INTRODUCTION TO THE SCRIPT "SDMData_30arc_MCLUST_for_NullModel_2016June21.R"

#
# A) What does this script do?
#
# B) What is needed to run this script?
#
#
############################################################################################################################
############################################################################################################################


############################################################################################################################
# 1) Load needed pakages
############################################################################################################################

setwd("G:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

library(sp)
library(raster)
library(geosphere)
library(ellipse)
library(mclust)


############################################################################################################################
# 2) 
############################################################################################################################


#read the file that has the brick with the species distribution models of all species
species.in.analysis <-  gsub(".tif$", "", list.files("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", pattern="*.tif$", full.names=F))
species.level.dir <- lapply(1:length(species.in.analysis), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species.in.analysis)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)

cell.side <- res(SDM.b)[1] #resolution, i.e., cell dimensions (length and width)
focal.sp.range <- raster(SDM.b, layer=1)
#create a vector with the cell number (or cell ID) of cells that are not-NA
no.na.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range)))]
s1.s2 <- coordinates(focal.sp.range)[no.na.cells,]




MCLUST.model.names <- paste("MCLUST", names(SDM.b), sep="_")


############################################################################################################################
# 3)
############################################################################################################################

dir.create("05_Wombling_Null_Models", showWarnings = F)
dir.create("05_Wombling_Null_Models/01_McClust", showWarnings = F)

for(h in 1:nlayers(SDM.b))
{
  if(!dir.exists(paste("05_Wombling_Null_Models/01_McClust/", species.in.analysis[[h]], sep=""))){
    dir.create(paste("05_Wombling_Null_Models/01_McClust/", species.in.analysis[[h]], sep=""), showWarnings = F)
    
    writeLines(paste(names(SDM.b)[h], ": ", h, " of ", length(names(SDM.b)), sep=""))

    focal.sp.range <- raster(SDM.b, layer=h)
    presence.cells.index <- extract(focal.sp.range, 1:ncell(focal.sp.range))>0
    presence.cells.index[is.na(presence.cells.index)] <- FALSE
    coo <- coordinates(raster(SDM.b, layer=1))[presence.cells.index,]
  
    if(length(coo[,1])>0){
    
      sample.coo <- coo[sample(1:nrow(coo), size=1000, replace=T),]
  
      #ptm <- proc.time()
      Mcluster.focal.sp <- Mclust(sample.coo, G=1:20)
      if(Mcluster.focal.sp$G == 20) Mcluster.focal.sp <- Mclust(sample.coo, G=10:30)
      if(Mcluster.focal.sp$G == 30) Mcluster.focal.sp <- Mclust(sample.coo, G=20:40)
      if(Mcluster.focal.sp$G == 40) Mcluster.focal.sp <- Mclust(sample.coo, G=3:50)
  
      assign(MCLUST.model.names[h],
         list(G=Mcluster.focal.sp$G,
              means=Mcluster.focal.sp$parameters$mean,
              sigmas=Mcluster.focal.sp$parameters$variance$sigma,
              pros=Mcluster.focal.sp$parameters$pro))
  
     save(list=MCLUST.model.names[h], file=paste("05_Wombling_Null_Models/01_McClust/", MCLUST.model.names[h], ".R", sep=""))
  
    }
    #plot(focal.sp.range, main=paste(names(SDM.b)[h], ": ", h, " of ", length(names(SDM.b)), sep=""))
    #points(t(Mcluster.focal.sp$parameters$mean), col="black", pch=19)
    #for (i in 1:Mcluster.focal.sp$G)
    #{
    #  points( 
    #    ellipse(Mcluster.focal.sp$parameters$variance$sigma[,,i], centre = Mcluster.focal.sp$parameters$mean[,i], level = 0.95, npoints = 100),
    #    type="l", col="black")
    #}
  }
}


