
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
#library(igraph)
library(geosphere)
library(ellipse)
library(mclust)


############################################################################################################################
# 2) 
############################################################################################################################

#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets") #from Ivan's laptop
#setwd("J:/Jimenez/Nicaragua_Biomes/Datasets") #from desktop 1ZTF at the Lehmann
#species.in.analysis <- read.table("SpeciesInAnalysis_2016May21.txt", header=T, sep=",")
species.in.analysis <- read.table("04_Wombling/Species_In_Analysis.txt", header=T, sep=",")
species.in.analysis[1:5,]
dim(species.in.analysis)

#species.genus.family <- 
#  paste(
#    unlist(lapply(strsplit(as.vector(species.in.analysis[,1]), split=c(" "), fixed = T), function(x) paste(x[1], x[2], sep="_"))),
#    species.in.analysis[,2], sep="_")
#class(species.genus.family)
#species.genus.family[1:5]

#read the file that has the brick with the species distribution models of all species
#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels") #from Ivan's laptop
#setwd("J:/Jimenez/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels") #from desktop 1ZTF at the Lehmann
#dir()
#SDM.b <- brick("SDMb_2015Oct6.grd")
SDM.b <- brick("04_Wombling/SDM_brick.grd")
SDM.b
nlayers(SDM.b)
names(SDM.b)
cell.side <- res(SDM.b)[1] #resolution, i.e., cell dimensions (length and width)
focal.sp.range <- raster(SDM.b, layer=1)
#create a vector with the cell number (or cell ID) of cells that are not-NA
no.na.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range)))]
#length(no.na.cells)
#dim(coordinates(focal.sp.range)[no.na.cells,])
#points(coordinates(focal.sp.range)[no.na.cells,])
s1.s2 <- coordinates(focal.sp.range)[no.na.cells,]

#AOO <- cellStats(SDM.b, 'sum')
#plot(AOO)
#summary(AOO)
#which(AOO<5100 & AOO>4900)
#hist(AOO, breaks=100)

#focal.sp.range <- raster(SDM.b, layer=3)
#focal.sp.range <- raster(SDM.b, layer=35)
#focal.sp.range <- raster(SDM.b, layer=43)
#h<-1

#sp.abbreviations <- unlist(lapply(strsplit(names(SDM.b), split=c("_"), fixed = T), function(x) x[1]))
MCLUST.model.names <- paste("MCLUST", names(SDM.b), sep="_")


############################################################################################################################
# 3)
############################################################################################################################

dir.create("05_Wombling_Null_Models", showWarnings = F)
dir.create("05_Wombling_Null_Models/01_McClust", showWarnings = F)

for(h in 1: nlayers(SDM.b))
{
  focal.sp.range <- raster(SDM.b, layer=h)
  #summary(extract(focal.sp.range, 1:ncell(focal.sp.range))>0)
  presence.cells.index <- extract(focal.sp.range, 1:ncell(focal.sp.range))>0
  #summary(presence.cells.index)
  presence.cells.index[is.na(presence.cells.index)] <- FALSE
  coo <- coordinates(raster(SDM.b, layer=1))[presence.cells.index,]
  #points(coo)
  #vertices <- chull(coo)
  #points(coo[vertices,], col="red")
  #polygon(coo[vertices,], bor="red", lwd=1)
  #areaPolygon(coo[vertices,], a=6378137, f=1/298.257223563)/1e6
  sample.coo <- coo[sample(1:nrow(coo), size=1000),]
  #points(sample.coo, col="red")
  
  #ptm <- proc.time()
  Mcluster.focal.sp <- Mclust(sample.coo, G=1:20)
  #proc.time() - ptm
  if(Mcluster.focal.sp$G == 20) Mcluster.focal.sp <- Mclust(sample.coo, G=10:30)
  if(Mcluster.focal.sp$G == 30) Mcluster.focal.sp <- Mclust(sample.coo, G=20:40)
  if(Mcluster.focal.sp$G == 40) Mcluster.focal.sp <- Mclust(sample.coo, G=3:50)
  
  #MCLUST.model.name <- paste("MCLUST", strsplit(names(SDM.b)[h], split=c("_"), fixed = T)[[1]][1], sep="_")
  assign(MCLUST.model.names[h],
         list(G=Mcluster.focal.sp$G,
              means=Mcluster.focal.sp$parameters$mean,
              sigmas=Mcluster.focal.sp$parameters$variance$sigma,
              pros=Mcluster.focal.sp$parameters$pro))
  
  #setwd("C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/MCLUST_models") #from Ivan's laptop
  #setwd("J:/Jimenez/Nicaragua_Biomes/NullDistributions/MCLUST_models") #from desktop 1ZTF at the Lehmann
  save(list=MCLUST.model.names[h], file=paste("05_Wombling_Null_Models/01_McClust/", MCLUST.model.names[h], ".R", sep=""))
  #MCLUST_acaAlo
  #rm(MCLUST_acaAlo)
  #load("MCLUST_acaAlo.R")
  
  #the line of code below asks for the elements in the R object that holds the results of the analysis,
  #note that further down the script examines several of these elements (e.g., BIC, classification, parameters).
  #attributes(Mcluster.focal.sp) 
  #summary(Mcluster.focal.sp)
  #examine Bayesian Information Criterion (BIC) values,
  #Mcluster.focal.sp$BIC
  #Mcluster.focal.sp$bic
  #dev.new()
  #plot(Mcluster.focal.sp, what=c("BIC"))
  #plot(Mcluster.focal.sp, what=c("BIC"), ylim=c(Mcluster.focal.sp$bic-100, Mcluster.focal.sp$bic))
  #abline(h=Mcluster.focal.sp$bic, lty=3)
  #abline(h=Mcluster.focal.sp$bic-2, lty=3)
  #abline(h=Mcluster.focal.sp$bic-10, lty=3)
  #abline(v=Mcluster.focal.sp$G)
  #plot(Mcluster.focal.sp, what=c("classification"), dimens=c(1,2))
  #Mcluster.focal.sp$parameters
  #Mcluster.focal.sp$parameters$mean
  #Mcluster.focal.sp$parameters$variance$sigma
  plot(focal.sp.range, main=paste(names(SDM.b)[h], ": ", h, " of ", length(names(SDM.b)), sep=""))
  points(t(Mcluster.focal.sp$parameters$mean), col="black", pch=19)
  for (i in 1:Mcluster.focal.sp$G)
  {
    points( 
      ellipse(Mcluster.focal.sp$parameters$variance$sigma[,,i], centre = Mcluster.focal.sp$parameters$mean[,i], level = 0.95, npoints = 100),
      type="l", col="black")
  }
}


############################################################################################################################
# 4)
############################################################################################################################

length(MCLUST.model.names)
length(unique(MCLUST.model.names))

which(is.na(
  match(MCLUST.model.names,
        unlist(lapply(strsplit(list.files("05_Wombling_Null_Models/01_McClust"), split=c("."), fixed = T), function(x) x[1])))
))




