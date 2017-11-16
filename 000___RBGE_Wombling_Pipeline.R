
###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
###############################################################################################################
###############################################################################################################
###############################################################################################################


############################################################################################################################
# 1) Load needed pakages/Prepare Study Area
############################################################################################################################

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

library(sp)
library(raster)
library(igraph)
library(ape)
library(RColorBrewer)

dir.create("04_Wombling", showWarnings=F)


############################################################################################################################
# 2) Define the study area and examine the relevant spatial links between adjacent grid cells. 
############################################################################################################################

############################################################################################################################
# 2.1) Create a raster mask with a resolution of 5*30 = 150 arc seconds resolution, which is 2.5 minutes
# or 0.04166667 X 0.04166667 degrees. You may skip this step (and go to 2.2) because the raster mask is
# already available.
############################################################################################################################


#Read in an earlier mask if created earlier
Nordeste.mask.0 <- raster("000_GIS_LAYERS/nordeste.tif")
Nordeste.mask.0[Nordeste.mask.0[1:length(Nordeste.mask.0)]==1] <- 0

#examine the total number of grid cells, the frequency of grid cells
#that are NA and not-NA, and determine cell number (or cell ID) of
#the cells that are not NA
ncell(Nordeste.mask.0) #total number of cells
#create a vector with the cell number (or cell ID) of cells that are not-NA
no.na.cells <- (1:ncell(Nordeste.mask.0))[!is.na(extract(Nordeste.mask.0, 1:ncell(Nordeste.mask.0)))]
length(no.na.cells)


############################################################################################################################
# 2.3) Define adjacent grid cells. You might want to skip to the end of this section, and use previously defined spatial
# links between grid cells, stored in the text file "cell_adj.txt".  
############################################################################################################################

# Define links between grid cells
source(paste(getwd(), "/01_Scripts/Wombling Modules/001___Produce_Cell_Adj_Object.R", sep=""))

#read text file with previously defined links between grid cells
cell.adj <- read.table("04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",", header=T)


#############################################################################
# 3) Create an R object of class "brick" with all species distribution models
#############################################################################

#list the modelled species with sufficiently high CBIs
species.in.analysis <-  gsub(".tif$", "", list.files("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", pattern="*.tif", full.names=F))

#save in a file the names of the species included in the analysis
write.table(species.in.analysis, file="04_Wombling/species_in_analysis.txt", quote=T, sep=",", row.names=F, col.names=F)

species.level.dir <- lapply(1:length(species.in.analysis), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species.in.analysis)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)


#save in a file the brick with the species distribution models of all species in the analysis that has 30 arc seconds resolution
writeRaster(SDM.b, filename="04_Wombling/SDM_brick.grd", bandorder='BIL')

#read the file that has the brick with the species distribution models of all species
SDM.b <- brick("04_Wombling/SDM_brick.grd")
brick.index.species.in.phylogeny <- names(SDM.b)

############################################################################################################################
# 4) Read and examine the phylogeny for the species in the analysis (the "Nordeste phylogeny")
############################################################################################################################
 
############################################################################################################################
# 5) Determine the subset of species included in the brick with the species distribution models that are
# represented in the phylogeny.
############################################################################################################################
  
  
 
  
############################################################################################################################
# 6) Calculate beta-diversity for the spatial links created in section 2 above.
############################################################################################################################
  
############################################################################################################################
# 6.1) Calculate taxonomic (i.e., species-based) beta-diversity using Simpson's and Sorensen's indices
############################################################################################################################

source(paste(getwd(), "/01_Scripts/Wombling Modules/002___Calculate_Beta_Diversity_Values.R", sep=""))

#read file with beta-diversity measured as taxonomic (i.e., species based) Sorensen's or Simpson's indices
obs.beta.sor <- read.table("04_Wombling/ObsBetaSor.txt", header=T, sep=",")
obs.beta.sim <- read.table("04_Wombling/ObsBetaSim.txt", header=T, sep=",")

  
############################################################################################################################
# 6.2) Calculate phylogenetic beta-diversity using Simpson's and Sorensen's indices
############################################################################################################################
  


############################################################################################################################
# 7) Define the beta-diversity metric to use in subsequent analyses. At this
# point there are four metrics available: taxonomic and phylogenetic versions of
# Sorensen's and Simpson's indexes.
############################################################################################################################

obs.beta <- obs.beta.sor
obs.beta <- obs.beta.sim
#obs.beta <- obs.phylo.beta.sor
#obs.beta <- obs.phylo.beta.sim


############################################################################################################################
# 8) Rank the values of the chosen beta-diversity metric. Note that the
# ranking involves random resolution of ties. Therefore, to ensure exact
# replication of results, you might want to skip to the end of this section
# and gather a previously created file with the ranks of beta-diversity values. 
############################################################################################################################

source(paste(getwd(), "/01_Scripts/Wombling Modules/003___Rank_Beta_Diversity_Values.R", sep=""))

#save ranks of observed beta values: useful because the raking has a random
#component to deal with ties, thus fully consistent results are unlikely
#using different iterations of the ranking
write.table(r.obs.beta, file="04_Wombling/Rank_ObsBetaSor.txt", row.names=F)
write.table(r.obs.beta, file="04_Wombling/Rank_ObsBetaSim.txt", row.names=F)
#write.table(r.obs.beta, file="04_Wombling/Rank_ObsPhyBetaSor.txt", row.names=F)
#write.table(r.obs.beta, file="04_Wombling/Rank_ObsPhyBetaSim.txt", row.names=F)

#read ranks of observed beta values

#read text file with previously defined links between grid cells
cell.adj <- read.table("04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",", header=T)

r.obs.beta <- read.table("04_Wombling/Rank_ObsBetaSim.txt", header=T, sep=",")[,1]
r.obs.beta <- read.table("04_Wombling/Rank_ObsBetaSor.txt", header=T, sep=",")[,1]
#r.obs.beta <- read.table("04_Wombling/Rank_ObsPhyBetaSor.txt", header=T, sep=",")
#r.obs.beta <- read.table("04_Wombling/Rank_ObsPhyBetaSim", header=T, sep=",")


############################################################################################################################
# 9) Map beta-diversity values for the spatial links created in section 2 (above).
# You may choose an "animated" loop (section 7.1 below) that displays increasingly lower values of beta-diversity,
# or "static" mapping (section 7.2 below) wherby a single set of high beta-diversity values is displayed.
############################################################################################################################


############################################################################################################################
# 9.1) Map beta-diversity values using an "animated" loop that displays increasingly lower values of beta-diversity.
############################################################################################################################

source(paste(getwd(), "/01_Scripts/Wombling Modules/004___Plot_Beta_Diversity_Values_In_Animated_Loop.R", sep=""))



############################################################################################################################
# 9.2) Map beta-diversity values using a "static" map to display a single set of high beta-diversity values.
############################################################################################################################

source(paste(getwd(), "/01_Scripts/Wombling Modules/005___Plot_Beta_Diversity_Values_In_Static_Form.R", sep=""))


############################################################################################################################
# 10) Sequentially remove spatial links (or deploy "candidate boundary elements") according to the rank created in section 8
# (above), and tally the number of resulting isolated regions (potential ecoregions, or "subgraphs" in the lingo of graphs
# theory). The removal of spatial links starts with the highest beta-diversity values, and proceeds to increasingly lower
# beta-diversity values. This step might take several hours or days, and you may skip to the end of this section to read a
# file with previously calculated number of isolated regions for a given number of spatial links removed, according to four
# beta-diversity metrics.
############################################################################################################################

#express the map of Nordeste as a "graph" (in the sense of graph theory),
#in which the "vertices" (or "vertexes") are grid cells, and the "edges"
#are spatial links between adjacent grid cells as defined in section 2.
#
#first create a matrix defining the edges:
cell.adj.char <- cbind(as.character(cell.adj[,1]), as.character(cell.adj[,2]))
dim(cell.adj.char)

#next define the graph and examine the result
Nordeste.graph <- graph_from_edgelist(cell.adj.char, directed = F)



#run a loop to sequentially remove spatial links and count the resulting number
#of subgraphs, which correspond to isolated regions (and potentially ecoregions).
number.of.subgraphs <- rep(NA, times=length(obs.beta[,1]))
plot(0, components(Nordeste.graph)$no, 
     xlim=c(0, sum(obs.beta[,1] >= quantile(obs.beta[,1], probs=0))),
     ylim=c(15,length(obs.beta[,1])),
     xlab="Candidate boundary elements deployed", ylab="Regions (or subgraphs)",
     pch=19, bty="n", cex.axis=1.5, cex.lab=1.5) 
#
start.time <- Sys.time()
# 
for(i in 1:length(obs.beta[,1]))
{
  print(i/143638*100)
  candidate.boundary.elements.to.deploy <- which(r.obs.beta > (length(r.obs.beta)-i))
  number.of.subgraphs[i] <- components(delete_edges(Nordeste.graph, candidate.boundary.elements.to.deploy))$no
  #points(i, number.of.subgraphs[i], pch=19)
}
difftime(Sys.time(), start.time, units="mins")
#this procedure might take about 1 or 2 minutes, depending on which computer is used
#

#save file with number of subgraphs, derived from the taxonomic (i.e, species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
write.table(number.of.subgraphs, file="04_Wombling/NumberSubgraphsSim.txt", sep=",", row.names=F)
write.table(number.of.subgraphs, file="04_Wombling/NumberSubgraphsSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSims.txt", sep=",", row.names=F)

#read file with number of subgraphs, derived from Sorensen's index
number.of.subgraphs <- read.table("04_Wombling/NumberSubgraphsSor.txt", header=T, sep=",")[,1]
number.of.subgraphs <- read.table("04_Wombling/NumberSubgraphsSim.txt", header=T, sep=",")[,1]
#number.of.subgraphs <- read.table("NumberSubgraphsPhyloSor.txt", header=T, sep=",")
#number.of.subgraphs <- read.table("NumberSubgraphsPhyloSim.txt", header=T, sep=",")



##################################################################################################
# 11) Calculatesuperfluidity, as defined by Oden et al. (1993, cited in the introduction).
# This step may take several hours or days, and you might skip to the end of this section
# to read a file with previously calculatedsuperfluidity values. Use code under 9.1 or 9.2 according
# to the version of package igraph and R you are using.
##################################################################################################


#define the percentiles of beta-diversity for whichsuperfluidity will be claculated,
#if you have not done so in section 8 (above)
beta.ranks.to.evaluate <- 143638 - round(143638*seq(0.05, 0.5, 0.05))

source(paste(getwd(), "/01_Scripts/Wombling Modules/007___Calculate_Superfluidity.R", sep=""))

#write files withsuperfluidity values, derived from the taxonomic (i.e., species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
write.table(cbind(evaluation.number.of.subgraphs,superfluidity), file="04_Wombling/superfluiditySor.txt", sep=",", row.names=F)
#write.table(cbind(evaluation.number.of.subgraphs,superfluidity), file="04_Wombling/superfluidityPhyloSor.txt", sep=",", row.names=F)
write.table(cbind(evaluation.number.of.subgraphs,superfluidity), file="04_Wombling/superfluiditySim.txt", sep=",", row.names=F)
#write.table(cbind(evaluation.number.of.subgraphs,superfluidity), file="04_Wombling/superfluidityPhyloSim.txt", sep=",", row.names=F)

#read file withsuperfluidity values, dderived from the taxonomic (i.e., species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
#superfluidity <- read.table("04_Wombling/superfluiditySor.txt", header=T, sep=",")
#superfluidity <- read.table("04_Wombling/superfluidityPhyloSor.txt", header=T, sep=",")
superfluidity <- read.table("04_Wombling/superfluiditySim.txt", header=T, sep=",")
#superfluidity <- read.table("04_Wombling/superfluidityPhyloSim.txt", header=T, sep=",")



##################################################################################################
# 12) Examine and map the regions (potential ecoregions) or "subgraphs" obtained in section 10
##################################################################################################

##################################################################################################
# 12.1) Obtain and examine data on the regions
##################################################################################################

#select the number of regions (potential ecoregions) or subgraphs to map
number.of.subgraphs[143638 - beta.ranks.to.evaluate]
pick.num.subgraphs <- 17908

source(paste(getwd(), "/01_Scripts/Wombling Modules/008a___Create_Raster_Map_SIM.R", sep=""))
source(paste(getwd(), "/01_Scripts/Wombling Modules/008b___Create_Raster_Map_SOR.R", sep=""))


Nordeste.mask.Regions <- raster(paste("04_Wombling/Nordeste_SIM_", pick.num.subgraphs,"_Regions.grd", sep=""))
Nordeste.mask.Regions <- raster(paste("04_Wombling/Nordeste_SOR_", pick.num.subgraphs,"_Regions.grd", sep=""))


##################################################################################################
# 12.2) Map the regions
##################################################################################################

require(rgdal)
BRA_ADM <- list.files("000_GIS_LAYERS/BRA_Adm_2/", pattern="[.]shp$", full.names=T)
BRA_ADM <- readOGR(BRA_ADM[1])


#define and examine colors to map regions
number.of.colors <- pick.num.subgraphs
subgraph.col <- colorRampPalette(brewer.pal(12,"Paired"))(number.of.colors)
subgraph.col <- sample(subgraph.col, length(subgraph.col))
subgraph.col[1:20]
plot(1:number.of.colors, col = subgraph.col, pch = 16, cex = 3)
length(subgraph.col)
subgraph.col <- rep(subgraph.col, length.out=pick.num.subgraphs)

#draw the map of regions at broad scale
plot(Nordeste.mask.Regions, col=subgraph.col)
plot(BRA_ADM, add=T)










#create a matrix with the species composition for each region
SpeciesRegions <- matrix(NA, nrow=pick.num.subgraphs, ncol=length(brick.index.species.in.phylogeny))
for(i in 1:length(brick.index.species.in.phylogeny))
{
  print(i/length(brick.index.species.in.phylogeny)*100)
  SpeciesRegions[,i] <- zonal(raster(SDM.b, layer=brick.index.species.in.phylogeny[i]), Nordeste.mask.Regions, fun='sum', digits=0, na.rm=TRUE)[,2] 
}
#examine resulting matrix
SpeciesRegions[1:5,1:5]
dim(SpeciesRegions)
#convert the matrix data to logical (true/false) data
SpeciesRegions <- SpeciesRegions > 0
colnames(SpeciesRegions) <- species.genus[]
SpeciesRegions[SpeciesRegions==T] <- 1
#examine the results
SpeciesRegions[1:5,1:5]
dim(SpeciesRegions)

#save files with species composition of for each region
write.table(SpeciesRegions, paste("04_Wombling/SpeciesRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
write.table(SpeciesRegions, paste("04_Wombling/SpeciesRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
#write.table(SizeRegions, paste("04_Wombling/SizeRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
#write.table(SizeRegions, paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")

#read files with species composition of for each region
SpeciesRegions <- read.table(paste("04_Wombling/SpeciesRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), header=T, sep=",")
SpeciesRegions <- read.table(paste("04_Wombling/SpeciesRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), header=T, sep=",")
head(SpeciesRegions)
dim(SpeciesRegions)
