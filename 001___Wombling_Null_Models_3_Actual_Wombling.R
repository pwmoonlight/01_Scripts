
############################################################################################################################
############################################################################################################################
#
#
############################################################################################################################
############################################################################################################################


############################################################################################################################
# 1) Load needed pakages/Prepare Study Area
############################################################################################################################

setwd("G:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

library(sp)
library(raster)
library(igraph)
library(ape)
library(RColorBrewer)


############################################################################################################################
# 2) Define the study area and examine the relevant spatial links between adjacent grid cells. 
############################################################################################################################

############################################################################################################################
# 2.1) Read in the mask of the study area
############################################################################################################################

#create and plot the Nordeste raster mask that indicates the grid cells that have climate data,
Nordeste.mask.0 <- raster("000_GIS_LAYERS/nordeste.tif")

############################################################################################################################
# 2.3) Read in the adjacent grid cells  
############################################################################################################################

cell.adj <- read.table("04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",", header=T)

#############################################################################
# 3) Create an R object of class "brick" with all species distribution models
#############################################################################

bricks <- gsub("Null_", "", list.files("05_Wombling_Null_Models/03_Null_Bricks", pattern="*.grd"))
bricks <- gsub("_SDMb.grd", "", bricks)
species.in.analysis <- read.csv("04_Wombling/species_in_analysis.txt", header=F)

rasterOptions(datatype="LOG1S") #specify datatype (true/false or presence/absence) to efficiently save the brick file  

# Create Folders in which to store the results
dir.create("05_Wombling_Null_Models/04_Null_Beta_Diversity", showWarnings = F)
lapply(bricks, function(x){dir.create(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", x, sep=""), showWarnings = F)})

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

brick.index.species.in.phylogeny <- as.character(species.in.analysis[,1])
SDM.b <- brick("04_wombling/SDM_brick.grd")
brick.index.species.in.phylogeny <- which(names(SDM.b) %in% gsub(" ", "_", brick.index.species.in.phylogeny))

for(x in 107:length(bricks)){
  SDM.b <- brick(paste("05_Wombling_Null_Models/03_Null_Bricks/Null_", bricks[[x]], "_SDMb.grd", sep=""))
  
  writeLines(paste("\nWorking on brick ", x, " of ", length(bricks), sep=""))
  
  writeLines("...Calculating beta diversity")
  source(paste(getwd(), "/01_Scripts/Wombling Modules/01_SIM_and_SOR_diversity.R", sep=""))

############################################################################################################################
# 6.2) Calculate phylogenetic beta-diversity using Simpson's and Sorensen's indices.
############################################################################################################################


############################################################################################################################
# 7) Rank the values of the beta-diversity metrics.
############################################################################################################################

  writeLines("...Ranking beta diversity values")
  source(paste(getwd(), "/01_Scripts/Wombling Modules/02_rank_SIM_and_SOR_diversity.R", sep=""))
  
}
  




############################################################################################################################
# 9) Map beta-diversity values for the spatial links created in section 2 (above).
# You may choose an "animated" loop (section 7.1 below) that displays increasingly lower values of beta-diversity,
# or "static" mapping (section 7.2 below) wherby a single set of high beta-diversity values is displayed.
############################################################################################################################

#If you haven't, go to section 2 and read (or create) the Nordeste raster mask that has a resolution of 5*30 = 150 arc seconds,
#which is 2.5 minutes or 0.04166667 X 0.04166667 degrees. That raster mask should be named "Nordeste.mask.0".

############################################################################################################################
# 9.1) Map beta-diversity values using an "animated" loop that displays increasingly lower values of beta-diversity.
############################################################################################################################

#plot the mask, broad scale
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F)

selected.quantile <- seq(0.95, 0.5, -0.05)
quantile(obs.beta, probs=selected.quantile[1])
values.at.above.quantile <- sum(obs.beta >= quantile(obs.beta, probs=selected.quantile[1]))
values.at.above.quantile

flag.candidate.boundary.elements <- r.obs.beta > (length(r.obs.beta)-values.at.above.quantile)
sum(flag.candidate.boundary.elements)

for(i in 1:length(selected.quantile))
{
  values.at.above.quantile <- sum(obs.beta >= quantile(obs.beta, probs=selected.quantile[i]))
  flag.candidate.boundary.elements <- r.obs.beta > (length(r.obs.beta)-values.at.above.quantile)

  from.coor <- xyFromCell(Nordeste.mask.0, cell.adj[flag.candidate.boundary.elements,1], spatial=FALSE)
  to.coor <- xyFromCell(Nordeste.mask.0, cell.adj[flag.candidate.boundary.elements,2], spatial=FALSE)
  arrows(from.coor[,1], from.coor[,2], to.coor[,1], to.coor[,2], length = 0, code = 2, col="red")
  legend("topright", paste("Quantile =", selected.quantile[i]), bty="o", bg="white", box.col="white")
  Sys.sleep(1)
}

#remove matrices with grid cell coordinates: they may be large and thus occupy significant space, are not
#needed after this section, and are easily re-created as needed.
rm(from.coor, to.coor, flag.candidate.boundary.elements)

############################################################################################################################
# 9.2) Map beta-diversity values using a "static" map to display a single set of high beta-diversity values.
############################################################################################################################

#plot the mask, broad scale
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F)
#plot the mask, narrower scale
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F, xlim=c(-87,-85), ylim=c(13,15))
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F, xlim=c(-87,-85), ylim=c(11,13))
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F, xlim=c(-85,-83), ylim=c(10.5,13))
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F, xlim=c(-85,-83), ylim=c(13,15.5))

selected.quantile <- 0.99
quantile(obs.beta, probs=selected.quantile)
values.at.above.quantile <- sum(obs.beta >= quantile(obs.beta, probs=selected.quantile))
values.at.above.quantile

flag.candidate.boundary.elements <- r.obs.beta > (length(r.obs.beta)-values.at.above.quantile)
sum(flag.candidate.boundary.elements)

from.coor <- xyFromCell(Nordeste.mask.0, cell.adj[flag.candidate.boundary.elements,1], spatial=FALSE)
to.coor <- xyFromCell(Nordeste.mask.0, cell.adj[flag.candidate.boundary.elements,2], spatial=FALSE)
arrows(from.coor[,1], from.coor[,2], to.coor[,1], to.coor[,2], length = 0, code = 2, col="red")
legend("topright", paste("Quantile =", selected.quantile), bty="o", bg="white", box.col="white")

#remove matrices with grid cell coordinates: they may be large and thus occupy significant space, are not
#needed after this section, and are easily re-created as needed.
rm(from.coor, to.coor, flag.candidate.boundary.elements)


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
as.matrix(cell.adj)[1:5,]
length(unique(c(cell.adj[,1], cell.adj[,2])))
length(no.na.cells)
cell.adj.char <- cbind(as.character(cell.adj[,1]), as.character(cell.adj[,2]))
cell.adj.char[1:5,]
dim(cell.adj.char)
#
#next define the graph and examine the result
Nordeste.graph <- graph_from_edgelist(cell.adj.char, directed = F)
Nordeste.graph
class(Nordeste.graph)
components(Nordeste.graph)$membership[1:100]
components(Nordeste.graph)$csize
components(Nordeste.graph)$no
length(components(Nordeste.graph)$membership)
length(V(Nordeste.graph))
length(unique(c(cell.adj[,1], cell.adj[,2])))
length(no.na.cells)
length(no.na.cells) - length(V(Nordeste.graph))
length(components(Nordeste.graph)$csize)
components(Nordeste.graph)$no
is.connected(Nordeste.graph)

#run a loop to sequentially remove spatial links and count the resulting number
#of subgraphs, which correspond to isolated regions (and potentially ecoregions).
number.of.subgraphs <- rep(NA, times=length(obs.beta))
plot(0, components(Nordeste.graph)$no, 
     xlim=c(0, sum(obs.beta >= quantile(obs.beta, probs=0))),
     ylim=c(15,50000),
     xlab="Candidate boundary elements deployed", ylab="Regions (or subgraphs)",
     pch=19, bty="n", cex.axis=1.5, cex.lab=1.5) 
#
start.time <- Sys.time()
# 
for(i in 1:length(obs.beta))
{
  print(i/35082*100)
  candidate.boundary.elements.to.deploy <- which(r.obs.beta > (length(r.obs.beta)-i))
  number.of.subgraphs[i] <- components(delete_edges(Nordeste.graph, candidate.boundary.elements.to.deploy))$no
  points(i, number.of.subgraphs[i], pch=19)
}
difftime(Sys.time(), start.time, units="mins")
#this procedure might take about 1 or 2 minutes, depending on which computer is used
#
#examine the results
class(number.of.subgraphs)
length(number.of.subgraphs)
summary(number.of.subgraphs)
number.of.subgraphs[1:100]

#plot the number of subgraphs (potential ecoregions) against the
#number of spatial links removed (or "candidate boundary elements" deployed).
plot(1:length(number.of.subgraphs), number.of.subgraphs,
     xlab="Candidate boundary elements deployed", ylab="Regions (or subgraphs)",
     pch=19, bty="n", cex.axis=1.5, cex.lab=1.5) 
abline(v=match(16, number.of.subgraphs), lty=3)
abline(h=16, lty=3)
abline(v=match(50000, number.of.subgraphs), lty=3)
abline(h=50000, lty=3)
abline(v=30829, lty=3)
abline(v=number.of.subgraphs[30829], col="red")
abline(v=35082 - beta.ranks.to.evaluate, lty=3)
abline(h=number.of.subgraphs[35082 - beta.ranks.to.evaluate], lty=3)

plot(1:length(number.of.subgraphs), number.of.subgraphs, xlim=c(0,6385), ylim=c(0,1600),
     xlab="Candidate boundary elements deployed", ylab="Regions (or subgraphs)",
     pch=19, bty="n", cex.axis=1.5, cex.lab=1.5) 
abline(v=35082 - beta.ranks.to.evaluate, lty=3)
abline(h=number.of.subgraphs[35082 - beta.ranks.to.evaluate], lty=3)

plot(1:length(number.of.subgraphs), number.of.subgraphs, xlim=c(400, 800), ylim=c(0,20),
     xlab="Candidate boundary elements deployed", ylab="Regions (or subgraphs)",
     pch=19, bty="n", cex.axis=1.5, cex.lab=1.5) 
abline(v=35082 - beta.ranks.to.evaluate, lty=3)
abline(v=match(pick.num.subgraphs, number.of.subgraphs), lty=3)
abline(h=number.of.subgraphs[35082 - beta.ranks.to.evaluate], lty=3)

#examine the number of regions (potential ecoregions or subgraphs) to evaluate against null models
number.of.subgraphs[35082 - beta.ranks.to.evaluate]

#save file with number of subgraphs, derived from the taxonomic (i.e, species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
write.table(number.of.subgraphs, file="04_Wombling/NumberSubgraphsSim.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="04_Wombling/NumberSubgraphsSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSims.txt", sep=",", row.names=F)

#read file with number of subgraphs, derived from Sorensen's index
#number.of.subgraphs <- read.table("NumberSubgraphsSor.txt", header=T, sep=",")
number.of.subgraphs <- read.table("04_Wombling/NumberSubgraphsSim.txt", header=T, sep=",")
#number.of.subgraphs <- read.table("NumberSubgraphsPhyloSor.txt", header=T, sep=",")
#number.of.subgraphs <- read.table("NumberSubgraphsPhyloSim.txt", header=T, sep=",")
head(number.of.subgraphs)
number.of.subgraphs <- number.of.subgraphs[,1]
number.of.subgraphs[1:5]
class(number.of.subgraphs)
length(number.of.subgraphs)
summary(number.of.subgraphs)
number.of.subgraphs[1:100]


##################################################################################################
# 11) Calculate superfluity, as defined by Oden et al. (1993, cited in the introduction).
# This step may take several hours or days, and you might skip to the end of this section
# to read a file with previously calculated superfluity values. Use code under 9.1 or 9.2 according
# to the version of package igraph and R you are using.
##################################################################################################

#define the percentiles of beta-diversity for which superfluity will be claculated,
#if you have not done so in section 8 (above)
beta.ranks.to.evaluate <- 35082 - round(35082*seq(0.05, 0.5, 0.05))

#define the number of regions (potentially ecoregions or "subgraphs") for which superfluity will be calculated
evaluation.number.of.subgraphs <- number.of.subgraphs[35082 - beta.ranks.to.evaluate]

#check that all evaluation points exist in the vector "number.of.subgraphs"
match(evaluation.number.of.subgraphs, number.of.subgraphs)  
sum(is.na(match(evaluation.number.of.subgraphs, number.of.subgraphs)))

#calculate superfluity
start.time <- Sys.time()
superfluity <- rep(NA, times=length(evaluation.number.of.subgraphs))
plot(min(evaluation.number.of.subgraphs), 0,
     xlim=c(min(evaluation.number.of.subgraphs), max(evaluation.number.of.subgraphs)),
     ylim=c(0, 200), 
     bty="n", xlab="Regions (or subgraphs)", ylab="Superfluity", cex.lab=1.5, cex.axis=1.5, type="n")
for (j in evaluation.number.of.subgraphs)
{
  print(paste(which(evaluation.number.of.subgraphs == j), "of", length(evaluation.number.of.subgraphs)))
  candidate.boundary.elements.to.delete <- which(r.obs.beta > (length(r.obs.beta) - match(j, number.of.subgraphs)))
  focal.graph <- delete_edges(Nordeste.graph, candidate.boundary.elements.to.delete)
  #focal.graph
  
  necessary.edge <- rep(NA, times=length(candidate.boundary.elements.to.delete))
  for(i in 1:length(candidate.boundary.elements.to.delete))
  {
    necessary.edge[i] <- ifelse(components(add_edges(focal.graph, cell.adj.char[candidate.boundary.elements.to.delete[i],]))$no < components(focal.graph)$no, 1, 0)
  }
  superfluity[which(evaluation.number.of.subgraphs==j)] <- sum(necessary.edge<1)/sum(necessary.edge>0)
  points(j, superfluity[which(evaluation.number.of.subgraphs==j)], pch=19)
}
difftime(Sys.time(), start.time, units="mins")
#this procedure might take about 14 minutes, depending on which computer is used
#
#examine the results
class(superfluity)
length(superfluity)
summary(superfluity)

#plot the results in linear scale
plot(evaluation.number.of.subgraphs, superfluity, pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Regions (or subgraphs)", ylab="Superfluity")

#plot the results using a semi-log graph
plot(evaluation.number.of.subgraphs, log(superfluity), pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Regions (or subgraphs)", ylab="Log (Superfluity)")

#plot the results using a log-log graph
plot(log(evaluation.number.of.subgraphs), log(superfluity), pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Log (Regions (or subgraphs))", ylab="Log (Superfluity)")

#write files with superfluity values, derived from the taxonomic (i.e., species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
#write.table(cbind(evaluation.number.of.subgraphs, superfluity), file="04_Wombling/SuperfluitySor.txt", sep=",", row.names=F)
#write.table(cbind(evaluation.number.of.subgraphs, superfluity), file="04_Wombling/SuperfluityPhyloSor.txt", sep=",", row.names=F)
write.table(cbind(evaluation.number.of.subgraphs, superfluity), file="04_Wombling/SuperfluitySim.txt", sep=",", row.names=F)
#write.table(cbind(evaluation.number.of.subgraphs, superfluity), file="04_Wombling/SuperfluityPhyloSim.txt", sep=",", row.names=F)

#read file with superfluity values, dderived from the taxonomic (i.e., species based) or phylogenetic
#versions of Sorensen's or Simpson's indices
#superfluity <- read.table("04_Wombling/SuperfluitySor.txt", header=T, sep=",")
#superfluity <- read.table("04_Wombling/SuperfluityPhyloSor.txt", header=T, sep=",")
#superfluity <- read.table("04_Wombling/SuperfluitySim.txt", header=T, sep=",")
#superfluity <- read.table("04_Wombling/SuperfluityPhyloSim.txt", header=T, sep=",")
head(superfluity)
class(superfluity)
length(superfluity)
summary(superfluity)
evaluation.number.of.subgraphs <- superfluity[,1]
superfluity <- superfluity[,2]


##################################################################################################
# 12) Examine and map the regions (potential ecoregions) or "subgraphs" obtained in section 10
##################################################################################################

##################################################################################################
# 12.1) Obtain and examine data on the regions
##################################################################################################

#select the number of regions (potential ecoregions) or subgraphs to map
number.of.subgraphs[35082 - beta.ranks.to.evaluate]
pick.num.subgraphs <- 4147
#make sure that the number you selected exists in the vector talling the number of subgraphs
match(pick.num.subgraphs, number.of.subgraphs)

#define the spatial links that should be removed to obtain the desired number of regions or subgraphs
candidate.boundary.elements.to.deploy <- which(r.obs.beta > (length(r.obs.beta) - match(pick.num.subgraphs, number.of.subgraphs)))

#remove the spatial links defined in the previous line of code and examine the results
modified.Nordeste.graph <- delete_edges(Nordeste.graph, candidate.boundary.elements.to.deploy)
components(modified.Nordeste.graph)$no
modified.Nordeste.graph
Nordeste.graph

#examine the resulting distribution of region size (or subgraph size)
SizeRegions <- components(modified.Nordeste.graph)$csize
histogram.subgraph.size <- hist(SizeRegions, breaks=seq(0.5,max(components(modified.Nordeste.graph)$csize)+0.5,1))
attributes(histogram.subgraph.size)
plot(histogram.subgraph.size$mids, histogram.subgraph.size$counts+1, 
     yaxt="n", log="xy", bty="n", type="h", pch=19, cex.axis=1.5, cex.lab=1.5,
     xlab="Region size (grid cells)", ylab="Regions", main=paste("Highest 50% CBE deployed (154188),", "total regions = ", pick.num.subgraphs))
axis(2, at=c(1,1+10^seq(0,4,1)), labels=c(0,10^seq(0,4,1)), cex.axis=1.5)
axis(2, at=c(2), labels=c(1), cex.axis=1.5)

#save file with region size
save(SizeRegions, file=paste("04_Wombling/SizeRegions", pick.num.subgraphs, ".R", sep=""))
#load file with region size
load(paste("04_Wombling/SizeRegions", pick.num.subgraphs, ".R", sep=""))

#obtain the coordinates of the center of the cells assigned to each region or subgraph
comp.coor <-  as.list (rep(NA, times=components(modified.Nordeste.graph)$no))
for(j in 1:components(modified.Nordeste.graph)$no)
{
  comp.coor[[j]] <- xyFromCell(Nordeste.mask.0, as.numeric(names(components(modified.Nordeste.graph)$membership))[components(modified.Nordeste.graph)$membership==j], spatial=FALSE)
}

#save the coordinates of the center of the cells in each region or subgraph
save(comp.coor, file=paste("04_Wombling/CompCoor_", pick.num.subgraphs, ".R", sep=""))

#load the coordinates of the center of the cells in each region or subgraph
load(paste("04_Wombling/CompCoor_", pick.num.subgraphs, ".R", sep=""))

#obtain cell number (or cell id) for each cell assigned to each region or subgraph
comp.cells <-  as.list (rep(NA, times=components(modified.Nordeste.graph)$no))
for(j in 1:components(modified.Nordeste.graph)$no)
{
  comp.cells[[j]] <- as.numeric(names(components(modified.Nordeste.graph)$membership))[components(modified.Nordeste.graph)$membership==j]
}

#create a raster of regions or subgraphs
Nordeste.mask.Regions <- Nordeste.mask.0
Nordeste.mask.Regions[] <- NA
Nordeste.mask.Regions[unlist(comp.cells)] <- rep(1:length(comp.cells), times=lapply(comp.cells, length))
Nordeste.mask.Regions
summary(Nordeste.mask.Regions)

#save the raster of regions or subgraphs
writeRaster(Nordeste.mask.Regions, paste("04_Wombling/Nordeste_, pick.num.subgraphs,"_Regions.grd", sep="))

##################################################################################################
# 12.2) Map the regions
##################################################################################################

#define and examine colors to map regions
number.of.colors <- 4147
subgraph.col <- colorRampPalette(brewer.pal(12,"Paired"))(number.of.colors)
subgraph.col <- sample(subgraph.col, length(subgraph.col))
subgraph.col[1:20]
plot(1:number.of.colors, col = subgraph.col, pch = 16, cex = 3)
length(subgraph.col)
subgraph.col <- rep(subgraph.col, length.out=pick.num.subgraphs)

#draw the map of regions at broad scale
plot(Nordeste.mask.Regions)
plot(Nordeste.mask.Regions, col=subgraph.col)
#plot regions, broad scale
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F)
#plot regions, narrower scale
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, xlim=c(-52,-43), ylim=c(-12,10))
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, xlim=c(-52,-43), ylim=c(-25,-12))
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, xlim=c(-43,-35), ylim=c(-12,10))
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, xlim=c(-43,-35), ylim=c(-25,-12))

#plot the mask, broad scale
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F)
#plot the mask, narrower scale
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, xlim=c(-52,-43), ylim=c(-12,10))
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, xlim=c(-52,-43), ylim=c(-25,-12))
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, xlim=c(-43,-35), ylim=c(-12,10))
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, xlim=c(-43,-35), ylim=c(-25,-12))
#plot the mask, even narrower scale

#create a matrix with the species composition for each region
SpeciesRegions <- matrix(NA, nrow=pick.num.subgraphs, ncol=length(brick.index.species.in.phylogeny))
for(i in 1:length(brick.index.species.in.phylogeny))
{
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
write.table(SpeciesRegions, "04_Wombling/SpeciesRegions", pick.num.subgraphs, ".txt", quote=T, sep=",")
write.table(SizeRegions, "04_Wombling/SizeRegions", pick.num.subgraphs, ".txt", quote=T, sep=",")

#read files with species composition of for each region
SpeciesRegions <- read.table("04_Wombling/SpeciesRegions", pick.num.subgraphs, ".txt", header=T, sep=",")
head(SpeciesRegions)
dim(SpeciesRegions)
