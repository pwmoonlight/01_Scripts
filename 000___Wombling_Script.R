
############################################################################################################################
############################################################################################################################
#
# INTRODUCTION TO THE SCRIPT "Wombling_SDMData_150arc_2017August04.R"
#
# A) What does this script do?
#
# This script identifies areas of Nordeste that could potentially be ecoregions, defined as as areas that are relatively
# homogeneous in plant species composition and encircled by boundaries that represent spatial discontinuities in plant species
# composition. To do so the script calculates observed beta-diversity between adjacent grid cells of 150 arc seconds across
# Nordeste, based on species distribution modes that estimate the geographic distributions of 891 species. The analytical
# method used in this script is known as "categorical wombling", described by:
# - Oden, N.L., Sokal, R.R., Fortin, M.J. and Hans Goebl. (1993) Categorical Wombling: Detecting regions of significant change
# in spatially located categorical variables. Geographical Analysis, Vol. 25, No. 4)
# - Fortin, M.J. and Drapeau, P. (1995) Delineation of Ecological Boundaries: Comparison of Approaches and Significance Tests.
# OIKOS 72: 323-332.
# See also:
# - Fortin, M. J & Dale, M. R. (2005). Spatial analysis: a guide for ecologists. Cambridge University Press.
#
# B) What is needed to run this script?
#
# If you want to run all processes in the script, you need the packages loaded in section 1 (below) and the
# following files, listed below their corresponding directory: 
#
# B.1)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Models)") #from Ivan's laptop
# setwd("J:/Jimenez/Nordeste_Biomes/Models") #from desktop 1ZTF at the Lehmann
# - raster files for all of the species distribution models, or a brick with all these files (see B.6 below)
#
# B.2)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Models_Summary_Best") #this is from Ivan's desktop
# setwd("J:/Jimenez/Nordeste_Biomes/Models_Summary_Best") #this is for desktop 1ZTF at the Lehmann
# - file summarizing the performance of species distribution models: "Summary_model_performance - 2015Sept30.csv"
#
# B.3)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Maps/Masks")# this is Ivan's directory
# setwd("J:/Jimenez/Nordeste_Biomes/Maps/Masks")#from desktop 1ZTF at the Lehmann
# - raster mask for Nordeste: "mask0_gadm_10kmBuffer_2015July23.tif"
#
# B.4) 
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Datasets") #from Ivan's laptop
# setwd("J:/Jimenez/Nordeste_Biomes/Datasets") #from desktop 1ZTF at the Lehmann
# - file with the phylogeny for the species included in the analysis:
# "intree.collapsed.dated.thorough.tre"
#
# B.5)
# List of species included in the analysis (all included in the Nordeste phylogeny)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Datasets") #from Ivan's laptop
# setwd("J:/Jimenez/Nordeste_Biomes/Datasets") #from desktop 1ZTF at the Lehmann
# "SpeciesInAnalysis_2016June21.txt"
#
# Running several of the processes in this script may take many hours or days. So you might want to have files with
# the results of several of those processes. If you decide to go that route, you would not need
# the raster files with the species distribution models, nor the file summarizing the performance of species distribution
# models. However, you would need the raster mask mentioned above (B.3), as well as at least some of the files listed below
# under their corresponding directory. In the heading of each script section you will find instructions as to when it is
# advisable to read each of these files.
#
# B.6)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Maps/Masks")# this is Ivan's directory
# setwd("J:/Jimenez/Nordeste_Biomes/Maps/Masks")#from desktop 1ZTF at the Lehmann
# - raster mask for Nordeste at a resolution of 150 arc seconds: "mask0_gadm_10kmBuffer_150arc_2015Oct27.tif"
#
# B.7)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Models/All_species_in_a_brick_BestModels") #from Ivan's laptop
# setwd("J:/Jimenez/Nordeste_Biomes/Models/All_species_in_a_brick_BestModels") #from desktop 1ZTF at the Lehmann
# - file that has the brick with the species distribution models of all species at a resolution of 150 arc seconds:
# "SDMb_150arc_2017June12.grd"
#
# B.8)
# setwd("C:/_transfer/Papers/Nordeste_Biomes/Datasets") #from Ivan's laptop
# setwd("J:/Jimenez/Nordeste_Biomes/Datasets") #from desktop 1ZTF at the Lehmann
# - definition of spatial links between grid cells, stored in the text file:
# "cell_adj_150arc.txt"
#
# - files with beta-diversity measured as taxonomic (i.e., species-based) and phylogenetic Sorensen's and Simpson's indexes:
# "ObsBetaSim_BestModels_150arc_2017June13", "ObsBetaSor_BestModels_150arc_2017June13",
# "ObsPhyloBetaSim_BestModels_150arc_2017June19.txt", "ObsPhyloBetaSor_BestModels_150arc_2017June19.txt"
# 
# - files with the ranks of beta-diversity values measured as taxonomic and phylogenetic Sorensen's and Simpson's indexes:
# "Rank_ObsBetaSim_BestModels_150arc_2017June13.txt", "Rank_ObsBetaSor_BestModels_150arc_2017June13.txt",
# "Rank_ObsPhyBetaSim_BestModels_150arc_2017June19.txt", "Rank_ObsPhyBetaSor_BestModels_150arc_2017June19.txt"
#
# - files with number of subgraphs, derived from taxonomic and phylogenetic Sorensen's and Simpson's indexes:
# "NumberSubgraphsSimBestModels_150arc_2017June13.txt", "NumberSubgraphsSorBestModels_150arc_2017June13.txt",
# "NumberSubgraphsPhyloSorBestModels_150arc_2017June13", "NumberSubgraphsPhyloSimBestModels_150arc_2017June13.txt", 
# 
# - file with superfluity values, derived from taxonomic and phylogenetic Sorensen's and Simpson's indexes:
# "SuperfluitySimBestModels_150arc_2017June13.txt", "SuperfluitySorBestModels_150arc_2017June13.txt", 
# "SuperfluityPhyloSorBestModels_150arc_2017June13.txt", "SuperfluityPhyloSimBestModels_150arc_2017June13.txt"
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
# 2.1) Create a raster mask with a resolution of 5*30 = 150 arc seconds resolution, which is 2.5 minutes
# or 0.04166667 X 0.04166667 degrees. You may skip this step (and go to 2.2) because the raster mask is
# already available.
############################################################################################################################

#create and plot the Nordeste raster mask that indicates the grid cells that have climate data,
Nordeste.mask.0 <- raster(list.files("03_Modelling/12_Thresholded_Models", pattern="*.tif", full.names = T)[[1]])
nordeste <- readOGR("000_GIS_LAYERS/nordeste.shp", layer="nordeste")
Nordeste.mask.0 <- crop(Nordeste.mask.0, nordeste)
Nordeste.mask.0 <- mask(Nordeste.mask.0, nordeste)
Nordeste.mask.0[Nordeste.mask.0[1:length(Nordeste.mask.0)]==1] <- 0
writeRaster(Nordeste.mask.0, "000_GIS_LAYERS/nordeste.tif", overwrite=T)
plot(Nordeste.mask.0)

#examine the properties of the mask
class(Nordeste.mask.0)
extent(Nordeste.mask.0)
res(Nordeste.mask.0) #this is the resolution of the mask
attributes(Nordeste.mask.0)
str(Nordeste.mask.0)

#Read in an earlier mask if created earlier
Nordeste.mask.0 <- raster("000_GIS_LAYERS/nordeste.tif")

#examine the total number of grid cells, the frequency of grid cells
#that are NA and not-NA, and determine cell number (or cell ID) of
#the cells that are not NA
ncell(Nordeste.mask.0) #total number of cells
length(Nordeste.mask.0[is.na(Nordeste.mask.0)]) # of NA cells
length(Nordeste.mask.0[!is.na(Nordeste.mask.0)]) #number of not-NA cells
#create a vector with the cell number (or cell ID) of cells that are not-NA
no.na.cells <- (1:ncell(Nordeste.mask.0))[!is.na(extract(Nordeste.mask.0, 1:ncell(Nordeste.mask.0)))]
length(no.na.cells)

############################################################################################################################
# 2.3) Define adjacent grid cells. You might want to skip to the end of this section, and use previously defined spatial
# links between grid cells, stored in the text file "cell_adj.txt".  
############################################################################################################################

#define adjacent grid cells, but only for cells that are not-NA in
#the raster with the species distribution model,
#the code below uses "rook's move neighbors" to define adjacent cells
cell.neig <- adjacent(Nordeste.mask.0, cells=no.na.cells, directions=4, pairs=T, target=no.na.cells, sorted=T)
cell.neig[1:10,]
dim(cell.neig)
#note that each pair of neighboring cells shows up twice in the
#matrix "cell.neig"; the following code removes duplicates 
flag.pairs.to.retain <- (cell.neig[,1] - cell.neig[,2])<0
length(flag.pairs.to.retain)
sum(flag.pairs.to.retain)
flag.pairs.to.retain[1:100]
cell.adj <- cell.neig[flag.pairs.to.retain,]
cell.adj[1:10,]
dim(cell.adj)

#plot some of the links between grid cells, but note that not all links need to appear
#in the plot, because an arbitrary part of the matrix cell.adj is selected
from.coor <- xyFromCell(Nordeste.mask.0, cell.adj[1999:4500,1], spatial=FALSE)
to.coor <- xyFromCell(Nordeste.mask.0, cell.adj[1999:4500,2], spatial=FALSE)
plot(Nordeste.mask.0, xlim=c(-51.05919,-32.36665), ylim=c(-10,-1.05))
#plot the center of grid cells
points(rbind(from.coor, to.coor), pch=19, cex=0.9, col="blue")
#plot the links
arrows(from.coor[,1], from.coor[,2], to.coor[,1], to.coor[,2], length = 0, code = 2, col="red")
#examine cells with no adjacent cells
no.na.cells[is.na(match(no.na.cells, unique(c(cell.adj[,1], cell.adj[,2]))))]
no.adj.coor <- xyFromCell(Nordeste.mask.0, no.na.cells[is.na(match(no.na.cells, unique(c(cell.adj[,1], cell.adj[,2]))))], spatial=FALSE)
plot(Nordeste.mask.0, xlim=c(-51.05919,-32.36665), ylim=c(-22.9,-1.05))
points(no.adj.coor, pch=19, cex=0.5, col="blue")
plot(Nordeste.mask.0, xlim=c(-51.05919,-32.36665), ylim=c(-22.9,-1.05))
points(no.adj.coor, pch=19, cex=0.5, col="blue")
plot(Nordeste.mask.0, xlim=c(-51.05919,-32.36665), ylim=c(-22.9,-1.05))
points(no.adj.coor, pch=19, cex=0.5, col="blue")

#save in a file the links between grid cells
dir.create("04_Wombling/Cell_Links", showWarnings=F)
write.table(cell.adj, file="04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",")

#read text file with previously defined links between grid cells
cell.adj <- read.table("04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",", header=T)


#############################################################################
# 3) Create an R object of class "brick" with all species distribution models
#############################################################################

#read and examine the file with information on the species to be included in the analysis,
#and the best performing SDM for each of those species:

Sum.SDM.Perf <- read.table("03_Modelling/CBI_results.csv", header=T, sep=",")
head(Sum.SDM.Perf)

#Remove any species form the analysis with CBI values < 0.5
Sum.SDM.Perf <- Sum.SDM.Perf[Sum.SDM.Perf[,7] > 0.5,]

#Remove any species form the analysis with AUC values < 0.75
Sum.SDM.Perf <- Sum.SDM.Perf[Sum.SDM.Perf[,13] > 0.7,]

species.in.analysis <- Sum.SDM.Perf[,1]
species.in.analysis

#save in a file the names of the species included in the analysis
write.table(species.in.analysis, file="04_Wombling/species_in_analysis.txt", quote=T, sep=",", row.names=F, col.names=F)

#read the raster files with SDMs for the species to be included in the analysis, make sure
#to use the correct directories within the loop

species.level.dir <- lapply(1:length(species.in.analysis), function(x){raster(paste("03_Modelling/12_Thresholded_Models/", as.character(species.in.analysis)[[x]], ".tif", sep=""))})

#create an R object of class "brick" with species distribution models for all species,
#first create a list of the names of raster objects (with species distribution models)
#that will be included in the "brick":
SDM.raster.names <- vector(mode = "list", length = length(species.level.dir))
for(i in 1:length(species.level.dir[]))
{
  SDM.raster.names[[i]]  <- species.level.dir[[i]]@file@name
}
#now create the "brick"
lapply(SDM.raster.names, eval) #checking this bit of code works
length(lapply(SDM.raster.names, eval)) #checking it is of the right length
SDM.b <- brick(lapply(SDM.raster.names, eval))
SDM.b <- crop(SDM.b, nordeste)
SDM.b <- mask(SDM.b, nordeste)
#check the result
class(SDM.b)
res(SDM.b)
dim(SDM.b)
SDM.b

#save in a file the brick with the species distribution models of all species in the analysis that has 30 arc seconds resolution
rasterOptions(datatype="LOG1S") #specify datatype (true/false or presence/absence) to efficiently save the brick file  
writeRaster(SDM.b, filename="04_Wombling/SDM_brick.grd", bandorder='BIL')

#save in a file the brick with the species distribution models of all species in the analysis
rasterOptions(datatype="LOG1S") #specify datatype (true/false or presence/absence) to efficiently save the brick file  

#read the file that has the brick with the species distribution models of all species
SDM.b <- brick("04_Wombling/SDM_brick.grd")
SDM.b
str(SDM.b)
nlayers(SDM.b)
names(SDM.b)


############################################################################################################################
# 4) Read and examine the phylogeny for the species in the analysis (the "Nordeste phylogeny")
############################################################################################################################

{
  
# read and examine attributes of the Nordeste phylogeny
#Nordeste.tree <- read.tree("intree.collapsed.dated.thorough.tre")
#str(Nordeste.tree)
#attributes(Nordeste.tree)
#hist(Nordeste.tree$edge.length, breaks=100)
#summary(Nordeste.tree$edge.length)
#Nordeste.tree$edge
#Nordeste.tree$tip.label
#edgelabels()
#Nordeste.tree$root.edge
#node.depth.edgelength(Nordeste.tree)
#max(node.depth.edgelength(Nordeste.tree))
#summary(node.depth.edgelength(Nordeste.tree))
#hist(node.depth.edgelength(Nordeste.tree), breaks=100)

# change two phylogeny tip labels, so that all names are the same in the phylogeny and the brick with the species distribution models
#Nordeste.tree$tip.label[Nordeste.tree$tip.label == "Antigonon_guatimalense"] <- "Antigonon_guatemalense"
#Nordeste.tree$tip.label[Nordeste.tree$tip.label == "Stemmadenia_donnellsmithii"] <- "Stemmadenia_donnell-smithii"

# plot the phylogeny
# first a simple option
#plot(Nordeste.tree, show.tip.label=FALSE)
# another version of the plot
#plot(Nordeste.tree, type="fan", cex=0.2, label.offset = 1,
#     rotate.tree=196, open.angle=15)
#axisPhylo(1, cex=1)
# yet another version of the plot
#par(mar= c(1, 1, 1, 1) + 0.1))
#win.graph(20, 20, 10)
#plot(Nordeste.tree, type="fan", cex=0.2, label.offset = 1,
#     rotate.tree=196, open.angle=15)
#axisPhylo(1, cex=1.5)
#}

############################################################################################################################
# 5) Determine the subset of species included in the brick with the species distribution models that are
# represented in the phylogeny.
############################################################################################################################


#species.in.analysis <- read.table("04_Wombling/species_in_analysis.txt", header=T, sep=",")
#class(species.in.analysis)
#class(species.in.analysis[,1])
#class(species.in.analysis[,2])
#dim(species.in.analysis)
#species.in.analysis[1:5,]

#species.genus <- 
#  unlist(lapply(strsplit(as.vector(species.in.analysis[,1]), split=c(" "), fixed = T), function(x) paste(x[1], x[2], sep="_")))
#class(species.genus)
#species.genus[1:5]
#length(species.genus)

#create an index to select the species included in the brick that are represented in the phylogeny
#brick.index.species.in.phylogeny <- which(!is.na(match(species.genus, Nordeste.tree$tip.label)))}


############################################################################################################################
# 6) Calculate beta-diversity for the spatial links created in section 2 above.
############################################################################################################################

############################################################################################################################
# 6.1) Calculate taxonomic (i.e., species-based) beta-diversity using Simpson's and Sorensen's indices
############################################################################################################################

brick.index.species.in.phylogeny <- names(SDM.b)

# Obtain, for all pairs of adjacent grid cells, terms to calculate beta-diversity.
# First, using a matrix indicate:
# species are unique to one of the grid cells, represented as "0",
# species unique to the other grid cell, represented as "Inf",
# and shared species, represented as a "1".
obs.beta.terms <-  SDM.b[cell.adj[,1]][,brick.index.species.in.phylogeny] / SDM.b[cell.adj[,2]][,brick.index.species.in.phylogeny]
dim(obs.beta.terms)

# Next, for all pairs of adjacent grid cells, calculate the following terms, see page 254 in
# Legendre and Legendre (2010, Numerical Ecology, Second Edition):
# the number of shared species 
a <- rowSums(obs.beta.terms>0.5 & obs.beta.terms<1.5, na.rm=T)
summary(a)
# the number of unique species in one of the grid cells 
b <- rowSums(obs.beta.terms<0.5, na.rm=T)
summary(b)
# the number of unique species in the other grid cell
c <- rowSums(is.infinite(obs.beta.terms))
summary(c)

# Calculate ecological distance using Sorensen's index
# (see page 286, equation 7.56 in Legendre and Legendre (2010, Numerical Ecology, Second Edition):
obs.beta.sor <- (b+c)/(2*a + b + c)
obs.beta.sor[which(is.na(obs.beta.sor))] <- 0
summary(obs.beta.sor)

# Calculate ecological distance using Simpson's index
# (see page 2230, equation 2 in Mouillot et al. (2013, Journal of Biogeography 40: 2228–2237):
obs.beta.sim <- pmin(b,c)/(a + pmin(b,c))
obs.beta.sim[which(is.na(obs.beta.sim))] <- 0
summary(obs.beta.sim)

# Graphically compare the Sorensen's and Simpson's indexes
plot(obs.beta.sor, obs.beta.sim, bty="n", cex.axis=1.5, cex.lab=1.5)
abline(0, 1, col="red")

write.table(obs.beta.sor, file="04_Wombling/ObsBetaSor.txt", sep=",")
write.table(obs.beta.sim, file="04_Wombling/ObsBetaSim.txt", sep=",")

#read file with beta-diversity measured as taxonomic (i.e., species based) Sorensen's or Simpson's indices
obs.beta.sor <- read.table("04_Wombling/ObsBetaSor.txt", header=T, sep=",")
obs.beta.sim <- read.table("04_Wombling/ObsBetaSim.txt", header=T, sep=",")
dim(obs.beta.sor)
head(obs.beta.sor)
obs.beta.sor[1:5,]
obs.beta.sor <- obs.beta.sor[,1]
obs.beta.sor[1:10]
#
dim(obs.beta.sim)
head(obs.beta.sim)
obs.beta.sim[1:5,]
obs.beta.sim <- obs.beta.sim[,1]
obs.beta.sim[1:10]


############################################################################################################################
# 6.2) Calculate phylogenetic beta-diversity using Simpson's and Sorensen's indices
############################################################################################################################

#{define the species list
#species.list <- species.genus[brick.index.species.in.phylogeny]

#create vectors needed in the loop to save results
#obs.phylo.beta.sim <- rep(NA, times=nrow(cell.adj))
#obs.phylo.beta.sor <- rep(NA, times=nrow(cell.adj))

#for (i in 1: nrow(cell.adj))
#{
#  #determine which of the species included in the Nordeste phylogeny occur in the first grid cell
#  species.1 <- species.list[as.logical(SDM.b.150arc[cell.adj[i,1]][,brick.index.species.in.phylogeny])]
#  #identify the edges of the Nordeste phylogeny that join the species occurring in the first grid cell to their most recent common ancestor 
#  edges.1 <- which.edge(Nordeste.tree, species.1)
#  #determine which of the species included in the Nordeste phylogeny occur in the second grid cell
#  species.2 <- species.list[as.logical(SDM.b.150arc[cell.adj[i,2]][,brick.index.species.in.phylogeny])]
#  #determine the edges of the Nordeste phylogeny that join the species occurring in the second grid cell to their most recent common ancestor 
#  edges.2 <- which.edge(Nordeste.tree, species.2)
#  #create a list of the species included in the Nordeste phylogeny that occur in one or both of the grid cells,
#  #avoiding duplication of species names
#  species.tot <- unique(c(species.1, species.2))
#  #determine the edges of the Nordeste phylogeny that join the species occurring in one or both of the grid cells to their most recent common ancestor
#  edges.tot <- which.edge(Nordeste.tree, species.tot)
#  #Calculate phylogenetic betadiversity, following equations 1-17 in Leprieur et al. (2012, Quantifying phylogenetic beta diversity: distinguishing
#  #between ‘true’ turnover of lineages and phylogenetic diversity gradients. PLoS ONE 7(8): e42760. doi:10.1371/journal.pone.0042760).
#  PDtot <- sum(Nordeste.tree$edge.length[unique(c(edges.1, edges.2))])
#  PD1 <- sum(Nordeste.tree$edge.length[edges.1])
#  PD2 <- sum(Nordeste.tree$edge.length[edges.2])
#  a <- PD1+PD2-PDtot
#  b <- PDtot-PD1
#  c <- PDtot-PD2
#  obs.phylo.beta.sim[i] <- min(b,c)/(a+min(b,c))
#  obs.phylo.beta.sor[i] <- (b+c)/(2*a + b + c)
#}

#examine results
#length(obs.phylo.beta.sim)
#head(obs.phylo.beta.sim)
#summary(obs.phylo.beta.sim)
#length(obs.phylo.beta.sor)
#head(obs.phylo.beta.sor)
#summary(obs.phylo.beta.sor)

#graphically compare the Sorensen's and Simpson's phylogenetic indexes
#plot(obs.phylo.beta.sor, obs.phylo.beta.sim, bty="n", cex.axis=1.5, cex.lab=1.5)
#abline(0, 1, col="red")

#plot(sort(obs.phylo.beta.sim))
#plot(sort(obs.phylo.beta.sor))

#save file with phylogenetic beta-diversity as Sorensen's index
#write.table(obs.phylo.beta.sor, file="ObsPhyloBetaSor_BestModels_150arc_2017June19.txt", sep=",")

#read file with phylogenetic beta-diversity as Sorensen's index
#obs.phylo.beta.sor <- read.table("ObsPhyloBetaSor_BestModels_150arc_2017June19.txt", header=T, sep=",")
#dim(obs.phylo.beta.sor)
#head(obs.phylo.beta.sor)
#obs.phylo.beta.sor[1:5,]
#obs.phylo.beta.sor <- obs.phylo.beta.sor[,1]
#obs.phylo.beta.sor[1:10]

#save file with phylogenetic beta-diversity as Simpson's index
#write.table(obs.phylo.beta.sim, file="ObsPhyloBetaSim_BestModels_150arc_2017June19.txt", sep=",")

#read file with phylogenetic beta-diversity as Simpson's index
#obs.phylo.beta.sim <- read.table("ObsPhyloBetaSim_BestModels_150arc_2017June19.txt", header=T, sep=",")
#dim(obs.phylo.beta.sim)
#head(obs.phylo.beta.sim)
#obs.phylo.beta.sim[1:5,]
#obs.phylo.beta.sim <- obs.phylo.beta.sim[,1]
#obs.phylo.beta.sim[1:10]
#}

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

r.obs.beta <- rank(obs.beta, ties.method="random")
class(r.obs.beta)
length(r.obs.beta)
summary(r.obs.beta, digits=6)
max(r.obs.beta)
length(unique(r.obs.beta))
r.obs.beta[1:100]

#define the percentiles of beta-diversity that will be evaluated against the null model
beta.ranks.to.evaluate <- 35082 - round(35082*seq(0.05, 0.5, 0.05))

plot(r.obs.beta, obs.beta, bty="n", cex.axis=1.5, cex.lab=1.5)
abline(v= beta.ranks.to.evaluate, lty=3)
abline(h= obs.beta[match(beta.ranks.to.evaluate, r.obs.beta)], lty=3)
#
plot(r.obs.beta, obs.beta, bty="n", cex.axis=1.5, cex.lab=1.5, xlim=c(6000,35082), ylim=c(0,0.1))
abline(v= beta.ranks.to.evaluate, lty=3)
abline(h= obs.beta[match(beta.ranks.to.evaluate, r.obs.beta)], lty=3)

#save ranks of observed beta values: useful because the raking has a random
#component to deal with ties, thus fully consistent results are unlikely
#using different iterations of the ranking
write.table(r.obs.beta, file="04_Wombling/Rank_ObsBetaSor.txt", row.names=F)
write.table(r.obs.beta, file="04_Wombling/Rank_ObsBetaSim.txt", row.names=F)
#write.table(r.obs.beta, file="04_Wombling/Rank_ObsPhyBetaSor.txt", row.names=F)
#write.table(r.obs.beta, file="04_Wombling/Rank_ObsPhyBetaSim.txt", row.names=F)

#read ranks of observed beta values
r.obs.beta <- read.table("04_Wombling/Rank_ObsBetaSor.txt", header=T, sep=",")
r.obs.beta <- read.table("04_Wombling/Rank_ObsBetaSim.txt", header=T, sep=",")
#r.obs.beta <- read.table("04_Wombling/Rank_ObsPhyBetaSor.txt", header=T, sep=",")
#r.obs.beta <- read.table("04_Wombling/Rank_ObsPhyBetaSim", header=T, sep=",")
class(r.obs.beta)
class(r.obs.beta[,1])
r.obs.beta[1:100,1]
dim(r.obs.beta)
head(r.obs.beta)
r.obs.beta <- as.numeric(r.obs.beta[,1])
class(r.obs.beta)
length(r.obs.beta)
summary(r.obs.beta, digits=6)
max(r.obs.beta)
length(unique(r.obs.beta))
r.obs.beta[1:100]


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
  #points(i, number.of.subgraphs[i], pch=19)
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
write.table(number.of.subgraphs, file="04_Wombling/NumberSubgraphsSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSor.txt", sep=",", row.names=F)
#write.table(number.of.subgraphs, file="NumberSubgraphsPhyloSims.txt", sep=",", row.names=F)

#read file with number of subgraphs, derived from Sorensen's index
number.of.subgraphs <- read.table("04_Wombling/NumberSubgraphsSor.txt", header=T, sep=",")
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
  print(paste(which(evaluation.number.of.subgraphs == j), "of", length(evaluation.number.of.subgraphs), "at", Sys.time()))
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
write.table(cbind(evaluation.number.of.subgraphs, superfluity), file="04_Wombling/SuperfluitySor.txt", sep=",", row.names=F)
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
save(SizeRegions, file=paste("04_Wombling/SizeRegions_SIM_", pick.num.subgraphs, ".R", sep=""))
save(SizeRegions, file=paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".R", sep=""))
#load file with region size
load(paste("04_Wombling/SizeRegions_SIM_", pick.num.subgraphs, ".R", sep=""))
load(paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".R", sep=""))

#obtain the coordinates of the center of the cells assigned to each region or subgraph
comp.coor <-  as.list (rep(NA, times=components(modified.Nordeste.graph)$no))
for(j in 1:components(modified.Nordeste.graph)$no)
{
  comp.coor[[j]] <- xyFromCell(Nordeste.mask.0, as.numeric(names(components(modified.Nordeste.graph)$membership))[components(modified.Nordeste.graph)$membership==j], spatial=FALSE)
}

#save the coordinates of the center of the cells in each region or subgraph
save(comp.coor, file=paste("04_Wombling/CompCoor_SIM_", pick.num.subgraphs, ".R", sep=""))
save(comp.coor, file=paste("04_Wombling/CompCoor_SOR_", pick.num.subgraphs, ".R", sep=""))

#load the coordinates of the center of the cells in each region or subgraph
load(paste("04_Wombling/CompCoor_SIM_", pick.num.subgraphs, ".R", sep=""))
load(paste("04_Wombling/CompCoor_SOR_", pick.num.subgraphs, ".R", sep=""))

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
writeRaster(Nordeste.mask.Regions, paste("04_Wombling/Nordeste_SIM_", pick.num.subgraphs,"_Regions.grd", sep=""))
writeRaster(Nordeste.mask.Regions, paste("04_Wombling/Nordeste_SOR_", pick.num.subgraphs,"_Regions.grd", sep=""))

##################################################################################################
# 12.2) Map the regions
##################################################################################################

#define and examine colors to map regions
number.of.colors <- pick.num.subgraphs
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
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, ylim=c(-12,10))
plot(Nordeste.mask.Regions, col=subgraph.col, useRaster=T, legend=F, ylim=c(-25,-12))

#plot the mask, broad scale
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F)
#plot the mask, narrower scale
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, ylim=c(-12,10))
plot(Nordeste.mask.0, col="white", useRaster=T, legend=F, ylim=c(-25,-12))
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
write.table(SpeciesRegions, paste("04_Wombling/SpeciesRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
write.table(SpeciesRegions, paste("04_Wombling/SpeciesRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
write.table(SizeRegions, paste("04_Wombling/SizeRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")
write.table(SizeRegions, paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), quote=T, sep=",")

#read files with species composition of for each region
SpeciesRegions <- read.table(paste("04_Wombling/SpeciesRegions_SIM_", pick.num.subgraphs, ".txt", sep=""), header=T, sep=",")
SpeciesRegions <- read.table(paste("04_Wombling/SpeciesRegions_SOR_", pick.num.subgraphs, ".txt", sep=""), header=T, sep=",")
head(SpeciesRegions)
dim(SpeciesRegions)
