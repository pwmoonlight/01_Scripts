
####################################################################################################################################
####################################################################################################################################
#
# INTRODUCTION TO THE SCRIPT "Nicaragua_clustering_150arc_2017August21.R"
#
# A) What does this script do?
# 
# The code in this script performs cluster analysis of the presence-absence matrix representing the species composition of Nicaragua
# regions defined by categorical wombling (Oden et al. 1993. Geographical Analysis, Vol. 25, No. 4) at a 150 arc sec scale, based on
# species distribution models that estimate the geographic distributions of 786 plant species. The Nicaragua regions identified by
# categorical wombling could potentially be ecoregions, defined as areas that are relatively homogeneous in plant species composition
# and encircled by boundaries that represent spatial discontinuities in plant species composition. Therefore, the cluster analysis
# performed in this script is primarily of interest because a given ecoregion might be composed of multiple spatially separated areas
# that, nonetheless, are characterized by similar species or phylogenetic composition. The clustering approach used here avoids bias
# stemming from the order of regions in the original presence–absence matrix by shuffling the matrix row order and constructing a
# consensus dendrogram (Dapporto et al. 2013. Ecography, 36: 1070–1075). This approach is applied to gauge support for different
# clusters using multiscale bootstrap (see Dapporto et al. 2013).
#
# B) What is needed to run this script?
#
# To run this script you need the packages loaded in section 1 (below) and the following files, listed below their corresponding
# directory: 
#
# B.1)
# setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets")
# setwd("C:/Users/rvisitor/Desktop/clustering")
# - presence-absence matrix representing the species composition of Nicaragua regions defined by categorical wombling; there are
#   several matrices for any particular wombling analysis, by example, the matrix for 1277 regions defined based on species 
#   beta-diversity as measured by Simpson's index is in a file named: "SpeciesRegions1277_TaxSim_2017June25.txt"
#
# B.2)
# setwd("C:/_transfer/Papers/Nicaragua_Biomes/Maps/Regions") #from Ivan's laptop
# - a raster mapping the regions defined by categorical wombling; again, there are several rasters for any particular wombling
#   analysis, by example, the raster file for 1277 regions defined according to species beta-diversity as measured by Simpson's index
#   is: "Nicaragua_1277Regions_150arc_SimBestModels_2017July14.grd"
#
# B.3)
# setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets")
# - a text file with information on region size, measured as number of grid cells. As before, there are several rasters for any
# particular wombling  analysis, by example, the file for 1277 regions defined according to species beta-diversity as measured
# by Simpson's index is: "SizeRegions1277_TaxSim_2017June25.txt"
#
# B.4)
# setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets") #from Ivan's laptop
# setwd("J:/Jimenez/Nicaragua_Biomes/Datasets") #from desktop 1ZTF at the Lehmann
# R object (a list) with the coordinates of the center of the cells assigned to each region. As an example, the file for 1277 regions
# defined according to species beta-diversity as measured by Simpson's index is: "CompCoor_1277_150arc_SimBestModels_2017June13.RData"
#
# B.5)
# Performing multi-scale bootstrap of the consensus dendrogram may take many hours or several days. So you might want to save files
# with consensus dendrograms and the respective multi-scale bootstrap analyses. Note that different runs of the script may result 
# in shomewhat different consensus dendrograms, and the multi-scale bootstrap analysis can only be meaningfully used in relation to
# a particular consensus dendrogram. As an example, below are the directories and names of the files with the 50% consensus tree and
# respecitve multi-scale bootstrap analyses for 1277 regions defined according to species beta-diversity as measured by Simpson's index:
# setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets")
# setwd("J:/Jimenez/Nicaragua_Biomes/Datasets")
# "SpeciesRegions1277_TaxSim_dendrogram50_2017July3.RData"
# "MultiBootDendrogram50_TaxSim_Regions1277_2017June24.txt"
#
####################################################################################################################################
####################################################################################################################################



### Prepare the working space
### -------------------------

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

library(raster)
library(sp)
library(recluster)
library(RColorBrewer)
library(phytools)

dir.create("06_Clustering", showWarnings = F)
dir.create("06_Clustering/Dendrograms", showWarnings = F)

####################################################################################################################################
# 2) Open and examine the presence-absence matrix representing species composition for each region identified by wombling,
# regions are rows and species columns, 0 = absence and 1 = presence.
####################################################################################################################################

pick.num.subgraphs <- 18318
pick.index <- "SOR"
 
sp.comp <- read.table(paste("04_Wombling/SpeciesRegions_", pick.index, "_", pick.num.subgraphs,".txt", sep=""), header=T, sep=",", row.names=1)
#convert data frame to matrix
sp.comp <- as.matrix(sp.comp)
dim(sp.comp)

load(paste("04_Wombling/SizeRegions_", pick.index,"_", pick.num.subgraphs, ".R", sep=""))

removals <- which(SizeRegions==1)
remainers <- 1:length(SizeRegions)
remainers <- remainers[-removals]
sp.comp <- sp.comp[-removals,]
dim(sp.comp)

####################################################################################################################################
# 3) Run cluster analysis. See: http://www.ecography.org/appendix/ecog-00444, the supplementary material for:
# Dapporto, L., Ramazzotti, M., Fattorini, S., Talavera, G., Vila, R. and Dennis, R. L. H. (2013), recluster:
# an unbiased clustering procedure for beta-diversity turnover. Ecography, 36: 1070–1075. )
####################################################################################################################################

sp.comp <- sp.comp[-which(rowSums(sp.comp)==0),]
dim(sp.comp)

#calculate dissimilarity matrix
#dismat <- recluster.dist(sp.comp, dist="sorensen")
dismat <- recluster.dist(sp.comp, dist="simpson")

#examine frequency o zero values and ties
recluster.hist(dismat)

#compute and plot the 100% consensus dendrogram (showing clusters occurring in all resampled dendrograms);
#dendrograms are resampled by reordering sites in the distance matrix ("dismat")
dendrogram.100 <- recluster.cons(dismat, tr = 100, p = 1, dist = "simpson") 
plot(dendrogram.100$cons)
plot(dendrogram.100$cons, type="fan", show.tip.label=F)
plot(dendrogram.100$cons, type="radial", show.tip.label=F)

#examine effect of row order bias on each node of the original dendrogram
original.dendrogram.node.strength <- recluster.node.strength(sp.comp, dist = "simpson", tr = 100, levels = 6)
plot.phylo(original.dendrogram.node.strength$tree) 
#plot.phylo(original.dendrogram.node.strength$tree, y.lim=c(1,100))
#plot.phylo(original.dendrogram.node.strength$tree, y.lim=c(1000,1200))
#plot.phylo(original.dendrogram.node.strength$tree, type="fan", show.node.label =T)
nodelabels(round(original.dendrogram.node.strength$result))

#calculate 50% consensus dendrogram, showing clusters occurring in at least half the resampled dendrograms (i.e., p=0.5);
#dendrograms are resampled by reordering sites in the distance matrix ("dismat") 
dendrogram.50 <- recluster.cons(dismat, tr = 100, p = 0.5, dist = "simpson")
str(dendrogram.50$cons) 
plot(dendrogram.50$cons)
plot(dendrogram.50$cons, type="fan", show.tip.label=F)
#plot(dendrogram.50$cons, type="radial", show.tip.label=F)
axisPhylo(1)
#to understand how branch lengths are calculated in the consensus dendrogram,
#see: http://blog.phytools.org/2011/03/for-fun-least-squares-phylogeny.html

#save the 50% consensus dendrogram
#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets")
#setwd("J:/Jimenez/Nicaragua_Biomes/Datasets")
save(dendrogram.50, file=paste("06_Clustering/Dendrograms/SpeciesRegions", pick.num.subgraphs,"_Tax",pick.index,"_dendrogram50_removals.RData", sep=""))
load(paste("06_Clustering/Dendrograms/SpeciesRegions", pick.num.subgraphs,"_Tax",pick.index,"_dendrogram50_removals.RData", sep=""))

#perform multi-scale bootstrap of the 50% consensus dendrogram, by constructing 1000 (bootstrapped) consensus trees,
#each made up of 100 resampled trees; the code below uses 10 bootstrap scales, in a sequence from x1 to x19 every x2;
#this may take a while, so it would be useful to save results
#multi.boot.dendrogram.50 <- recluster.multi(dendrogram.50$cons, sp.comp, tr=20, p=0.5, boot=100, levels=10, step=2)

ptm <- Sys.time()
multi.boot.dendrogram.50 <- recluster.multi(dendrogram.50$cons, sp.comp, tr=100, p=0.5, boot=1000, levels=10, step=2)
class(multi.boot.dendrogram.50)
Sys.time() - ptm

write.table(multi.boot.dendrogram.50, file=paste("06_Clustering/Dendrograms/SpeciesRegions", pick.num.subgraphs,"_Tax", pick.index,"_MultiBootDendrogram50_removals.txt", sep=""), sep=",")
multi.boot.dendrogram.50 <- read.table(paste("06_Clustering/Dendrograms/SpeciesRegions", pick.num.subgraphs,"_Tax", pick.index,"_MultiBootDendrogram50_removals.txt", sep=""), header=T, sep=",")

class(multi.boot.dendrogram.50)
dim(multi.boot.dendrogram.50)
head(multi.boot.dendrogram.50)

#
max.boot.support <- apply(multi.boot.dendrogram.50, 1, max)
#
rows.90.boot.support <- as.vector(which(max.boot.support>=90))
nodes.90.boot.support <- rows.90.boot.support + length(dendrogram.50$cons$tip.label)
length(nodes.90.boot.support) #number of strongly supported nodes
dendrogram.50$cons$Nnode #total number of nodes
#
scale.nodes.90.boot.support <- apply(multi.boot.dendrogram.50[rows.90.boot.support,], 1, function(x) min(which(x>90)))
#
height.nodes.dendrogram50 <- round(max(node.depth.edgelength(dendrogram.50$cons)) - node.depth.edgelength(dendrogram.50$cons), digits=6)
height.nodes.dendrogram50
height.nodes.90.boot.support <- height.nodes.dendrogram50[nodes.90.boot.support]


hist(height.nodes.dendrogram50, breaks=seq(0,max(height.nodes.dendrogram50)+0.01,0.01),
	ylim=c(0, 500), xlab="Node height", ylab="Node frequency", bty="n", cex.lab=1.5, cex.axis=1.5)
par(new=T)
hist(height.nodes.90.boot.support, breaks=seq(0,max(height.nodes.dendrogram50)+0.01,0.01),
	ylim=c(0, 500), col="gray70", xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
legend("topright", c("All nodes", "Strongly supported nodes"), fill=c("white", "gray70"), cex=1.5)
#abline(v=quantile(height.nodes.90.boot.support, probs = seq(0.25, 1, 0.25)), lty=3)

hist(height.nodes.dendrogram50, breaks=seq(0,max(height.nodes.dendrogram50)+0.01,0.01),
	ylim=c(0, 2000), xlab="Node height", ylab="Node frequency", bty="n", cex.lab=1.5, cex.axis=1.5)
par(new=T)
hist(height.nodes.90.boot.support, breaks=seq(0,max(height.nodes.dendrogram50)+0.01,0.01),
	ylim=c(0, 2000), col="gray70", xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
legend("topright", c("All nodes", "Strongly supported nodes"), fill=c("white", "gray70"), cex=1.5)
#abline(v=quantile(height.nodes.90.boot.support, probs = seq(0.25, 1, 0.25)), lty=3)



boxplot(height.nodes.dendrogram50, height.nodes.90.boot.support,
	ylab="Node height", names=c("All nodes","Strongly supported nodes"),
	cex.axis=1.5, cex.lab=1.5, frame.plot=F)

plot(height.nodes.90.boot.support, scale.nodes.90.boot.support,
	xlab="Height of strongly supported nodes", ylab="Bootstrap scale at which node is strongly suppoted",
	pch=19, bty="n", cex.lab=1.5, cex.axis=1.5)



#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Figures")
pdf(file = paste("06_Clustering/Dendrograms/Dendrogram50_Tax", pick.index,"_Regions", pick.num.subgraphs,"_removals.pdf", sep=""), width=7, height=0.15*length(dendrogram.50$cons$tip.label))
par(mar=c(1,1,1,1))
plot.phylo(dendrogram.50$cons, show.tip.label=T, label.offset=0.001)
axisPhylo(backward = T, line=-6)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "lightblue", col = "black", cex=0.5)
dev.off()

#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Datasets")
#get region size
load(paste("04_Wombling/SizeRegions_", pick.index, "_", pick.num.subgraphs, ".R", sep=""))
#get region coordinate extent
load(paste("04_Wombling/CompCoor_", pick.index, "_", pick.num.subgraphs, ".R", sep=""))


comp.coor <- comp.coor[-removals]
SizeRegions <- SizeRegions[-removals]


extent.regions <- lapply(comp.coor, function(x) c(min(x[,1]), max(x[,1]), min(x[,2]), max(x[,2])))
extent.regions <- unlist(lapply(extent.regions, function(x) paste(round(x[1:4],1), collapse=",")))
class(extent.regions)
length(extent.regions)
extent.regions[12]

alternative.tip.label <- paste(1:length(extent.regions), SizeRegions, extent.regions, sep="|")
alternative.tip.label <- alternative.tip.label[match(as.numeric(dendrogram.50$cons$tip.label), remainers)]
#check that labels match
#cbind(alternative.tip.label, dendrogram.50$cons$tip.label)

plot.phylo(dendrogram.50$cons, show.tip.label=F)
axisPhylo(backward = T, line=0)
axis(2)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "lightblue", col = "black", cex=0.5)
tiplabels(alternative.tip.label, tip=1:length(dendrogram.50$cons$tip.label), adj = c(0, 0.5), frame = "none", col = "black", bg = "transparent")

#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Figures")
pdf(file = paste("06_Clustering/Dendrograms/Dendrogram50_Tax", pick.index,"_Regions", pick.num.subgraphs, "b_removals.pdf", sep=""), width=7, height=0.15*length(dendrogram.50$cons$tip.label))
par(mar=c(1,1,1,1))
plot.phylo(dendrogram.50$cons, show.tip.label=T, tip.color="transparent", label.offset=0.07)
axisPhylo(backward = T, line=-6)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "lightblue", col = "black", cex=0.5)
tiplabels(alternative.tip.label, tip=1:length(dendrogram.50$cons$tip.label), adj = c(0, 0.5), frame = "none", col = "black", bg = "transparent")
dev.off()



#################

#set node height threshold
node.height.threshold <- 0.2

#graphically examine threshold on dendrogram
plot.phylo(dendrogram.50$cons, show.tip.label=F)
axisPhylo(backward = T, line=0)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "lightblue", col = "black", cex=0.5)
axis(2)
abline(v=max(node.depth.edgelength(dendrogram.50$cons))-node.height.threshold, lty=1, col="red")


#determine how many and which nodes meet the threshold
sum(height.nodes.90.boot.support<node.height.threshold)
nodes.threshold <- nodes.90.boot.support[height.nodes.90.boot.support<node.height.threshold]
scale.nodes.threshold <- scale.nodes.90.boot.support[height.nodes.90.boot.support<node.height.threshold]
height.nodes.threshold <- height.nodes.90.boot.support[height.nodes.90.boot.support<node.height.threshold]

#create a vector with the highest (strongly supported) nodes that meet the height threshold,
#excluding (strongly suported) nodes that descend from other (strongly suported) nodes that
#meet the height threshold 
descendant.nodes <- vector()
for(i in 1:length(nodes.threshold))
{
	#all tips and nodes descending from a given node
	node.tips.des <- getDescendants(dendrogram.50$cons, nodes.threshold[i])
	#only the nodes descending from a a given node
	descendant.nodes <- append(descendant.nodes, node.tips.des[node.tips.des>length(dendrogram.50$cons$tip.label)])
	#length(dendrogram.50$cons$tip.label)
	#dendrogram.50$cons$Nnode
}

highest.nodes.threshold <- nodes.threshold[is.na(match(nodes.threshold, descendant.nodes))]
length(highest.nodes.threshold)
scale.highest.nodes.threshold <- scale.nodes.threshold[is.na(match(nodes.threshold, descendant.nodes))]
length(scale.highest.nodes.threshold)
height.highest.nodes.threshold <- height.nodes.threshold[is.na(match(nodes.threshold, descendant.nodes))]
length(height.highest.nodes.threshold)

#create a list with the tips descending from highest (strongly supported) nodes that meet
#the height threshold
tips.from.highest.nodes.threshold <- vector("list", length(highest.nodes.threshold))
for(i in 1:length(highest.nodes.threshold))
{
	#all tips and nodes descending from a given node
	node.tips.des <- getDescendants(dendrogram.50$cons, highest.nodes.threshold[i])
	#only the tips descending from a a given node
	tips.from.highest.nodes.threshold[[i]] <- node.tips.des[node.tips.des<=length(dendrogram.50$cons$tip.label)]
	#length(dendrogram.50$cons$tip.label)
	#dendrogram.50$cons$Nnode
}

length(unlist(tips.from.highest.nodes.threshold))


#define and examine colors to show regions and region clusters in dendrogram and raster map,
#for information about using color in thematic maps see:
#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
#https://www.r-bloggers.com/r-using-rcolorbrewer-to-colour-your-figures-in-r/

number.of.colors <- length(dendrogram.50$cons$tip.label) - length(unlist(tips.from.highest.nodes.threshold)) + length(tips.from.highest.nodes.threshold) 
regions.color <- colorRampPalette(brewer.pal(12,"Paired"))(number.of.colors)
set.seed(25) #works ok for 408 regions, 0.001 threshold
regions.color <- sample(regions.color, length(regions.color))
regions.color[1:20]
plot(1:number.of.colors, col = regions.color, pch = 16, cex = 3)

dendrogram.tips.not.clustered <- (1:length(dendrogram.50$cons$tip.label))[is.na(match(1:length(dendrogram.50$cons$tip.label), unlist(tips.from.highest.nodes.threshold)))]
length(dendrogram.tips.not.clustered)

dendrogram.tip.color <- rep(NA, times=length(dendrogram.50$cons$tip.label))
dendrogram.tip.color[dendrogram.tips.not.clustered] <- regions.color[1:length(dendrogram.tips.not.clustered)]
dendrogram.tip.color
sum(is.na(dendrogram.tip.color))
length(unlist(tips.from.highest.nodes.threshold))
#
length(dendrogram.tips.not.clustered)
length(highest.nodes.threshold)
length(tips.from.highest.nodes.threshold)
number.of.colors

for(i in 1:length(tips.from.highest.nodes.threshold))
{
	dendrogram.tip.color[tips.from.highest.nodes.threshold[[i]]] <- regions.color[i+length(dendrogram.tips.not.clustered)]
}
sum(is.na(dendrogram.tip.color))

#assign color to the largest region
dendrogram.tip.color[dendrogram.50$cons$tip.label=="1"] <- "gray70"


dendrogram.edge.color <- rep("black", times=nrow(dendrogram.50$cons$edge))
for(i in 1:length(tips.from.highest.nodes.threshold))
{
	targe.edges <- which.edge(dendrogram.50$cons, tips.from.highest.nodes.threshold[[i]])
	dendrogram.edge.color[targe.edges] <- regions.color[i+length(dendrogram.tips.not.clustered)]
}

dendrogram.edge.width <- rep(1, times=nrow(dendrogram.50$cons$edge))
for(i in 1:length(tips.from.highest.nodes.threshold))
{
	targe.edges <- which.edge(dendrogram.50$cons, tips.from.highest.nodes.threshold[[i]])
	dendrogram.edge.width[targe.edges] <- 4
}


#graphically examine threshold on dendrogram
plot.phylo(dendrogram.50$cons, show.tip.label=T, tip.color=dendrogram.tip.color, edge.color=dendrogram.edge.color, edge.width=dendrogram.edge.width)
axisPhylo(backward = T, line=0)
#nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "transparent", col = "black", cex=1)
axis(2)
abline(v=max(node.depth.edgelength(dendrogram.50$cons))-node.height.threshold, lty=1, col="red")
#nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "transparent", col = "black", cex=1)



pdf(file = paste("06_Clustering/Dendrograms/Dendrogram50_Tax", pick.index,"_Regions", pick.num.subgraphs,"_", node.height.threshold, "_threshold_removals.pdf", sep=""), width=7, height=0.15*length(dendrogram.50$cons$tip.label))
par(mar=c(1,1,1,1))
plot.phylo(dendrogram.50$cons, show.tip.label=T, tip.color=dendrogram.tip.color, edge.color=dendrogram.edge.color, edge.width=dendrogram.edge.width, label.offset=0.001)
axisPhylo(backward = T, line=-6)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "transparent", col = "black", cex=0.5)
abline(v=max(node.depth.edgelength(dendrogram.50$cons))-node.height.threshold, lty=1, col="black")
dev.off()

pdf(file = paste("06_Clustering/Dendrograms/Dendrogram50_Tax", pick.index,"_Regions", pick.num.subgraphs,"_", node.height.threshold, "_threshold_b_removals.pdf", sep=""), width=7, height=0.15*length(dendrogram.50$cons$tip.label))
par(mar=c(1,1,1,1))
plot.phylo(dendrogram.50$cons, show.tip.label=T, tip.color="transparent", edge.color=dendrogram.edge.color, edge.width=dendrogram.edge.width, label.offset=0.07)
axisPhylo(backward = T, line=-6)
nodelabels(scale.nodes.90.boot.support, nodes.90.boot.support, frame = "circle", bg = "transparent", col = "black", cex=0.5)
tiplabels(alternative.tip.label, tip=1:length(dendrogram.50$cons$tip.label), adj = c(0, 0.5), frame = "none", col = dendrogram.tip.color, bg = "transparent")
abline(v=max(node.depth.edgelength(dendrogram.50$cons))-node.height.threshold, lty=1, col="black")
dev.off()


####################################################################################################################################
# 4) Produce a raster map showing the regions and region clusters using dendrogram colors 
####################################################################################################################################

#read the raster of regions
Nordeste.mask.Regions <- raster(paste("04_Wombling/Nordeste_", pick.index,"_", pick.num.subgraphs,"_Regions.grd", sep=""))


raster.region.color <- rep(NA, times=pick.num.subgraphs)
raster.region.color[as.numeric(dendrogram.50$cons$tip.label)] <- dendrogram.tip.color

# NA regions (i.e. transition zones) can be coloured here
# Black
# raster.region.color[is.na(raster.region.color)] <- "#000000"
# White
# raster.region.color[is.na(raster.region.color)] <- "#FFFFFF"
# Grey
 raster.region.color[is.na(raster.region.color)] <- "#CCCCCC"




#plot regions, broad scale
plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F)
plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F,
xlab="Longitude", ylab="Longitude", cex.axis=1.4, cex.lab=1.5)


dir.create("06_Clustering/Maps", showWarnings = F)

pdf(file = paste("06_Clustering/Maps/Map_Dendrogram50_Tax", pick.index, "_Regions", pick.num.subgraphs, "_", node.height.threshold, "_threshold_removals.pdf", sep=""))
plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F)
grid()
dev.off()




####################################################################################################################################
# 5) Explore the spatial distribution of clusters
####################################################################################################################################


topo.Nordeste.mask.0 <- raster("000_GIS_LAYERS/Brazil_Masked_GIS_Layers/alt/alt_at_MODIS_CHIRPS_res.tif")
topo.Nordeste.mask.0 <- crop(topo.Nordeste.mask.0, Nordeste.mask.Regions)
topo.Nordeste.mask.0 <- mask(topo.Nordeste.mask.0, Nordeste.mask.Regions)

plot(topo.Nordeste.mask.0, colNA="black")


#close all graphical devices, and then open two graphical devices
dev.new()
dev.new()
available.graphical.windows <- dev.list()
available.graphical.windows

for(i in 1:length(tips.from.highest.nodes.threshold))
{
	dev.set(available.graphical.windows[3])
	zoom(dendrogram.50$cons, tips.from.highest.nodes.threshold[[i]], subtree = FALSE, 
	col = regions.color[1+length(dendrogram.tips.not.clustered):length(regions.color)][i],
	tip.color=regions.color[1+length(dendrogram.tips.not.clustered):length(regions.color)][i],
	edge.width=1)
	dev.set(available.graphical.windows[4])
	raster.region.color <- rep("gray70", times=length(dendrogram.50$cons$tip.label))
	raster.region.color[tips.from.highest.nodes.threshold[[i]]] <- dendrogram.tip.color[tips.from.highest.nodes.threshold[[i]]]
	plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F)
	Sys.sleep(5)
}




for(i in 1:length(tips.from.highest.nodes.threshold))
{
	dev.set(available.graphical.windows[3])
	zoom(dendrogram.50$cons, tips.from.highest.nodes.threshold[[i]], subtree = FALSE, 
	col = "red", tip.color="red",	edge.width=1) #, tip)
	dev.set(available.graphical.windows[4])
	raster.region.color <- rep("gray70", times=length(dendrogram.50$cons$tip.label))
	raster.region.color[tips.from.highest.nodes.threshold[[i]]] <- "red"
	plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F,
	main=paste("Cluster ", i, " of ", length(tips.from.highest.nodes.threshold), ", 90% bootstrap scale = ", scale.highest.nodes.threshold[i], ", node height = ", round(height.highest.nodes.threshold[i],3), sep=""))
	Sys.sleep(3)
}


dendrogram.50.alternative <- dendrogram.50$cons
dendrogram.50.alternative$tip.label <- unlist(lapply(strsplit(alternative.tip.label, split="|", fixed = T), function(x) paste(x[1], x[2], sep=" | ")))


for(i in 1:length(tips.from.highest.nodes.threshold))
{
	dev.set(available.graphical.windows[3])
	zoom(dendrogram.50.alternative, tips.from.highest.nodes.threshold[[i]], subtree = F, 
	col = "red", tip.color="red",	edge.width=1)
	par(new=T)
	par(mfcol=c(1,1))
	title(paste("Cluster ", i, " of ", length(tips.from.highest.nodes.threshold), ", 90% bootstrap scale = ", scale.highest.nodes.threshold[i], ", node height = ", round(height.highest.nodes.threshold[i],3), sep=""))
	dev.set(available.graphical.windows[4])
	raster.region.color <- rep("gray70", times=length(dendrogram.50$cons$tip.label))
	raster.region.color[as.numeric(dendrogram.50$cons$tip.label[tips.from.highest.nodes.threshold[[i]]])] <- "red"
	plot(Nordeste.mask.Regions, col=raster.region.color, useRaster=T, legend=F)
	Sys.sleep(9)
}





for(i in 1:length(tips.from.highest.nodes.threshold))
{
	dev.set(available.graphical.windows[3])
	par(mar=c(1,1,1,1))
	zoom(dendrogram.50.alternative, tips.from.highest.nodes.threshold[[i]], subtree = F, 
		col = "red", tip.color="red",	edge.width=1)
	par(new=T)
	par(mfcol=c(1,1))
	title(paste("Cluster ", i, " of ", length(tips.from.highest.nodes.threshold), ", 90% bootstrap scale = ", scale.highest.nodes.threshold[i], ", node height = ", round(height.highest.nodes.threshold[i],3), sep=""))
	dev.set(available.graphical.windows[4])
	raster.region.color <- rep("transparent", times=length(dendrogram.50$cons$tip.label))
	raster.region.color[as.numeric(dendrogram.50$cons$tip.label[tips.from.highest.nodes.threshold[[i]]])] <- "red"
	plot(topo.Nordeste.mask.0, colNA="black", box=T, axes=T, useRaster=T, col=gray.colors(100, start = 0.7, end = 1, gamma = 2.2, alpha = NULL),
		xlab="Longitude", ylab="Latitude", horizontal=T,
		legend.lab="Altitude (m)", legend.mar=-1,
		smallplot= c(0.12,0.9,0.07,0.095), bigplot=c(0.12,0.9,0.33,0.99))
	plot(Nordeste.mask.Regions, col=raster.region.color, colNA="black", useRaster=T, legend=F, add=T)
	Sys.sleep(5)
}







##################################################
# Notes
##################################################


class(a.consensus0.5dendrogram)
str(a.consensus0.5dendrogram)
str(a.consensus0.5dendrogram$cons)
plot(a.consensus0.5dendrogram$cons, direction="downwards", show.node.label = TRUE)
plot.phylo(a.consensus0.5dendrogram$cons, use.edge.length = TRUE)
axisPhylo()
a.consensus0.5dendrogram$cons$edge.length[5:7]
a.consensus0.5dendrogram$cons$edge[1:20,]
a.consensus0.5dendrogram$cons$tip.label[c(57,84)]



plot(1:10, col = colorRampPalette(c('red','blue','green'))(10), pch = 16, cex = 3)
plot(1:9, col = brewer.pal(9,"Blues"), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(9,"Blues"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(12,"Paired"))(100), pch = 16, cex = 3)
plot(1:12, col = brewer.pal(12,"Paired"), pch = 16, cex = 3)
plot(1:12, col = brewer.pal(8,"Accent"), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(8,"Accent"))(100), pch = 16, cex = 3)
plot(1:8, col = brewer.pal(8,"Dark2"), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(8,"Dark2"))(100), pch = 16, cex = 3)

plot(1:100, col = colorRampPalette(brewer.pal(11,"Spectral"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(11,"BrBG"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(11,"PiYG"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(11,"PRGn"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(11,"RdYlBu"))(100), pch = 16, cex = 3)
plot(1:100, col = colorRampPalette(brewer.pal(11,"RdYlGn"))(100), pch = 16, cex = 3)

plot(1:100, col = rainbow(6), pch = 16, cex = 3)

qual.col.pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector <- unlist(mapply(brewer.pal, qual.col.pals$maxcolors, rownames(qual_col_pals)))

plot(1:length(col.vector), col = col.vector, pch = 16, cex = 3)

plot(1:500, col=colorRampPalette(col.vector)(500), pch = 16, cex = 3)





map.color <-  grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
plot(1:length(map.color), col = map.color, pch = 16, cex = 3)
plot(1:20, col = map.color, pch = 16, cex = 3)

