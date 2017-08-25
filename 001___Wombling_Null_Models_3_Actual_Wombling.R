
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
no.na.cells <- (1:ncell(Nordeste.mask.0))[!is.na(extract(Nordeste.mask.0, 1:ncell(Nordeste.mask.0)))]

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

for(x in 1:length(bricks)){
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

for(x in 1:length(bricks)){
  SDM.b <- brick(paste("05_Wombling_Null_Models/03_Null_Bricks/Null_", bricks[[x]], "_SDMb.grd", sep=""))
  writeLines(paste("\nWorking on brick ", x, " of ", length(bricks), sep=""))
  
#read in the beta diversity measures
  obs.beta.sim <- read.table(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim.csv", sep=""), sep=",", header=T)[,2]
  obs.beta.sor <- read.table(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor.csv", sep=""), sep=",", header=T)[,2]
  r.obs.beta.sim <- read.table(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim_ranked.csv", sep=""), sep=",", header=T)[,2]
  r.obs.beta.sor <- read.table(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor_ranked.csv", sep=""), sep=",", header=T)[,2]
  
#first create a matrix defining the edges:
  cell.adj.char <- cbind(as.character(cell.adj[,1]), as.character(cell.adj[,2]))
  #
#next define the graph
  Nordeste.graph.sim <- graph_from_edgelist(cell.adj.char, directed = F)
  Nordeste.graph.sor <- graph_from_edgelist(cell.adj.char, directed = F)
  
#run a loop to sequentially remove spatial links and count the resulting number
#of subgraphs, which correspond to isolated regions (and potentially ecoregions).
  number.of.subgraphs.sim <- rep(NA, times=length(obs.beta.sim))
  number.of.subgraphs.sor <- rep(NA, times=length(obs.beta.sor))
  #
  start.time <- Sys.time()
  # 
  for(i in 1:length(obs.beta.sim))
  {
    percents <- round(length(obs.beta.sim)*seq(1, 0, -0.05))
    if(i %in% percents){
      writeLines(paste("... ", round(i/length(obs.beta.sim)*100, 0), "% at ", Sys.time(), sep=""))
    }
    candidate.boundary.elements.to.deploy.sim <- which(r.obs.beta.sim > (length(r.obs.beta.sim)-i))
    candidate.boundary.elements.to.deploy.sor <- which(r.obs.beta.sor > (length(r.obs.beta.sor)-i))
    number.of.subgraphs.sim[i] <- components(delete_edges(Nordeste.graph.sim, candidate.boundary.elements.to.deploy.sim))$no
    number.of.subgraphs.sor[i] <- components(delete_edges(Nordeste.graph.sor, candidate.boundary.elements.to.deploy.sor))$no
  }
  
  write.csv(number.of.subgraphs.sim, file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim_number_of_subgraphs.csv", sep=""))
  write.csv(number.of.subgraphs.sor, file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor_number_of_subgraphs.csv", sep=""))
}



##################################################################################################
# 11) Calculate superfluidity, as defined by Oden et al. (1993, cited in the introduction).
# This step may take several hours or days, and you might skip to the end of this section
# to read a file with previously calculated superfluidity values. Use code under 9.1 or 9.2 according
# to the version of package igraph and R you are using.
##################################################################################################

#define the percentiles of beta-diversity for which superfluidity will be claculated,
#if you have not done so in section 8 (above)
beta.ranks.to.evaluate <- 35082 - round(35082*seq(0.05, 0.5, 0.05))

number.of.subgraphs.sim <- read.csv("04_Wombling/NumberSubgraphsSim.txt")
number.of.subgraphs.sor <- read.csv("04_Wombling/NumberSubgraphsSor.txt")
evaluation.number.of.subgraphs.sim <- read.csv("04_Wombling/superfluiditySim.txt")[,1]
evaluation.number.of.subgraphs.sor <- read.csv("04_Wombling/superfluiditySor.txt")[,1]

cell.adj.char <- cbind(as.character(cell.adj[,1]), as.character(cell.adj[,2]))
Nordeste.graph <- graph_from_edgelist(cell.adj.char, directed = F)

for(x in 121:112){
  writeLines(paste("brick ", x, " of ", length(bricks), sep=""))
  null.number.of.subgraphs.sim <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim_number_of_subgraphs.csv", sep=""))
  null.number.of.subgraphs.sor <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor_number_of_subgraphs.csv", sep=""))

  null.r.obs.beta.sim <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim_ranked.csv", sep=""))[,2]
  null.r.obs.beta.sor <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor_ranked.csv", sep=""))[,2]
  
  null.superfluidity.sim <- rep(NA, times=length(evaluation.number.of.subgraphs.sim))
  null.superfluidity.sor <- rep(NA, times=length(evaluation.number.of.subgraphs.sor))
  
  Nordeste.graph.sim <- Nordeste.graph
  Nordeste.graph.sor <- Nordeste.graph
  
  for(j in 1:length(evaluation.number.of.subgraphs.sim)){
    writeLines(paste("...", j, "of", length(evaluation.number.of.subgraphs.sim), "at", Sys.time()))
    candidate.boundary.elements.to.delete.sim <- which(null.r.obs.beta.sim > (length(null.r.obs.beta.sim) - match(evaluation.number.of.subgraphs.sim[[j]], null.number.of.subgraphs.sim[,2])))
    candidate.boundary.elements.to.delete.sor <- which(null.r.obs.beta.sor > (length(null.r.obs.beta.sor) - match(evaluation.number.of.subgraphs.sor[[j]], null.number.of.subgraphs.sor[,2])))
    
    focal.graph.sim <- delete_edges(Nordeste.graph, candidate.boundary.elements.to.delete.sim)
    focal.graph.sor <- delete_edges(Nordeste.graph, candidate.boundary.elements.to.delete.sor)
    
    necessary.edge.sim <- rep(NA, times=length(candidate.boundary.elements.to.delete.sim))
    writeLines("   ...Working on Simpson")
    for(i in 1:length(candidate.boundary.elements.to.delete.sim))
    {
      necessary.edge.sim[i] <- ifelse(components(add_edges(focal.graph.sim, cell.adj.char[candidate.boundary.elements.to.delete.sim[i],]))$no < components(focal.graph.sim)$no, 1, 0)
    }
    writeLines("   ...Working on Sorenson")
    necessary.edge.sor <- rep(NA, times=length(candidate.boundary.elements.to.delete.sor))
    for(i in 1:length(candidate.boundary.elements.to.delete.sor))
    {
      necessary.edge.sor[i] <- ifelse(components(add_edges(focal.graph.sor, cell.adj.char[candidate.boundary.elements.to.delete.sor[i],]))$no < components(focal.graph.sor)$no, 1, 0)
    }
    null.superfluidity.sim[j] <- sum(necessary.edge.sim<1)/sum(necessary.edge.sim>0)
    null.superfluidity.sor[j] <- sum(necessary.edge.sor<1)/sum(necessary.edge.sor>0)
  }
  write.csv(cbind(evaluation.number.of.subgraphs.sim, null.superfluidity.sim), file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/superfluidity_sim.csv", sep=""))
  write.csv(cbind(evaluation.number.of.subgraphs.sor, null.superfluidity.sor), file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/superfluidity_sor.csv", sep=""))
}



############################################################################################################################
# 12) Test the significance of the superfluidity
#     PART 1 :- COllate the data
############################################################################################################################

subgraphs.NULL.sim <- matrix(NA, nrow=length(bricks), ncol=35082)
subgraphs.NULL.sor <- matrix(NA, nrow=length(bricks), ncol=35082)

# Find those null bricks for which superfluidity has been calculated
sims_w_superfluidity <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks, "/superfluidity_sim.csv", sep="")))
sors_w_superfluidity <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks, "/superfluidity_sor.csv", sep="")))

# Collate the superfluidity values into matrices
superfluidity.NULL.sim <- matrix(NA, nrow=length(sims_w_superfluidity), ncol=10)
superfluidity.NULL.sor <- matrix(NA, nrow=length(sors_w_superfluidity), ncol=10)

for(x in 1:length(sims_w_superfluidity)){
  superfluidity.NULL.sim[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[sims_w_superfluidity][[x]], "/superfluidity_sim.csv", sep=""))[,3]
}
for(x in 1:length(sors_w_superfluidity)){
  superfluidity.NULL.sor[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[sors_w_superfluidity][[x]], "/superfluidity_sor.csv", sep=""))[,3]
}


# Very occasionally superfluidity will be infinite (i.e. when ALL boundary elements contribute)
# Replace these with 200 for ease of plotting
superfluidity.NULL.sim[superfluidity.NULL.sim==Inf] <- 200
superfluidity.NULL.sor[superfluidity.NULL.sor==Inf] <- 200

# Read in the emperical superfluidity values
superfluidity.OBS.sim <- read.csv("04_Wombling/superfluiditySim.txt")
evaulation.number.of.subgraphs.sim <- superfluidity.OBS.sim[,1]
superfluidity.OBS.sim <- superfluidity.OBS.sim[,2]
superfluidity.OBS.sor <- read.csv("04_Wombling/superfluiditySor.txt")
evaulation.number.of.subgraphs.sor <- superfluidity.OBS.sor[,1]
superfluidity.OBS.sor <- superfluidity.OBS.sor[,2]

############################################################################################################################
# 13) Test the significance of the superfluidity
#     PART 2 :- Plot the data
############################################################################################################################

plot(evaluation.number.of.subgraphs.sim, superfluidity.NULL.sim[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(evaluation.number.of.subgraphs.sim)), ylim=c(0, max(superfluidity.NULL.sim)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(superfluidity.NULL.sim))
{
  points(evaluation.number.of.subgraphs.sim, superfluidity.NULL.sim[i,], type="l", col="gray70")
}
points(evaluation.number.of.subgraphs.sim, superfluidity.OBS.sim, type="o", col="blue", lwd=2)

plot(evaluation.number.of.subgraphs.sor, superfluidity.NULL.sor[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(evaluation.number.of.subgraphs.sor)), ylim=c(0, max(superfluidity.NULL.sor)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(superfluidity.NULL.sor))
{
  points(evaluation.number.of.subgraphs.sor, superfluidity.NULL.sor[i,], type="l", col="gray70")
}
points(evaluation.number.of.subgraphs.sor, superfluidity.OBS.sor, type="o", col="blue", lwd=2)


############################################################################################################################
# 14) Test the significance of the superfluidity
#     PART 2 :- Test for Significance
############################################################################################################################


p.superfluidity.sim <- rep(NA, times=ncol(superfluidity.NULL.sim))
p.superfluidity.sor <- rep(NA, times=ncol(superfluidity.NULL.sor))
#calculate p-value
for(i in 1:ncol(superfluidity.NULL.sim))
{
  p.superfluidity.sim[i] <- sum(superfluidity.OBS.sim[i] >= superfluidity.NULL.sim[,i])/nrow(superfluidity.NULL.sim)
}	

p.superfluidity.sim

#calculate p-value
for(i in 1:ncol(superfluidity.NULL.sim))
{
  p.superfluidity.sor[i] <- sum(superfluidity.OBS.sor[i] >= superfluidity.NULL.sor[,i])/nrow(superfluidity.NULL.sor)
}	

p.superfluidity.sor


############################################################################################################################
# 15) Test the significance of the number of subgraphs
#     PART 1 :- COllate the data
############################################################################################################################

subgraphs.NULL.sim <- matrix(NA, nrow=length(bricks), ncol=35082)
subgraphs.NULL.sor <- matrix(NA, nrow=length(bricks), ncol=35082)

# Find those null bricks for which the number of subgraphs has been calculated
sims_w_no_o_subgraphs <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks, "/sim_number_of_subgraphs.csv", sep="")))
sors_w_no_o_subgraphs <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks, "/sor_number_of_subgraphs.csv", sep="")))

# Collate the superfluidity values into matrices
no.o.subgraphs.NULL.sim <- matrix(NA, nrow=length(sims_w_no_o_subgraphs), ncol=35082)
no.o.subgraphs.NULL.sor <- matrix(NA, nrow=length(sors_w_no_o_subgraphs), ncol=35082)

for(x in 1:length(sims_w_no_o_subgraphs)){
  no.o.subgraphs.NULL.sim[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[sims_w_no_o_subgraphs][[x]], "/sim_number_of_subgraphs.csv", sep=""))[,2]
}
for(x in 1:length(sors_w_no_o_subgraphs)){
  no.o.subgraphs.NULL.sor[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[sors_w_no_o_subgraphs][[x]], "/sor_number_of_subgraphs.csv", sep=""))[,2]
}

# Read in the emperical superfluidity values
no.o.subgraphs.OBS.sim <- read.csv("04_Wombling/NumberSubgraphsSim.txt")
no.o.subgraphs.OBS.sim <- no.o.subgraphs.OBS.sim[,1]
no.o.subgraphs.OBS.sor <- read.csv("04_Wombling/NumberSubgraphsSor.txt")
no.o.subgraphs.OBS.sor <- no.o.subgraphs.OBS.sor[,1]


############################################################################################################################
# 15) Test the significance of the number of subgraphs
#     PART 2 :- Plot the data
############################################################################################################################

plot(1:35082, no.o.subgraphs.NULL.sim[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(35082)), ylim=c(0, max(no.o.subgraphs.NULL.sim)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(no.o.subgraphs.NULL.sim))
{
  points(1:35082, no.o.subgraphs.NULL.sim[i,], type="l", col="gray70")
}
points(1:35082, no.o.subgraphs.OBS.sim, type="o", col="blue", lwd=1)

plot(1:35082, no.o.subgraphs.NULL.sor[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(35082)), ylim=c(0, max(no.o.subgraphs.NULL.sor)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(no.o.subgraphs.NULL.sor))
{
  points(1:35082, no.o.subgraphs.NULL.sor[i,], type="l", col="gray70")
}
points(1:35082, no.o.subgraphs.OBS.sor, type="o", col="blue", lwd=1)


############################################################################################################################
# 16) Test the significance of the number of subgraphs
#     PART 2 :- Test for Significance
############################################################################################################################


p.no.o.subgraphs.sim <- rep(NA, times=ncol(no.o.subgraphs.NULL.sim))
p.no.o.subgraphs.sor <- rep(NA, times=ncol(no.o.subgraphs.NULL.sor))
#calculate p-value
for(i in 1:ncol(no.o.subgraphs.NULL.sim))
{
  p.no.o.subgraphs.sim[i] <- sum(no.o.subgraphs.OBS.sim[i] <= no.o.subgraphs.NULL.sim[,i])/nrow(no.o.subgraphs.NULL.sim)
}	

p.no.o.subgraphs.sim
plot(p.no.o.subgraphs.sim)

#calculate p-value
for(i in 1:ncol(no.o.subgraphs.NULL.sor))
{
  p.no.o.subgraphs.sor[i] <- sum(no.o.subgraphs.OBS.sor[i] <= no.o.subgraphs.NULL.sor[,i])/nrow(no.o.subgraphs.NULL.sor)
}	

p.no.o.subgraphs.sor
plot(p.no.o.subgraphs.sor)
