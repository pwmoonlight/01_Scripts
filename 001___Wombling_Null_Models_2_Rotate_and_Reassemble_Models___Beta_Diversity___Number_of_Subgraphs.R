
############################################################################################################################
# 1) Load packages, define working directories and the iterations to perform
############################################################################################################################

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

#load pakages
library(sp)
library(raster)
library(igraph)
library(ape)
library(RColorBrewer)

#define the iterations to perform
iterations.to.perform <- 1000

dir.create("05_Wombling_Null_Models/03_Null_Beta_Diversity")
dir.create("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs")


############################################################################################################################
# 2)  Read in Required Objects
############################################################################################################################



############################################################################################################################
# 2.1) Read in the thresholded distributions
############################################################################################################################

#read the file that has the brick with the species distribution models of all species
species.in.analysis <-  gsub("\\.", "-", gsub("_", " ", gsub("MCLUST_", "", gsub(".R$", "", list.files("05_Wombling_Null_Models/01_McClust/", pattern="*.R$", full.names=F)))))
species.level.dir <- lapply(1:length(species.in.analysis), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species.in.analysis)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)

nlayers(SDM.b)
cell.side <- res(SDM.b)[1] #resolution, i.e., cell dimensions (length and width)
focal.sp.range <- raster(SDM.b, layer=1)
no.na.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range)))]
s1.s2 <- coordinates(focal.sp.range)[no.na.cells,]
dim(s1.s2)



############################################################################################################################
# 2.2) Read in the mask of the study area
############################################################################################################################

#create and plot the Nordeste raster mask that indicates the grid cells that have climate data,
Nordeste.mask.0 <- raster("000_GIS_LAYERS/nordeste.tif")
no.na.cells <- (1:ncell(Nordeste.mask.0))[!is.na(extract(Nordeste.mask.0, 1:ncell(Nordeste.mask.0)))]



############################################################################################################################
# 2.3) Read in the adjacent grid cells  
############################################################################################################################

cell.adj <- read.table("04_Wombling/Cell_Links/cell_adj_150arc.txt", sep=",", header=T)
cell.adj.char <- cbind(as.character(cell.adj[,1]), as.character(cell.adj[,2]))



############################################################################################################################
# 2.4) Read in the file names 
############################################################################################################################

#AOO <- rep(NA, times=length(species.level.dir))
#for(i in 1:length(species.level.dir)){
#  print(i/length(species.level.dir)*100)
#  AOO[i] <- cellStats(species.level.dir[[i]], 'sum')
#}

#save(AOO, file=paste("04_Wombling/AOO.R", sep=""))
load(file=paste("04_Wombling/AOO.R", sep=""))

species.in.analysis <- read.csv(file="04_Wombling/species_in_analysis.txt", sep=",", header=F)[,1]
brick.index.species.in.phylogeny <- gsub("-", ".", gsub(" ", "_", as.character(species.in.analysis)))


############################################################################################################################
# 2.5) Read in the AOOs of each species  
############################################################################################################################

MCLUST.model.names <- gsub(" ", "_", paste("MCLUST", species.in.analysis, sep="_"))

### Find any files that will have been renamed by this process
### These will need to be manually renamed before the next (very long!) stage
for(x in MCLUST.model.names){
  if(grepl("-", x)){
    print(x)
  }
}

############################################################################################################################
# 3) Do absolutely everything
############################################################################################################################



ptm <- Sys.time()
dir.create("05_Wombling_Null_Models/02_Null_Distributions", showWarnings = F)

for(null.iterations in 1:iterations.to.perform)
{
  brick <- null.iterations
  if(!file.exists(paste("05_Wombling_Null_Models/02_Null_Distributions/", null.iterations, sep="")))
  {
    writeLines(paste("\nSimulating Distributions for brick ", brick, sep=""))

    dir.create(paste("05_Wombling_Null_Models/02_Null_Distributions/", null.iterations, sep=""), showWarnings = F)
    #dir.create(paste("05_Wombling_Null_Models/03_Null_Bricks/", null.iterations, sep=""), showWarnings = F)
    #the vector "final.j" will store the number of attempts needed to rotate the centroids of each species
    final.j <- rep(NA, times=length(MCLUST.model.names)) 
    writeLines(paste("...0% at ", Sys.time(), sep=""))
    for(h in 1:length(MCLUST.model.names))
    {
      percents <- round(length(MCLUST.model.names)*seq(1, 0, -0.02))
      if(h %in% percents){
        writeLines(paste("...", round(h/length(MCLUST.model.names)*100, 0), "% at ", Sys.time(), sep=""))
      }

      #load Mclust model and assign to:
      load(paste("05_Wombling_Null_Models/01_McClust/", MCLUST.model.names[[h]], ".R", sep=""))
      Mcluster.focal.sp <- eval(as.name(gsub("-", ".", MCLUST.model.names[h])))
      rm(list=c(MCLUST.model.names[h]))
      
      focal.sp.range <- raster(SDM.b, layer=h)
      occurrence.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range))) & extract(focal.sp.range, 1:ncell(focal.sp.range))>0]
      occu.coor <- coordinates(focal.sp.range)[occurrence.cells,]
      rm(occurrence.cells)
      
      ###############
      # Shift and rotate centroids of range
      ###############
      
      #define the number of (shifted and rotated) simulated ranges to be obtained "k"
      k <- 1
      #define the number of iterations "m" to be performed
      #This is a very large number. This number of iterations will not actually be performed.
      #Instead, the first suitable iteration will be kept.
      m <- 100000
      
      #set the value of i and j equal to one
      i <- 1
      j <- 1
      
      #find rotation and shift for means
      while(k>=i & m>=j)
      {
        
        #randomly sample a point within the study area that will be the centroid of the shifted rotated geographic range;
        #a raster cell within the study region is randomly picked, then the coordinates of the center of that cell are
        #slightly modified by adding a number randomly sampled from a uniform distribution between -cell.side/2 and
        #+cell.side/2, that way the centroid of the shifted rotated geographic range is not constrained to be at the
        #center of any raster cell
        centroid.null <- s1.s2[sample(1:nrow(s1.s2), size=1),] + runif(2, min=-cell.side/2, max=cell.side/2)
        
        #center means of the focal species
        centered.means <- scale(t(Mcluster.focal.sp$means), center = TRUE, scale = FALSE)
        
        #create rotation matrix for an angle "theta.means" in radians, the angle is randomly sampled from a uniform distribution between 0 and 2*pi
        theta.means <- runif(1, min=0, max=2*pi)
        rotation.matrix.means <- matrix(c(cos(theta.means), sin(theta.means), -sin(theta.means), cos(theta.means)), nrow = 2, ncol = 2)
        
        #shift and rotate means
        shifted.rotated.means <- t(rotation.matrix.means %*% t(centered.means) + centroid.null)
        
        #update the count of j
        j <- j+1
        
        #if not 3/4 shifted and rotated localities fall within the study region go to the next iteration,
        #otherwise place the shifted and rotated localities in the position i of the list "shifted.rotated.ranges"   
        if(sum(is.na(extract(focal.sp.range,  shifted.rotated.means)))>length(shifted.rotated.means[,1])/4) next
        i <- 2
      }
      
      final.j[h] <- j 
      
      ###############
      # Shift and rotate range, based on null centroid and random angle obtained above
      ###############
      
      centered.occu.coor <-  occu.coor - matrix(rep(colMeans(t(Mcluster.focal.sp$means)), times=nrow(occu.coor)), ncol=2, byrow=T)
      rm(occu.coor)
      shifted.rotated.range <- t(rotation.matrix.means %*% t(centered.occu.coor) + centroid.null)
      shifted.rotated.range <- shifted.rotated.range[!is.na(extract(focal.sp.range,  shifted.rotated.range)),]
      occurrences.r <- focal.sp.range
      occurrences.r[no.na.cells] <- 0
      
      if(class(shifted.rotated.range) == "numeric"){
        shifted.rotated.range <- t(as.data.frame(shifted.rotated.range))
      }
      
      occurrences.r <- rasterize(shifted.rotated.range, occurrences.r, rep(1, times=nrow(shifted.rotated.range)), update=T)
      rm(focal.sp.range, centroid.null, Mcluster.focal.sp, shifted.rotated.range)

      absence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))<1
      absence.cells.index.r[is.na(absence.cells.index.r)] <- FALSE
      absence.cells.r <- (1:ncell(occurrences.r))[absence.cells.index.r]
      
      presence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))>0
      presence.cells.index.r[is.na(presence.cells.index.r)] <- FALSE
      presence.cells.r <- (1:ncell(occurrences.r))[presence.cells.index.r]
      
      AOO.fraction.to.add <- as.vector(AOO[h] - length(presence.cells.r))
      
      # If the whole range falls within the study area, proceed
      if(AOO.fraction.to.add < 1) writeRaster(occurrences.r, filename=paste("05_Wombling_Null_Models/02_Null_Distributions/", null.iterations, "/", names(SDM.b)[h], "_NULL.grd", sep=""))
      if(AOO.fraction.to.add < 1) next 
      
      # If not, add cells around the edge of the range until the null AOO equals the AOO
      while(AOO.fraction.to.add > 0)
      {
        colo.cells <- adjacent(occurrences.r, presence.cells.r, directions=4, pairs=F, target=absence.cells.r, sorted=FALSE, include=FALSE, id=FALSE) 
        coo.colo.cells <- xyFromCell(occurrences.r, colo.cells, spatial=FALSE)
        if(length(colo.cells) <= AOO.fraction.to.add) occurrences.r <- rasterize(matrix(coo.colo.cells, nrow=length(colo.cells), ncol=2), occurrences.r, rep(1, times=length(colo.cells)), update=T)
        if(length(colo.cells) > AOO.fraction.to.add) additional.occurrences <- sample(1:nrow(coo.colo.cells), size= AOO.fraction.to.add)
        if(length(colo.cells) > AOO.fraction.to.add) occurrences.r <- rasterize(matrix(coo.colo.cells[additional.occurrences,], nrow=length(additional.occurrences), ncol=2), occurrences.r, rep(1, times=length(additional.occurrences)), update=T)
        
        absence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))<1
        absence.cells.index.r[is.na(absence.cells.index.r)] <- FALSE
        absence.cells.r <- (1:ncell(occurrences.r))[absence.cells.index.r]
        
        presence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))>0
        presence.cells.index.r[is.na(presence.cells.index.r)] <- FALSE
        presence.cells.r <- (1:ncell(occurrences.r))[presence.cells.index.r]
        AOO.fraction.to.add <- as.vector(AOO[h] - length(presence.cells.r))
      }
      
      rm(colo.cells, absence.cells.index.r, presence.cells.index.r, presence.cells.r, AOO.fraction.to.add)

      # Save the NULL species raster
      writeRaster(occurrences.r, filename=paste("05_Wombling_Null_Models/02_Null_Distributions/", null.iterations, "/", names(SDM.b)[h], "_NULL.grd", sep=""), overwrite=T)

      #save file with vector "final.j", which stores the number of attempts needed to rotate the centroids of each species
      write.table(final.j, file = paste("05_Wombling_Null_Models/Null_", null.iterations, "_final_j.csv", sep=""), quote = TRUE, sep = ",")
    }

    ###############
    # Calculate the Beta Diversity
    ###############
    
    writeLines(paste("\nCalculating the Beta Diversity of brick ", brick, sep=""))
    dir.create(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", brick, sep=""))

    writeLines(paste("...Creating Brick", sep=""))
    null.SDM.b <- lapply(1:length(species.in.analysis), function(x){raster(paste("05_Wombling_Null_Models/02_Null_Distributions/", brick, "/", brick.index.species.in.phylogeny[[x]], "_NULL.gri", sep=""))})
    SDM.raster.names <- vector(mode = "list", length = length(null.SDM.b))
    for(i in 1:length(null.SDM.b[]))
    {
      null.SDM.b[[i]]@data@names  <- brick.index.species.in.phylogeny[[i]]
    }
    for(x in SDM.raster.names){eval(x)}
    
    writeLines(paste("...Stacking Brick", sep=""))
    null.SDM.b <- stack(null.SDM.b)
    
    writeLines("...Calculating beta diversity")
    source(paste(getwd(), "/01_Scripts/Wombling Modules/01_SIM_and_SOR_diversity.R", sep=""))
    
    writeLines("...Ranking beta diversity values")
    source(paste(getwd(), "/01_Scripts/Wombling Modules/02_rank_SIM_and_SOR_diversity.R", sep=""))
    
    rm(SDM.raster.names)
    
    
    ###############
    # Calculate the Number of Subgraphs
    ###############
        
    writeLines(paste("\nCalculating the Number of Subgraphs of brick ", brick, sep=""))
    
    #read in the beta diversity measures
    obs.beta.sim <- read.table(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", brick, "/sim.csv", sep=""), sep=",", header=T)[,2]
    obs.beta.sor <- read.table(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", brick, "/sor.csv", sep=""), sep=",", header=T)[,2]
    r.obs.beta.sim <- read.table(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", brick, "/sim_ranked.csv", sep=""), sep=",", header=T)[,2]
    r.obs.beta.sor <- read.table(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", brick, "/sor_ranked.csv", sep=""), sep=",", header=T)[,2]
    
    #next define the graph
    Nordeste.graph.sim <- graph_from_edgelist(cell.adj.char, directed = F)
    Nordeste.graph.sor <- Nordeste.graph.sim
    
    #run a loop to sequentially remove spatial links and count the resulting number
    #of subgraphs, which correspond to isolated regions (and potentially ecoregions).
    
    number.of.subgraphs.sim <- rep(NA, times=length(obs.beta.sim))
    number.of.subgraphs.sor <- rep(NA, times=length(obs.beta.sor))
    
    writeLines(paste("...0% at ", Sys.time(), sep=""))
    for(i in 1:length(obs.beta.sim))
    {
      percents <- round(length(obs.beta.sim)*seq(1, 0, -0.05))
      if(i %in% percents){
        writeLines(paste("...", round(i/length(obs.beta.sim)*100, 0), "% at ", Sys.time(), sep=""))
      }
      candidate.boundary.elements.to.deploy.sim <- which(r.obs.beta.sim > (length(r.obs.beta.sim)-i))
      candidate.boundary.elements.to.deploy.sor <- which(r.obs.beta.sor > (length(r.obs.beta.sor)-i))
      number.of.subgraphs.sim[i] <- components(delete_edges(Nordeste.graph.sim, candidate.boundary.elements.to.deploy.sim))$no
      number.of.subgraphs.sor[i] <- components(delete_edges(Nordeste.graph.sor, candidate.boundary.elements.to.deploy.sor))$no
    }
    
    dir.create(paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", brick, sep=""))
    write.csv(number.of.subgraphs.sim, file=paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", brick, "/sim_number_of_subgraphs.csv", sep=""))
    write.csv(number.of.subgraphs.sor, file=paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", brick, "/sor_number_of_subgraphs.csv", sep=""))
    
    rm(obs.beta.sim, obs.beta.sor, percents, candidate.boundary.elements.to.deploy.sim, candidate.boundary.elements.to.deploy.sor, null.SDM.b)
  }
}

Sys.time() - ptm











############################################################################################################################
# 12) Test the significance of the superfluidity
#     PART 1 :- COllate the data
############################################################################################################################

subgraphs.NULL.sim <- matrix(NA, nrow=length(bricks), ncol=35082)
subgraphs.NULL.sor <- matrix(NA, nrow=length(bricks), ncol=35082)

# Find those null bricks for which superfluidity has been calculated
sims_w_superfluidity <- which(file.exists(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", bricks, "/superfluidity_sim.csv", sep="")))
sors_w_superfluidity <- which(file.exists(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", bricks, "/superfluidity_sor.csv", sep="")))

# Collate the superfluidity values into matrices
superfluidity.NULL.sim <- matrix(NA, nrow=length(sims_w_superfluidity), ncol=10)
superfluidity.NULL.sor <- matrix(NA, nrow=length(sors_w_superfluidity), ncol=10)

for(x in 1:length(sims_w_superfluidity)){
  superfluidity.NULL.sim[x,] <- read.csv(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", bricks[sims_w_superfluidity][[x]], "/superfluidity_sim.csv", sep=""))[,3]
}
for(x in 1:length(sors_w_superfluidity)){
  superfluidity.NULL.sor[x,] <- read.csv(paste("05_Wombling_Null_Models/03_Null_Beta_Diversity/", bricks[sors_w_superfluidity][[x]], "/superfluidity_sor.csv", sep=""))[,3]
}


# Very occasionally superfluidity will be infinite (i.e. when ALL boundary elements contribute)
# Replace these with 200 for ease of plotting
superfluidity.NULL.sim[superfluidity.NULL.sim==Inf] <- 200
superfluidity.NULL.sor[superfluidity.NULL.sor==Inf] <- 200

# Read in the emperical superfluidity values
superfluidity.OBS.sim <- read.csv("04_Wombling/SuperfluiditySim.txt")
evaulation.number.of.subgraphs.sim <- superfluidity.OBS.sim[,1]
superfluidity.OBS.sim <- superfluidity.OBS.sim[,2]
superfluidity.OBS.sor <- read.csv("04_Wombling/SuperfluiditySor.txt")
evaulation.number.of.subgraphs.sor <- superfluidity.OBS.sor[,1]
superfluidity.OBS.sor <- superfluidity.OBS.sor[,2]

############################################################################################################################
# 13) Test the significance of the superfluidity
#     PART 2 :- Plot the data
############################################################################################################################

plot(evaulation.number.of.subgraphs.sim, superfluidity.NULL.sim[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(evaulation.number.of.subgraphs.sim)), ylim=c(0, max(superfluidity.NULL.sim)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(superfluidity.NULL.sim))
{
  points(evaulation.number.of.subgraphs.sor, superfluidity.NULL.sim[i,], type="l", col="gray70")
}
points(evaulation.number.of.subgraphs.sim, superfluidity.OBS.sim, type="o", col="blue", lwd=2)

plot(evaulation.number.of.subgraphs.sor, superfluidity.NULL.sor[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(evaulation.number.of.subgraphs.sor)), ylim=c(0, max(superfluidity.NULL.sor)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(superfluidity.NULL.sor))
{
  points(evaulation.number.of.subgraphs.sor, superfluidity.NULL.sor[i,], type="l", col="gray70")
}
points(evaulation.number.of.subgraphs.sor, superfluidity.OBS.sor, type="o", col="blue", lwd=2)


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

bricks <- 1:123

# Read in the emperical superfluidity values
no.o.subgraphs.OBS.sim <- read.csv("04_Wombling/NumberSubgraphsSim.txt")
no.o.subgraphs.OBS.sim <- no.o.subgraphs.OBS.sim[,1]
no.o.subgraphs.OBS.sor <- read.csv("04_Wombling/NumberSubgraphsSor.txt")
no.o.subgraphs.OBS.sor <- no.o.subgraphs.OBS.sor[,1]

subgraphs.NULL.sim <- matrix(NA, nrow=length(bricks), ncol=length(no.o.subgraphs.OBS.sim))
subgraphs.NULL.sor <- matrix(NA, nrow=length(bricks), ncol=length(no.o.subgraphs.OBS.sim))

# Find those null bricks for which the number of subgraphs has been calculated
sims_w_no_o_subgraphs <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", bricks, "/sim_number_of_subgraphs.csv", sep="")))
sors_w_no_o_subgraphs <- which(file.exists(paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", bricks, "/sor_number_of_subgraphs.csv", sep="")))

# Collate the superfluidity values into matrices
no.o.subgraphs.NULL.sim <- matrix(NA, nrow=length(sims_w_no_o_subgraphs), ncol=length(no.o.subgraphs.OBS.sim))
no.o.subgraphs.NULL.sor <- matrix(NA, nrow=length(sors_w_no_o_subgraphs), ncol=length(no.o.subgraphs.OBS.sim))

for(x in 1:length(sims_w_no_o_subgraphs)){
  no.o.subgraphs.NULL.sim[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", bricks[sims_w_no_o_subgraphs][[x]], "/sim_number_of_subgraphs.csv", sep=""))[,2]
}
for(x in 1:length(sors_w_no_o_subgraphs)){
  no.o.subgraphs.NULL.sor[x,] <- read.csv(paste("05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/", bricks[sors_w_no_o_subgraphs][[x]], "/sor_number_of_subgraphs.csv", sep=""))[,2]
}



############################################################################################################################
# 15) Test the significance of the number of subgraphs
#     PART 2 :- Plot the data
############################################################################################################################

plot(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.NULL.sim[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(length(no.o.subgraphs.OBS.sim))), ylim=c(0, max(no.o.subgraphs.NULL.sim)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(no.o.subgraphs.NULL.sim))
{
  points(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.NULL.sim[i,], type="l", col="gray70")
}
points(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.OBS.sim, type="o", col="blue", lwd=1)

plot(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.NULL.sor[1,], type="l", col="gray70", 
     bty="n", xlim=c(0, max(length(no.o.subgraphs.OBS.sim))), ylim=c(0, max(no.o.subgraphs.NULL.sor)),
     cex.axis=1.5, cex.lab=1.5, xlab="Evaluation number of subgraphs", ylab="superfluidity")
for(i in 2:nrow(no.o.subgraphs.NULL.sor))
{
  points(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.NULL.sor[i,], type="l", col="gray70")
}
points(1:length(no.o.subgraphs.OBS.sim), no.o.subgraphs.OBS.sor, type="o", col="blue", lwd=1)


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

dir.create("05_Wombling_Null_Models/05_Significance", showWarnings = F)
write.csv(p.no.o.subgraphs.sor, file="05_Wombling_Null_Models/05_Significance/sor_no_o_subgraphs_signif.csv")
write.csv(p.no.o.subgraphs.sim, file="05_Wombling_Null_Models/05_Significance/sim_no_o_subgraphs_signif.csv")

write.csv(no.o.subgraphs.NULL.sor, file="05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/collated_sor.csv")
write.csv(no.o.subgraphs.NULL.sim, file="05_Wombling_Null_Models/04_Null_Number_of_Subgraphs/collated_sim.csv")
