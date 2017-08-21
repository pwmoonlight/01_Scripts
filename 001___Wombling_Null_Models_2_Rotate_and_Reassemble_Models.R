
############################################################################################################################
############################################################################################################################
#
# INTRODUCTION TO THE SCRIPT "SDMData_30arc_rotate_range_for_NullModel_2016June21.R"
#
# A) What does this script do?
#
# B) What is needed to run this script?
#
#
############################################################################################################################
############################################################################################################################


############################################################################################################################
# 1) Load pakages, define working directories and the iterations to perform
############################################################################################################################

setwd("G:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

#load pakages
library(sp)
library(raster)
#library(igraph)
#library(geosphere)
#library(ellipse)
#library(mclust)

#define directory to load the brick with SDMs
#path.for.brick <- "C:/_transfer/Papers/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels" #from Ivan's laptop
#path.for.brick <- "J:/Jimenez/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels" #from desktop 1ZTF at the Lehmann
#path.for.brick <- "H:/Jimenez/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels" #from desktop 1ZTF at the Lehmann

#define directory to datasets (used to load list of species in the analysis)
#path.for.datasets <- "C:/_transfer/Papers/Nicaragua_Biomes/Datasets" #from Ivan's laptop
#path.for.datasets <- "J:/Jimenez/Nicaragua_Biomes/Datasets" #from desktop 1ZTF at the Lehmann
#path.for.datasets <- "H:/Jimenez/Nicaragua_Biomes/Datasets" #from desktop 1ZTF at the Lehmann

#define directory to load results from MCLUST models
#path.for.MCLUST <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/MCLUST_models" #from Ivan's laptop
#path.for.MCLUST <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/MCLUST_models" #from desktop 1ZTF at the Lehmann
#path.for.MCLUST <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/MCLUST_models" #from desktop 1ZTF at the Lehmann

#define directory to null distributions
# A
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/A_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/A_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/A_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# B
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/B_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/B_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/B_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# C
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/C_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/C_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/C_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# D
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/D_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/D_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/D_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# E
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/E_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/E_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/E_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# F
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/F_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/F_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/F_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# G
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/G_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/G_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/G_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# H
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/H_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/H_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/H_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# I
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/I_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/I_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/I_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# J
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/J_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/J_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/J_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# K
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/K_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/K_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/K_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# L
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/L_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/L_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/L_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# M
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/M_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/M_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/M_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# N
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/N_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/N_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/N_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# O
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/O_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/O_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/O_NULL_Distributions" #from desktop 1ZTF at the Lehmann
# P
#path.for.NullDistributions <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/P_NULL_Distributions" #from Ivan's laptop
#path.for.NullDistributions <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/P_NULL_Distributions" #from desktop 1ZTF at the Lehmann
#path.for.NullDistributions <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/P_NULL_Distributions" #from desktop 1ZTF at the Lehmann


#define directory to null bricks
#path.for.NullBricks <- "C:/_transfer/Papers/Nicaragua_Biomes/NullDistributions/NULL_bricks" #from Ivan's laptop
#path.for.NullBricks <- "J:/Jimenez/Nicaragua_Biomes/NullDistributions/NULL_bricks" #from desktop 1ZTF at the Lehmann
#path.for.NullBricks <- "H:/Jimenez/Nicaragua_Biomes/NullDistributions/NULL_bricks" #from desktop 1ZTF at the Lehmann

#define the iterations to perform
iterations.to.perform <- 100


############################################################################################################################
# 2) 
############################################################################################################################

setwd(path.for.brick)
#dir()
#SDM.b <- brick("SDMb_2015Oct6.grd")
SDM.b <- brick("04_Wombling/SDM_brick.grd")
SDM.b
nlayers(SDM.b)
#names(SDM.b)
cell.side <- res(SDM.b)[1] #resolution, i.e., cell dimensions (length and width)
focal.sp.range <- raster(SDM.b, layer=1)
#create a vector with the cell number (or cell ID) of cells that are not-NA
no.na.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range)))]
#length(no.na.cells)
#dim(coordinates(focal.sp.range)[no.na.cells,])
#points(coordinates(focal.sp.range)[no.na.cells,])
s1.s2 <- coordinates(focal.sp.range)[no.na.cells,]
dim(s1.s2)
AOO <- cellStats(SDM.b, 'sum')
#plot(AOO)
#summary(AOO)
#which(AOO<5100 & AOO>4900)
#hist(AOO, breaks=100)

#setwd(path.for.datasets)
#species.in.analysis <- read.table("SpeciesInAnalysis_2016May21.txt", header=T, sep=",")
species.in.analysis <- read.table("04_Wombling/Species_In_Analysis.txt", sep=",")
species.in.analysis[1:5,]
dim(species.in.analysis)

#species.genus.family <- 
#paste(
#unlist(lapply(strsplit(as.vector(species.in.analysis[,1]), split=c(" "), fixed = T), function(x) paste(x[1], x[2], sep="_"))),
#species.in.analysis[,2], sep="_")
#class(species.genus.family)
#species.genus.family[1:5]
#length(species.genus.family)

MCLUST.model.names <- gsub(" ", "_", paste("MCLUST", species.in.analysis[,1], sep="_"))


############################################################################################################################
# 3)
############################################################################################################################

#h <- 1
#plot(focal.sp.range, colNA="gray70", main=species.genus.family[h], sub="Observed geographic range", cex.sub=1.5)
#dev.new()
#dev.set(2)

plot(c(1,200), c(1,100), type="n", bty="n", xlab="Time (hours)", ylab= "Null model iterations", cex.lab=1.5, cex.axis=1.5)

#ptm <- proc.time()
ptm <- Sys.time()

dir.create("05_Wombling_Null_Models/02_Null_Distributions", showWarnings = F)
dir.create("05_Wombling_Null_Models/03_Null_Bricks", showWarnings = F)

for(null.iterations in 1:iterations.to.perform)
{
	#the vector "final.j" will store the number of attempts needed to rotate the centroids of each species
	final.j <- rep(NA, times=length(MCLUST.model.names)) 
	#h <- 1
	for(h in 1: length(MCLUST.model.names))
	{
		#load Mclust model and assign to:
		#setwd(path.for.MCLUST)
		load(paste("05_Wombling_Null_Models/01_McClust/", MCLUST.model.names[[h]], ".R", sep=""))
 		Mcluster.focal.sp <- eval(as.name(MCLUST.model.names[h]))
		#rm(Mcluster.focal.sp)
		rm(list=c(MCLUST.model.names[h]))
		
		#dev.set(2)
		focal.sp.range <- raster(SDM.b, layer=h)
		#plot(focal.sp.range, colNA="gray70", main=paste(species.genus.family[h], h, "of", length(MCLUST.model.names)), sub="Observed geographic range, centroids and vcv", cex.sub=1.5)
		occurrence.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range))) & extract(focal.sp.range, 1:ncell(focal.sp.range))>0]
		occu.coor <- coordinates(focal.sp.range)[occurrence.cells,]
		rm(occurrence.cells)
		#length(occurrence.cells)
		#dim(occu.coor)
		#points(t(Mcluster.focal.sp$means), col="black", pch=19)
		#for (i in 1:Mcluster.focal.sp$G)
		#	{
		#		points( 
		#		ellipse(Mcluster.focal.sp$sigmas[,,i], centre = Mcluster.focal.sp$means[,i], level = 0.95, npoints = 100),
		#		type="l", col="black")
		#	}
		#dev.new()
		#dev.set(3)

		###############
		# Shift and rotate centroids of range
		###############

		#define the number of (shifted and rotated) simulated ranges to be obtained "k"
		k <- 1
		#define the number of iterations "m" to be performed
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
			#points(centroid.null[1], centroid.null[2], pch=19, cex=1)
			#legend("topright", paste("searching range =", i), bg="white")

			#center means of the focal species
			centered.means <- scale(t(Mcluster.focal.sp$means), center = TRUE, scale = FALSE)
			#points(centered.means, pch=19, cex=0.5, col="red")

			#create rotation matrix for an angle "theta.means" in radians, the angle is randomly sampled from a uniform distribution between 0 and 2*pi
			theta.means <- runif(1, min=0, max=2*pi)
			rotation.matrix.means <- matrix(c(cos(theta.means), sin(theta.means), -sin(theta.means), cos(theta.means)), nrow = 2, ncol = 2)

			#shift and rotate means
			shifted.rotated.means <- t(rotation.matrix.means %*% t(centered.means) + centroid.null)
			#points(shifted.rotated.means, pch=19, cex=0.5, col="red")
			
			#update the count of j
			j <- j+1

			#if not all shifted and rotated localities fall within the study region go to the next iteration,
			#otherwise place the shifted and rotated localities in the position i of the list "shifted.rotated.ranges"   
			if(sum(is.na(extract(focal.sp.range,  shifted.rotated.means)))>0) next
			i <- 2
		}
		
		final.j[h] <- j 

		#dev.set(3)
		#plot(focal.sp.range, col="transparent", colNA="gray70", main=species.genus.family[h], sub="Shifted rotated centroids", cex.sub=1.5)
		#points(shifted.rotated.means, pch=19, cex=0.5, col="red")

		###############
		# Shift and rotate range, based on null centroid and random angle obtained above
		###############
		
		centered.occu.coor <-  occu.coor - matrix(rep(colMeans(t(Mcluster.focal.sp$means)), times=nrow(occu.coor)), ncol=2, byrow=T)
		rm(occu.coor)
		shifted.rotated.range <- t(rotation.matrix.means %*% t(centered.occu.coor) + centroid.null)
		shifted.rotated.range <- shifted.rotated.range[!is.na(extract(focal.sp.range,  shifted.rotated.range)),]
		occurrences.r <- focal.sp.range
		occurrences.r[no.na.cells] <- 0
		occurrences.r <- rasterize(shifted.rotated.range, occurrences.r, rep(1, times=nrow(shifted.rotated.range)), update=T)
		rm(focal.sp.range)
		rm(centroid.null)
		rm(Mcluster.focal.sp)
		rm(shifted.rotated.range)		

		absence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))<1
		#summary(absence.cells.index.r)
		absence.cells.index.r[is.na(absence.cells.index.r)] <- FALSE
		absence.cells.r <- (1:ncell(occurrences.r))[absence.cells.index.r]
		#length(absence.cells.r)
		#str(absence.cells.r)
		#coo.absence.cells.r <- coordinates(occurrences.r)[absence.cells.index.r,]
		#dim(coo.absence.cells.r)
		#plot(occurrences.r)
		#points(coo.absence.cells.r, cex=0.5, pch=19)

		presence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))>0
		#summary(presence.cells.index.r)
		presence.cells.index.r[is.na(presence.cells.index.r)] <- FALSE
		presence.cells.r <- (1:ncell(occurrences.r))[presence.cells.index.r]
		#length(presence.cells.r)
		#str(presence.cells.r)
		#coo.presence.cells.r <- coordinates(occurrences.r)[presence.cells.index.r,]
		#dim(coo.presence.cells.r)
		#plot(occurrences.r)
		#points(coo.presence.cells.r, cex=0.5, pch=19)

		#plot(occurrences.r, add=T)
		#dim(occu.coor)
		#dim(centered.occu.coor)
		#length(shifted.rotated.range.cells)
		#dim(shifted.rotated.range)
		#points(shifted.rotated.range)
		#plot(occurrences.r, colNA="gray70", main=species.genus.family[h], sub="Initial rotated range and centroids", cex.sub=1.5)
		#points(shifted.rotated.means, pch=19, cex=1, col="red")
		AOO.fraction.to.add <- as.vector(AOO[h] - length(presence.cells.r))
		#if(AOO.fraction.to.add < 1) setwd(path.for.NullDistributions)
		if(AOO.fraction.to.add < 1) writeRaster(occurrences.r, filename=paste(names(SDM.b)[h], "_NULL.grd", sep=""))
		if(AOO.fraction.to.add < 1) next 

		while(AOO.fraction.to.add > 0)
		{
			colo.cells <- adjacent(occurrences.r, presence.cells.r, directions=4, pairs=F, target=absence.cells.r, sorted=FALSE, include=FALSE, id=FALSE) 
			#length(colo.cells)
			#str(colo.cells)
			coo.colo.cells <- xyFromCell(occurrences.r, colo.cells, spatial=FALSE)
			if(length(colo.cells) <= AOO.fraction.to.add) occurrences.r <- rasterize(matrix(coo.colo.cells, nrow=length(colo.cells), ncol=2), occurrences.r, rep(1, times=length(colo.cells)), update=T)
			if(length(colo.cells) > AOO.fraction.to.add) additional.occurrences <- sample(1:nrow(coo.colo.cells), size= AOO.fraction.to.add)
			if(length(colo.cells) > AOO.fraction.to.add) occurrences.r <- rasterize(matrix(coo.colo.cells[additional.occurrences,], nrow=length(additional.occurrences), ncol=2), occurrences.r, rep(1, times=length(additional.occurrences)), update=T)

			absence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))<1
			#summary(absence.cells.index.r)
			absence.cells.index.r[is.na(absence.cells.index.r)] <- FALSE
			absence.cells.r <- (1:ncell(occurrences.r))[absence.cells.index.r]
			#length(absence.cells.r)
			#str(absence.cells.r)
			#coo.absence.cells.r <- coordinates(occurrences.r)[absence.cells.index.r,]
			#dim(coo.absence.cells.r)
			#plot(occurrences.r)
			#points(coo.absence.cells.r, cex=0.5, pch=19)

			presence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))>0
			#summary(presence.cells.index.r)
			presence.cells.index.r[is.na(presence.cells.index.r)] <- FALSE
			presence.cells.r <- (1:ncell(occurrences.r))[presence.cells.index.r]
			#length(presence.cells.r)
			#str(presence.cells.r)
			#coo.presence.cells.r <- coordinates(occurrences.r)[presence.cells.index.r,]
			#dim(coo.presence.cells.r)
			#plot(occurrences.r)
			#points(coo.presence.cells.r, cex=0.5, pch=19)

			AOO.fraction.to.add <- as.vector(AOO[h] - length(presence.cells.r))
		}

		rm(colo.cells)
 		rm(absence.cells.index.r)
		rm(presence.cells.index.r)
		rm(presence.cells.r)
		rm(AOO.fraction.to.add)


		#plot(occurrences.r, colNA="gray70", main=species.genus.family[h], sub="Final rotated range and centroids", cex.sub=1.5)
		#points(shifted.rotated.means, pch=19, cex=1, col="red")
		#cellStats(occurrences.r, 'sum')
	
		#setwd(path.for.NullDistributions)
		writeRaster(occurrences.r, filename=paste("05_Wombling_Null_Models/02_Null_Distributions/", names(SDM.b)[h], "_NULL.grd", sep=""))
		rm(occurrences.r)
	}

	#create an R object of class "brick" with species distribution models for all species,
	#first create a list of the names of raster objects (with species distribution models)
	#that will be included in the "brick":
	#setwd(path.for.NullDistributions)
	SDM.raster.names <- list.files("05_Wombling_Null_Models/02_Null_Distributions", pattern="*.grd", full.names = T)


	#create the null "brick"
	#lapply(SDM.raster.names, eval) #checking this bit of code works
	#length(lapply(SDM.raster.names, eval)) #checking it is of the right length
	Null_SDM.b <- brick(lapply(SDM.raster.names, eval))
	#check the result
	#class(Null_SDM.b)
	#res(Null_SDM.b)
	#dim(Null_SDM.b)
	#Null_SDM.b

	#save in a file the brick with the species distribution models of all species in the analysis
	#setwd(path.for.NullBricks)
	random.number <- as.numeric(paste(sample(1:9, 10, replace=T), collapse=""))
	writeRaster(Null_SDM.b,	filename=paste("05_Wombling_Null_Models/03_Null_Bricks/Null_", random.number, "_", format(Sys.time(), "%d%b%Y_%H%M%S"), "_SDMb.grd", sep=""), bandorder='BIL')
	#remove null brick
	rm(Null_SDM.b)
	removeTmpFiles(0) #remove temporary files
	#showTmpFiles()
	#remove null distribution files
	#setwd(path.for.NullDistributions)
	file.remove(list.files("05_Wombling_Null_Models/02_Null_Distributions/", full.names = T))
	#save file with vector "final.j", which stores the number of attempts needed to rotate the centroids of each species
	write.table(final.j, file = paste("05_Wombling_Null_Models/Null_", random.number, "_", format(Sys.time(), "%d%b%Y_%H%M%S"), "_final_j", sep=""), quote = TRUE, sep = ",")
	#plot iterations number
	points(difftime(Sys.time(), ptm, units="hours"), null.iterations, pch=19, cex=1)
}

#proc.time() - ptm
Sys.time() - ptm

############################################################################################################################
# 4)
############################################################################################################################

#read the file that has the brick with the species distribution models of all species
#setwd("C:/_transfer/Papers/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels") #from Ivan's laptop
#setwd("J:/Jimenez/Nicaragua_Biomes/Models/All_species_in_a_brick_BestModels") #from desktop 1ZTF at the Lehmann
#dir()
#SDM.b <- brick("SDMb_2015Oct6.grd")
#SDM.b
#nlayers(SDM.b)
#names(SDM.b)

