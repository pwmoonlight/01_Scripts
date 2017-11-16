
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


#define the iterations to perform
iterations.to.perform <- 1000


############################################################################################################################
# 2) 
############################################################################################################################



#read the file that has the brick with the species distribution models of all species
species.in.analysis <-  gsub("\\.", "-", gsub("_", " ", gsub("MCLUST_", "", gsub(".R$", "", list.files("05_Wombling_Null_Models/01_McClust/", pattern="*.R$", full.names=F)))))
species.level.dir <- lapply(1:length(species.in.analysis), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species.in.analysis)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)


AOO <- rep(NA, times=length(species.level.dir))
for(i in 1:length(species.level.dir)){
  print(i/length(species.level.dir)*100)
  AOO[i] <- cellStats(species.level.dir[[i]], 'sum')
}

save(AOO, file=paste("04_Wombling/AOO.R", sep=""))
load(file=paste("04_Wombling/AOO.R", sep=""))



nlayers(SDM.b)
cell.side <- res(SDM.b)[1] #resolution, i.e., cell dimensions (length and width)
focal.sp.range <- raster(SDM.b, layer=1)
no.na.cells <- (1:ncell(focal.sp.range))[!is.na(extract(focal.sp.range, 1:ncell(focal.sp.range)))]
s1.s2 <- coordinates(focal.sp.range)[no.na.cells,]
dim(s1.s2)

load(file=paste("04_Wombling/AOO.R", sep=""))
#plot(AOO)
#summary(AOO)
#which(AOO<5100 & AOO>4900)
#hist(AOO, breaks=100)


MCLUST.model.names <- gsub(" ", "_", paste("MCLUST", species.in.analysis, sep="_"))


############################################################################################################################
# 3)
############################################################################################################################

plot(c(1,200), c(1,100), type="n", bty="n", xlab="Time (hours)", ylab= "Null model iterations", cex.lab=1.5, cex.axis=1.5)

ptm <- Sys.time()

dir.create("05_Wombling_Null_Models/02b_Null_Distributions", showWarnings = F)
dir.create("05_Wombling_Null_Models/03_Null_Bricks", showWarnings = F)

for(null.iterations in iterations.to.perform:1)
  {
  if(!file.exists(paste("05_Wombling_Null_Models/02b_Null_Distributions/", null.iterations, sep="")))
    {
    dir.create(paste("05_Wombling_Null_Models/02b_Null_Distributions/", null.iterations, sep=""), showWarnings = F)
    dir.create(paste("05_Wombling_Null_Models/03_Null_Bricks/", null.iterations, sep=""), showWarnings = F)
  	#the vector "final.j" will store the number of attempts needed to rotate the centroids of each species
  	final.j <- rep(NA, times=length(MCLUST.model.names)) 
  	for(h in length(MCLUST.model.names):1)
  	  {
	    writeLines(paste("Iteration ", null.iterations, " of ", iterations.to.perform, ": Species ", h, " of ", length(MCLUST.model.names), ": ", MCLUST.model.names[[h]], sep=""))
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
    	rm(focal.sp.range)
	    rm(centroid.null)
	  	rm(Mcluster.focal.sp)
    	rm(shifted.rotated.range)		

		  absence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))<1
  	  absence.cells.index.r[is.na(absence.cells.index.r)] <- FALSE
	    absence.cells.r <- (1:ncell(occurrences.r))[absence.cells.index.r]

		  presence.cells.index.r <- extract(occurrences.r, 1:ncell(occurrences.r))>0
  	  presence.cells.index.r[is.na(presence.cells.index.r)] <- FALSE
	    presence.cells.r <- (1:ncell(occurrences.r))[presence.cells.index.r]

	    AOO.fraction.to.add <- as.vector(AOO[h] - length(presence.cells.r))
		
  		# If the whole range falls within the study area, proceed
  		if(AOO.fraction.to.add < 1) writeRaster(occurrences.r, filename=paste("05_Wombling_Null_Models/02b_Null_Distributions/", null.iterations, "/", names(SDM.b)[h], "_NULL.grd", sep=""))
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

  		rm(colo.cells)
  		rm(absence.cells.index.r)
	  	rm(presence.cells.index.r)
  		rm(presence.cells.r)
  		rm(AOO.fraction.to.add)

      # Save the NULL species raster
	  	writeRaster(occurrences.r, filename=paste("05_Wombling_Null_Models/02b_Null_Distributions/", null.iterations, "/", names(SDM.b)[h], "_NULL.grd", sep=""), overwrite=T)
	  	
	  	
	  	#save file with vector "final.j", which stores the number of attempts needed to rotate the centroids of each species
	  	write.table(final.j, file = paste("05_Wombling_Null_Models/Null_", null.iterations, "_final_j.csv", sep=""), quote = TRUE, sep = ",")
	  	
     	}

  	#create an R object of class "brick" with species distribution models for all species,
  	#first create a list of the names of raster objects (with species distribution models)
  	#that will be included in the "brick":
  	#SDM.raster.names <- list.files("05_Wombling_Null_Models/02b_Null_Distributions", pattern="*.grd", full.names = T)

  	#create the null "brick"
    #Null_SDM.b <- brick(lapply(SDM.raster.names, eval))

	  #save in a file the brick with the species distribution models of all species in the analysis
    #random.number <- as.numeric(paste(sample(1:9, 10, replace=T), collapse=""))
    #writeRaster(Null_SDM.b,	filename=paste("05_Wombling_Null_Models/03_Null_Bricks/Null_", random.number, "_", format(Sys.time(), "%d%b%Y_%H%M%S"), "_SDMb.grd", sep=""), bandorder='BIL')
	  #remove null brick
	  #rm(Null_SDM.b)
	  #removeTmpFiles(0) #remove temporary files

    # Delete the individual species rasters to save disk space
  	#file.remove(list.files("05_Wombling_Null_Models/02b_Null_Distributions/", full.names = T))


  	#plot iterations number as a measure of progress
    #points(difftime(Sys.time(), ptm, units="hours"), null.iterations, pch=19, cex=1)
    }
  }

Sys.time() - ptm

