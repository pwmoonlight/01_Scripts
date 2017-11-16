###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
###############################################################################################################
###############################################################################################################
###############################################################################################################

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
save(SizeRegions, file=paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".R", sep=""))
#load file with region size
load(paste("04_Wombling/SizeRegions_SOR_", pick.num.subgraphs, ".R", sep=""))

#obtain the coordinates of the center of the cells assigned to each region or subgraph
comp.coor <-  as.list (rep(NA, times=components(modified.Nordeste.graph)$no))
for(j in 1:components(modified.Nordeste.graph)$no)
{
  print(j/components(modified.Nordeste.graph)$no*100)
  comp.coor[[j]] <- xyFromCell(Nordeste.mask.0, as.numeric(names(components(modified.Nordeste.graph)$membership))[components(modified.Nordeste.graph)$membership==j], spatial=FALSE)
}

#save the coordinates of the center of the cells in each region or subgraph
save(comp.coor, file=paste("04_Wombling/CompCoor_SOR_", pick.num.subgraphs, ".R", sep=""))

#load the coordinates of the center of the cells in each region or subgraph
load(paste("04_Wombling/CompCoor_SOR_", pick.num.subgraphs, ".R", sep=""))

#obtain cell number (or cell id) for each cell assigned to each region or subgraph
comp.cells <-  as.list (rep(NA, times=components(modified.Nordeste.graph)$no))
for(j in 1:components(modified.Nordeste.graph)$no)
{
  print(j/components(modified.Nordeste.graph)$no*100)
  comp.cells[[j]] <- as.numeric(names(components(modified.Nordeste.graph)$membership))[components(modified.Nordeste.graph)$membership==j]
}

#create a raster of regions or subgraphs
Nordeste.mask.Regions <- Nordeste.mask.0
Nordeste.mask.Regions[] <- NA
Nordeste.mask.Regions[unlist(comp.cells)] <- rep(1:length(comp.cells), times=lapply(comp.cells, length))

#save the raster of regions or subgraphs
writeRaster(Nordeste.mask.Regions, paste("04_Wombling/Nordeste_SOR_", pick.num.subgraphs,"_Regions.grd", sep=""))
