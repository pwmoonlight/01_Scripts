

###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
###############################################################################################################
###############################################################################################################
###############################################################################################################



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

rm(cell.neig, flag.pairs.to.retain, from.coor, to.coor, no.adj.coor)
