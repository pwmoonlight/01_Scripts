#plot the mask, broad scale
plot(Nordeste.mask.0, col="gray70", useRaster=T, legend=F)

selected.quantile <- seq(0.95, 0.5, -0.05)
quantile(obs.beta[,1], probs=selected.quantile[1])
values.at.above.quantile <- sum(obs.beta[,1] >= quantile(obs.beta[,1], probs=selected.quantile[1]))
values.at.above.quantile

flag.candidate.boundary.elements <- r.obs.beta > (length(r.obs.beta)-values.at.above.quantile)
sum(flag.candidate.boundary.elements)

for(i in 1:length(selected.quantile))
{
  values.at.above.quantile <- sum(obs.beta[,1] >= quantile(obs.beta[,1], probs=selected.quantile[i]))
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