

#define the number of regions (potentially ecoregions or "subgraphs") for whichsuperfluidity will be calculated
evaluation.number.of.subgraphs <- number.of.subgraphs[143638 - beta.ranks.to.evaluate]

#check that all evaluation points exist in the vector "number.of.subgraphs"
match(evaluation.number.of.subgraphs, number.of.subgraphs)  
sum(is.na(match(evaluation.number.of.subgraphs, number.of.subgraphs)))

#calculatesuperfluidity
start.time <- Sys.time()
superfluidity <- rep(NA, times=length(evaluation.number.of.subgraphs))
plot(min(evaluation.number.of.subgraphs), 0,
     xlim=c(min(evaluation.number.of.subgraphs), max(evaluation.number.of.subgraphs)),
     ylim=c(0, 200), 
     bty="n", xlab="Regions (or subgraphs)", ylab="superfluidity", cex.lab=1.5, cex.axis=1.5, type="n")
for (j in evaluation.number.of.subgraphs)
{
  writeLines(paste(which(evaluation.number.of.subgraphs == j), "of", length(evaluation.number.of.subgraphs), "at", Sys.time()))
  candidate.boundary.elements.to.delete <- which(r.obs.beta > (length(r.obs.beta) - match(j, number.of.subgraphs)))
  focal.graph <- delete_edges(Nordeste.graph, candidate.boundary.elements.to.delete)
  #focal.graph
  
  necessary.edge <- rep(NA, times=length(candidate.boundary.elements.to.delete))
  for(i in 1:length(candidate.boundary.elements.to.delete))
  {
    necessary.edge[i] <- ifelse(components(add_edges(focal.graph, cell.adj.char[candidate.boundary.elements.to.delete[i],]))$no < components(focal.graph)$no, 1, 0)
  }
  superfluidity[which(evaluation.number.of.subgraphs==j)] <- sum(necessary.edge<1)/sum(necessary.edge>0)
  points(j,superfluidity[which(evaluation.number.of.subgraphs==j)], pch=19)
}
difftime(Sys.time(), start.time, units="mins")
#this procedure might take about 14 minutes, depending on which computer is used
#
#examine the results
class(superfluidity)
length(superfluidity)
summary(superfluidity)

#plot the results in linear scale
plot(evaluation.number.of.subgraphs,superfluidity, pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Regions (or subgraphs)", ylab="superfluidity")

#plot the results using a semi-log graph
plot(evaluation.number.of.subgraphs, log(superfluidity), pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Regions (or subgraphs)", ylab="Log (superfluidity)")

#plot the results using a log-log graph
plot(log(evaluation.number.of.subgraphs), log(superfluidity), pch=19,
     bty="n", cex.axis=1.5, cex.lab=1.5, type="o",
     xlab="Log (Regions (or subgraphs))", ylab="Log (superfluidity)")
