###############################################################################################################
 ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

# Obtain, for all pairs of adjacent grid cells, terms to calculate beta-diversity.
# First, using a matrix indicate:
# species are unique to one of the grid cells, represented as "0",
# species unique to the other grid cell, represented as "Inf",
# and shared species, represented as a "1".
obs.beta.terms <-  SDM.b[cell.adj[,1]][,brick.index.species.in.phylogeny] / SDM.b[cell.adj[,2]][,brick.index.species.in.phylogeny]

# Next, for all pairs of adjacent grid cells, calculate the following terms, see page 254 in
# Legendre and Legendre (2010, Numerical Ecology, Second Edition):
# the number of shared species 
a <- rowSums(obs.beta.terms>0.5 & obs.beta.terms<1.5, na.rm=T)
# the number of unique species in one of the grid cells 
b <- rowSums(obs.beta.terms<0.5, na.rm=T)
# the number of unique species in the other grid cell
c <- rowSums(is.infinite(obs.beta.terms))

# Calculate ecological distance using Sorensen's index
# (see page 286, equation 7.56 in Legendre and Legendre (2010, Numerical Ecology, Second Edition):
obs.beta.sor <- (b+c)/(2*a + b + c)
obs.beta.sor[which(is.na(obs.beta.sor))] <- 0

# Calculate ecological distance using Simpson's index
# (see page 2230, equation 2 in Mouillot et al. (2013, Journal of Biogeography 40: 2228–2237):
obs.beta.sim <- pmin(b,c)/(a + pmin(b,c))
obs.beta.sim[which(is.na(obs.beta.sim))] <- 0

write.csv(obs.beta.sim, file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sim.csv", sep=""))
write.csv(obs.beta.sor, file=paste("05_Wombling_Null_Models/04_Null_Beta_Diversity/", bricks[[x]], "/sor.csv", sep=""))

rm(obs.beta.terms, a, b, c)