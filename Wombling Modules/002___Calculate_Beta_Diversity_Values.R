brick.index.species.in.phylogeny <- names(SDM.b)

# Obtain, for all pairs of adjacent grid cells, terms to calculate beta-diversity.
# First, using a matrix indicate:
# species are unique to one of the grid cells, represented as "0",
# species unique to the other grid cell, represented as "Inf",
# and shared species, represented as a "1".
obs.beta.terms <-  SDM.b[cell.adj[,1]][,brick.index.species.in.phylogeny] / SDM.b[cell.adj[,2]][,brick.index.species.in.phylogeny]
dim(obs.beta.terms)

# Next, for all pairs of adjacent grid cells, calculate the following terms, see page 254 in
# Legendre and Legendre (2010, Numerical Ecology, Second Edition):
# the number of shared species 
a <- rowSums(obs.beta.terms>0.5 & obs.beta.terms<1.5, na.rm=T)
summary(a)
# the number of unique species in one of the grid cells 
b <- rowSums(obs.beta.terms<0.5, na.rm=T)
summary(b)
# the number of unique species in the other grid cell
c <- rowSums(is.infinite(obs.beta.terms))
summary(c)

# Calculate ecological distance using Sorensen's index
# (see page 286, equation 7.56 in Legendre and Legendre (2010, Numerical Ecology, Second Edition):
obs.beta.sor <- (b+c)/(2*a + b + c)
obs.beta.sor[which(is.na(obs.beta.sor))] <- 0
summary(obs.beta.sor)

# Calculate ecological distance using Simpson's index
# (see page 2230, equation 2 in Mouillot et al. (2013, Journal of Biogeography 40: 2228–2237):
obs.beta.sim <- pmin(b,c)/(a + pmin(b,c))
obs.beta.sim[which(is.na(obs.beta.sim))] <- 0
summary(obs.beta.sim)

# Graphically compare the Sorensen's and Simpson's indexes
plot(obs.beta.sor, obs.beta.sim, bty="n", cex.axis=1.5, cex.lab=1.5)
abline(0, 1, col="red")

write.table(obs.beta.sor, file="04_Wombling/ObsBetaSor.txt", sep=",")
write.table(obs.beta.sim, file="04_Wombling/ObsBetaSim.txt", sep=",")

