
r.obs.beta <- rank(obs.beta[,1], ties.method="random")
class(r.obs.beta)
length(r.obs.beta)
summary(r.obs.beta, digits=6)
max(r.obs.beta)
length(unique(r.obs.beta))
r.obs.beta[1:100]

#define the percentiles of beta-diversity that will be evaluated against the null model
beta.ranks.to.evaluate <- 143638 - round(143638*seq(0.05, 0.5, 0.05))

plot(r.obs.beta, obs.beta[,1], bty="n", cex.axis=1.5, cex.lab=1.5)
abline(v= beta.ranks.to.evaluate, lty=3)
abline(h= obs.beta[match(beta.ranks.to.evaluate, r.obs.beta),1], lty=3)
#
plot(r.obs.beta, obs.beta[,1], bty="n", cex.axis=1.5, cex.lab=1.5, xlim=c(6000,143638), ylim=c(0,0.1))
abline(v= beta.ranks.to.evaluate, lty=3)
abline(h= obs.beta[match(beta.ranks.to.evaluate, r.obs.beta),1], lty=3)

