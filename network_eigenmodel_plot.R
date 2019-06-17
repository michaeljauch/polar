# This script creates Figure 1. 

# Load the lambda draws for each of the methods 

load('~/Desktop/PPE/stan_network_eigenmodel_draws.RData')
stan_lambdas <- draws$lam
rm(draws)

geodesic_lambdas <- 
  read.csv('~/Desktop/PPE/geodesic_network_eigenmodel_draws.csv',
           header=FALSE)

gibbs_lambdas <- 
  read.csv('~/Desktop/PPE/gibbs_network_eigenmodel_draws.csv', 
           header=TRUE)


# Rearrange so that the lambdas are decreasing order

stan_lambdas <- stan_lambdas[,c(1,3,2)]

geodesic_lambdas <- geodesic_lambdas[,c(1,3,2)]

gibbs_lambdas <- gibbs_lambdas[,c(1,2,3)]

col_stan <- 'black'
col_geodesic <- 'blue'
col_gibbs <- 'red'
lwd_stan <- .5
lwd_geodesic <- 2
lwd_gibbs <- 3

pdf(file='~/eigenmodel_trace.pdf', width=12, height=6)
par(mfrow=c(3,1))

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(0, 5, 4.5, 3))
plot(stan_lambdas[1:1000,3], type='l', ylim=c(120,160),
     ylab=expression(lambda[1]), lty=1, lwd=lwd_stan, col=col_stan, xaxt='n')
points(geodesic_lambdas[1:1000,3], type='l', lty=2, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,3], type='l', lty=3, lwd=lwd_gibbs, col=col_gibbs)

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(2.25, 5, 2.25, 3))
plot(stan_lambdas[1:1000,2], type='l', ylim=c(70,115), 
     ylab=expression(lambda[2]), lty=1, lwd=lwd_stan, col=col_stan, xaxt='n')
points(geodesic_lambdas[1:1000,2], type='l', lty=2, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,2], type='l', lty=3, lwd=lwd_gibbs, col=col_gibbs)

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(4.5, 5, 0, 3))
plot(stan_lambdas[1:1000,1], type='l', ylim=c(-130,-82), 
     ylab=expression(lambda[3]), lty=1, lwd=lwd_stan, xlab='Iteration', col=col_stan)
points(geodesic_lambdas[1:1000,1], type='l', lty=2, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,1], type='l', lty=3, lwd=lwd_gibbs, col=col_gibbs)
dev.off()
