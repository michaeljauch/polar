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

col_stan <- 'grey'
col_geodesic <- 'black'
col_gibbs <- 'grey35'
lwd_stan <- 1.5
lwd_geodesic <- 1.5
lwd_gibbs <- 5
lty_stan <- 1
lty_geodesic <- 1
lty_gibbs <- 3

pdf(file='~/Desktop/eigenmodel_trace.pdf', width=12, height=8)
par(mfrow=c(3,1))

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(0, 5, 4.5, 3))
plot(stan_lambdas[1:1000,3], type='l', ylim=c(120,160),
     ylab=expression(lambda[1]), lty=lty_stan, lwd=lwd_stan, col=col_stan, xaxt='n')
points(geodesic_lambdas[1:1000,3], type='l', lty=lty_geodesic, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,3], type='l', lty=lty_gibbs, lwd=lwd_gibbs, col=col_gibbs)

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(2.25, 5, 2.25, 3))
plot(stan_lambdas[1:1000,2], type='l', ylim=c(70,115), 
     ylab=expression(lambda[2]), lty=lty_stan, lwd=lwd_stan, col=col_stan, xaxt='n')
points(geodesic_lambdas[1:1000,2], type='l', lty=lty_geodesic, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,2], type='l', lty=lty_gibbs, lwd=lwd_gibbs, col=col_gibbs)

par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(4.5, 5, 0, 3))
plot(stan_lambdas[1:1000,1], type='l', ylim=c(-130,-82), 
     ylab=expression(lambda[3]), lty=lty_stan, lwd=lwd_stan, xlab='Iteration', col=col_stan)
points(geodesic_lambdas[1:1000,1], type='l', lty=lty_geodesic, lwd=lwd_geodesic, col=col_geodesic)
points(gibbs_lambdas[1:1000,1], type='l', lty=lty_gibbs, lwd=lwd_gibbs, col=col_gibbs)
dev.off()
