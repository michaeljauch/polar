# This script creates Figures 2-4. 

library(fda)
library(rstan)
library(invgamma)
load('~/fPCA_weather_fit.RDATA')

draws <- rstan::extract(fit)
num_iter <- length(draws$phi)

# We compute the posterior mean M of UDV' 
M <- matrix(rep(0, n*p), nrow=n)
for(t in 1:num_iter){
  if(k==1){
    M <- draws$D[t]* draws$U[t,,] %*% t(draws$V[t,,])/num_iter
  } else{
    M <- draws$U[t,,] %*% diag(draws$D[t,]) %*% t(draws$V[t,,])/num_iter
  }
}

# Our points estimates of U, D, and V are computed via as follows from M: 
svdM <- svd(M, nu=k, nv=k)

############################
# Figure 3: MLE comparison #
############################

pdf(file='~/MLEcomparison.pdf', width=12, height=4)
par(mfrow=c(1,k), cex.axis=1.5, cex.lab=2, cex.main=2, cex.sub=2)
for(i in 1:k){
  plot(days, svdY$v[,i], type='l', xlab='Day', ylab='', lwd=1, 
       col=gray.colors(100)[40])
  if(i == 3){
    points(days, -svdM$v[,i], type='l', lwd=2)
  } else{
    points(days, svdM$v[,i], type='l', lwd=2)
  }
}
dev.off()

###################################################
# Figure 4: histogram of rho and posterior curves #
###################################################

pdf(file='~/miscfig.pdf', width=12, height=6)

par(mfrow=c(1,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)

# Left panel 
hist(draws$rho, breaks=20, main='', xlab=expression(rho), freq=FALSE,
     xlim=c(20,46))
x = seq(20,46,.01)
alpha <- 39.35719
beta <- 1172.206
points(x, dinvgamma(x=x, shape=alpha, rate=beta), type='l', lwd=2)

# Right panel 
plot(days, svdM$v[,3], type='l', xlab='Day', lwd=2, col='white', ylab='')
idx <- seq(1,101, 10)
for(i in idx){
  points(days, draws$V[i,,3], type='l', col=gray.colors(100)[60])
  
}
points(days, svdM$v[,3], type='l', xlab='Day', lwd=2)

dev.off()

##############################################
# Figure 2: raw data and interpretable fPCAs #
##############################################

rawY <- t(CanadianWeather$dailyAv[days,,1]) # Load data 

pdf(file='~/Desktop/PPE/interpfig.pdf', width=12, height=12)

par(mfrow=c(2,2))
ylim_vec <- c(-38, 30)

# Plot raw data in the first panel 
par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
par(mar = c(4.5, 4.5, 4.5, 4.5))
plot(days, rawY[1,], xlab='Day', ylab=expression(paste(degree,"C")), type='l',
     ylim=ylim_vec, main='Raw data')
for(i in 2:n){
  points(days, rawY[i,], type='l')
}

# When plotting +'s and -'s, we only plot every 5 points because it's clearer.
idx <- seq(1, length(days), 5)
reduced_days <- days[idx]

# Plot fPCAs in other 3 panels 
for(i in 1:k){
  par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
  par(mar = c(4.5, 4.5, 4.5, 4.5))
  plot(days, cmeans, xlab='Day', ylab=expression(paste(degree,"C")), 
       type='l', ylim=ylim_vec, main=paste('PC', i, sep=' '))
  points(reduced_days, cmeans[idx] + 50*svdM$v[idx,i], pch='+')
  points(reduced_days, cmeans[idx] - 50*svdM$v[idx,i], pch='-')
}

dev.off()
