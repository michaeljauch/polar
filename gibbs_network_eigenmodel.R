# This script implements a Gibbs sampler for the network eigenmodel of 
# Hoff (2009) applied to the protein interaction data of Butland. It is based on
# code appearing in the Vignette 'Bayesian analysis of matrix data with 
# rstiefel' associated with the rstiefel package. The eigenmodel in the vignette
# uses different variable names compared to the eigenmodel in Hoff (2009). 
# In particular, 'c' in Hoff (2009) becomes 'theta' in the vignette. The main 
# changes I've made to the code from the vignette are the replacement of the 
# function rZ_fc, the addition of comments, and the absense of thinning.  

library(rstiefel)
library(truncnorm)

Y <- as.matrix(read.table('~/hoff.dat', sep=' '))

# R is the rank, t2.lambda is the prior variance of the lambdas, t2.theta is the 
# prior variance of theta. 
R<-3 ; t2.lambda<-dim(Y)[1] ; t2.theta<-100

## Starting values
theta<-qnorm(mean(c(Y),na.rm=TRUE))
L<-diag(0,R)
set.seed(1)
U<-rustiefel(dim(Y)[1],R)

## MCMC
# LPS will store posterior samples of the lambdas, sorted 
# TPS will store posterior sampes of theta
# MPS will store the sum of the posterior samples of U%*%L%*%t(U) 
LPS <- TPS <- NULL; MPS <- matrix(0,dim(Y),dim(Y))
burnin <- 5000
nsamples <- 5000 
total <- burnin + nsamples
for(s in 1:total)
{
 print(s)
 if(s==burnin + 1){
   start_time <- proc.time()
 }
  
  # Update Z 
  #Z<-rZ_fc(Y,theta+U%*%L%*%t(U)) 
  Z <- matrix(rep(0,dim(Y)[1]^2), nrow=dim(Y)[1])
  mean_mat <- theta+U%*%L%*%t(U)
  for(i in 2:dim(Y)[1]){
    for(j in 1:(i-1)){
      if(Y[i,j]==1){
        Z[i,j] <- rtruncnorm(n=1, mean=mean_mat[i,j], sd=1, a=0, b=Inf) 
      } else{
        Z[i,j] <- rtruncnorm(n=1, mean=mean_mat[i,j], sd=1, a=-Inf, b=0)
      }
    }
  }
  Z <- Z + t(Z) + diag(rnorm(n=dim(Y)[1], mean=diag(mean_mat), sd=sqrt(2))) 

  # Update theta
  E<-Z-U%*%L%*%t(U)
  v.theta<-1/(1/t2.theta + choose(dim(Y)[1],2))
  e.theta<-v.theta*sum(E[upper.tri(E)])
  theta<-rnorm(1,e.theta,sqrt(v.theta))

  # Update L 
  E<-Z-theta
  v.lambda<-2*t2.lambda/(2+t2.lambda)
  e.lambda<-v.lambda*diag(t(U)%*%E%*%U/2)
  L<-diag(rnorm(R,e.lambda,sqrt(v.lambda)))

  # Update U 
  U<-rbing.matrix.gibbs(E/2,L,U)
  
## Output
  if(s>burnin)  #& s%%10==0
  {
   LPS<-rbind(LPS,sort(diag(L))) ; TPS<-c(TPS,theta) ; MPS<-MPS+U%*%L%*%t(U)
  }
}

end_time <- proc.time()
time_after_burn <- end_time - start_time

write.csv(LPS,'~/gibbs_network_eigenmodel_draws.csv', 
          row.names=FALSE)
