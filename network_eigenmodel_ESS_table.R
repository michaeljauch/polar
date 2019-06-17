# This script computes the effective sample sizes in Table 1. 

library(mcmcse)

# Load the lambda draws for each of the methods 

load('~/stan_network_eigenmodel_draws.RData')
stan_lambdas <- draws$lam
rm(draws)

geodesic_lambdas <- 
  read.csv('~/geodesic_network_eigenmodel_draws.csv',
           header=FALSE)

gibbs_lambdas <- 
  read.csv('~/gibbs_network_eigenmodel_draws.csv', 
           header=TRUE)

# Rearrange so that the lambdas are decreasing order

stan_lambdas <- stan_lambdas[,c(1,3,2)]

geodesic_lambdas <- geodesic_lambdas[,c(1,3,2)]

gibbs_lambdas <- gibbs_lambdas[,c(1,2,3)]

# Compute ess for each of the methods

ess(stan_lambdas)/nrow(stan_lambdas)

ess(geodesic_lambdas)/nrow(geodesic_lambdas)

ess(gibbs_lambdas)/nrow(gibbs_lambdas)
