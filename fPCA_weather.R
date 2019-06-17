set.seed(1)

# This script calls stan to perform posterior simulation for the Bayesian 
# functional principal component analysis of the Canadian Weather data set, 
# as described in the paper. 

output_file <- '~/fPCA_weather_fit.RDATA'
stan_file <- '~/fPCA_weather.stan'
days <- seq(1,365,1) #Days we include in the analysis.
k <- 3 # Number of principal components
chains <- 4 # Number of Markov chains run
warmup <- 300 # Number of warmup iterations per chain 
num_iter <- 1250 # Number of post warmup iterations per chain 

# The matrices U,D, and V (and other) parameters in the PCA model are only 
# identifiable up to simultaneous permutations and sign changes. When running 
# multiple chains, chains may end up in distinct but symmetric modes. In this case, 
# convergence diagnostics for these parameters will indicate a 
# lack of convergence. The matrix M = U D V' is identifiable, as are sig2, phi,
# and rho. 

library(fda) # Includes CanadianWeather data set 
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Y <- t(CanadianWeather$dailyAv[days,,1]) # Load data 
n <- nrow(Y) # Number of rows (weather stations)
p <- ncol(Y) # Number of columns (days) 

# Subtract column means
cmeans <- colMeans(Y)
Y <- Y - rep(1,n) %*% t(cmeans)

# Subtract row means 
rmeans <- rowMeans(Y)
Y <- Y - rmeans %*% t(rep(1,p)) 

# The priors for sig2 and the diagonal elements of D are data dependent, as 
# described in the paper, so we need to compute their MLEs (under the model 
# without the autoregressive process) and other quantities.

# Compute squared Frobenius norm of Y
Ynorm2 <- norm(Y, type="F")^2

# Compute MLE of sig2
svdY <- svd(Y, nu=k, nv=k)
Yhat <- svdY$u %*% diag(svdY$d[1:k]) %*% t(svdY$v)
sig2hat <- sum((Y - Yhat)^2)/(n*p)

# Set prior parameters for sig2 and D: 
nu <- 1
s2 <- (nu/2 +1)/(nu/2)*sig2hat # (nu, s2) are parameters of inv gamma for sig2
tau <- sqrt((Ynorm2 - n*p*sig2hat)/k) # sd of prior for diag elements of D 

# alpha and beta are the parameters of the inverse gamma prior for rho. 
# They are chosen so that the inverse gamma prior has mean 365/(2*pi*2) and 
# standard deviation 5. 

alpha <- 39.35719
beta <- 1172.206

# The following list stores all the data we need to pass to the Stan sampler. 
dat <- list('n'=n, 'p'=p, 'Y'=Y, 'k'=k, 'days'=days, 'nu'=nu, 's2'=s2, 
            'tau'=tau, 'alpha'=alpha, 'beta'= beta)

fit <- stan(file = stan_file, data = dat, iter = num_iter + warmup, 
            chains = chains, refresh=1, thin=1, warmup=warmup, 
            control=list(max_treedepth=13, adapt_delta = 0.99))

save(list=c('fit', 'days', 'n', 'p', 'k', 'Y', 'rmeans', 'cmeans', 'svdY'),
     file=output_file)




