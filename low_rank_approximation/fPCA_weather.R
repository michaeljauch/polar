# install.packages("fda")
library(cmdstanr)
library(fda)
library(posterior)
library(invgamma)
library(loo)
data("CanadianWeather")

days <- seq(1,365,1) #Days we include in the analysis.
k <- 3 # Number of principal components
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
            'tau'=tau, 'alpha'=alpha, 'beta'= beta,
            L = 5/2 * max(days),
            B = 30)

mod_fpca_approx <- cmdstan_model("../polar/fPCA_approx.stan", force_recompile = T)

fit_mod_approx <- mod_fpca_approx$sample(data = dat,
                           iter_warmup = 400,
                           iter_sampling = 400,
                           adapt_delta = 0.999,
                           max_treedepth = 13,
                           parallel_chains = 4)

mod_var <- mod_fpca_approx$variational(data = dat,
                                       iter = 10000,
                                       init = 0.1,
                                       adapt_engaged = T,
                                       algorithm = "lowrank",
                                       rank = 6,
                                       mcse_cut = 0.02,
                                       rhat_cut = 1.01,
                                       eval_window = 75,
                                       ess_cut = 50,
                                       window_size = 0.8,
                                       seed = 14892
)


# defaults
# ess_cut = 20
# eval_window = 100
# mcse_cut = 0.02
# rhat_cut = 1.2
# window_size = 0.5

fit <- mod$variational(
  data = data_list,
  algorithm = "lowrank",
  rank = 2,
  eval_window = 200,
  window_size = 0.7,
  rhat_cut = 1.2,
  mcse_cut = 0.01,
  ess_cut = 50,
  iter = 50000
)


odraws <- mod_var$draws(c("z_U", "z_V", "xi", "sig2", "rho"))
init_approx <- sapply(c("z_U", "z_V", "xi", "sig2", "rho"),
                function(variable) {as.numeric(subset(odraws, variable=variable))})


printed_info <- capture.output(
  mod_fpca_approx$variational(data = dat, iter = 100000,
                              seed = 12489123
  )
)


mod_var$cmdstan_diagnose()
mod_var_opt <- mod_fpca_approx$optimize(data = dat,
                                    init = function() { init_approx },
                                    iter = 100000,
                                   # seed = 2343,
                                    algorithm = "bfgs"
                                   # tol_rel_grad = 1e5
                                    # tol_obj = 1e-13,
                                    # tol_grad = 1e-9,
                                    # tol_param = 1e-9
                                    
)
mod_var$metadata()

mod_var$cmdstan_diagnose()

log_lik <- mod_var$draws("log_lik", format = "draws_array")
loo_out <- loo(log_lik, r_eff = relative_eff(log_lik))
print(loo_out)

# We compute the posterior mean M of UDV' 
M <- mod_var$draws("M_out", format = "list")
# M_mean <- summarise_draws(M)
M_mat <- matrix(apply(as_draws_matrix(M), 2, mean), n, p)

matplot(t(M_mat), type='l', xlab='Years', ylab='rate')
legend('bottomright', inset=.05, legend=colnames(m1), 
       pch=1, horiz=TRUE)

svdM <- svd(M_mat, nu=k, nv=k)

par(mfrow=c(1,k), cex.axis=1.5, cex.lab=2, cex.main=2, cex.sub=2)
for(i in 1:k){
  plot(days, svdY$v[,i], type='l', xlab='Day', ylab='', lwd=1, 
       col=gray.colors(100)[40])
  if(i == 3){
    points(days, -svdM$v[,i], type='l', lwd=2)
    # points(days, V_mat[, 3], type='l', lwd=2)
  } else{
    points(days, svdM$v[,i], type='l', lwd=2)
  }
}

draws_rho <- mod_var$draws("rho")
draws_V <- mod_var$draws("V", format = "list")
V_mat <- matrix(apply(as_draws_matrix(draws_V), 2, mean), p, k)

par(mfrow=c(1,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)

# Left panel 
hist(draws_rho, breaks=20, main='', xlab=expression(rho), freq=FALSE,
     xlim=c(20,46))
x = seq(20,46,.01)
alpha <- 39.35719
beta <- 1172.206
points(x, dinvgamma(x=x, shape=alpha, rate=beta), type='l', lwd=2)

# Right panel 
plot(days, -svdM$v[,3], type='l', xlab='Day', lwd=2, col='white', ylab='')

points(days, -V_mat[, 3], type='l')
points(days, -svdM$v[,3], type='l', xlab='Day', lwd=2)



file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod <- cmdstan_model(file)
data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))
fit <- mod$variational(
  data = data_list
)
