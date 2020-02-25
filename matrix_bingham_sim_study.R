set.seed(1)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan_file_default <- '~/matrix_bingham_default.stan'
stan_file_alt <- '~/matrix_bingham_alt.stan'

chains <- 4 # Number of Markov chains run in parallel during a repetition
warmup <- 1500 # Number of warmup iterations per chain 
num_iter <- 2500 # Number of post warmup iterations per chain 
nrep <- 10 # The number of repetitions for each combination of dimensions 

# The combinations of dimensions which we will try
dim_list <- 
  list(
    list(p=3, k=1, rho=1),
    list(p=3, k=2, rho=1),
    list(p=3, k=3, rho=1),
    list(p=3, k=1, rho=2),
    list(p=3, k=2, rho=2),
    list(p=3, k=3, rho=2),
    list(p=10, k=1, rho=1),
    list(p=10, k=5, rho=1),
    list(p=10, k=10, rho=1),
    list(p=10, k=1, rho=2),
    list(p=10, k=5, rho=2),
    list(p=10, k=10, rho=2)
    )

# Functions to construct correlation matrix
sqexpcorr <- function(dist, rho){
  return(exp(-dist^2/(2*rho^2)))
}

constructSig <- function(p, rho){
  Sig <- matrix(nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sig[i,j] <- sqexpcorr(abs(i-j), rho)
    }
  }
  return(Sig)
}

# Initialize a data frame to store the results
results_df <- 
  data.frame('p'=integer(), 
             'k'=integer(), 
             'rho'=double(),
             'conditional'=character(), 
             'Rhat'=double(),
             'n_eff'=double(), 
             'min_eBFMI'=double(),
             'divergences'= integer(),
             stringsAsFactors = FALSE)


# Loop through the combination of dimensions 
for(i in 1:length(dim_list)){
  
  
  p <- dim_list[[i]]$p
  k <- dim_list[[i]]$k
  rho <- dim_list[[i]]$rho
  
  # Construct Sig & SigInv matrices 
  
  Sig <- constructSig(p, rho)
  Sigchol <- chol(Sig)
  SigInv <- chol2inv(Sig)
  
  # Simulate nrep times for both choices of conditional dist. 
  dat <- list('p'=p, 'k'=k, 'SigInv'=SigInv)
  
  for(rep in 1:nrep){
    
    print(c(p,k,rho,rep))
    
    fit_default <- stan(file = stan_file_default, data = dat, 
                        iter = num_iter + warmup, chains = chains, refresh=250, 
                        thin=1, warmup=warmup, seed=rep, 
                        control=list(max_treedepth=13, adapt_delta = 0.999))
    
    divergences_default <- 0
    for(i in 1:chains){
      divergences_defaults <- divergences_default + 
        sum(get_sampler_params(fit_default, inc_warmup=FALSE)[[i]][,'divergent__'])
    }
    
    results_df <- 
      rbind(results_df, 
            data.frame('p'=p, 'k'=k, 'rho'=rho, conditional='default', 
                       'Rhat'=summary(fit_default)$summary['lp__','Rhat'], 
                       'n_eff'=summary(fit_default)$summary['lp__','n_eff'],
                       'min_eBFMI'=min(get_bfmi(fit_default)),
                       'divergences'=divergences_default))

    fit_alt <- stan(file = stan_file_alt, data = dat, iter = num_iter + warmup, 
                    chains = chains, refresh=250, thin=1, warmup=warmup, seed=rep,
                    control=list(max_treedepth=13, adapt_delta = 0.999))
    
    divergences_alt <- 0
      for(i in 1:chains){
        divergences_alt <- divergences_alt + 
          sum(get_sampler_params(fit_alt, inc_warmup=FALSE)[[i]][,'divergent__'])
      }
    
    results_df <- 
      rbind(results_df, 
            data.frame('p'=p, 'k'=k, 'rho'=rho, conditional='alt', 
                       'Rhat'=summary(fit_alt)$summary['lp__','Rhat'], 
                       'n_eff'=summary(fit_alt)$summary['lp__','n_eff'],
                       'min_eBFMI'=min(get_bfmi(fit_alt)),
                       'divergences'=divergences_alt))

  }
}

#save(results_df, file='~/results_df_second_run.RData')

library(sqldf)
query <- 
"
SELECT
  p,
  k,
  rho,
  conditional,
	min(n_eff/10000),
	median(n_eff/10000),
	max(Rhat),
	min(min_eBFMI),
	max(divergences)
FROM 
	results_df
GROUP BY
	p,
	k,
	rho,
	conditional
ORDER BY
	p, 
	k, 
	rho,
	conditional
"

query_results <- sqldf(query)

#save(query_results, file='~/query_results_second_run.RData')
#load('~/query_results_second_run.RData')

# Round query results 
query_results$`min(n_eff/10000)` <- round(query_results$`min(n_eff/10000)`, digits=2)
query_results$`median(n_eff/10000)` <- round(query_results$`median(n_eff/10000)`, digits=2)
query_results$`max(Rhat)` <- round(query_results$`max(Rhat)`, digits=3)
query_results$`min(min_eBFMI)` <- round(query_results$`min(min_eBFMI)`, digits=2)

# Low dim. moderate correlation table: 
query_results[which(query_results$p==3 & query_results$rho==1),]

# Low dim. high correlation table: 
query_results[which(query_results$p==3 & query_results$rho==2),]

# Higher dim. moderate correlation table: 
query_results[which(query_results$p==10 & query_results$rho==1),]

# Higher dim. high correlation table: 
query_results[which(query_results$p==10 & query_results$rho==2),]

