functions {
  matrix polar (matrix X) {
    int k = cols(X);
    matrix[k, k] X_cross = crossprod(X);
    vector[k] eval = eigenvalues_sym( X_cross);
    vector[k] eval_trans = inv_sqrt(eval);
    matrix[k, k] evec = eigenvectors_sym( X_cross ); 
 
   return X * tcrossprod(diag_pre_multiply(sqrt(eval_trans), evec)); 
  }
  
  real lambda(real L, int m) {
		real lam;
		lam = (m*pi()/(2*L))^2;
				
		return lam;
	}
	real spd(real alpha, real rho, real w) {
		real S;
		S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
				
		return S;
	}
	vector phi(real L, int m, vector x) {
		vector[rows(x)] fi;
		fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
				
		return fi;
	}
}

data{
  int<lower=1> n; // Number of rows
  int<lower=1> p; // Number of columns
  matrix[n,p] Y; // Data matrix
  int<lower=1> k; // Number of principal components
  vector[p] days; // Locations of data in time
  real<lower=0> nu; // Hyperparameter of inverse gamma prior for sig2 
  real<lower=0> s2; // Hyperparameter of inverse gamma prior for sig2
  real<lower=0> tau; // SD parameter of prior for D
  real<lower=0> alpha; // Hyperparameter of inverse gamma prior for rho
  real<lower=0> beta; // Hyperparameter of inverse gamma prior for rho
  
  real L;						//boundary condition factor
	int<lower=1> B;				//nº of basis functions		

}
transformed data {
  matrix[p, B] PHI;
  
	for (b in 1:B) 
	  PHI[, b] = phi(L, b, days); 
}

parameters{
  vector[n * k] z_U; 
  vector[B * k] z_V; 
  vector<lower=0>[k] D; 
  real<lower=-1, upper=1> xi; // Autoregressive parameter 
  real<lower=0> sig2; // Marginal variance of autoregressive process
  real<lower=0> rho; // Length scale parameter 
}

transformed parameters{
  matrix[n, k] U = polar(to_matrix(z_U, n, k));
  matrix[p, k] V;
  
  vector[B] diagSPD;
	matrix[B, k] SPD_beta;
	
	for(b in 1:B)
		diagSPD[b] =  sqrt(spd(1, rho, sqrt(lambda(L, b)))); 

  SPD_beta = diag_pre_multiply(diagSPD, to_matrix(z_V, B, k));
	  
  V = polar(PHI * SPD_beta); 
}

model{ 
  matrix[n, p] M = diag_post_multiply(U, D) * V'; 
  
  // Prior specification
  z_V ~ std_normal(); 
  z_U ~ std_normal(); 
  D ~ normal(0, tau);
  sig2 ~ inv_gamma(nu / 2.0, nu * s2 / 2.0);
  rho ~ inv_gamma(alpha, beta); 
  target += -.5 * log1m(square(xi)); // Arc-sine prior for phi (see Fosdick et al. 2012)
  
  // Likelihood
  for (i in 1:n)
    Y[i, 2:p] ~ normal(M[i, 2:p] + xi * Y[i, 1:(p - 1)], sqrt(sig2 * (1 - square(xi))));

}
generated quantities {
  matrix[n, p] M_out = diag_post_multiply(U, D) * V'; 
  real log_lik[n];
  
  for (i in 1:n) {
    log_lik[i] = 0;
    for (j in 2:p) {
     log_lik[i] += normal_lpdf(Y[i, j] | M_out[i, j] + xi * Y[i, j - 1], sqrt(sig2 * (1 - square(xi))));
    }
  }
}