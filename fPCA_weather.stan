data{
  int<lower=1> n; // Number of rows
  int<lower=1> p; // Number of columns
  matrix[n,p] Y; // Data matrix
  int<lower=1> k; // Number of principal components
  real days[p]; // Locations of data in time
  real<lower=0> nu; // Hyperparameter of inverse gamma prior for sig2 
  real<lower=0> s2; // Hyperparameter of inverse gamma prior for sig2
  real<lower=0> tau; // SD parameter of prior for D
  real<lower=0> alpha; // Hyperparameter of inverse gamma prior for rho
  real<lower=0> beta; // Hyperparameter of inverse gamma prior for rho
}

parameters{
  vector[n*k] z_U; 
  vector[p*k] z_V; 
  vector<lower=0>[k] D; 
  real<lower=-1, upper=1> phi; // Autoregressive parameter 
  real<lower=0> sig2; // Marginal variance of autoregressive process
  real<lower=0> rho; // Length scale parameter 
}

transformed parameters{
  matrix[n,k] U; 
  matrix[p,k] V;
    {
    matrix[p,p] K; // Parameter of MACG prior for V
    matrix[p,p] L_K; // Cholesky decomposition of K  
    matrix[n,k] X_U; // U is the orth. component of the polar decomp. of X_U
    matrix[p,k] X_V; // V is the orth. component of the polar decomp. of X_U
    vector[k] eval_U; // Eigenvalues of X_U'*X_U
    vector[k] eval_trans_U; // Transformation of eigenvalues for polar decomp.
    matrix[k,k] evec_U; // Eigenvectors of X_U'*X_U
    vector[k] eval_V; // Eigenvalues of X_V'*X_V
    vector[k] eval_trans_V;  // Transformation of eigenvalues for polar decomp.
    matrix[k,k] evec_V; // Eigenvectors of X_V'*X_V
    
    // Constructing U as orthogonal component of the polar decomposition of X_U 
    X_U = to_matrix(z_U, n, k); 
    eval_U = eigenvalues_sym(X_U'*X_U);
    for(i in 1:k){
      eval_trans_U[i] = 1.0/sqrt(eval_U[i]);
    }
    evec_U = eigenvectors_sym(X_U'*X_U); 
    U = X_U*evec_U*diag_matrix(eval_trans_U)*evec_U'; 
    
    // Constructing V as orthogonal component of the polar decomposition of X_V 
    K = cov_exp_quad(days, 1, rho) + diag_matrix(rep_vector(1e-7, p));
    L_K = cholesky_decompose(K);
    X_V = L_K*to_matrix(z_V, p, k);
    eval_V = eigenvalues_sym(X_V'*X_V);
    for(i in 1:k){
      eval_trans_V[i] = 1.0/sqrt(eval_V[i]);
    }
    evec_V = eigenvectors_sym(X_V'*X_V); 
    V = X_V*evec_V*diag_matrix(eval_trans_V)*evec_V'; 
  }
}

model{ 
  matrix[n,p] M;
  M = U*diag_matrix(D)*V'; 
  
  // Prior specification
  z_V ~ normal(0, 1.0); 
  z_U ~ normal(0, 1.0); 
  D ~ normal(0, tau);
  sig2 ~ inv_gamma(nu/2.0, nu*s2/2.0);
  rho ~ inv_gamma(alpha, beta); 
  target += -.5*log(1-phi^2); // Arc-sine prior for phi (see Fosdick et al. 2012)
  
  // Likelihood
  for(i in 1:n){
    Y[i,2:p] ~ normal(M[i,2:p] + phi*Y[i,1:(p-1)], sqrt(sig2*(1-phi^2)));
  }
  
}
