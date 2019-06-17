functions {
  // Define function which returns a vector w(A) of the 
  // lower triangular elements of A obtained as in Magnus.
  vector lower_tri(matrix A) {
    int m;
    int counter; 
    vector[dims(A)[1]*(dims(A)[1]-1)/2] w_A;
    m = dims(A)[1];
    w_A = rep_vector(0, m*(m-1)/2);
    counter = 1; 
    for (j in 1:(m-1)) {
      for(i in (j+1):m) {
        w_A[counter] = A[i,j];
        counter = counter + 1; 
      }
    }
    return w_A;
  }
  
  vector std_norm_cdf(vector x) {
    int m;
    vector[dims(x)[1]] cdf_vec;
    m = dims(x)[1];
    cdf_vec = rep_vector(0, m);
    for (i in 1:m) {
      cdf_vec[i] = Phi_approx(x[i]);
    }
    return cdf_vec;
  }
}
  
data{
  int<lower=0> p; // Number of vertices in network
  int<lower=0, upper=1> y[p*(p-1)/2]; // Binary network vector
  int<lower=1> k; // Rank of latent matrix
  real<lower=0> sd_lam; // sd of trunc normal prior on lam
  real<lower=0> sd_c; // sd of normal prior on c
}

parameters{
  vector[p*k] z; 
  vector[k] lam; 
  real c; 
}

transformed parameters{
  matrix[p,k] Q; 
  {
  matrix[p,k] X; 
  vector[k] eval;
  vector[k] eval_trans;
  matrix[k,k] evec; 
  X = to_matrix(z, p, k); 
  eval = eigenvalues_sym(X'*X);
  for(l in 1:k){
    eval_trans[l] = 1/sqrt(eval[l]);
  }
  evec = eigenvectors_sym(X'*X); 
  Q = X*evec*diag_matrix(eval_trans)*evec'; 
  }
}

model{
  real ltP[p*(p-1)/2]; // lower triangular entries of probability matrix
  
  // Priors 
  lam ~ normal(0, sd_lam); 
  c ~ normal(0,sd_c);
  z ~ normal(0,1); 

  // Likelihood 
  ltP = to_array_1d(std_norm_cdf(c + lower_tri(quad_form(diag_matrix(lam),Q')))); 
  y ~ bernoulli(ltP);
}
