// This is the Stan code we use to simulate from the matrix Bingham 
// via polar expansion with the default conditional distribution. 
// This code is called in the script matrix_bingham_sim_study.R 

data{
  int<lower=0> p; // Number of rows of the matrix 
  int<lower=1> k; // Number of columns of the matrix 
  cov_matrix[p] SigInv; // Inverse of parameter of matrix Bingham (SPD)
}

parameters{
  matrix[p,k] X; // X is the p x k real matix. 
}

transformed parameters{
  // We calculate the orthogonal matrix Q from X using the polar decomp.
  matrix[p,k] Q; 
  {
  vector[k] eval;
  vector[k] eval_trans;
  matrix[k,k] evec; 
  eval = eigenvalues_sym(X'*X);
  for(l in 1:k){
    eval_trans[l] = 1/sqrt(eval[l]);
  }
  evec = eigenvectors_sym(X'*X); 
  Q = X*evec*diag_matrix(eval_trans)*evec'; 
  }
}

model{
  // Here we specify the log of the density f_X^alt. 
  target += - trace(X'*X/2) - trace(Q'*SigInv*Q/2); 
}
