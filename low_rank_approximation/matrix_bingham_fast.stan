// This is the Stan code we use to simulate from the matrix Bingham 
// via polar expansion with the alternative conditional distribution. 
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
  matrix[p, k] Q; 
  matrix[k, k] X_cross = crossprod(X);
  
  {
  vector[k] eval = eigenvalues_sym( X_cross);
  vector[k] eval_trans = inv_sqrt(eval);
  matrix[k, k] evec = eigenvectors_sym( X_cross ); 
 
   Q = X * tcrossprod(diag_pre_multiply(sqrt(eval_trans), evec)); 
  }
}

model{
  target += -0.5 * (trace(X_cross) + trace_quad_form(SigInv, Q) ); 
}
