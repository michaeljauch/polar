data{
  int<lower=0> p; // Number of rows
  int<lower=1> k; // Number of columns
  cov_matrix[p] SigInv; // Inverse of parameter of matrix Bingham (SPD)
}

parameters{
  matrix[p,k] X; 
}

transformed parameters{
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
    target += - trace(X'*X/2) - trace(Q'*SigInv*Q/2); 
}
