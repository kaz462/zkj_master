functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}


data {
  int<lower=1> n;
  int<lower=1> N;
  int y[N];
  vector[N] pop_tn;
  vector[N] TTime;
  int<lower=1> K; //number of predictors for poisson part
  matrix[N, 3] x_lambda; 
  matrix[N, 2] x_theta; 
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs
  real gamma1;
  real gamma2;
}


transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[n] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[n] lambda;       // eigenvalues of invsqrtD * W * invsqrtD

  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:n) D_sparse[i] = sum(W[i]);
  {
    vector[n] invsqrtD;  
    for (i in 1:n) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }  
}


parameters {
  vector[K] beta_m;     // coefficients in Poisson part
  vector[2] beta_z;     // coefficients in zero-inflated part
  //real a;               // coefficent for measurement error modeling
  vector[n] phi;
  real<lower = 0> sd_tau;
  real<lower = 0, upper = 1> alpha;  
}

transformed parameters {
  real <lower = 0> tau;
  tau = inv(square(sd_tau));
 }

model{
  vector[N] theta;
  vector[N] m1_log;   
  vector[N] m;  
 // vector[N] R; 
  theta = inv_logit(x_theta * beta_z);
  //R = Phi(TTime*a+  to_vector(rep_matrix(phi, 14)));
  m1_log = x_lambda * beta_m; 
  m=exp(m1_log+log(pop_tn));//+log(R)
  
    for (i in 1:N) { 
        if (y[i] == 0)
          target += log_sum_exp(bernoulli_lpmf(1 | theta[i]),
                        bernoulli_lpmf(0 | theta[i])
                        + poisson_lpmf(y[i] | m[i]));
        else
          target += bernoulli_lpmf(0 | theta[i]) + poisson_lpmf(y[i] | m[i]);
      } 

  alpha ~ beta(gamma1,gamma2);
  sd_tau ~ uniform(0,10);
  phi ~ sparse_car(tau, alpha, W_sparse, D_sparse, lambda, n, W_n);
  beta_m ~ multi_normal(rep_vector(0, 3),diag_matrix(rep_vector(100, 3)));
  beta_z ~ multi_normal(rep_vector(0, 2),diag_matrix(rep_vector(100, 2)));
 // a ~ normal(0,100);
}

