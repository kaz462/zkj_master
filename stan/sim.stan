data {
  int<lower=1> n;
  int<lower=1> N;
  vector[N] pop_tn;
  vector[N] TTime;
  int<lower=1> K; // number of predictors for poisson part
  matrix[N, K] x_lambda; 
  matrix[N, 2] x_theta; 
  matrix<lower = 0, upper = 1>[n, n] W; // adjacency matrix
}

transformed data {
  vector[K] mu = rep_vector(0, K);
  vector[83] zeros = rep_vector(0,83);
  matrix<lower = 0>[n, n] D;
  vector[n] W_rowsums;
  for (i in 1:n) {
    W_rowsums[i] = sum(W[i, ]);
    }
  D = diag_matrix(W_rowsums);
}

parameters {
  vector[2] beta_z;           // coefficients for zero inflated part
  vector[K] beta_m;           // coefficients for poisson part
  real<lower = 0> tau;
  real<lower = 0, upper = 1> alpha;  
  real a;
}

transformed parameters{
  matrix[83,83] Omega;
  cov_matrix[83] Sigma;
  Omega=tau*(D-alpha*W);
  Sigma=inverse_spd(Omega);
}

generated quantities {
  vector[N] sim_theta;
  vector[N] sim_m; 
  vector[N] sim_R;
  int sim_y[N]; 
  vector[n] sim_phi;
  
  sim_phi=multi_normal_rng(zeros,Sigma);
  sim_theta = inv_logit(x_theta * beta_z);
  //sim_R = Phi(TTime*a); 
  sim_R = inv_logit(TTime*a);
  sim_m = exp(x_lambda * beta_m + to_vector(rep_matrix(sim_phi, 14))+log(sim_R)+log(pop_tn));
  
  for (i in 1:N) { 
        if (binomial_rng(1,sim_theta[i])>0)
          sim_y[i] = 0;
        else
          sim_y[i] = poisson_rng(sim_m[i]);
      } 
}
