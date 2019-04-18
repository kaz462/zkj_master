load("data/collec.RData")
x_lambda <- as.matrix(cbind(rep(1,83*14),scale(xx1),scale(xx2)))
x_theta <- as.matrix(cbind(rep(1,83*14),scale(xx1)))
input_data <- list("n"=83,"N"=1162,"pop_tn"=pop_tn,"K"=3, "x_lambda"=x_lambda,
                   "TTime"=TTime,"x_theta"=x_theta,"W"=W,"W_n"=W_n)

params <- list(
  "beta_z" = c(-1.98, 0.57),
  "beta_m" = c(-8.44, 0.17, -1.06),
  "tau" = 0.25,
  "alpha" = 0.5,
  "a" = -0.5
)

##sim_ME.stan
sim_data <- simulate_data(
  file = "stan/sim.stan",
  data_name = "sim",
  input_data = input_data,
  param_values = params,
  vars = c("sim_theta","sim_m","sim_y","sim_phi","sim_R"),
  nsim = 220,
  path = "simulated_a_5_alpha5/",
  seed = 1234,
  use_cores = 6
)