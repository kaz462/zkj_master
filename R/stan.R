#March 10, 2019
#library(spdep)
#library(raster)

library(rstan) 
library(rstansim)
library(bayesplot)
library(xtable)
library(parallel)

#cl = makeCluster(4, type = "FORK", outfile = "")
#options(mc.cores = 4)
rstan_options(auto_write = TRUE)


#############################################################################
#### read data 
## uncomment the following when running in terminal
## setwd("..")

#source("R/data.R")
load("data/collec.RData")
x_lambda <- as.matrix(cbind(rep(1,83*14),scale(xx1),scale(xx2)))
x_theta <- as.matrix(cbind(rep(1,83*14),scale(xx1)))
#############################################################################





#############################################################################
#### simulation
## library(glmmTMB)
# dat <- as.data.frame(cbind(yy,scale(xx1),scale(xx2),county))
# tmb.fit <-
#   glmmTMB::glmmTMB(
#     yy ~ scale(xx1) + scale(xx2) + (1 | county),
#     zi =  ~ scale(xx1),
#     family = poisson,
#     data = dat,
#     offset = log(0.55 * pop_tn)
#   )
# confint(tmb.fit)

summary_ysim <- function(ysim) {
  # range
  range(ysim)
  # proportion
  sum(ysim == 0) / length(ysim)
  sum(yy == 0) / length(yy)
  # mean of simulated y (non-zeros)
  mean(ysim[ysim != 0])
  mean(yy[yy != 0])
  # table
  xtable(t(table(ysim)),caption='Table of ysim')
  xtable(t(table(yy)),caption='Table of yy')
  # hist
  par(mfrow=c(1,2))
  hist(
    ysim,
    xlab = "Simulated counts" ,
    xlim = c(0, max(ysim) + 10),
    ylim = c(0, 1000),
    breaks = seq(0, max(ysim), by = 1)
  )
  hist(
    yy,
    xlab = "Real counts",
    xlim = c(0, max(ysim) + 10),
    ylim = c(0, 1000),
    breaks = seq(0, max(ysim), by = 1)
  )
}


#############################################################################





#############################################################################
#### model
## parameters for beta prior
beta2 <- function(para) {
  para[3] <- 0.025
  para[4] <- 0.975
  sum.sqr <- function(x) {
    sum1 <- (pbeta(para[1], shape1 = x[1], shape2 = x[2]) - para[3])^2
    sum2 <- (pbeta(para[2], shape1 = x[1], shape2 = x[2]) - para[4])^2
    ans <- sum1 + sum2
  }
  opt <- optim(c(1, 1), sum.sqr, lower = c(0.01, 0.01), method = "L-BFGS-B")
  #opt<-optim(c(1,1),sum.sqr)
  list(a = opt$par[1], b = opt$par[2])
}

#aa<-0.6;bb<-0.8
#sbeta <- beta2(c(0.08,0.12))
#sbeta <- beta2(c(0.6,0.8))
sbeta <- beta2(c(0.7,0.9))
#sbeta
gamma1 <- sbeta$a
gamma2 <- sbeta$b


## stan model
inits1 <-
  list(
    beta_m = rep(1, 3),
    beta_z = rep(0, 2),
    phi = rep(0, 83),
    tau = 0.5,
    alpha = 0.2,
    a = 0.1
  )
inits2 <-
  list(
    beta_m = rep(-0.1, 3),
    beta_z = rep(0.5, 2),
    phi = rep(1, 83),
    tau = 1,
    alpha = 0.6,
    a = 1
  )
inits <- list(inits1,inits2)

system ("if [ ! -d out_a8_alpha1 ]; then mkdir out_a8_alpha1; fi") 

for (i in 1L:200) {
  test_data <- readRDS(sprintf("simulated_a8_alpha1/sim_%d.rds", i))
  ysim <- test_data$y
  stanfit <-
    stan(
      "stan/model.stan",
      data = list(
        n = 83, N = 1162, y = ysim, pop_tn = pop_tn,
        K = 3, x_theta = x_theta, TTime = TTime,
        x_lambda = x_lambda, W = W, W_n = W_n,
        gamma1 = gamma1, gamma2 = gamma2
      ),
      chains = 2, warmup = 1000, iter = 2000,
      save_warmup = FALSE, init = inits,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
  
  saveRDS(stanfit, file = sprintf("out_a8_alpha1/stanfit_%d.rds", i))
  gc(verbose = F)
}

## output
#print(stanfit_ME, digits=2,pars=c("beta_z","beta_m","tau","alpha"),
#       probs=c(0.025, 0.5, 0.975)) 
#print(stanfit, digits=2,pars=c("beta_z","beta_m","tau","alpha","a"),
#      probs=c(0.025, 0.5, 0.975)) 
#trace <- traceplot(stanfit,pars=c("beta_z","beta_m","tau","alpha","a"));trace
#pairs(stanfit_ME,pars=c("beta_z","beta_m"));
#pairs(stanfit,pars=c("beta_z","beta_m","a"))
#############################################################################
#fit_all <- lapply(1:144, function(i) {extract(readRDS(sprintf("out/stanfit_%d.rds", i)))})
#fit_sample <- lapply(fit_all, function(x) {list(beta_z1=x$beta_z[,1], beta_z2=x$beta_z[,2], 
#                                                beta_m1=x$beta_m[,1], beta_m2=x$beta_m[,2], 
#                                                beta_m3=x$beta_m[,3], a=x$a, alpha=x$alpha, 
#                                                tau=x$tau)})
#est <- sapply(fit_sample, function(x) {as.vector(sapply(x, mean))})
#apply(t(est), 2, mean)
