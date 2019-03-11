#March 10, 2019
library(RColorBrewer)
library(maps)
library(e1071)
library(classInt)
library(coda)
#library("R2WinBUGS")
library("sp")
library(maptools)
library(StanHeaders) 
library(ggplot2) 
library(rstan) 
library(spdep)
library(rstan)  
library(glmmBUGS)
library(classInt)
library(maptools)
library(maps)
library(RColorBrewer)
library(xtable)
library(parallel)

cl = makeCluster(4, type = "FORK", outfile = "")
options(mc.cores = 4)  

#setwd("/home/kaz/Jan14")
#############################################################################
#### read data 
## change the directory in data.R to the directory of rawdata file
source("R/data.R")
x_lambda <- as.matrix(cbind(rep(1,83*14),scale(xx1),scale(xx2)))
x_theta <- as.matrix(cbind(rep(1,83*14),scale(xx1)))
#############################################################################





#############################################################################
#### simulation
## glmmTMB
dt <- as.data.frame(cbind(yy,scale(xx1),scale(xx2),county))
library(glmmTMB)
obj <- glmmTMB(yy~scale(xx1)+scale(xx2)+(1|county),zi=~scale(xx1),
               family=poisson,data=dt,offset=log(0.55*pop_tn))
out <- confint(obj)[,3]
##  simulate y with true x
library(rstansim)
library(raster)
input_data <- list("n"=83,"N"=1162,"pop_tn"=pop_tn,"K"=3, "x_lambda"=x_lambda,
                   "TTime"=TTime,"x_theta"=x_theta,"W"=W,"W_n"=W_n)
#params <- list("beta_z"=c(-2,0.6),"beta_m"=c(-8.5,0.06,-1),
#               "tau"=0.3,"alpha"=0.7,"a"=0.7)

#params <- list("beta_z"=c(-2,0.6),"beta_m"=c(-8.5,0.06,-0.8),
#               "tau"=0.3,"alpha"=0.7,"a"=0.75)

#params <- list("beta_z"=c(-1.98,0.57),"beta_m"=c(-8.45,0.06,-1.06),
#               "tau"=0.3,"alpha"=0.7,"a"=0.76)




#params <- list("beta_z"=c(-1.98,0.57),"beta_m"=c(-8.45,0.17,-1.06),
#               "tau"=1,"alpha"=0.7,"a"=0.75)
params <- list("beta_z"=c(-1.98,0.57),"beta_m"=c(-8.45,0.17,-1.06),
               "tau"=0.75,"alpha"=0.7,"a"=0.75)


##sim_ME.stan
sim_data <- simulate_data(
  file = "stan/sim_phi.stan",
  data_name = "sim",
  input_data = input_data,
  param_values = params,
  vars = c("sim_theta","sim_m","sim_y","sim_phi","sim_R"),
  nsim = 200,
  path = "simulated/",
  seed = 1234
)

test_data <- readRDS("simulated/sim_2.rds")


sprintf("simulated/sim_%d.rds", 2)


# range
ysim <- test_data$y; range(ysim)
# proportion
sum(ysim==0)/length(ysim);sum(yy==0)/length(yy)
# mean of simulated y (non-zeros)
mean(ysim[ysim!=0]);mean(yy[yy!=0])
# table
xtable(t(table(ysim)),caption='Table of ysim')
xtable(t(table(yy)),caption='Table of yy')
# hist
par(mfrow=c(1,2))
hist(ysim,xlab="Simulated counts" ,xlim=range(ysim),ylim=c(0,1000),breaks=seq(0,max(ysim),by=1))
hist(yy,xlab="Real counts",xlim=range(ysim),ylim=c(0,1000),breaks=seq(0,max(ysim),by=1))
#############################################################################





#############################################################################
#### model
## parameters for beta prior
beta2 <- function(para)
{ para[3]<-0.025;para[4]<-0.975
sum.sqr<-function(x){
  sum1<-(pbeta(para[1],shape1=x[1],shape2=x[2])-para[3])^2
  sum2<-(pbeta(para[2],shape1=x[1],shape2=x[2])-para[4])^2
  ans<-sum1+sum2}
opt<-optim(c(1,1),sum.sqr,lower=c(0.01,0.01),method="L-BFGS-B")
#opt<-optim(c(1,1),sum.sqr)
out<-list(a=opt$par[1],b=opt$par[2])
print(out)
}

#aa<-0.6;bb<-0.8
sbeta<-beta2(c(0.6,0.8))
gamma1 <- sbeta$a
gamma2 <- sbeta$b


stan_i <- function(i) {
  test_data <- readRDS(sprintf("simulated/sim_%d.rds", i))
  ysim <- test_data$y
  stan("stan/model_phi.stan", data=list(n=83,N=1162,y=ysim,pop_tn=pop_tn,K=3,
                                                  x_theta=x_theta,TTime=TTime,
                                                  x_lambda=x_lambda,W=W,W_n=W_n,
                                                  gamma1=gamma1,gamma2=gamma2),
                 chains=2, warmup=500, iter=800, save_warmup=FALSE,init = inits,
                 control = list(adapt_delta = 0.95,max_treedepth = 15));
}
fit.all <- lapply(seq(5), stan_i)
fit.all <- parLapply(cl = cl, X = seq(4), stan_i)





## stan model
inits1 <- list(beta_m=rep(1,3),beta_z=rep(0,2),phi=rep(0,83),tau=0.5,alpha=0.2,a=0.1)
inits2 <- list(beta_m=rep(-0.1,3),beta_z=rep(0.5,2),phi=rep(1,83),tau=1,alpha=0.6,a=1)
inits <- list(inits1,inits2)



stanfit = stan("stan/model_phi.stan", data=list(n=83,N=1162,y=ysim,pop_tn=pop_tn,K=3,
                                        x_theta=x_theta,TTime=TTime,
                                        x_lambda=x_lambda,W=W,W_n=W_n,
                                        gamma1=gamma1,gamma2=gamma2),
                chains=2, warmup=2000, iter=6000, save_warmup=FALSE,init = inits,
                control = list(adapt_delta = 0.95,max_treedepth = 15));


## output
print(stanfit, digits=2,pars=c("beta_z","beta_m","tau","alpha","a"),
      probs=c(0.025, 0.5, 0.975)) 
trace<-traceplot(stanfit,pars=c("beta_z","beta_m","tau","alpha","a"));trace
pairs(stanfit,pars=c("beta_z","beta_m","a"))
#############################################################################

stanfit.p <- summary(stanfit)

print(stanfit, digits=2,pars=c("beta_z","beta_m","tau","alpha"),
      probs=c(0.025, 0.5, 0.975)) 
trace<-traceplot(stanfit,pars=c("beta_z","beta_m","tau","alpha"))
pairs(stanfit,pars=c("beta_z","beta_m"))

parse_exp <- function(s) {
  parse(text = as.expression(s))
}

var_names <- list(
  "beta_z[1]" = parse_exp("beta[z1]"), 
  "beta_z[2]" = parse_exp("beta[z2]"), 
  "beta_m[1]" = parse_exp("beta[m1]"), 
  "beta_m[2]"= parse_exp("beta[m2]"), 
  "beta_m[3]" = parse_exp("beta[m3]"), 
  "tau" = parse_exp("tau"), 
  "alpha" = parse_exp("alpha"), 
  "a" = parse_exp("a")
)

var_labeller <- function(variable,value){
  return(var_names[value])
}

mcmc_trace(as.array(stanfit), 
           pars = c("beta_z[1]", "beta_z[2]", "beta_m[1]", "beta_m[2]", "beta_m[3]", "tau", "alpha", "a"), 
           facet_args = list(labeller = var_labeller))

