library(rstan)
library(cowplot)

###########################################################Function###########################################################
## results ME model
results_ME <- function(file){
  fit_all <- lapply(1:200, function(i) {extract(readRDS(sprintf("%s/stanfit_%d.rds", file, i)))})
  fit_sample <- lapply(fit_all, function(x) {list(beta_z1=x$beta_z[,1], beta_z2=x$beta_z[,2], 
                                                  beta_m1=x$beta_m[,1], beta_m2=x$beta_m[,2], 
                                                  beta_m3=x$beta_m[,3], a=x$a, alpha=x$alpha, 
                                                  tau=x$tau)})
  est <- sapply(fit_sample, function(x) {as.vector(sapply(x, mean))})
  est <- data.frame(apply(t(est), 2, mean))
  rownames(est) <- c("beta_z1", "beta_z2", "beta_m1", "beta_m2", "beta_m3", "a", "alpha", "tau")
  colnames(est) <- "Estimation"
  est
}

## results Naive model
results_Naive <- function(file){
  fit_all <- lapply(1:200, function(i) {extract(readRDS(sprintf("%s/stanfit_%d.rds", file, i)))})
  fit_sample <- lapply(fit_all, function(x) {list(beta_z1=x$beta_z[,1], beta_z2=x$beta_z[,2], 
                                                  beta_m1=x$beta_m[,1], beta_m2=x$beta_m[,2], 
                                                  beta_m3=x$beta_m[,3], alpha=x$alpha, 
                                                  tau=x$tau)})
  est <- sapply(fit_sample, function(x) {as.vector(sapply(x, mean))})
  est <- data.frame(apply(t(est), 2, mean))
  rownames(est) <- c("beta_z1", "beta_z2", "beta_m1", "beta_m2", "beta_m3", "alpha", "tau")
  colnames(est) <- "Estimation"
  est
}

## bias plot 
plot_bias <- function(Difference, alpha, a){
  model <-c(rep("ME Model", 7), rep("Naive Model", 7))
  Parameters <- rep(c("a", "b", "c", "d", "e", "f", "g"), 2)
  mydata <- data.frame(Parameters, Difference, model)
  p <-ggplot(mydata, aes(Parameters, Difference, fill=model))
  p +geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette="Paired") +
    theme_bw()+ scale_x_discrete(labels = c(expression(beta[z1]),
                                            expression(beta[z2]),
                                            expression(beta[m1]),
                                            expression(beta[m2]),
                                            expression(beta[m3]),
                                            expression(alpha),
                                            expression(tau))) +
    ylim(0,0.6) + ylab("Relative Bias (in absolute value)") +  
    xlab("") + ggtitle(bquote(list(a==.(a), alpha==.(alpha)))) +
    theme(legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5))
}
#############################################################################################################################





###########################################################Results###########################################################
a0_alpha7 <- results(file="out_a0_alpha7")
a0_alpha7_ME <- results_ME(file="out_a0_alpha7_ME")
a2_alpha7 <- results(file="out_a2_alpha7")
a2_alpha7_ME <- results_ME(file="out_a2_alpha7_ME")
a5_alpha7 <- results(file="out_a5_alpha7")
a5_alpha7_ME <- results_ME(file="out_a5_alpha7_ME")

a0_alpha4 <- results(file="out_a0_alpha4")
a0_alpha4_ME <- results_ME(file="out_a0_alpha4_ME")
a2_alpha4 <- results(file="out_a2_alpha4")
a2_alpha4_ME <- results_ME(file="out_a2_alpha4_ME")
a5_alpha4 <- results(file="out_a5_alpha4")
a5_alpha4_ME <- results_ME(file="out_a5_alpha4_ME")

a0_alpha1 <- results(file="out_a0_alpha1")
a0_alpha1_ME <- results_ME(file="out_a0_alpha1_ME")
a2_alpha1 <- results(file="out_a2_alpha1")
a2_alpha1_ME <- results_ME(file="out_a2_alpha1_ME")
a5_alpha1 <- results(file="out_a5_alpha1")
a5_alpha1_ME <- results_ME(file="out_a5_alpha1_ME")

true <- c(-1.98, 0.57, -8.44, 0.17, -1.06, 0.5, 0.25)
Difference0_7 <- abs(c((true-t(as.vector(a0_alpha7))[-6])/true, (true-t(as.vector(a0_alpha7_ME))[1:7])/true))
Difference2_7 <- abs(c((true-t(as.vector(a2_alpha7))[-6])/true, (true-t(as.vector(a2_alpha7_ME))[1:7])/true))
Difference5_7 <- abs(c((true-t(as.vector(a5_alpha7))[-6])/true, (true-t(as.vector(a5_alpha7_ME))[1:7])/true))
Difference7_5 <- abs(c((true-t(as.vector(a7_alpha5))[-6])/true, (true-t(as.vector(a7_alpha5_ME))[1:7])/true))

true <- c(-1.98, 0.57, -8.44, 0.17, -1.06, 0.5, 0.25)
Difference0_4 <- abs(c((true-t(as.vector(a0_alpha4))[-6])/true, (true-t(as.vector(a0_alpha4_ME))[1:7])/true))
Difference2_4 <- abs(c((true-t(as.vector(a2_alpha4))[-6])/true, (true-t(as.vector(a2_alpha4_ME))[1:7])/true))
Difference5_4 <- abs(c((true-t(as.vector(a5_alpha4))[-6])/true, (true-t(as.vector(a5_alpha4_ME))[1:7])/true))
Difference1_5 <- abs(c((true-t(as.vector(a1_alpha5))[-6])/true, (true-t(as.vector(a1_alpha5_ME))[1:7])/true))

true <- c(-1.98, 0.57, -8.44, 0.17, -1.06, 0.5, 0.25)
Difference0_1 <- abs(c((true-t(as.vector(a0_alpha1))[-6])/true, (true-t(as.vector(a0_alpha1_ME))[1:7])/true))
Difference2_1 <- abs(c((true-t(as.vector(a2_alpha1))[-6])/true, (true-t(as.vector(a2_alpha1_ME))[1:7])/true))
Difference5_1 <- abs(c((true-t(as.vector(a5_alpha1))[-6])/true, (true-t(as.vector(a5_alpha1_ME))[1:7])/true))
Difference_5_5 <- abs(c((true-t(as.vector(a_5_alpha5))[-6])/true, (true-t(as.vector(a_5_alpha5_ME))[1:7])/true))

bias0_7 <- plot_bias(Difference0_7, alpha = 0.7, a=0)
bias2_7 <- plot_bias(Difference2_7, alpha = 0.7, a=0.2)
bias5_7 <- plot_bias(Difference5_7, alpha = 0.7, a=0.5)

bias0_4 <- plot_bias(Difference0_4, alpha = 0.4, a=0)
bias2_4 <- plot_bias(Difference2_4, alpha = 0.4, a=0.2)
bias5_4 <- plot_bias(Difference5_4, alpha = 0.4, a=0.5)

bias0_1 <- plot_bias(Difference0_1, alpha = 0.1, a=0)
bias2_1 <- plot_bias(Difference2_1, alpha = 0.1, a=0.2)
bias5_1 <- plot_bias(Difference5_1, alpha = 0.1, a=0.5)


plot_grid(bias0_7, bias2_7, bias5_7, 
          bias0_4, bias2_4, bias5_4, 
          bias0_1, bias2_1, bias5_1, 
          get_legend(bias2_1 + theme(legend.position=c(2.7, 0.85))),
          ncol = 3)
######################################################################################################################