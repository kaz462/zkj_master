library(rstan)
library(cowplot)
## read in results "stanfit" and "stanfit_ME" from real data analysis
  post <- extract(stanfit)
  post_ME <- extract(stanfit_ME)
  beta_z1 <- post$beta_z[,1]
  beta_z2 <- post$beta_z[,2]
  beta_m1 <- post$beta_m[,1]
  beta_m2 <- post$beta_m[,2]
  beta_m3 <- post$beta_m[,3]
  
  beta_z1_ME <- post_ME$beta_z[,1]
  beta_z2_ME <- post_ME$beta_z[,2]
  beta_m1_ME <- post_ME$beta_m[,1]
  beta_m2_ME <- post_ME$beta_m[,2]
  beta_m3_ME <- post_ME$beta_m[,3]  
  
## plot posterior distributions with and without adjustment of ME
####################################################### m1
  ME_model <- data.frame(beta_m1=beta_m1)
  Naive_model <- data.frame(beta_m1=beta_m1_ME)
  ME_model$model <- "ME_model"
  Naive_model$model <- "Naive_model"
  data <- rbind(ME_model,Naive_model)
  
  m1<-ggplot(data, aes(beta_m1, fill = model)) + 
    geom_histogram(data=ME_model, aes(y = ..density..),
                   bins = 30,
                   colour = "grey",
                   fill = "white")+
    geom_density(alpha = 0.2, adjust = 2)+
    xlab(expression(beta[m1])) +
    theme_bw() +
    theme(legend.position = "none")
####################################################### m2
  ME_model <- data.frame(beta_m2=beta_m2)
  Naive_model <- data.frame(beta_m2=beta_m2_ME)
  ME_model$model <- "ME_model"
  Naive_model$model <- "Naive_model"
  data <- rbind(ME_model,Naive_model)

  m2<-ggplot(data, aes(beta_m2, fill = model)) + 
    geom_histogram(data=ME_model, aes(y = ..density..),
                   bins = 30,
                   colour = "grey",
                   fill = "white")+
    geom_density(alpha = 0.2, adjust = 2)+
    xlab(expression(beta[m2])) +
    theme_bw() +
    theme(legend.position = "none")  
####################################################### m3
  ME_model <- data.frame(beta_m3=beta_m3)
  Naive_model <- data.frame(beta_m3=beta_m3_ME)
  ME_model$model <- "ME_model"
  Naive_model$model <- "Naive_model"
  data <- rbind(ME_model,Naive_model)  
  
  m3<-ggplot(data, aes(beta_m3, fill = model)) + 
    geom_histogram(data=ME_model, aes(y = ..density..),
                   bins = 30,
                   colour = "grey",
                   fill = "white")+
    geom_density(alpha = 0.2, adjust = 2)+
    xlab(expression(beta[m3])) +
    theme_bw() +
    theme(legend.position = "none")  
####################################################### z1
  ME_model <- data.frame(beta_z1=beta_z1)
  Naive_model <- data.frame(beta_z1=beta_z1_ME)
  ME_model$model <- "ME_model"
  Naive_model$model <- "Naive_model"
  data <- rbind(ME_model,Naive_model)  
  
  z1 <- ggplot(data, aes(beta_z1, fill = model)) + 
    geom_histogram(data=ME_model, aes(y = ..density..),
                   bins = 30,
                   colour = "grey",
                   fill = "white")+
    geom_density(alpha = 0.2, adjust = 2)+
    xlab(expression(beta[z1])) +
    theme_bw() +
    theme(legend.position = "none")  
####################################################### z2
  ME_model <- data.frame(beta_z2=beta_z2)
  Naive_model <- data.frame(beta_z2=beta_z2_ME)
  ME_model$model <- "ME_model"
  Naive_model$model <- "Naive_model"
  data <- rbind(ME_model,Naive_model)  
  
  z2 <- ggplot(data, aes(beta_z2, fill = model)) + 
    geom_histogram(data=ME_model, aes(y = ..density..),
                   bins = 30,
                   colour = "grey",
                   fill = "white")+
    geom_density(alpha = 0.2, adjust = 2)+
    xlab(expression(beta[z2])) +
    theme_bw() +
    theme(legend.position = "none")    

## arrange plots
  plot_grid(m1, m2, m3,
            z1, z2,
            get_legend(m1 + theme(legend.position=c(0.6, 0.5))),
            ncol = 3)

