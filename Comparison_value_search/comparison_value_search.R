###########################################################################
# Description:  This R script allows reproducing the results and plots in 
#               Section 5 of Supplemental_simulations.pdf
#
# User instructions: The folder 'Reproduce_simulations_results' should be
#               saved on the user's computer. The variable mainpath should
#               be replaced by the location of that folder on the user's
#               computer.
#
#
# Last updated: March 11, 2019
###########################################################################

remove(list = ls())
library(ggplot2)
library(ggpubr)
library(hydroGOF)
library(DTRreg)

# mainpath refers to the location of the folder 'Reproduce_simulations_results'
mainpath <- "C:/Users/Gabrielle/Google Drive/McGill - PhD/Thesis/dw-SM/JASA_updatedSERA/Supplementary Material/Reproduce_simulations_results/"
secpath <- "Comparison_value_search/"
setwd(paste(mainpath, secpath, sep = ""))

expit <- function(x) exp(x) / (1 + exp(x))
largen <- 10000


#### ---- Comparison DWSurv and dynamic MSM with true blip (-0.8,0.9) ---- ####
theta1 <- c(3.7, 1.5, 0.8, -0.8, 0.9)
dmsm1 <- dsurv1 <- matrix(NA, nrow = 1000, ncol = 2)
for(n in c(500,1000,5000)){
  print(n)
  dmsm1 <- dsurv1 <- matrix(NA, nrow = 1000, ncol = 2)
  for(i in 1:1000){
    # simulate data
    X11 <- runif(n, 0.5, 1.5)
    X12 <- rbinom(n, 1, 0.6)
    A1 <- rbinom(n, 1, expit(2*X11 - 1))
    delta <- rbinom(n, 1, expit(X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X11[delta == 1] + theta1[3]*X12[delta == 1] + theta1[4]*A1[delta == 1] + theta1[5]*A1[delta == 1]*X11[delta == 1] + rnorm(sum(delta), sd = 0.3)
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    mydata <- data.frame(X11, X12, A1, delta, Y)
    
    # dynamic msm
    # estimate IPCW and IPTW
    tm <- glm(A1 ~ X11, family = binomial)
    pt <- predict(tm, type = "response", data = mydata)
    to <- glm(delta ~ X12, family = binomial)
    mydata$po <- po <- predict(to, type = "response", data = mydata)
    mydata$w <- w <- 1/(A1*pt + (1-A1)*(1-pt))*1/((po*delta) + (1-po)*(1-delta))
    
    # estimate theta for a sequence of value theta I[X1 > theta]
    eyopt <- rep(NA, length(seq(0.5,1.5,by=0.01)))
    ltheta <- seq(0.5,1.5,by=0.01)
    count <- 1
    for(theta in seq(0.5,1.5,by=0.01)){
      tempdat <- mydata
      tempdat$Aopt <- ifelse(tempdat$X11 > theta, 1, 0)
      tempdat1 <- tempdat[which(tempdat$A1 == tempdat$Aopt & tempdat$delta == 1),]
      tempdat1$nw <- tempdat1$w/sum(tempdat1$w)
      eyopt[count] <- sum(tempdat1$Y*tempdat1$nw)
      count <- count + 1
    }
    dmsm1[i,2] <- mytheta <- ltheta[which(eyopt == max(eyopt))][1]
    
    # dwsurv
    # fit model and extract estimates
    model <- dWSurv(time = list(~Y), blip.mod = list(~ X11), treat.mod = list(A1 ~ X11), tf.mod = list( ~ X11 + X12), cens.mod = list(delta ~ X12), var.estim = "none", data = mydata)
    blip <- as.numeric(coef(model)$`stage1`)
    dsurv1[i,2] <- -blip[1]/blip[2]
    
    # simulate large dataset with optimal treatment
    X11 <- runif(largen, 0.5, 1.5) 
    X12 <- rbinom(largen, 1, 0.6)
    A1dwsurv <- ifelse(X11 > (dsurv1[i,2]), 1, 0)
    A1dmsm <- ifelse(X11 > (dmsm1[i,2]), 1, 0)
    error <- rnorm(largen, sd = 0.3)
    logTdwsurv <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1dwsurv + theta1[5]*A1dwsurv*X11 + error
    logTdmsm <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1dmsm + theta1[5]*A1dmsm*X11 + error
    dsurv1[i,1] <- mean(logTdwsurv)
    dmsm1[i,1] <- mean(logTdmsm)
  }
  if(n==500){
    fplot1 <- as.data.frame(rbind(dmsm1, dsurv1))
  }else{
    fplot1 <- as.data.frame(rbind(fplot1,dmsm1, dsurv1))
  }
}


#### ---- Comparison DWSurv and dynamic MSM with true blip (-0.1,0.2) ---- ####
theta1 <- c(3.7, 1.5, 0.8, -0.15, 0.2)
dmsm2 <- dsurv2 <- matrix(NA, nrow = 1000, ncol = 2)
for(n in c(500,1000,5000)){
  print(n)
  dmsm2 <- dsurv2 <- matrix(NA, nrow = 1000, ncol = 2)
  for(i in 1:1000){
    # simulate data
    X11 <- runif(n, 0.5, 1.5)
    X12 <- rbinom(n, 1, 0.6)
    A1 <- rbinom(n, 1, expit(2*X11 - 1))
    delta <- rbinom(n, 1, expit(X12 + 0.1))
    logT <- theta1[1] + theta1[2]*X11[delta == 1] + theta1[3]*X12[delta == 1] + theta1[4]*A1[delta == 1] + theta1[5]*A1[delta == 1]*X11[delta == 1] + rnorm(sum(delta), sd = 0.3)
    C <- rexp(n - sum(delta), rate = 1/300)
    Y <- rep(NA, n)
    Y[delta == 1] <- exp(logT)
    Y[delta == 0] <- C
    mydata <- data.frame(X11, X12, A1, delta, Y)
    
    # dynamic msm
    # estimate IPCW and IPTW
    tm <- glm(A1 ~ X11, family = binomial)
    pt <- predict(tm, type = "response", data = mydata)
    to <- glm(delta ~ X12, family = binomial)
    mydata$po <- po <- predict(to, type = "response", data = mydata)
    mydata$w <- w <- 1/(A1*pt + (1-A1)*(1-pt))*1/((po*delta) + (1-po)*(1-delta))
    
    # estimate theta for a sequence of value theta I[X1 > theta]
    eyopt <- rep(NA, length(seq(0.5,1.5,by=0.01)))
    ltheta <- seq(0.5,1.5,by=0.01)
    count <- 1
    for(theta in seq(0.5,1.5,by=0.01)){
      tempdat <- mydata
      tempdat$Aopt <- ifelse(tempdat$X11 > theta, 1, 0)
      tempdat1 <- tempdat[which(tempdat$A1 == tempdat$Aopt & tempdat$delta == 1),]
      tempdat1$nw <- tempdat1$w/sum(tempdat1$w)
      eyopt[count] <- sum(tempdat1$Y*tempdat1$nw)
      count <- count + 1
    }
    dmsm2[i,2] <- mytheta <- ltheta[which(eyopt == max(eyopt))][1]
    
    # dwsurv
    # fit model and extract estimates
    model <- dWSurv(time = list(~Y), blip.mod = list(~ X11), treat.mod = list(A1 ~ X11), tf.mod = list( ~ X11 + X12), cens.mod = list(delta ~ X12), var.estim = "none", data = mydata)
    blip <- as.numeric(coef(model)$`stage1`)
    dsurv2[i,2] <- -blip[1]/blip[2]
    
    # simulate large dataset with optimal treatment
    X11 <- runif(largen, 0.5, 1.5) 
    X12 <- rbinom(largen, 1, 0.6)
    A1dwsurv <- ifelse(X11 > (dsurv2[i,2]), 1, 0)
    A1dmsm <- ifelse(X11 > (dmsm2[i,2]), 1, 0)
    error <- rnorm(largen, sd = 0.3)
    logTdwsurv <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1dwsurv + theta1[5]*A1dwsurv*X11 + error
    logTdmsm <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1dmsm + theta1[5]*A1dmsm*X11 + error
    dsurv2[i,1] <- mean(logTdwsurv)
    dmsm2[i,1] <- mean(logTdmsm)
  }
  if(n==500){
    fplot2 <- as.data.frame(rbind(dmsm2, dsurv2))
  }else{
    fplot2 <- as.data.frame(rbind(fplot2,dmsm2, dsurv2))
  }
}


#### ---- plots ---- ####
# X11 <- runif(largen, 0.5, 1.5) 
# X12 <- rbinom(largen, 1, 0.6)
# A1true <- ifelse(X11 > (-theta1[4]/theta1[5]), 1, 0)
# error <- 0 #rnorm(largen, sd = 0.3)
# logTtrue <- theta1[1] + theta1[2]*X11 + theta1[3]*X12 + theta1[4]*A1true + theta1[5]*A1true*X11
# true_ElogYopt <- mean(logTtrue)
# true_ElogYopt

colnames(fplot1) <- c("ElogYopt", "thetahat")
fplot1$Method <- rep(c(rep("DynamicMSM", 1000), rep("DWSurv", 1000)),3)
fplot1$n <- factor(rep(c("n=500","n=1000","n=5000"), each = 2000), levels = c("n=500","n=1000","n=5000"))
fplot1$et <- exp(fplot1$ElogYopt)

p1 <- ggplot(fplot1, aes(x = Method, y = ElogYopt)) + geom_boxplot() + scale_y_continuous(name = "Expected log-survival time") + scale_x_discrete(name = "Method") + ggtitle("Boxplot of mean log-survival time under optimal treatment") + theme_bw() + facet_grid(.~n)
p2 <- ggplot(fplot1, aes(x = Method, y = thetahat)) + geom_boxplot() + geom_hline(aes(yintercept = 0.8/0.9)) + scale_y_continuous(name = expression(hat(theta))) + scale_x_discrete(name = "Method") + ggtitle(substitute(paste("Boxplot of ",hat(theta)))) + theme_bw() + facet_grid(.~n)
#png(paste(mainpath, secpath, "dwsurv_dmsm1.png", sep = ""), width = 1000, height = 300)
png("C:/Users/Gabrielle/Google Drive/McGill - PhD/Thesis/dw-SM/JASA_updatedSERA/Response to reviewers/Extra_simulations/March 11/dwsurv_dmsm1.png", width = 800, height = 600)
ggarrange(p1,p2,nrow=2,common.legend = TRUE, labels = c("A)", "B)"))
dev.off()

colnames(fplot2) <- c("ElogYopt", "thetahat")
fplot2$Method <- rep(c(rep("DynamicMSM", 1000), rep("DWSurv", 1000)),3)
fplot2$n <- factor(rep(c("n=500","n=1000","n=5000"), each = 2000), levels = c("n=500","n=1000","n=5000"))
fplot2$et <- exp(fplot2$ElogYopt)

p1 <- ggplot(fplot2, aes(x = Method, y = ElogYopt)) + geom_boxplot() + scale_y_continuous(name = "Expected log-survival time") + scale_x_discrete(name = "Method") + ggtitle("Boxplot of mean log-survival time under optimal treatment") + theme_bw() + facet_grid(.~n)
p2 <- ggplot(fplot2, aes(x = Method, y = thetahat)) + geom_boxplot() + geom_hline(aes(yintercept = 0.15/0.2)) + scale_y_continuous(name = expression(hat(theta)), limits = c(0.5,1.5)) + scale_x_discrete(name = "Method") + ggtitle(substitute(paste("Boxplot of ",hat(theta)))) + theme_bw() + facet_grid(.~n)
#png(paste(mainpath, secpath, "dwsurv_dmsm2.png", sep = ""), width = 1000, height = 300)
png("C:/Users/Gabrielle/Google Drive/McGill - PhD/Thesis/dw-SM/JASA_updatedSERA/Response to reviewers/Extra_simulations/March 11/dwsurv_dmsm2.png", width = 800, height = 600)
ggarrange(p1,p2,nrow=2,common.legend = TRUE, labels = c("A)", "B)"))
dev.off()
save.image("C:/Users/Gabrielle/Google Drive/McGill - PhD/Thesis/dw-SM/JASA_updatedSERA/Response to reviewers/Extra_simulations/March 11/comparison_dmsm_dwsurv.RData")

