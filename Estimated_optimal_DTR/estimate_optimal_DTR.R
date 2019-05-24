library(survival)
library(DTRreg)
expit <- function(x) exp(x) / (1 + exp(x))
# mainpath: location of the folder "Reproduce_simulations_results"
mainpath <- ""

# true parameters
theta1 <- c(6.31578947, 1.5, -0.8, 0.1, 0.1)
theta2 <- c(4, 1.1052632, -0.2105263, -0.9, 0.6, -0.1052632)
lambda <- 1/300
p <- 0.9
beta <- 2

#### independent censoring, 30% censoring, scenario 1, n=500,1000,10000 ####
set.seed(1020)

for(n in c(500, 1000, 10000)){
  # store proportions
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  Huang <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("I30",n))
  for(i in 1:1000){
    #### simulate survival times with 30% independent censoring ####
    X1 <- runif(n, 0.1, 1.29)
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2)
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    delta <- rbinom(n, size = 1, prob = 0.6)
    eta2 <- rbinom(n, 1, prob = 0.8)
    delta2 <- delta[eta2 == 1]

    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    A2opt <- ifelse(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (A2opt - A2[eta2 == 1 & delta == 1])*(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3)
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)

    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)

    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)

    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), data = mydata)

    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)

    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)

    #### fit Huang and Ning ####
    logY2 <- log(Y2[eta2 == 1 & delta == 1])
    logY <- log(Y1 + Y2)
    delta_C <- abs(delta - 1)
    Y_C <- exp(logY)
    delta22 <- rep(NA, n)
    delta22[eta2 == 1] <- delta2
    logY22 <- rep(NA, n)
    logY22[eta2 == 1 & delta == 1] <- logY2
    logY11 <- rep(NA, n)
    logY11[eta2 == 1 & delta == 1] <- log(T1)

    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, delta22, delta_C, logY, logY11, logY22, Y_C, eta2, trueA1opt, trueA2opt))
    mydata <- mydata[order(mydata$Y_C), ]

    eta2 <- mydata$eta2
    X1 <- mydata$X1
    X14 <- mydata$X14
    A1 <- mydata$A1
    X2 <- mydata$X2
    X23 <- mydata$X23
    A2 <- mydata$A2
    delta <- mydata$delta
    delta2 <- mydata$delta22[eta2 == 1]
    logY <- mydata$logY
    logY1 <- mydata$logY11[eta2 == 1 & delta == 1]
    logY2 <- mydata$logY22[eta2 == 1 & delta == 1]
    trueA1opt <- mydata$trueA1opt
    trueA2opt <- mydata$trueA2opt

    km <- survfit(Surv(Y_C, delta_C) ~ 1, data = mydata)
    mydata$ipcw <- summary(km, times = mydata$Y_C)$surv

    # stage 2
    data2 <- as.data.frame(cbind(int = rep(1, sum(eta2)), X1 = X1[eta2 == 1], X2 = X2[eta2 == 1], X23 = X23[eta2 == 1], A2 = A2[eta2 == 1], A2X2 = A2[eta2 == 1]*X2[eta2 == 1]))
    H2 <- cbind(data2$int, data2$X2, data2$X23, data2$X1, data2$A2, data2$A2X2)
    w <- 1/mydata$ipcw
    H22 <- H2[delta2 == 1, ]
    w22 <- w[delta == 1 & eta2 == 1]
    psiHuang2 <- as.vector(solve(t(H22) %*% (w22*H22)) %*% t(H22) %*% (w22*logY2))[(ncol(H2)-1):ncol(H2)]

    # pseudo-outcome
    logYtilde <- logY
    A2opt <- ifelse(psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logYtilde[eta2 == 1 & delta == 1] <- log(exp(logY1) + exp(logY2 + (A2opt - A2[eta2 == 1 & delta == 1]) * (psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1])))

    # first stage estimation
    data1 <- as.data.frame(cbind(int = rep(1, n), X1, X14, A1, A1X1 = A1*X1))
    H1 <- cbind(data1$int, data1$X1, data1$X14, data1$A1, data1$A1X1)
    H11 <- H1[delta == 1, ]
    logYtilde <- logYtilde[delta == 1]
    w11 <- w[delta == 1]
    psiHuang1 <- as.vector(solve(t(H11) %*% (w11 * H11)) %*% t(H11) %*% (w11 * logYtilde))[(ncol(H1)-1):ncol(H1)]


    estA1optHuang <- ifelse(psiHuang1[1] + psiHuang1[2]*X1 > 0, 1, 0)
    estA2optHuang <- ifelse(psiHuang2[1] + psiHuang2[2]*X2 > 0, 1, 0)

    Huang[i,1] <- sum(trueA1opt == estA1optHuang)/n
    Huang[i,2] <- sum(trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1])/sum(eta2)
    Huang[i,3] <- sum((trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1optHuang[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvIndep30_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_HuangIndep30_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
  write.table(Huang, file = name2)
}

#### Dependent censoring, 30% censoring, scenario 1, n=500,1000,10000 ####
set.seed(345)

for(n in c(500, 1000, 10000)){
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  Huang <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("D30",n))
  for(i in 1:1000){
    #### simulate survival times with 30% dependent censoring (baseline) ####
    X1 <- runif(n, 0.1, 1.29)
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2)
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    delta <- rbinom(n, size = 1, prob = expit(2*X1 - 0.4))
    eta2 <- rbinom(n, 1, prob = 0.8)
    delta2 <- delta[eta2 == 1]
    
    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1 & delta == 1])*(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3)
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)
    
    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)
    
    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)
    
    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ X1, delta ~ X1), data = mydata)
    
    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)
    
    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)
    
    
    #### fit Huang and Ning ####
    logY2 <- log(Y2[eta2 == 1 & delta == 1])
    logY <- log(Y1 + Y2)
    delta_C <- abs(delta - 1)
    Y_C <- exp(logY)
    coxm <- coxph(Surv(Y_C, delta_C) ~ X1)
    bh <- basehaz(coxm)[,1]
    ipcw <- exp(-bh*exp(X1*coef(coxm)))
    
    data2 <- as.data.frame(cbind(int = rep(1, sum(eta2)), X1 = X1[eta2 == 1], X2 = X2[eta2 == 1], X23 = X23[eta2 == 1], A2 = A2[eta2 == 1], A2X2 = A2[eta2 == 1]*X2[eta2 == 1]))
    H2 <- cbind(data2$int, data2$X2, data2$X23, data2$X1, data2$A2, data2$A2X2)
    w <- 1/ipcw
    H22 <- H2[delta2 == 1, ]
    w22 <- w[delta == 1 & eta2 == 1]
    psiHuang2 <- as.vector(solve(t(H22) %*% (w22*H22)) %*% t(H22) %*% (w22*logY2))[(ncol(H2)-1):ncol(H2)]
    
    # pseudo-outcome
    logYtilde <- logY
    A2opt <- ifelse(psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logYtilde[eta2 == 1 & delta == 1] <- log(T1 + exp(logY2 + (A2opt - A2[eta2 == 1 & delta == 1]) * (psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1])))
    
    # first stage estimation
    data1 <- as.data.frame(cbind(int = rep(1, n), X1, X14, A1, A1X1 = A1*X1))
    H1 <- cbind(data1$int, data1$X1, data1$X14, data1$A1, data1$A1X1)
    H11 <- H1[delta == 1, ]
    logYtilde <- logYtilde[delta == 1]
    w11 <- w[delta == 1]
    psiHuang1 <- as.vector(solve(t(H11) %*% (w11 * H11)) %*% t(H11) %*% (w11 * logYtilde))[(ncol(H1)-1):ncol(H1)]
    
    estA1optHuang <- ifelse(psiHuang1[1] + psiHuang1[2]*X1 > 0, 1, 0)
    estA2optHuang <- ifelse(psiHuang2[1] + psiHuang2[2]*X2 > 0, 1, 0)
    
    Huang[i,1] <- sum(trueA1opt == estA1optHuang)/n
    Huang[i,2] <- sum(trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1])/sum(eta2)
    Huang[i,3] <- sum((trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1optHuang[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvDep30_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_HuangDep30_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
  write.table(Huang, file = name2)
}

#### independent censoring, 60% censoring, scenario 1, n=500,1000,10000 ####
set.seed(4409)

for(n in c(500, 1000, 10000)){
  # store proportions
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  Huang <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("I60",n))
  for(i in 1:1000){
    #### simulate survival times with 60% independent censoring ####
    X1 <- runif(n, 0.1, 1.29)
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2)
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    delta <- rbinom(n, size = 1, prob = 0.4)
    eta2 <- rbinom(n, 1, prob = 0.8)
    delta2 <- delta[eta2 == 1]

    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    A2opt <- ifelse(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (A2opt - A2[eta2 == 1 & delta == 1])*(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3)
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)

    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)

    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)

    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), data = mydata)

    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)

    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)

    #### fit Huang and Ning ####
    logY2 <- log(Y2[eta2 == 1 & delta == 1])
    logY <- log(Y1 + Y2)
    delta_C <- abs(delta - 1)
    Y_C <- exp(logY)
    delta22 <- rep(NA, n)
    delta22[eta2 == 1] <- delta2
    logY22 <- rep(NA, n)
    logY22[eta2 == 1 & delta == 1] <- logY2
    logY11 <- rep(NA, n)
    logY11[eta2 == 1 & delta == 1] <- log(T1)

    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, delta22, delta_C, logY, logY11, logY22, Y_C, eta2, trueA1opt, trueA2opt))
    mydata <- mydata[order(mydata$Y_C), ]

    eta2 <- mydata$eta2
    X1 <- mydata$X1
    X14 <- mydata$X14
    A1 <- mydata$A1
    X2 <- mydata$X2
    X23 <- mydata$X23
    A2 <- mydata$A2
    delta <- mydata$delta
    delta2 <- mydata$delta22[eta2 == 1]
    logY <- mydata$logY
    logY1 <- mydata$logY11[eta2 == 1 & delta == 1]
    logY2 <- mydata$logY22[eta2 == 1 & delta == 1]
    trueA1opt <- mydata$trueA1opt
    trueA2opt <- mydata$trueA2opt

    km <- survfit(Surv(Y_C, delta_C) ~ 1, data = mydata)
    mydata$ipcw <- summary(km, times = mydata$Y_C)$surv

    # stage 2
    data2 <- as.data.frame(cbind(int = rep(1, sum(eta2)), X1 = X1[eta2 == 1], X2 = X2[eta2 == 1], X23 = X23[eta2 == 1], A2 = A2[eta2 == 1], A2X2 = A2[eta2 == 1]*X2[eta2 == 1]))
    H2 <- cbind(data2$int, data2$X2, data2$X23, data2$X1, data2$A2, data2$A2X2)
    w <- 1/mydata$ipcw
    H22 <- H2[delta2 == 1, ]
    w22 <- w[delta == 1 & eta2 == 1]
    psiHuang2 <- as.vector(solve(t(H22) %*% (w22*H22)) %*% t(H22) %*% (w22*logY2))[(ncol(H2)-1):ncol(H2)]

    # pseudo-outcome
    logYtilde <- logY
    A2opt <- ifelse(psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logYtilde[eta2 == 1 & delta == 1] <- log(exp(logY1) + exp(logY2 + (A2opt - A2[eta2 == 1 & delta == 1]) * (psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1])))

    # first stage estimation
    data1 <- as.data.frame(cbind(int = rep(1, n), X1, X14, A1, A1X1 = A1*X1))
    H1 <- cbind(data1$int, data1$X1, data1$X14, data1$A1, data1$A1X1)
    H11 <- H1[delta == 1, ]
    logYtilde <- logYtilde[delta == 1]
    w11 <- w[delta == 1]
    psiHuang1 <- as.vector(solve(t(H11) %*% (w11 * H11)) %*% t(H11) %*% (w11 * logYtilde))[(ncol(H1)-1):ncol(H1)]


    estA1optHuang <- ifelse(psiHuang1[1] + psiHuang1[2]*X1 > 0, 1, 0)
    estA2optHuang <- ifelse(psiHuang2[1] + psiHuang2[2]*X2 > 0, 1, 0)

    Huang[i,1] <- sum(trueA1opt == estA1optHuang)/n
    Huang[i,2] <- sum(trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1])/sum(eta2)
    Huang[i,3] <- sum((trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1optHuang[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvIndep60_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_HuangIndep60_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
  write.table(Huang, file = name2)
}

#### Dependent censoring, 60% censoring, scenario 1, n=500,1000,10000 ####
set.seed(20)

for(n in c(500, 1000, 10000)){
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  Huang <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("D60", n))
  for(i in 1:1000){
    #### simulate survival times with 60% dependent censoring (baseline) ####
    X1 <- runif(n, 0.1, 1.29)
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2)
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    delta <- rbinom(n, size = 1, prob = expit(2*X1 - 1.8))
    eta2 <- rbinom(n, 1, prob = 0.8)
    delta2 <- delta[eta2 == 1]

    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1 & delta == 1])*(theta2[4] + theta2[5]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3)
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)

    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)

    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)

    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ X1, delta ~ X1), data = mydata)

    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)

    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)


    #### fit Huang and Ning ####
    logY2 <- log(Y2[eta2 == 1 & delta == 1])
    logY <- log(Y1 + Y2)
    delta_C <- abs(delta - 1)
    Y_C <- exp(logY)
    coxm <- coxph(Surv(Y_C, delta_C) ~ X1)
    bh <- basehaz(coxm)[,1]
    ipcw <- exp(-bh*exp(X1*coef(coxm)))

    data2 <- as.data.frame(cbind(int = rep(1, sum(eta2)), X1 = X1[eta2 == 1], X2 = X2[eta2 == 1], X23 = X23[eta2 == 1], A2 = A2[eta2 == 1], A2X2 = A2[eta2 == 1]*X2[eta2 == 1]))
    H2 <- cbind(data2$int, data2$X2, data2$X23, data2$X1, data2$A2, data2$A2X2)
    w <- 1/ipcw
    H22 <- H2[delta2 == 1, ]
    w22 <- w[delta == 1 & eta2 == 1]
    psiHuang2 <- as.vector(solve(t(H22) %*% (w22*H22)) %*% t(H22) %*% (w22*logY2))[(ncol(H2)-1):ncol(H2)]

    # pseudo-outcome
    logYtilde <- logY
    A2opt <- ifelse(psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logYtilde[eta2 == 1 & delta == 1] <- log(T1 + exp(logY2 + (A2opt - A2[eta2 == 1 & delta == 1]) * (psiHuang2[1] + psiHuang2[2]*X2[eta2 == 1 & delta == 1])))

    # first stage estimation
    data1 <- as.data.frame(cbind(int = rep(1, n), X1, X14, A1, A1X1 = A1*X1))
    H1 <- cbind(data1$int, data1$X1, data1$X14, data1$A1, data1$A1X1)
    H11 <- H1[delta == 1, ]
    logYtilde <- logYtilde[delta == 1]
    w11 <- w[delta == 1]
    psiHuang1 <- as.vector(solve(t(H11) %*% (w11 * H11)) %*% t(H11) %*% (w11 * logYtilde))[(ncol(H1)-1):ncol(H1)]

    estA1optHuang <- ifelse(psiHuang1[1] + psiHuang1[2]*X1 > 0, 1, 0)
    estA2optHuang <- ifelse(psiHuang2[1] + psiHuang2[2]*X2 > 0, 1, 0)

    Huang[i,1] <- sum(trueA1opt == estA1optHuang)/n
    Huang[i,2] <- sum(trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1])/sum(eta2)
    Huang[i,3] <- sum((trueA2opt[eta2 == 1] == estA2optHuang[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1optHuang[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvDep60_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_HuangDep60_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
  write.table(Huang, file = name2)
}


#### Censoring dependent on time-varying covariates (dWSurv only), 30% censoring, scenarios 1-4, n=500,1000,10000 ####
set.seed(994)

for(n in c(500, 1000, 10000)){
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("D30", n))
  for(i in 1:1000){
    #### simulate survival times with 30% dependent censoring (time-varying) ####
    X1 <- runif(n, 0.1, 1.29) 
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2) 
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    
    pC1 <- expit(2*X1 - 0.2)/(1 - 0.0805*0.8)
    C1 <- rbinom(n, 1, pC1)
    eta2 <- rep(0, n)
    eta2[C1 == 1] <- rbinom(sum(C1), 1, 0.8)
    delta2 <- rbinom(sum(eta2), 1, 1 - expit(0.3 - 2*X2[eta2 == 1]))
    delta <- rep(1, n)
    delta[C1 == 0] <- 0
    delta[eta2 == 1] <- delta2
    
    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    trueA2opt <- ifelse(theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1 & delta == 1])*(theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3) 
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)
    
    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)
    
    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)
    
    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ X1, delta ~ X1), data = mydata)
    
    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)
    
    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvX1X2Dep30_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
}

#### Censoring dependent on time-varying covariates (dWSurv only), 60% censoring, scenarios 1-4, n=500,1000,10000 ####
set.seed(276)

for(n in c(500, 1000, 10000)){
  dWSurv <- matrix(NA, ncol = 3, nrow = 1000)
  print(c("D60", n))
  for(i in 1:1000){
    #### simulate survival times with 60% dependent censoring (time-varying) ####
    X1 <- runif(n, 0.1, 1.29) 
    X14 <- X1^4
    A1 <- rbinom(n, size = 1, prob = expit(2*X1 - 1))
    X2 <- runif(n, 0.9, 2) 
    X23 <- X2^3
    A2 <- rbinom(n, size = 1, prob = expit(-2*X2 + 2.8))
    
    pC1 <- expit(2*X1 - 1.8)/(1 - 0.1466*0.8)
    C1 <- rbinom(n, 1, pC1)
    eta2 <- rep(0, n)
    eta2[C1 == 1] <- rbinom(sum(C1), 1, 0.8)
    delta2 <- rbinom(sum(eta2), 1, 1 - expit(1 - 2*X2[eta2 == 1]))
    delta <- rep(1, n)
    delta[C1 == 0] <- 0
    delta[eta2 == 1] <- delta2
    
    logT2 <- theta2[1] + theta2[2]*X2[eta2 == 1 & delta == 1] + theta2[3]*X23[eta2 == 1 & delta == 1] + theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] + theta2[6]*X1[eta2 == 1 & delta == 1] + rnorm(sum(eta2*delta), sd = 0.3)
    trueA2opt <- ifelse(theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1] > 0, 1, 0)
    logT2opt <- logT2 + (trueA2opt - A2[eta2 == 1 & delta == 1])*(theta2[4]*A2[eta2 == 1 & delta == 1] + theta2[5]*A2[eta2 == 1 & delta == 1]*X2[eta2 == 1 & delta == 1])
    logT <- theta1[1] + theta1[2]*X1 + theta1[3]*X14 + theta1[4]*A1 + theta1[5]*A1*X1 + rnorm(n,sd = 0.3) 
    T1 <- exp(logT[eta2 == 1 & delta == 1]) - exp(logT2opt)
    if(min(T1) <= 0){
      next
    }
    C <- rexp(n - sum(delta), rate = 1/300)
    
    Y2 <- rep(NA, n)
    Y1 <- rep(NA, n)
    eta2d0 <- eta2[delta == 0]
    C1 <- rep(NA, length(C))
    C2 <- rep(NA, length(C))
    for(j in 1:length(C))
    {
      if(eta2d0[j] == 0){
        C1[j] <-  C[j]
        C2[j] <- 0
      }else{
        C1[j] <- runif(1, 0, C[j])
        C2[j] <- C[j] - C1[j]
      }
    }
    Y2[delta == 0] <- C2
    Y1[delta == 0] <- C1
    Y1[delta == 1 & eta2 == 1] <- T1
    Y1[delta == 1 & eta2 == 0] <- exp(logT[delta == 1 & eta2 == 0])
    Y2[delta == 1 & eta2 == 0] <- 0
    Y2[delta == 1 & eta2 == 1] <- exp(logT2)
    
    trueA1opt <- ifelse(theta1[4] + theta1[5]*X1 > 0, 1, 0)
    trueA2opt <- ifelse(theta2[4] + theta2[5]*X2 > 0, 1, 0)
    
    #### fit dWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- dWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ X1, delta ~ X1), data = mydata)
    
    estA1opt <- ifelse(mod$psi[[1]][1] + mod$psi[[1]][2]*X1 > 0, 1, 0)
    estA2opt <- ifelse(mod$psi[[2]][1] + mod$psi[[2]][2]*X2 > 0, 1, 0)
    
    dWSurv[i,1] <- sum(trueA1opt == estA1opt)/n
    dWSurv[i,2] <- sum(trueA2opt[eta2 == 1] == estA2opt[eta2 == 1])/sum(eta2)
    dWSurv[i,3] <- sum((trueA2opt[eta2 == 1] == estA2opt[eta2 == 1]) & (trueA1opt[eta2 == 1] == estA1opt[eta2 == 1]))/sum(eta2)
  }
  name1 <- paste(mainpath, "Estimated_optimal_DTR/Results/1_dWSurvX1X2Dep60_", n, ".txt", sep = "")
  write.table(dWSurv, file = name1)
}
