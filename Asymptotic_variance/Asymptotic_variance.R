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


#### independent censoring, 30% censoring, n=500,1000,10000 ####
set.seed(1020)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
  print(c("I30", n))
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_indep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_indep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_indep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_indep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}

#### Dependent censoring, 30% censoring, n=500,1000,10000 ####
set.seed(345)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
  print(c("D30", n))
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
    C <- (-log(runif(n - sum(delta), 0, 1))/(lambda*exp(beta*X1[delta == 0])))^(1/p)
    
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(~Y1, ~Y2), status = delta, blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list( ~ 1, ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_dep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_dep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_dep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_dep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}

#### independent censoring, 60% censoring, n=500,1000,10000 ####
set.seed(4409)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
  print(c("I60", n))
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(Y1, Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_indep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_indep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_indep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_indep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}

#### Dependent censoring, 60% censoring, n=500,1000,10000 ####
set.seed(20)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
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
    C <- (-log(runif(n - sum(delta), 0, 1))/(lambda*exp(beta*X1[delta == 0])))^(1/p)
    
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_dep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_dep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_dep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_dep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}

#### Censoring dependent on time-varying covariates, 30% censoring, n=500,1000,10000 ####
set.seed(994)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
  print(c("X1X230", n))
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_X1X2dep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_X1X2dep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_X1X2dep30_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_X1X2dep30_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}

#### Censoring dependent on time-varying covariates, 60% censoring, n=500,1000,10000 ####
set.seed(276)

for(n in c(500, 1000, 10000)){
  # store estimates
  est_psidWSurv1 <- est_psidWSurv2 <- se_psidWSurv1 <- se_psidWSurv2 <- matrix(NA, ncol = 2, nrow = 1000)
  print(c("X1X260", n))
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
    
    #### fit DWSurv ####
    mydata <- as.data.frame(cbind(X1, X14, A1, X2, X23, A2, delta, Y1, Y2))
    mod <- DWSurv(time = list(~Y1, ~Y2), blip.mod = list(~ X1, ~ X2), treat.mod = list(A1 ~ X1, A2 ~ X2), tf.mod = list( ~ X1 + X14, ~ X2 + X23 + X1), cens.mod = list(delta ~ 1, delta ~ 1), var.estim = "asymptotic", data = mydata)
    
    est_psidWSurv1[i,] <- mod$psi[[1]]
    est_psidWSurv2[i,] <- mod$psi[[2]]
    se_psidWSurv1[i,] <- c(sqrt(mod$covmat[[1]][4,4]), sqrt(mod$covmat[[1]][5,5]))
    se_psidWSurv2[i,] <- c(sqrt(mod$covmat[[2]][5,5]), sqrt(mod$covmat[[2]][6,6]))
  }
  # export results
  name1 <- paste(mainpath, "Asymptotic_variance/Results/est_X1X2dep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/est_X1X2dep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(est_psidWSurv1, file = name1)
  write.table(est_psidWSurv2, file = name2)
  
  name1 <- paste(mainpath, "Asymptotic_variance/Results/se_X1X2dep60_psidWSurv1_", n, ".txt", sep = "")
  name2 <- paste(mainpath, "Asymptotic_variance/Results/se_X1X2dep60_psidWSurv2_", n, ".txt", sep = "")
  write.table(se_psidWSurv1, file = name1)
  write.table(se_psidWSurv2, file = name2)
}