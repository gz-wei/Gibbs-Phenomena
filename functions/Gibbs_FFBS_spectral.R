###Name: Spectral Gibbs Forward Filtering Backward Sampling for dynamical model estimation
###Author: Guanzhou Wei
###Last revision: 02/08/2022
###Required packages: MCMCpack, Rfast
###Required functions: KF_Non_Missing.R


library(expm)
library(MCMCpack)
library(Rfast)
library(rasterImage)

rm(list = ls())

source(here::here("functions", "KF_Non_Missing.R"))
load(here::here("data", "N=100", "y.tilde.lp.RData"))
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "G_nit.R"))
K = 4
Omega <- Function_Omega(K)
v <- c(0.01, 0)
G <- G_nit(Omega,v)
step <- 1
G.step <- expm(step*G)
G.A <- rbind(
  cbind(G.step, diag(K^2)),
  cbind(matrix(0,K^2,K^2), diag(K^2))
)
F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))
Nr = 100
F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
if (file.exists(here::here("data", "simulation", F.name))){
  load(here::here("data", "simulation", F.name))
} else{
  F <- Function_F(Nr, K, Omega)
  save(F, file = here::here("data", "simulation", F.name))
}
F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))


Gibbs_FFBS_spectral <- function(N.sample, y.tilde.lp, G.A, F.tilde.A){
  T <- length(y.tilde.lp)
  m0 = rnorm(2*K^2, 0, 1)
  C0 = diag(0.1, 2*K^2)
  
  ### IG prior 
  sigma2.v <- numeric(N.sample)
  alpha.v =  1
  beta.v = 0.1
  # sigma2.v[1] <- rinvgamma(1, alpha.v, beta.v) ### set the initial value
  sigma2.v[1] <- 3e-05

  sigma2.w1 <- numeric(N.sample)
  alpha.w1 = 1
  beta.w1 = 0.1
  ### set the initial value 
  #sigma2.w1[1] <- rinvgamma(1, alpha.w1, beta.w1) ### set the initial value 
  sigma2.w1[1] <- 3e-05
  sigma2.w2 <- numeric(N.sample)
  alpha.w2 = 1
  beta.w2 = 0.1
  #sigma2.w2[1] <- rinvgamma(1, alpha.w2, beta.w2) ### set the initial value 
  sigma2.w1[1] <- 2e-06
  
  V <- diag(sigma2.v[1], nrow(F.tilde.A))
  W <- diag(c(rep(sigma2.w1[1], 1/2*nrow(C0)), rep(sigma2.w2[1], 1/2*nrow(C0))))
  
  ### IW prior 
  # W1 <- list()
  # nu.w1 <- 1/2*nrow(C0)
  # Phi.w1 <- diag(0.01, 1/2*nrow(C0)) 
  # W1[[1]] <- rwish(nu.w1, Phi.w1)
  # W2 <- list()
  # nu.w2 <- 1/2*nrow(C0)
  # Phi.w2 <- diag(0.01, 1/2*nrow(C0)) 
  # W2[[1]] <- rwish(nu.w2, Phi.w2)
  # W <- list()
  # W[[1]] <- rbind(
  #   cbind(W1[[1]], matrix(0, nrow(W1[[1]]), ncol(W1[[1]]))),
  #   cbind(matrix(0, nrow(W2[[1]]), ncol(W2[[1]])), W2[[1]])
  # )
  
  #k = 2
  for(k in 2:N.sample) {
    ## Forward Filtering
    fit.KF <- KF_Non_Missing(y.tilde.lp, F.tilde.A, G.A, V, W, m0, C0)
    ## plot the FF performance
    # for (t in 1:T){
    #   image(matrix(F.A%*%fit.KF$m.flt[[t]], Nr, Nr))
    #   Sys.sleep(0.3)
    # }
    
    ## Backward Sampling
    theta <- matrix(0, nrow(G.A), T)
    theta[, T] <- t(rmvnorm(1, fit.KF$m.flt[[T+1]], fit.KF$C.flt[[T+1]]))
    for (t in (T-1):1){
      ht <- fit.KF$m.flt[[t+1]] + fit.KF$C.flt[[t+1]]%*%t(G.A)%*%solve(fit.KF$C.prd[[t+1]])%*%(theta[,t+1] - fit.KF$m.prd[[t+1]])
      Ht <- fit.KF$C.flt[[t+1]] - fit.KF$C.flt[[t+1]]%*%t(G.A)%*%solve(fit.KF$C.prd[[t+1]])%*%G.A%*%fit.KF$C.flt[[t+1]]
      #theta[, t] <- t(rmvnorm(1, ht, Ht))
      theta[, t] <- t(rmvnorm(1, ht, as.matrix(nearPD(Ht)$mat)))
    }
    ## plot the BS performance
    # for (t in 1:T){
    #   image(matrix(F.A%*%theta[, t], Nr, Nr))
    #   Sys.sleep(0.3)
    # }
    
    ## update parameters
    # residual calculation
    y.rsd <- matrix(0, nrow(F.tilde.A), T)
    state.rsd <- matrix(0, nrow(G.A), T-1)
    for (t in 1:T){
      y.rsd[,t] <- y.tilde.lp[[t]] - c(F.tilde.A%*%theta[, t])
    }
    for (t in 1:(T-1)){
      state.rsd[,t] <- theta[,t+1] - G.A%*%theta[,t]
    }
    #### closed-form posterior updating 
    ## update v
    alpha.v <- alpha.v + nrow(F.tilde.A)*T/2
    beta.v <- beta.v + sum(y.rsd^2)/2
    sigma2.v[k] <- rinvgamma(1, alpha.v, beta.v)
    ## update w1
    alpha.w1 <- alpha.w1 + 1/2*nrow(G.A)*(T-1)/2
    row.w1 <- 1:(1/2*nrow(G.A))
    beta.w1 <- beta.w1 + sum(state.rsd[row.w1,]^2)/2
    sigma2.w1[k] <- rinvgamma(1, alpha.w1, beta.w1)
    ## update w2
    alpha.w2 <- alpha.w2 + 1/2*nrow(G.A)*(T-1)/2
    row.w2 <- (1/2*nrow(G.A)+1):nrow(G.A)
    beta.w2 <- beta.w2 + sum(state.rsd[row.w2,]^2)/2
    sigma2.w2[k] <- rinvgamma(1, alpha.w2, beta.w2)

    V <- diag(sigma2.v[k], nrow(F.tilde.A))
    W <- diag(c(rep(sigma2.w1[k], 1/2*nrow(C0)), rep(sigma2.w2[k], 1/2*nrow(C0))))
    ## IW prior update 
    # nu.w1 <- nu.w1 + T-1
    # row.w1 <- 1:(1/2*nrow(G.A))
    # Phi.w1 <- Phi.w1 + state.rsd[row.w1,]%*%t(state.rsd[row.w1,])
    # W1[[k]] <- riwish(nu.w1, Phi.w1)
    # nu.w2 <- nu.w2 + T-1
    # row.w2 <- (1/2*nrow(G.A)+1):nrow(G.A)
    # Phi.w2 <- Phi.w2 + state.rsd[row.w2,]%*%t(state.rsd[row.w2,])
    # W2[[k]] <- riwish(nu.w2, Phi.w2)
    # W[[k]] <- rbind(
    #   cbind(W1[[k]], matrix(0, nrow(W1[[k]]), ncol(W1[[k]]))),
    #   cbind(matrix(0, nrow(W2[[k]]), ncol(W2[[k]])), W2[[k]])
    # )
  }
  return(list(sigma2.w1=sigma2.w1, sigma2.w2=sigma2.w2, sigma2.v=sigma2.v))
}

fit.test = Gibbs_FFBS_spectral(30000, y.tilde.lp, G.A, F.tilde.A)
save(fit.test, file = here::here("data", "simulation", "MCMC.RData"))
load(here::here("data", "simulation", "MCMC.RData"))
plot(fit.test$sigma2.v)
plot(fit.test$sigma2.w1)
plot(fit.test$sigma2.w2)
plot(fit.test$sigma2.v[25000:30000])
plot(fit.test$sigma2.w1[25000:30000])
plot(fit.test$sigma2.w2[25000:30000])

fit.test$sigma2.v[30000]
fit.test$sigma2.w1[30000]
fit.test$sigma2.w2[30000]