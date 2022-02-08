library(expm)
library(MCMCpack)
library(Rfast)

rm(list = ls())

source(here::here("functions", "KF_Non_Missing.R"))
load(here::here("data", "N=100", "y.tilde.lp.RData"))
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "G_nit.R"))
K = 10 
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
N.sample = 20

Gibbs_FFBS_spectral <- function(N.sample, y.tilde.lp, G.A, F.tilde.A){
  T <- length(y.tilde.lp)
  sigma2.v <- numeric(N.sample)
  sigma2.w <- numeric(N.sample)
  alpha.v =  1
  beta.v = 0.1
  alpha.w = 1
  beta.w = 0.1
  sigma2.v[1] <- rinvgamma(1, alpha.v, beta.v)
  #sigma2.w[1] <- rinvgamma(1, alpha.w, beta.w)
  sigma2.w[1] <- 2
  m0 = rnorm(2*K^2, 0, 1)
  C0 = diag(0.1, 2*K^2)
  V <- diag(sigma2.v[1], nrow(F.tilde.A))
  W <- diag(sigma2.w[1], 2*K^2)
  
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
      theta[, t] <- t(rmvnorm(1, ht, Ht))
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
    # closed-form posterior updating 
    alpha.v <- alpha.v + nrow(F.tilde.A)*T/2
    beta.v <- beta.v + sum(y.rsd^2)/2
    #alpha.w <- alpha.w + nrow(F.tilde.A)*(T-1)/2
    #beta.w <- beta.w + sum(state.rsd^2)/2
    sigma2.v[k] <- rinvgamma(1, alpha.v, beta.v)
    #sigma2.w[k] <- rinvgamma(1, alpha.w, beta.w)
    V <- diag(sigma2.v[k], nrow(F.tilde.A))
    #W <- diag(sigma2.w[k], 2*K^2)
  }
}

plot(sigma2.v)
plot(sigma2.w)