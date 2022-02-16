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


l <- function(x){
  v <- x[1] 
  w1 <- x[2] 
  w2 <- x[3]
  T <- length(y.tilde.lp)
  V <- diag(v, length(y.tilde.lp[[1]]))
  K = 4
  W <- diag(c(rep(w1,K^2), rep(w2,K^2)))
  m0 = rnorm(2*K^2, 0, 1)
  C0 = diag(0.1, 2*K^2)
  F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))
  fit.KF <- KF_Non_Missing(y.tilde.lp, F.tilde.A, G.A, V, W, m0, C0)
  sum.t <- numeric(T)
  for (t in 1:T){
    tempt1 <- fit.KF$C.prd[[t]][1:K^2, 1:K^2] + V
    tempt2 <- y.tilde.lp[[t]] - fit.KF$m.prd[[t]][1:K^2]
    sum.t[t] <- -1/2*(log(det(tempt1)) + t(tempt2)%*%solve(tempt1)%*%tempt2)
  }
  return(-sum(sum.t))
}

l(c(0.005^2, 0.005^2, 0.009^2))
optim(c(0.005^2,0.005^2, 0.001^2), l, method = "BFGS")
