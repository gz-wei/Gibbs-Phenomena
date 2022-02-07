library(expm)
library(MASS)

rm(list = ls()) 
### KF for the original low-pass data stream 
# load the simulated data stream 
load(here::here("data", "N=100", "y.tilde.lp.RData"))
# compute the transition matrix G.A and F.tilde.A
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
K = 10 
Omega <- Function_Omega(K)
v <- c(0.02, 0)
G <- G_nit(Omega,v)
step <- 1
G.step <- expm(step*G)
G.A <- rbind(
  cbind(G.step, diag(K^2)),
  cbind(matrix(0,K^2,K^2), diag(K^2))
)
F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))
### KF with the specified V and W 
source(here::here("functions", "KF_Non_Missing.R"))
V <- diag(0.01, length(y.tilde.lp[[1]]))
W <- diag(0.01, 2*K^2)
m0 = rnorm(2*K^2, 0, 1)
C0 = diag(0.1, 2*K^2)
start_time <- Sys.time()
fit.KF <- KF_Non_Missing(y.tilde.lp, F.tilde.A, G.A, V, W, m0, C0)
end_time <- Sys.time()
end_time - start_time
### plot the filtered data stream
Nr = 100
F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
if (file.exists(here::here("data", "simulation", F.name))){
  load(here::here("data", "simulation", F.name))
} else{
  F <- Function_F(Nr, K, Omega) 
  save(F, file = here::here("data", "simulation", F.name))
}
F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
N.step <- length(y.tilde.lp)
for (i in 1:N.step){
  tempt <- matrix(F.A%*%fit.KF$m.flt[[i+1]], Nr, Nr)
  image(tempt)
  Sys.sleep(0.3)
}
### test the physical meaning 
image(matrix(F.A%*%fit.KF$m.flt[[10]], Nr, Nr))
image(matrix(F.A%*%G.A%*%fit.KF$m.flt[[10]], Nr, Nr))
image(matrix(F.A%*%fit.KF$m.flt[[11]], Nr, Nr))
### plot the source term 
for (i in 1:N.step){
  beta <- fit.KF$m.flt[[i+1]][(K^2+1):(2*K^2)]
  tempt <- matrix(F%*%beta, Nr, Nr)
  image(tempt)
  Sys.sleep(0.3)
}


### KF for the flipped low-pass data stream 
# load the simulated data stream 
load(here::here("data", "N=100", "y.tilde.flp.RData"))
# compite G*
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
K.f = 2*K
Omega.f <- Function_Omega(K.f)
F.f.name = paste(c("F.f", "Nr", as.character(2*Nr), "K.f", as.character(K.f), "RData"), collapse = ".")
if(file.exists(here::here("data", "simulation", F.f.name))){
  load(here::here("data", "simulation", F.f.name))
}else{
  F.f <- Function_F(2*Nr, K.f, Omega.f) 
  save(F.f, file = here::here("data", "simulation", F.f.name))
}
dim <- Nr
I <- diag(1, dim)
J <- I[dim:1,]
P <- rbind(I,J)
if(file.exists(here::here("data", "simulation", "H.RData"))){
  load(here::here("data", "simulation", "H.RData"))
}else{
  H <- ginv(F.f)%*%kronecker(P,P)%*%F
  save(H, file = here::here("data", "simulation", "H.RData"))
}
G.f <- H%*%G%*%ginv(H)
G.step.f <- expm(step*G.f)
G.f.A <- rbind(
  cbind(G.step.f, diag(K.f^2)),
  cbind(matrix(0,K.f^2, K.f^2), diag(K.f^2))
)
# compute F*
F.f.tilde.A <- cbind(diag(K.f^2), matrix(0, K.f^2, K.f^2))
### KF with the specified V and W 
V.f <- diag(0.01, length(y.tilde.flp[[1]]))
dim(V.f)
W.f <- diag(0.01, 2*K.f^2)
m0.f = rnorm(2*K.f^2, 0, 1)
C0.f = diag(0.1, 2*K.f^2)
if(file.exists(here::here("data", "simulation", "Fit.KF.f.RData"))){
  load(here::here("data", "simulation", "Fit.KF.f.RData"))
} else {
  start_time <- Sys.time()
  fit.KF.f <- KF_Non_Missing(y.tilde.flp, F.f.tilde.A, G.f.A, V.f, W.f, m0.f, C0.f)
  end_time <- Sys.time()
  print(end_time - start_time)
  save(fit.KF.f, file = here::here("data", "simulation", "Fit.KF.f.RData"))
}

### plot the filtered data stream
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
N.step <- length(y.tilde.flp)
for (i in 1:N.step){
  tempt <- matrix(F.f.A%*%fit.KF.f$m.flt[[i+1]], 2*Nr, 2*Nr)[1:Nr, 1:Nr]
  image(tempt)
  Sys.sleep(0.3)
}
### test the physical meaning 
image(matrix(F.f.A%*%fit.KF.f$m.flt[[10]], 2*Nr, 2*Nr))
image(matrix(F.f.A%*%G.f.A%*%fit.KF.f$m.flt[[10]], 2*Nr, 2*Nr))
image(matrix(F.f.A%*%fit.KF.f$m.flt[[11]], 2*Nr, 2*Nr))
### plot the source term 
for (i in 1:N.step){
  beta.f <- fit.KF.f$m.flt[[i+1]][(K.f^2+1):(2*K.f^2)]
  tempt <- matrix(F.f%*%beta.f, 2*Nr, 2*Nr)[1:Nr, 1:Nr]
  image(tempt)
  Sys.sleep(0.3)
}


