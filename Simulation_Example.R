library(expm)
library(MASS)
library(matrixcalc)

rm(list = ls()) 
#***********************************************************************************
### KF for the original low-pass data stream 
#***********************************************************************************
# load the simulated data stream 
K = 10
y.tilde.lp.name = paste(c("y.tilde.lp", as.character(K), "RData"), collapse = ".")
load(here::here("data", "data_simulation", y.tilde.lp.name))
# compute the transition matrix G.A and F.tilde.A
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
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
### KF with the specified V and W 
source(here::here("functions", "KF_Non_Missing.R"))
V <- diag(0.005, nrow(F.tilde.A))
W1 <- diag(0.005, K^2)
W2 <- diag(0.001, K^2)
W <- rbind(
  cbind(W1, matrix(0, nrow(W1), ncol(W1))),
  cbind(matrix(0, nrow(W2), ncol(W2)), W2)
)
m0 = rnorm(2*K^2, 0, 1)
C0 = diag(0.1, 2*K^2)
fit0.name = paste(c("fit0", "K", as.character(K), "RData"), collapse = ".")
if(file.exists(here::here("data", "simulation_example", fit0.name))){
  load(here::here("data", "simulation_example", fit0.name))
} else{
  start_time <- Sys.time()
  Fit0.KF <- KF_Non_Missing(y.tilde.lp, F.tilde.A, G.A, V, W, m0, C0)
  save(Fit0.KF, file = here::here("data", "simulation_example", fit0.name))
  end_time <- Sys.time()
  end_time - start_time
}
### plot the filtered data stream
Nr = 100
F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
if (file.exists(here::here("data", "data_simulation", F.name))){
  load(here::here("data", "data_simulation", F.name))
} else{
  F <- Function_F(Nr, K, Omega) 
  save(F, file = here::here("data", "data_simulation", F.name))
}
F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
N.step <- length(y.tilde.lp)
for (i in 1:N.step){
  tempt <- matrix(F.A%*%Fit0.KF$m.flt[[i+1]], Nr, Nr)
  image(tempt)
  Sys.sleep(0.3)
}

#***********************************************************************************
### KF for the flipped low-pass data stream 
#***********************************************************************************
# load the simulated data stream 
K.f = 2*K
y.tilde.flp.lp.name = paste(c("y.tilde.flp.lp", as.character(K.f), "RData"), collapse = ".")
load(here::here("data", "data_simulation", y.tilde.flp.lp.name))
# compute G*
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
Omega.f <- Function_Omega(K.f)
F.f.name = paste(c("F.f", "Nr", as.character(2*Nr), "K.f", as.character(K.f), "RData"), collapse = ".")
if(file.exists(here::here("data", "data_simulation", F.f.name))){
  load(here::here("data", "data_simulation", F.f.name))
}else{
  F.f <- Function_F(2*Nr, K.f, Omega.f) 
  save(F.f, file = here::here("data", "data_simulation", F.f.name))
}
dim <- Nr
I <- diag(1, dim)
J <- I[dim:1,]
P <- rbind(I,J)
if(file.exists(here::here("data", "simulation_example", "H.RData"))){
  load(here::here("data", "simulation_example", "H.RData"))
}else{
  H <- ginv(F.f)%*%kronecker(P,P)%*%F
  save(H, file = here::here("data", "simulation_example", "H.RData"))
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
V.f <- H%*%V%*%t(H)
W1.f <- H%*%W1%*%t(H)
W2.f <- H%*%W2%*%t(H)
W.f <- rbind(
  cbind(W1.f, matrix(0, nrow(W1.f), ncol(W1.f))),
  cbind(matrix(0, nrow(W2.f), ncol(W2.f)), W2.f)
)
W.f <- as.matrix(nearPD(W.f)$mat)
m0.f = rnorm(2*K.f^2, 0, 1)
C0.f = diag(0.1, 2*K.f^2)
if(file.exists(here::here("data", "simulation_example", "Fit.KF.f.RData"))){
  load(here::here("data", "simulation_example", "Fit.KF.f.RData"))
} else {
  start_time <- Sys.time()
  fit.KF.f <- KF_Non_Missing(y.tilde.flp, F.f.tilde.A, G.f.A, V.f, W.f, m0.f, C0.f)
  end_time <- Sys.time()
  print(end_time - start_time)
  save(fit.KF.f, file = here::here("data", "simulation_example", "Fit.KF.f.RData"))
}

### plot the filtered data stream
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
N.step <- length(y.tilde.flp)
for (i in 1:N.step){
  tempt <- matrix(F.f.A%*%fit.KF.f$m.flt[[i+1]], 2*Nr, 2*Nr)[1:Nr, 1:Nr]
  image(tempt)
  Sys.sleep(0.3)
}
