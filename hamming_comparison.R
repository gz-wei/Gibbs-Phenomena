

rm(list = ls()) 

source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
source(here::here("functions", "KF_Non_Missing.R"))
source(here::here('functions', 'Hamming_window.R'))
load(here::here("data", "data_simulation", "y.sim.RData"))

Nr = 100 
h = Func_hamming(Nr, Nr)

### run existing approach with different K and save the fitted model
for (K in c(10, 14, 20)) {
  ## obtain y.tilde with the Hamming window 
  y.hw.tilde.lp <- list()
  Omega <- Function_Omega(K)
  F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
  if (file.exists(here::here("data", "simulation_example", F.name))){
    load(here::here("data", "simulation_example", F.name))
  } else {
    start_time <- Sys.time()
    F <- Function_F(Nr, K, Omega)
    save(F, file = here::here("data", "simulation_example", F.name))
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  for (i in 1:20) {
    dat <-  matrix(c(h)*c(y.sim[[i]]), Nr, Nr)
    coef <- Function_coef_FFT(dat, Omega, Nr)$rv
    y.hw.tilde.lp[[i]] <- coef
  }
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
  V <- diag(0.005, nrow(F.tilde.A))
  W1 <- diag(0.005, K^2)
  W2 <- diag(0.001, K^2)
  W <- rbind(
    cbind(W1, matrix(0, nrow(W1), ncol(W1))),
    cbind(matrix(0, nrow(W2), ncol(W2)), W2)
  )
  m0 = rnorm(2*K^2, 0, 1)
  C0 = diag(0.1, 2*K^2)
  fitHW.name = paste(c("fitHW", "K", as.character(K), "RData"), collapse = ".")
  if(file.exists(here::here("data", "simulation_example", fitHW.name))){
    load(here::here("data", "simulation_example", fitHW.name))
  } else{
    start_time <- Sys.time()
    Fit0.KF <- KF_Non_Missing(y.hw.tilde.lp, F.tilde.A, G.A, V, W, m0, C0)
    save(Fit0.KF, file = here::here("data", "simulation_example", fitHW.name))
    end_time <- Sys.time()
    end_time - start_time
  }
}

### compare the MAEs
MAE1 <- numeric(5)
MAE2 <- numeric(5)
MAE3 <- numeric(5)

for (K in c(10, 14, 20)){
  model.name = paste(c("fitHW.K", as.character(K), "RData"), collapse =".")
  load(here::here("data", "simulation_example", model.name))
  F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
  load(here::here("data", "simulation_example", F.name))
  F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
  for (i in 11:20){
    tempt = matrix(F.A%*%Fit0.KF$m.flt[[i+1]], Nr, Nr) - matrix(y.sim[[i]], Nr, Nr)
    if (K == 10){
      MAE1[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    } else if(K == 14) {
      MAE2[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    } else if(K == 20) {
      MAE3[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
    }
  }
}

MAE4 <- numeric(5)
load(here::here("data", "simulation_example", "Fit.KF.f.RData"))
load(here::here("data", "simulation_example", "F.f.Nr.200.K.f.20.RData"))
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
for (i in 11:20) {
  tempt = matrix(F.f.A%*%fit.KF.f$m.flt[[i+1]], 2*Nr, 2*Nr)[1:100, 1:100] - matrix(y.sim[[i]], Nr, Nr)
  MAE4[i-10] =  sum(abs(tempt))/(nrow(tempt)*ncol(tempt))
}

MAE1
MAE2
MAE3
MAE4

















