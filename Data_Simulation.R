library(expm)
rm(list = ls()) 

## simultate a source
source(here::here("functions", "Source2.R"))
Nr <- 100
s <- expand.grid(
  x = seq(0, 1-1/Nr, length.out = Nr),
  y = seq(0, 1-1/Nr, length.out = Nr)
)
dat <- matrix(Source(c(0.1, 0.0), NA, NA, NA, s), Nr, Nr)
#range(dat)
image(dat)

## compute the spectral coefficients using my hand-writing Fourier transformation
# source(here::here("functions", "Function_Omega.R"))
# source(here::here("functions", "Function_coef.R"))
# source(here::here("functions", "ft.R"))
# start_time <- Sys.time()
# N <- 100
# Omega <- Function_Omega(N)
# coef <- Function_coef(dat, N, Omega)
# #save(coef, file = here::here("data", "N=100", "coef.RData"))
# end_time <- Sys.time()
# end_time - start_time

## compute the spectral coefficients using FFT 
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
N <- 100
Omega <- Function_Omega(N)
coef <- Function_coef_FFT(dat, Omega, Nr)$rv
#rasterImage2(z = Re(Function_coef_FFT(dat, Omega, Nr)$my.s))
#rasterImage2(z = Im(Function_coef_FFT(dat, Omega, Nr)$my.s))


## recover the source image 
# compute the Fourier bases F
if (file.exists(here::here("data", "N=100", "F.RData"))){
  load(here::here("data", "N=100", "F.RData"))
} else{
  start_time <- Sys.time()
  source(here::here("functions", "Function_F.R"))
  F <- Function_F(Nr, N, Omega)
  save(F, file = here::here("data", "N=100", "F.RData"))
  end_time <- Sys.time()
  print(end_time - start_time)
}
dat_rc <- matrix(F%*%coef, Nr, Nr)
image(dat_rc)
#range(dat_rc)

## compute the transition matrix
if (file.exists(here::here("data", "N=100", "G.step.0.01.RData"))){
  load(here::here("data", "N=100", "G.step.0.01.RData"))
} else{
  source(here::here("functions", "G_nit.R"))
  start_time <- Sys.time()
  v <- c(0.01, 0)
  G <- G_nit(Omega, v)
  step <- 1
  G.step <- expm(step*G)
  save(G.step, file = here::here("data", "N=100", "G.step.0.01.RData"))
  end_time <- Sys.time()
  print(end_time - start_time)
}

## visualization
N.step <- 30
alpha <- coef
noise.y.tilde <- 0.005
noise.coef <- 0.01
noise.source <- 0.005
y.sim <- list()
y.sim[[1]] <- pmax(F%*%alpha,0) + F%*%rnorm(nrow(F), 0, noise.y.tilde)
for (i in 2:N.step){
  coef <- coef + rnorm(length(coef), 0, noise.source)
  alpha <- G.step%*%alpha + coef + rnorm(length(alpha), 0, noise.coef)
  y.sim[[i]] <- pmax(F%*%alpha,0) + F%*%rnorm(nrow(F), 0, noise.y.tilde)
}
for (i in 1:N.step){
  tempt <- matrix(y.sim[[i]], Nr, Nr)
  image(tempt)
  print(range(tempt))
  Sys.sleep(0.3)
}
save(y.sim, file = here::here("data", "N=100", "y.sim.RData"))

##### test the edge 
# dat.test <- matrix(y.sim[[20]], Nr, Nr)
# image(dat.test)
# range(dat.test)
# source(here::here("functions", "Function_Omega.R"))
# source(here::here("functions", "Function_coef.R"))
# source(here::here("functions", "ft.R"))
# source(here::here("functions", "Function_F.R"))
# N.t <- 60
# Omega.t <- Function_Omega(N.t)
# coef.t <- Function_coef(dat.test, N.t, Omega.t)
# F.t <- Function_F(Nr, N.t, Omega.t)
# dat.rc.test <- matrix(F.t%*%coef.t, Nr, Nr)
# image(dat.rc.test)
# range(dat.rc.test)

## comparasion with the flipping
# dim <- Nr
# I11 <- diag(1, dim)
# I12 <- I11[dim:1,]
# I1 <- rbind(I11, I12)
# I2 <- cbind(I11, I12)
# dat.flip <- I1%*%dat.test%*%I2
# image(dat.flip)
# 
# N.f <- 20
# Omega.f <- Function_Omega(N.f)
# coef.f <- Function_coef(dat.flip, N.f, Omega.f)
# F.f <- Function_F(2*Nr, N.f, Omega.f)
# dat.rc.f <- matrix(F.f%*%coef.f, 2*Nr, 2*Nr)
# image(dat.rc.f)
# range(dat.rc.test)

## low-pass filtering for the original data with N=10
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef.R"))
source(here::here("functions", "ft.R"))
source(here::here("functions", "Function_F.R"))
y.tilde.lp <- list()
N = 8
Omega <- Function_Omega(N)
F <- Function_F(Nr, N, Omega)
for (i in 1:N.step) {
  dat <-  matrix(y.sim[[i]], Nr, Nr)
  coef <- Function_coef_FFT(dat, Omega, Nr)$rv
  y.tilde.lp[[i]] <- coef
}
for (i in 1:N.step){
  tempt <- matrix(F%*%y.tilde.lp[[i]], Nr, Nr)
  image(tempt)
  Sys.sleep(0.3)
}
save(y.tilde.lp, file = here::here("data", "N=100", "y.tilde.lp.RData"))


## low-pass filtering for the flipped data 
dim <- Nr
I11 <- diag(1, dim)
I12 <- I11[dim:1,]
I1 <- rbind(I11, I12)
I2 <- cbind(I11, I12)
y.f <- list()
for (i in 1:N.step){
  dat <-  matrix(y.sim[[i]], Nr, Nr)
  y.f[[i]] <- I1%*%dat%*%I2
}
for (i in 1:N.step){
  tempt <- matrix(y.f[[i]], 2*Nr, 2*Nr)
  image(tempt)
  Sys.sleep(0.3)
}
# K.f = 2K
N.f <- 20
Omega.f <- Function_Omega(N.f)
if (file.exists(here::here("data", "N=100", "F.f.20.RData"))){
  load(here::here("data", "N=100", "F.f.20.RData"))
} else {
  start_time <- Sys.time()
  F.f <- Function_F(2*Nr, N.f, Omega.f)
  save(F.f, here::here("data", "N=100", "F.f.20.RData"))
  end_time <- Sys.time()
  print(end_time - start_time)
}
y.tilde.flp <- list()
for (i in 1:N.step) {
  dat <-  matrix(y.f[[i]], 2*Nr, 2*Nr)
  coef.f <- Function_coef_FFT(dat, Omega.f, 2*Nr)$rv
  y.tilde.flp[[i]] <- coef.f
}
for (i in 1:N.step){
  tempt <- matrix(F.f%*%y.tilde.flp[[i]], 2*Nr, 2*Nr)
  image(tempt)
  Sys.sleep(0.3)
}
#save(y.tilde.flp, file = here::here("data", "N=100", "y.tilde.flp.RData"))

