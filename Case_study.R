library(rasterImage)
library(SpatialVx) # for OF function
library(tidyverse)
library(expm)
library(MASS)


rm(list = ls()) 

### plot the dataset
load(here::here("data", "case", "data.rf.RData"))
T = length(data.rf)
Nr <- 100
s <- expand.grid(
  x = seq(0, 1-1/Nr, length.out = Nr),
  y = seq(0, 1-1/Nr, length.out = Nr)
)
for (t in 1:T){
  rasterImage2(
    z = data.rf[[t]], 
    zlim = c(0, 120),
    x = seq(0, 1-1/Nr, length.out = Nr),
    y = seq(0, 1-1/Nr, length.out = Nr)
  )
  Sys.sleep(0.5)
}

### low-pass filtering for the original data with K=10
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
source(here::here("functions", "Function_F.R"))
K = 10
Omega <- Function_Omega(K)
F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
if (file.exists(here::here("data", "case", F.name))){
  load(here::here("data", "case", F.name))
} else{
  F <- Function_F(Nr, K, Omega)
  save(F, file = here::here("data", "case", F.name))
} 
### plot the low-pass data stream 
data.rf.lp <- list()
for (t in 1:T) {
  coef <- Function_coef_FFT(data.rf[[t]], Omega, Nr)$rv
  data.rf.lp[[t]] <- coef
  tempt <- matrix(F%*%coef, Nr, Nr)
  rasterImage2(z=tempt, zlim = c(0,60))
  Sys.sleep(0.5)
}

### estimation of the velocity field uisng OF method (package issue)
# data(hump)
# initial <- hump$initial
# final <- hump$final
# look <- OF(final, xhat=initial, W=9, verbose=TRUE)
# plot(look, full=TRUE)
# Nr <- nrow(initial)
# s <- expand.grid(
#   x = seq(0, 1-1/Nr, length.out = Nr),
#   y = seq(0, 1-1/Nr, length.out = Nr)
# )
# speed <- matrix(look$err.mag.lin, Nr, Nr)
# angle <- matrix(look$err.ang.lin, Nr, Nr)
# v <- data.frame(s, speed=c(speed), angle=c(angle/180*pi))
# ggplot(v, aes(x, y, fill=speed,angle=angle,radius=scales::rescale(speed, c(0.01, 0.02)))) +
#   geom_raster() +
#   geom_spoke(arrow = arrow(length = unit(0.015, 'inches'))) +
#   scale_fill_viridis_c(option = "plasma", limits = c(NA,NA)
# )
# if (file.exists(here::here("data", "case", "OF.v.RData"))){
#   load(here::here("data", "case", "OF.v.RData"))
# } else{
#   OF.v <- list()
#   for (i in 1:(T-1)){
#     initial <- data.rf.lp[[i]]
#     final <- data.rf.lp[[i+1]]
#     OF.v[[i]] <- OF(final, xhat=initial, W=20, verbose = TRUE)
#   }
#   save(OF.v, file=here::here("data", "case", "OF.v.Rdata"))
# } 
# v <- list()
# for (i in 1:(T-1)){
#   speed <- matrix(OF.v[[i]]$err.mag.lin, Nr, Nr)
#   angle <- matrix(OF.v[[i]]$err.ang.lin, Nr, Nr)
#   v[[i]] <- data.frame(s, speed=c(speed), angle=c(angle/180*pi))
# }
# source(here::here("functions", "velocity.R"))
# v.smt <- velocity(v, 3, 5)
# v.smt.y <- v.smt$speed*sin(v.smt$angle)
# v.smt.x <- v.smt$speed*cos(v.smt$angle)
# ggplot(v.smt, aes(x, y, fill=speed,angle=angle,radius=scales::rescale(speed, c(0.005, 0.01)))) +
#   geom_raster() +
#   geom_spoke(arrow = arrow(length = unit(0.015, 'inches'))) +
#   scale_fill_viridis_c(option = "plasma", limits = c(NA,NA)
# )

### calculate the transition matrix 
source(here::here("functions", "G_nit.R"))
v = c(1/100, -1/200) ## for a given velocity field (need to be modified later)
G <- G_nit(Omega, v)
step <- 1
G.step <- expm(step*G)
G.A <- rbind(
  cbind(G.step, diag(K^2)),
  cbind(matrix(0,K^2,K^2), diag(K^2))
)
F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))

### estimate parameters 


### KF with the specified V and W 
source(here::here("functions", "KF_Non_Missing.R"))
V <- diag(0.01, length(data.rf.lp[[1]]))
W <- diag(0.01, 2*K^2)
m0 = rnorm(2*K^2, 0, 1)
C0 = diag(0.1, 2*K^2)
if (file.exists(here::here("data", "case", "fit.KF.RData"))) {
  load(here::here("data", "case", "fit.KF.RData"))
} else{
  start_time <- Sys.time()
  fit.KF <- KF_Non_Missing(data.rf.lp, F.tilde.A, G.A, V, W, m0, C0)
  end_time <- Sys.time()
  print(end_time - start_time)
  save(fit.KF, file = here::here("data", "case", "fit.KF.RData"))
}
### plot the filtered data stream
F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
for (i in 1:T){
  tempt <- matrix(F.A%*%fit.KF$m.flt[[i+1]], Nr, Nr)
  rasterImage2(
    z=tempt, zlim = c(0,60), 
    x = seq(0, 1-1/Nr, length.out = Nr),
    y = seq(0, 1-1/Nr, length.out = Nr)
  )
  Sys.sleep(0.5)
}

### for the flipping method 
### plot the flipping datastream
dim <- Nr
I11 <- diag(1, dim)
I12 <- I11[dim:1,]
I1 <- rbind(I11, I12)
I2 <- cbind(I11, I12)
data.rf.flp <- list()
for (t in 1:T){
  dat <-  data.rf[[t]]
  data.rf.flp[[t]] <- I1%*%dat%*%I2
  rasterImage2(
    z = data.rf.flp[[t]], 
    zlim = c(0, 120),
    x = seq(0, 2-1/(2*Nr), length.out = 2*Nr),
    y = seq(0, 2-1/(2*Nr), length.out = 2*Nr)
  )
  Sys.sleep(0.5)
}
### low-pass filtering for the flipped data with K.f=20
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
source(here::here("functions", "Function_F.R"))
K.f = 2*K
Omega.f <- Function_Omega(K.f)
F.f.name <- paste(c("F.f", "Nr", as.character(2*Nr), "K.f", as.character(2*K), "RData"), collapse = ".")
if (file.exists(here::here("data", "case", F.f.name))){
  load(here::here("data", "case", F.f.name))
} else{
  F.f <- Function_F(2*Nr, K.f, Omega.f)
  save(F.f, file = here::here("data", "case", F.f.name))
} 
### plot the low-pass data stream
data.rf.flp.lp <- list()
for (t in 1:T) {
  coef.f <- Function_coef_FFT(data.rf.flp[[t]], Omega.f, 2*Nr)$rv
  data.rf.flp.lp[[t]] <- coef.f
  tempt <- matrix(F.f%*%coef.f, 2*Nr, 2*Nr)
  rasterImage2(z=tempt, zlim = c(0,120))
  Sys.sleep(0.3)
}

### calculate the transition matrix G*
dim <- Nr
I <- diag(1, dim)
J <- I[dim:1,]
P <- rbind(I,J)
if (file.exists(here::here("data", "case", "H.RData"))){
  load(here::here("data", "case", "H.RData"))
} else{
  H <- ginv(F.f)%*%kronecker(P,P)%*%F
  save(H, file = here::here("data", "case", "H.RData"))
}
G.f <- H%*%G%*%ginv(H)
G.step.f <- expm(step*G.f)
G.f.A <- rbind(
  cbind(G.step.f, diag(K.f^2)),
  cbind(matrix(0,K.f^2, K.f^2), diag(K.f^2))
)
# compute F*
F.f.tilde.A <- cbind(diag(K.f^2), matrix(0, K.f^2, K.f^2))
# dim(F.f.tilde.A)

### KF with the specified V and W 
#V.f <- diag(0.01, length(data.rf.flp.lp[[1]]))
#W.f <- diag(0.01, 2*K.f^2)
V.f <- I1%*%V%*%I2
W.f <- G.f%*%W
m0.f = rnorm(2*K.f^2, 0, 1)
C0.f = diag(0.1, 2*K.f^2)
if (file.exists(here::here("data", "case", "fit.KF.f.RData"))){
  load(here::here("data", "case", "fit.KF.f.RData"))
} else{
  start_time <- Sys.time()
  fit.KF.f <- KF_Non_Missing(data.rf.flp.lp, F.f.tilde.A, G.f.A, V.f, W.f, m0.f, C0.f)
  end_time <- Sys.time()
  print(end_time - start_time)
  save(fit.KF.f, file = here::here("data", "case", "fit.KF.f.RData"))
}
### plot the filtered data stream
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
for (i in 1:T){
  tempt <- matrix(F.f.A%*%fit.KF.f$m.flt[[i+1]], 2*Nr, 2*Nr)[1:Nr, 1:Nr]
  rasterImage2(z=tempt, zlim = c(0,120))
  Sys.sleep(0.3)
}
