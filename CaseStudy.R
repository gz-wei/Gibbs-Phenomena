library(rasterImage)
library(tidyverse)
library(expm)
library(MASS)
library(mgcv)
library(matrixStats)
library(fields)
library(MARSS)


rm(list = ls()) 
#*******************************************************************************************
### plot the dataset
#*******************************************************************************************
load(here::here("data", "case", "data.rf.RData"))
T = length(data.rf)
Nr <- 100
s <- expand.grid(
  x = seq(0, 1, length.out = Nr),
  y = seq(0, 1, length.out = Nr)
)
for (t in 1:T){
  rasterImage2(
    z = data.rf[[t]], 
    zlim = c(0, 120),
    x = seq(0, 1, length.out = Nr),
    y = seq(0, 1, length.out = Nr)
  )
  Sys.sleep(0.6)
}

#*******************************************************************************************
### low-pass filtering for the original data with K=10
#*******************************************************************************************
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
source(here::here("functions", "Function_F.R"))
K = 10
Omega <- Function_Omega(K)
F.name <- paste(c("F", "Nr", as.character(Nr), "K", as.character(K), "RData"), collapse = ".")
if (file.exists(here::here("data", "case", F.name))){
  load(here::here("data", "case", F.name))
} else{
  start_time <- Sys.time()
  F <- Function_F(Nr, K, Omega)
  save(F, file = here::here("data", "case", F.name))
  end_time <- Sys.time()
  print(end_time - start_time)
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

#*******************************************************************************************
### estimate the wind field uing the radar images 
#*******************************************************************************************
## track the velocity field using my hand-coded COTREC function 
load(here::here("data", "case", "data.rd.RData"))
if(file.exists(here::here("data", "case", "fit.v.RData"))){
  load(here::here("data", "case", "fit.v.RData"))
} else{
  start_time <- Sys.time()
  source(here::here("functions", "COTREC.R"))
  fit.v = list()
  for (i in 1:(length(data.rd)-1)){
    fit.v[[i]] = COTREC(data.rd[[i]], data.rd[[i+1]], 5, 4, 4)
  }
  save(fit.v, file = here::here("data", "case", "fit.v.RData"))
  end_time <- Sys.time()
  print(end_time - start_time)
}
## smooth the velocity feild 
tempt.vx = tempt.vy = matrix(0, Nr, Nr)
for (i in 1:(length(data.rd)-1)){
  tempt.vx = fit.v[[i]]$vs.x + c(tempt.vx)
  tempt.vy = fit.v[[i]]$vs.y + c(tempt.vy)
}
vs.x = tempt.vx/(length(data.rd)-1)
vs.y = tempt.vy/(length(data.rd)-1)
site = expand.grid(
  x = seq(0,1,length.out = nrow(data.rd[[1]])),
  y = seq(0,1,length.out = ncol(data.rd[[1]]))
)
v.s = data.frame(
  x = site$x,
  y = site$y,
  v.x = vs.x,
  v.y = vs.y
)
## plot the smoothed velocity field
ggplot(v.s, aes(x = x, y = y, fill = v.x)) +
  geom_segment(aes(xend = x+v.x, yend = y+v.y),
               arrow = arrow(length = unit(0.05, "cm")), size = 0.15)
## plot the low-resolution velocity field
every_n <- function(x, by = by) {
  x <- sort(x)
  x[seq(1, length(x), by = by)]
}
keepx <- every_n(unique(v.s$x), by = 6)
keepy <- every_n(unique(v.s$y), by = 6)
vsub <- filter(v.s, x %in% keepx  &  y %in% keepy)
ggplot(vsub, aes(x = x, y = y)) +
  geom_segment(aes(xend = x+v.x, yend = y+v.y),
               arrow = arrow(length = unit(0.05, "cm")), size = 0.15)

#*******************************************************************************************
### calculate the transition matrix in the dynamical mode 
#*******************************************************************************************
### compute the diffusion parameters 
vs.x = matrix(vs.x, Nr, Nr)
vs.y = matrix(vs.y, Nr, Nr)
col.dif.x <- rowDiffs(vs.x)/(1/Nr)
row.dif.x <- colDiffs(vs.x)/(1/Nr)
col.dif.y <- rowDiffs(vs.y)/(1/Nr)
row.dif.y <- colDiffs(vs.y)/(1/Nr)
p1 <- cbind(col.dif.x, rep(NA,Nr))
p2 <- rbind(row.dif.y, rep(NA,Nr))
p3 <- rbind(colDiffs(vs.x), rep(NA,Nr))
p4 <- cbind(rowDiffs(vs.y), rep(NA,Nr))
D <- 0.28*1/Nr*1/Nr*sqrt((p1-p2)^2+(p3+p4)^2)
D.smooth <- image.smooth(D)$z
D = matrix(D.smooth, Nr, Nr)
col.dif <- rowDiffs(D)
row.dif <- colDiffs(D)
D.D.x <- cbind(col.dif, col.dif[,2])/(1/Nr)
D.D.y <- rbind(row.dif, row.dif[2,])/(1/Nr)
D.ifm <- data.frame(D=c(D), D.D.x=c(D.D.x), D.D.y=c(D.D.y))

G.name <- paste(c("G", "K", as.character(K), "RData"), collapse = ".")
if(file.exists(here::here("data", "case", G.name))){
  load(here::here("data", "case", G.name))
} else{
  start_time <- Sys.time()
  source(here::here("functions", "G_ad.R"))
  source(here::here("functions", "Phi.R"))
  v = data.frame(v.x = vs.x, v.y = vs.y)
  G <- G_ad(1/Nr, v, K.ifm, Omega)
  save(G, file = here::here("data", "case", G.name))
  end_time <- Sys.time()
  print(end_time-start_time)
}

step <- 1
G.step <- expm(step*G)
G.A <- rbind(
  cbind(G.step, diag(K^2)),
  cbind(matrix(0,K^2,K^2), diag(K^2))
)

#*******************************************************************************************
### KF with estimated parameters for the original data stream
#*******************************************************************************************
source(here::here("functions", "KF_Non_Missing.R"))
F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))
V <- diag(0.001, length(data.rf.lp[[1]]))
W1 <- diag(0.0049, K^2)
W2 <- diag(0.0001, K^2)
W <- rbind(
  cbind(W1, matrix(0, nrow(W1), ncol(W1))),
  cbind(matrix(0, nrow(W2), ncol(W2)), W2)
)
m0 = rnorm(2*K^2, 0, 1)
C0 = diag(0.1, 2*K^2)
fit0.KF <- KF_Non_Missing(data.rf.lp, F.tilde.A, G.A, V, W, m0, C0)
### plot the filtered data stream
F.A <- cbind(F, matrix(0, nrow(F), ncol(F)))
for (i in 1:T){
  tempt <- matrix(F.A%*%fit0.KF$m.flt[[i+1]], Nr, Nr)
  rasterImage2(
    z=tempt, zlim = c(0,60), 
    x = seq(0, 1, length.out = Nr),
    y = seq(0, 1, length.out = Nr)
  )
  Sys.sleep(0.5)
}
### test the physical meaning 
image(matrix(F.A%*%fit0.KF$m.flt[[10]], Nr, Nr))
image(matrix(F.A%*%G.A%*%fit0.KF$m.flt[[10]], Nr, Nr))
image(matrix(F.A%*%fit0.KF$m.flt[[11]], Nr, Nr))


#*******************************************************************************************
### for the flipping method 
#*******************************************************************************************
### plot the flipping datastream
dat_flp = function(data, method){
  dim = nrow(data)
  if(method =="rightup"){
    I11 <- diag(1, dim)
    I12 <- I11[dim:1,]
    I1 <- rbind(I11, I12)
    I2 <- cbind(I11, I12)
    output = I1%*%dat%*%I2
  } else if(method =="rightdown"){
    I11 <- diag(1, dim)
    I12 <- I11[dim:1,]
    I1 <- rbind(I11, I12)
    I2 <- cbind(I12, I11)
    output = I1%*%dat%*%I2
  } else if(method =="leftup"){
    I11 <- diag(1, dim)
    I12 <- I11[dim:1,]
    I1 <- rbind(I12, I11)
    I2 <- cbind(I11, I12)
    output = I1%*%dat%*%I2
  } else if(method =="leftdown"){
    I11 <- diag(1, dim)
    I12 <- I11[dim:1,]
    I1 <- rbind(I12, I11)
    I2 <- cbind(I12, I11)
    output = I1%*%dat%*%I2
  } else{
    stop("method can only be one of rightup, rightdown, leftup, leftdown")
  }
  return(list(data.flp = output, I11 = I11, I12 = I12))
}
method = "rightup"
data.rf.flp <- list()
for (t in 1:T){
  dat <-  data.rf[[t]]
  data.rf.flp[[t]] <- dat_flp(dat, method)$data.flp
  rasterImage2(
    z = data.rf.flp[[t]], 
    zlim = c(0, 120),
    x = seq(0, 2, length.out = 2*Nr),
    y = seq(0, 2, length.out = 2*Nr)
  )
  Sys.sleep(0.5)
}

#*******************************************************************************************
### low-pass filtering for the flipped data with K.f=20
#*******************************************************************************************
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

#*******************************************************************************************
### calculate the transition matrix G* and F*
#*******************************************************************************************
dim <- Nr
I <- dat_flp(data.rf[[1]], method)$I11
J <- dat_flp(data.rf[[1]], method)$I12
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
F.f.tilde.A <- cbind(diag(K.f^2), matrix(0, K.f^2, K.f^2))


#*******************************************************************************************
### KF with estimated parameters for the flipped data stream
#*******************************************************************************************
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
if (file.exists(here::here("data", "case", "fit.f.KF.RData"))){
  load(here::here("data", "case", "fit.f.KF.RData"))
} else{
  start_time <- Sys.time()
  fit.f.KF <- KF_Non_Missing(data.rf.flp.lp, F.f.tilde.A, G.f.A, V.f, W.f, m0.f, C0.f)
  end_time <- Sys.time()
  print(end_time - start_time)
  save(fit.f.KF, file = here::here("data", "case", "fit.KF.f.RData"))
}
### plot the filtered data stream
F.f.A <- cbind(F.f, matrix(0, nrow(F.f), ncol(F.f)))
for (i in 1:T){
  tempt <- matrix(F.f.A%*%fit.f.KF$m.flt[[i+1]], 2*Nr, 2*Nr)[1:Nr, 1:Nr]
  rasterImage2(z=tempt, zlim = c(0,60))
  Sys.sleep(0.3)
}
### test the physical meaning 
image(matrix(F.f.A%*%fit.f.KF$m.flt[[10]], 2*Nr, 2*Nr)[1:Nr, 1:Nr])
image(matrix(F.f.A%*%G.f.A%*%fit.f.KF$m.flt[[10]], 2*Nr, 2*Nr)[1:Nr, 1:Nr])
image(matrix(F.f.A%*%fit.f.KF$m.flt[[11]], 2*Nr, 2*Nr)[1:Nr, 1:Nr])
