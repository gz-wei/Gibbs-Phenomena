library(MARSS)

m(list = ls()) 
#*************************************************************
### plot the dataset
#*************************************************************
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

#*************************************************************
### estimate the wind field uing the radar images 
#*************************************************************
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

#*************************************************************
### low-pass filtering for the original data with K=6
#*************************************************************
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_coef_FFT.R"))
source(here::here("functions", "Function_F.R"))
K = 6
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

#*************************************************************
### calculate the transition matrix (K=6) in the dynamical mode
#*************************************************************
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
  G <- G_ad(1/Nr, v, D.ifm, Omega)
  save(G, file = here::here("data", "case", G.name))
  end_time <- Sys.time()
  print(end_time-start_time)
}




source(here::here("functions", "KF_Non_Missing.R"))
load(here::here("data", "N=100", "y.tilde.lp.RData"))
source(here::here("functions", "Function_Omega.R"))
source(here::here("functions", "Function_F.R"))
source(here::here("functions", "G_nit.R"))
K = 6
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