###Name: Kalman filter for the dataset without missing value
###Author: Guanzhou Wei
###Last revision: 12/04/2021

# rm(list = ls()) 
# load(here::here("data", "G_step.6.RData"))
# load(here::here("data", "F.6.RData"))
# load(here::here("data", "y_sim.RData"))
# N = 6
# y = y_sim
# m0 = rnorm(2*N^2, 0, 1)
# C0 = diag(0.1, 2*N^2)
# F = cbind(F, matrix(0, nrow(F), ncol(F)))
# G = G_full <- rbind(
#   cbind(G_step, diag(N^2)),
#   cbind(matrix(0,N^2,N^2), diag(N^2))
# )
# V <- diag(0.1, length(y[[1]]))
# W <- diag(0.1, 2*N^2)


KF_Non_Missing <- function(y, F, G, V, W, m0, C0){
  T <- length(y)
  m.flt = m.prd = C.flt = C.prd = vector("list")
  m.flt[[1]] <- m0
  C.flt[[1]] <- C0
  for (t in 1:T){
    m.prd[[t]] <- G%*%m.flt[[t]]
    C.prd[[t]] <- G%*%C.flt[[t]]%*%t(G) + W
    Q <-  F%*%C.prd[[t]]%*%t(F) + V
    Q.inv <- solve(Q)
    m.flt[[t+1]] <- m.prd[[t]]+C.prd[[t]]%*%t(F)%*%Q.inv%*%(y[[t]]-F%*%m.prd[[t]])
    C.flt[[t+1]] <- C.prd[[t]]-C.prd[[t]]%*%t(F)%*%Q.inv%*%F%*%C.prd[[t]]
  }
  m.prd[[T+1]] <- G%*%m.flt[[t]]
  C.prd[[T+1]] <- G%*%C.flt[[t]]%*%t(G) + W
  return(list(m.flt=m.flt, m.prd=m.prd, C.flt=C.flt, C.prd=C.prd))
}


