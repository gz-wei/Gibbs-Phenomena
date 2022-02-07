# source(here::here("functions", "KF_Non_Missing.R"))
# y.tilde <- data.rf.lp
# v = 0.01
# w = 0.01
# K = 20
l <- function(G.A, C0, y.tilde, v, w, K){
  V <- diag(v, length(y.tilde[[1]]))
  W <- diag(w, 2*K^2)
  m0 = rnorm(2*K^2, 0, 1)
  C0 = diag(0.1, 2*K^2)
  F.tilde.A <- cbind(diag(K^2), matrix(0, K^2, K^2))
  fit.KF <- KF_Non_Missing(data.rf.lp, F.tilde.A, G.A, V, W, m0, C0)
  sum.t <- numeric(T)
  for (t in 1:T){
    tempt1 <- fit.KF$C.prd[[t]][1:K^2, 1:K^2] + V
    tempt2 <- y.tilde[[t]] - fit.KF$m.prd[[t]][1:K^2]
    sum.t[t] <- -1/2*(log(det(tempt1)) + t(tempt2)%*%solve(tempt1)%*%tempt2)
  }
  return(sum(sum.t))
}