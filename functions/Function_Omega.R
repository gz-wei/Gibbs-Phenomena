Function_Omega <- function(N){
  Omega1 <- data.frame(
    k1 = c(0, 0, N/2, N/2),
    k2 = c(0, N/2, 0, N/2)
  )
  Omega2_p1 <- data.frame(
    k1 = seq(1, N/2-1, 1),
    k2 = N/2
  )
  Omega2_p2 <- expand.grid(
    k1 = seq(0, N/2, 1),
    k2 = seq(1, N/2-1, 1)
  ) 
  Omega2_p3 <- expand.grid(
    k1 = seq(1, N/2-1, 1),
    k2 = seq(-N/2+1, 0, 1)
  ) 
  Omega2 <- rbind(Omega2_p1, Omega2_p2, Omega2_p3)
  return(list(Omega1=Omega1, Omega2=Omega2))
}