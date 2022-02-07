G_nit <- function(Omega, v){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  N_Gr <- 2*nrow(Omega2) + nrow(Omega1)
  G <- matrix(0, N_Gr, N_Gr)
  k1 <- nrow(Omega1)
  k2 <- nrow(Omega2)
  #for real part in Omega2 
  for (m in 1:k2){
    for (k in (k1+k2+1):(k1+2*k2)){
      if (k-(k1+k2) == m){
        G[m+k1, k] <- -2*pi*(v[1]*Omega2[k-(k1+k2),1]+v[2]*Omega2[k-(k1+k2),2])
      }
      else{
        G[m+k1, k] <- 0
      }
    }
  }
  ##for imaginary part in Omega2
  for (m in 1:k2){
    for (k in (k1+1):(k1+k2)){
      if (k-k1 == m){
        G[m+k1+k2, k] <- 2*pi*(v[1]*Omega2[(k-k1),1]+v[2]*Omega2[(k-k1),2])
      }
      else{
        G[m+k1+k2, k] <- 0
      }
    }
  }
  return(G)
}
