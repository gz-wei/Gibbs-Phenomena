Function_F <- function(N_rcr, N, Omega){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  sites_rc <- expand.grid(
    seq(0, 1-1/N_rcr, length.out = N_rcr),
    seq(0, 1-1/N_rcr, length.out = N_rcr)
  )
  F1 <- matrix(0, nrow(sites_rc), nrow(Omega1))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega1)){
      F1[i,j] <- cos(2*pi*sites_rc[i,1]*Omega1[j,1]+2*pi*sites_rc[i,2]*Omega1[j,2])
    }
  }
  F2c <- matrix(0, nrow(sites_rc), nrow(Omega2))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega2)){
      F2c[i,j] <- 2*cos(2*pi*sites_rc[i,1]*Omega2[j,1]+2*pi*sites_rc[i,2]*Omega2[j,2])
    }
  }
  F2s <- matrix(0, nrow(sites_rc), nrow(Omega2))
  for (i in 1:nrow(sites_rc)){
    for(j in 1:nrow(Omega2)){
      F2s[i,j] <- 2*sin(2*pi*sites_rc[i,1]*Omega2[j,1]+2*pi*sites_rc[i,2]*Omega2[j,2])
    }
  }
  rv <- cbind(F1, F2c, F2s)
  return(rv)
}
