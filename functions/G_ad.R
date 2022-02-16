###Name: Trainsition matrix calculation with physical parameters
###Author: Guanzhou Wei
###Last revision (dd/mm/yyyy): 06/06/2021
###Required functions: "Phi1.int.R","Phi2.int.R","Phi3.int.R","Phi5.int.R","Phi6.int.R","Phi7.int.R","Phi8.int.R"


G_ad <- function(delta,v,K.ifm,Omega){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  N_Gr <- 2*nrow(Omega2) + nrow(Omega1)
  G <- matrix(0, N_Gr, N_Gr)
  k1 <- nrow(Omega1)
  k2 <- nrow(Omega2)
  ##for real part in Omega1 
  ##for special case, m = 1
  for (k in 1:k1){
    G[1, k] <- Phi1.int(Omega1[k,1],Omega1[k,2],Omega1[1,1],Omega1[1,2],delta,v) +
      Phi5.int(Omega1[k,1],Omega1[k,2],Omega1[1,1],Omega1[1,2],delta,K.ifm)
  }
  for (k in (k1+1):(k1+k2)){
    G[1, k] <- 2*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[1,1],Omega1[1,2],delta,v) + 
      2*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[1,1],Omega1[1,2],delta,K.ifm)
  }
  for (k in (k1+k2+1):(k1+2*k2)){
    G[1, k] <- 2*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[1,1], Omega1[1,2],delta,v) + 
      2*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[1,1], Omega1[1,2],delta,K.ifm)
  }
  ##for m!=0 in Omega1 
  for (m in 2:k1){
    for (k in 1:k1){
      G[m, k] <- 2*Phi1.int(Omega1[k,1],Omega1[k,2],Omega1[m,1],Omega1[m,2],delta,v) + 
        2*Phi5.int(Omega1[k,1],Omega1[k,2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m, k] <- 4*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[m,1],Omega1[m,2],delta,v) + 
        4*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m, k] <- 4*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[m,1],Omega1[m,2],delta,v) + 
        4*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega1[m,1],Omega1[m,2],delta,K.ifm)
    }
  }
  ##for real part in Omega2 
  for (m in 1:k2){
    for (k in 1:k1){
      G[m+k1, k] <- Phi1.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        Phi5.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m+k1, k] <- 2*Phi1.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        2*Phi5.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m+k1, k] <- 2*Phi2.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega2[m,1],Omega2[m,2],delta,v) +
        2*Phi6.int(Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
  }
  ##for imaginary part in Omega2
  for (m in 1:k2){
    for (k in 1:k1){
      G[m+k1+k2, k] <- Phi3.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        Phi7.int(Omega1[k,1],Omega1[k,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+1):(k1+k2)){
      G[m+k1+k2, k] <- 2*Phi3.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,v) + 
        2*Phi7.int(Omega2[k-k1,1],Omega2[k-k1,2],Omega2[m,1],Omega2[m,2],delta,K.ifm)
    }
    for (k in (k1+k2+1):(k1+2*k2)){
      G[m+k1+k2, k] <- -2*Phi1.int(Omega2[m,1],Omega2[m,2],Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],delta,v) +
        2*Phi8.int(Omega2[m,1],Omega2[m,2],Omega2[k-(k1+k2),1],Omega2[k-(k1+k2),2],delta,K.ifm)
    }
  }
  return(G)
}
