###Name: Coefficients calculation under Fourier transform
###Author: Guanzhou Wei
###Last revision (dd/mm/yyyy): 06/06/2021
###Required functions: "ft.R"


Function_coef <- function(dat, N, Omega){
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  coef1 <- numeric(nrow(Omega1))
  for(i in 1:nrow(Omega1)){
    coef1[i] <- ft(dat,Omega1[i,1],Omega1[i,2])$real
  }
  coef2c <- numeric(nrow(Omega2))
  for(i in 1:nrow(Omega2)){
    coef2c[i] <- ft(dat,Omega2[i,1],Omega2[i,2])$real
  }
  coef2s <- numeric(nrow(Omega2))
  for(i in 1:nrow(Omega2)){
    coef2s[i] <- ft(dat,Omega2[i,1],Omega2[i,2])$imaginary
  }
  rv <- c(coef1, coef2c, coef2s)
  return(rv)
}
