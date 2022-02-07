Function_coef_FFT <- function(dat, Omega, Nr) {
  n.x = n.y = Nr
  tmp <- t(mvfft(t(dat)))
  FT = mvfft(tmp)/n.x/n.y
  mvfft.x = seq(0,n.x-1,1)
  mvfft.y = seq(0,n.y-1,1)
  mvfft.z =  FT
  my.x = seq(-n.x/2+1, n.x/2, 1)
  my.y = seq(-n.y/2+1, n.y/2, 1)
  my.s = array(0/0, dim=c(n.x,n.y))
  for (i in 1:n.x){
    for (j in 1:n.y){
      target.x = my.x[i]
      target.y = my.y[j]
      
      case1 = which(mvfft.x == target.x)
      case2 = which(mvfft.y == target.y) 
      
      if (length(case1)+length(case2)==2){
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case1,case2]), 
                             imaginary = Im(FT[case1,case2]))
      }
      
      if (length(case1)-length(case2)==1){
        case3 = which(mvfft.x == target.x)
        case4 = which(mvfft.y == target.y+n.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)-length(case1)==1){
        case3 = which(mvfft.x == target.x+n.x)
        case4 = which(mvfft.y == target.y ) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = Im(FT[case3,case4]))
      }
      
      if (length(case2)+length(case1)==0){
        case3 = which(mvfft.x == -target.x)
        case4 = which(mvfft.y == -target.y) 
        my.s[i,j] =  complex(length.out = 0, real = Re(FT[case3,case4]), 
                             imaginary = -Im(FT[case3,case4]))
      }
    }
  }
  Omega1 <- Omega$Omega1
  Omega2 <- Omega$Omega2
  coef1 <- numeric(nrow(Omega1))
  for(i in 1:nrow(Omega1)){
    id.x <- which(my.x == Omega1[i,1])
    id.y <- which(my.y == Omega1[i,2])
    coef1[i] <- Re(my.s[id.x, id.y])
  }
  coef2c <- numeric(nrow(Omega2))
  for(i in 1:nrow(Omega2)){
    id.x <- which(my.x == Omega2[i,1])
    id.y <- which(my.y == Omega2[i,2])
    coef2c[i] <- Re(my.s[id.x, id.y])
  }
  coef2s <- numeric(nrow(Omega2))
  for(i in 1:nrow(Omega2)){
    id.x <- which(my.x == Omega2[i,1])
    id.y <- which(my.y == Omega2[i,2])
    coef2s[i] <- -Im(my.s[id.x, id.y])
  }
  rv <- c(coef1, coef2c, coef2s)
  return(list(rv = rv, my.s = my.s))
}

