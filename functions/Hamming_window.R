Func_hamming = function(N1, N2){
  rv = matrix(0, N1, N2)
  for (i in 1:N1){
    for (j in 1:N2){
      rv[i,j] = (0.54 - 0.46*cos(2*pi*(i-1)/(N1-1)))*(0.54 - 0.46*cos(2*pi*(j-1)/(N2-1)))
    }
  }
  return(rv)
}