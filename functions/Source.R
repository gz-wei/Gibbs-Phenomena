Source <- function(x, x_src1, intensity, boundray, range) {
  output <- numeric(nrow(x))
  for(i in 1:nrow(x)) {
    output[i] = 
      ifelse(norm(x[i,]-x_src1,type="2") > boundray, 0, intensity/(2*pi*range^2)*exp(-norm(x_src1-x[i, ],type = "2")/(2*range^2)))[1] 
  }
  return(output)
}

