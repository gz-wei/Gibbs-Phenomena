Source <- function(x_src1, x_src2, x_src3, x_src4, x) {
  s1 = 3
  h1 = 0.18
  s2 = 2
  h2 = 0.15
  s3 = 2
  h3 = 0.15
  s4 = 3
  h4 = 0.25
  output <- numeric(nrow(x))
  for(i in 1:nrow(x)) {
    output[i] = s1/(2*pi*h1^2)*exp(-norm(x_src1-x[i, ],type="2")/(2*h1^2)) + 
      ifelse(is.na(x_src2),0,s2/(2*pi*h2^2)*exp(-norm(x_src2-x[i, ],type = "2")/(2*h2^2)))[1] +
      ifelse(is.na(x_src3),0,s3/(2*pi*h3^2)*exp(-norm(x_src3-x[i, ],type = "2")/(2*h3^2)))[1] +
      ifelse(is.na(x_src4),0,s4/(2*pi*h4^2)*exp(-norm(x_src4-x[i, ],type = "2")/(2*h4^2)))[1]
  }
  return(output)
}