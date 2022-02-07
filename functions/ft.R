# function for fourier transform 

ft <- function(dat, k1, k2){
  f_real <- function(n) {dat[n[1], n[2]]*cos(2*pi/N*(n[1]-1)*k1+2*pi/N*(n[2]-1)*k2)}
  f_imaginary <- function(n) {dat[n[1], n[2]]*sin(2*pi/N*(n[1]-1)*k1+2*pi/N*(n[2]-1)*k2)}
  N <- nrow(dat)
  n1 <- seq(1, N, 1)
  n2 <- seq(1, N, 1)
  z <- as.matrix(expand.grid(n1, n2))
  rv_real <- sum(apply(z, 1, f_real))/N^2
  rv_imaginary <- sum(apply(z, 1, f_imaginary))/N^2
  return(list(real = rv_real, imaginary = rv_imaginary))
}