Phi1.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {2*pi*(s[3]*k1+s[4]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi2.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {-2*pi*(s[3]*k1+s[4]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row) 
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi3.int <- function(k1, k2, m1, m2, delta, v) {
  f <- function(s) {2*pi*(s[3]*k1+s[4]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), v))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi5.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2))*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi6.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2))*cos(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi7.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*cos(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*sin(2*pi*s[1]*k1+2*pi*s[2]*k2))*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}

Phi8.int <- function(k1, k2, m1, m2, delta, K.ifm) {
  f <- function(s) {(-4*pi^2*(k1^2+k2^2)*s[3]*sin(2*pi*s[1]*k1+2*pi*s[2]*k2)-2*pi*(s[4]*k1+s[5]*k2)*cos(2*pi*s[1]*k1+2*pi*s[2]*k2))*sin(2*pi*s[1]*m1+2*pi*s[2]*m2)}
  N.row <- 1/delta
  x <- seq(0, 1-1/N.row, length.out = N.row)
  y <- seq(0, 1-1/N.row, length.out = N.row)
  z <- as.matrix(cbind(expand.grid(x, y), K.ifm))
  rv <- sum(apply(z, 1, f))*delta^2
  return(rv)
}