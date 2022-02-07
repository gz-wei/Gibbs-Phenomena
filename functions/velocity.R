###Name: Average velocity using the optical flow method
###Author: Guanzhou Wei
###Last revision (dd/mm/yyyy): 08/06/2021
###Required packages: "fields"

# rm(list = ls()) 
# library(fields)
# load(here::here("data", "dat.velocity.Rdata"))
# N.s <- 19
# N.e <- 19

velocity <- function(dat.velocity, N.s, N.e){
  N_r <- sqrt(nrow(dat.velocity[[1]]))
  speed <- matrix(0, N_r^2, (N.e-N.s+1))
  s <- expand.grid(
    x = seq(0, 1-1/N_r, length.out = N_r),
    y = seq(0, 1-1/N_r, length.out = N_r)
  )
  for (i in N.s:N.e){
    speed[ ,(i-N.s+1)] <- dat.velocity[[i]]$speed
  }
  angle <- matrix(0, N_r^2, (N.e-N.s+1))
  for (i in N.s:N.e){
    angle[ ,(i-N.s+1)] <- dat.velocity[[i]]$angle
  }
  speed.avg <- matrix(apply(speed, 1, mean), N_r, N_r)
  angle.avg <- matrix(apply(angle, 1, mean), N_r, N_r)
  speed.avg.smooth <- image.smooth(speed.avg)$z/N_r
  angle.avg.smooth <- image.smooth(angle.avg)$z
  dat.velocity.avg.smooth <- data.frame(s, speed=c(speed.avg.smooth), angle=c(angle.avg.smooth))
  return(dat.velocity.avg.smooth)
}

