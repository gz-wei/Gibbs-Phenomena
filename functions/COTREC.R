# library(selextR)
# library(rasterImage)
# library(tidyverse)
# library(grid)
# library(fields)
# library(bnstruct)
# 
# rm(list = ls())
# 
# ### load the data for testing function
# SVY21 <- '+proj=tmerc +lat_0=1.366666666666667 +lon_0=103.8333333333333 +k=1 +x_0=28001.642 +y_0=38744.572 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
# scan <- readRDS(here::here("data", "radar", "scan_list_19"))
# radar.grd = SpatialPoints(coordinates(scan[[1]]$data),proj4string=CRS(SVY21)) # grid points in SVY21 projection system
# Grid = coordinates(radar.grd)
# n.radar.grd = length(radar.grd)
# radar.z1 = scan[[1]]$data@data[,1]
# grd.x = unique(coordinates(Grid)[,1])
# grd.y = rev(unique(coordinates(Grid)[,2]))
# n.x = length(grd.x)
# n.y = length(grd.y)
# Image1 = matrix(radar.z1,nrow=n.y)[, n.y:1]
# rasterImage2(z = Image1)
# radar.z2 = scan[[6]]$data@data[,1]
# Image2 = matrix(radar.z2,nrow=n.y)[, n.y:1]
# rasterImage2(z = Image2)

COTREC <- function(Image1, Image2, WindowSize, overlap, SearchSize) {
  n.x = nrow(Image1)
  n.y = ncol(Image1)
  # WindowSize = 17
  # overlap = 15
  # SearchSize = 4
  if (WindowSize <= overlap|WindowSize%%2 == 0){
    stop("WindowSize must be an odd number and greater than overlap")
  } else {
    ### determine the maximum moving step 
    MovingSize = WindowSize - overlap 
    MovingStep = floor((n.x - WindowSize)/MovingSize)
    ### define the velocity grid
    gr.v.x = (WindowSize + 1)/2 + MovingSize*(0:MovingStep)
    gr.v.y = (WindowSize + 1)/2 + MovingSize*(0:MovingStep)
    ### using correlation to track the velocity field
    ## defind a function to extract the image value within the moving window 
    MovingWindow = function(site, Image){
      output = matrix(0, WindowSize^2, nrow(site))
      for (i in 1:nrow(site)){
        x.s = site[i, 1] - (WindowSize-1)/2
        x.e = site[i, 1] + (WindowSize-1)/2
        y.s = site[i, 2] - (WindowSize-1)/2
        y.e = site[i, 2] + (WindowSize-1)/2
        output[,i] = c(Image[x.s:x.e, y.s:y.e])
      }
      return(output)
    }
    v.x = matrix(0, nrow(Image1), ncol(Image1))
    v.y = matrix(0, nrow(Image1), ncol(Image1))
    for (x in gr.v.x){
      for (y in gr.v.y){
        MovingWindow1 = MovingWindow(t(c(x,y)), Image1)
        if (sum((MovingWindow1 == rep(-32, WindowSize^2))*1) == WindowSize^2){
          v.x[x,y] = 0
          v.y[x,y] = 0
        } else{
          search.x.id = which( abs(gr.v.x-x) <= SearchSize )
          search.y.id = which( abs(gr.v.y-y) <= SearchSize )
          search.gr = expand.grid(x = gr.v.x[search.x.id], y = gr.v.y[search.y.id])
          MovingWindow2 = MovingWindow(search.gr, Image2)
          options(warn = - 1)  
          cor = cor( MovingWindow1, MovingWindow2)
          if (sum(is.na(cor)*1) == ncol(cor)){
            v.x[x,y] = 0
            v.y[x,y] = 0
          } else {
            v.id = which.max(cor)
            v.x[x,y] = (search.gr[v.id, 1] - x)/nrow(Image1)
            v.y[x,y] = (search.gr[v.id, 2] - y)/ncol(Image2) 
          }
        }
      }
    }
  }
  v.x[v.x == 0] = NA
  v.y[v.y == 0] = NA
  site = expand.grid(
    x = seq(0,1,length.out = nrow(Image1)),
    y = seq(0,1,length.out = ncol(Image1))
  )
  v = data.frame(
    x = site$x,
    y = site$y,
    v.x = c(v.x),
    v.y = c(v.y)
  )
  fit.vx = gam(v.x ~ s(x)+s(y), data = v, method = "REML")
  fit.vy = gam(v.y ~ s(x)+s(y), data = v, method = "REML")
  site.pred = data.frame(x = site$x, y = site$y)
  vs.x = predict(fit.vx, newdata = site.pred)
  vs.y = predict(fit.vy, newdata = site.pred)
  return(list(v.x = v.x, v.y = v.y, vs.x = vs.x, vs.y = vs.y))
}


# fit.test = COTREC(Image1, Image2, 17, 15, 4)
# site = expand.grid(
#   x = seq(0,1,length.out = nrow(Image1)),
#   y = seq(0,1,length.out = ncol(Image1))
# )
# v = data.frame(
#   x = site$x,
#   y = site$y, 
#   v.x = c(fit.test$v.x),
#   v.y = c(fit.test$v.y)
# )
# every_n <- function(x, by = by) {
#   x <- sort(x)
#   x[seq(1, length(x), by = by)]
# }
# keepx <- every_n(unique(v$x), by = 5)
# keepy <- every_n(unique(v$y), by = 5)
# vsub <- filter(v, x %in% keepx  &  y %in% keepy)
# ggplot(vsub, aes(x = x, y = y, fill = v.x)) +
#   geom_segment(aes(xend = x+v.x, yend = y+v.y), 
#                arrow = arrow(length = unit(0.05, "cm")), size = 0.15)

