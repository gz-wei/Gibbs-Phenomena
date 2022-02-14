library(selextR)
library(rasterImage)

rm(list = ls()) 

SVY21 <- '+proj=tmerc +lat_0=1.366666666666667 +lon_0=103.8333333333333 +k=1 +x_0=28001.642 +y_0=38744.572 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
scan <- readRDS(here::here("data", "radar", "scan_list_19"))

radar.grd = SpatialPoints(coordinates(scan[[1]]$data),proj4string=CRS(SVY21)) # grid points in SVY21 projection system
Grid = coordinates(radar.grd) 
n.radar.grd = length(radar.grd) 
radar.z1 = scan[[1]]$data@data[,1]
grd.x = unique(coordinates(Grid)[,1])
grd.y = rev(unique(coordinates(Grid)[,2])) 
n.x = length(grd.x)
n.y = length(grd.y)
Image1 = matrix(radar.z1,nrow=n.y)[, n.y:1]
rasterImage2(z = Image1)
radar.z2 = scan[[6]]$data@data[,1]
Image2 = matrix(radar.z2,nrow=n.y)[, n.y:1]
rasterImage2(z = Image2)

TREC <- function(Image1, Image2, WindowSize, overlap) {
  WindowSize = 17
  overlap = 16
  if (WindowSize <= overlap |WindowSize%%2 == 0){
    stop("WindowSize must be an odd number and greater than overlap")
  } else {
    gr.v.s = (WindowSize + 1)/2
    gr.v.e = 1
    
  }
}