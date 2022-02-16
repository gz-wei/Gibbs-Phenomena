#library(sp)
#library(rgdal)
#library(rgeos)
#library(gstat)
library(RColorBrewer)
#library(spacetime)
library(selextR)
library(tidyverse)
#library(timeSeries)
#library(fields)
#library(MASS)
#library(mvtnorm)
#library(matrixcalc)
library(rasterImage)


rm(list = ls()) 

SVY21 <- '+proj=tmerc +lat_0=1.366666666666667 +lon_0=103.8333333333333 +k=1 +x_0=28001.642 +y_0=38744.572 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
frequency <- 5 * 60 # Radar scan frequency

### Plotting information
singList <- list("sp.polygons", sgBd)
myList <- list("sp.polygons", myBd)
inList <- list("sp.polygons", inBd)
sp.layout.list <- list(singList, myList, inList)
# png(filename = here::here("data", "Singapore.png"),
#     width = 4, height =4, unit="in", pointsize = 12,
#     bg = "white", res = 600)
par(mar=c(0.1,0.1,0.1,0.1))
plot(sgBd, axes=FALSE, xlim=c(-20000,80000), ylim=c(-10000,90000), col="blue")
plot(myBd, add=TRUE, col="grey")
plot(inBd, add=TRUE, col="grey")
x_coord <- c(2000,  58000, 58000,2000, 2000)
y_coord <- c(20000, 20000,51000, 51000, 20000)
xym <- cbind(x_coord, y_coord)
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps,add=TRUE,lwd=3,border="red")
# dev.off()

### sketch the plot
colPalette <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'YlOrRd'))(100), .85)
colPalette2 <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Blues'))(100), 0.99)
colPalette1 <- rev(rainbow(32, start=5/6, end=3/6, alpha=0.8))
at.points <- seq(from=-35, to=75, length=33)
scan <- readRDS(here::here("data", "radar", "scan_list_19"))
plot(scan[[1]], boundaries=TRUE)

# extract the information of radar images
scan <- readRDS(here::here("data", "radar", "scan_list_19"))
radar.grd = SpatialPoints(coordinates(scan[[1]]$data),proj4string=CRS(SVY21)) # grid points in SVY21 projection system
Grid = coordinates(radar.grd) 
n.radar.grd = length(radar.grd) 
radar.z = scan[[1]]$data@data[,1]
rainfall = (10^(radar.z/10)/200)^(5/8)
grd.x = unique(coordinates(Grid)[,1])
grd.y = rev(unique(coordinates(Grid)[,2])) 
n.x = length(grd.x)
n.y = length(grd.y)
rasterImage2(
  z = matrix(radar.z,nrow=n.y)[, n.y:1],
  #zlim = c(0, 120), 
  x = grd.x,
  y = grd.y
)

# animation id = 10, 19, 22, 30, 34, 43, 71, #318, 327, 334
scan <- readRDS(here::here("data", "radar", "scan_list_34"))
time.id <- seq(1, length(scan), 5)
for (i in time.id){
  radar.z = scan[[i]]$data@data[,1]
  rainfall = (10^(radar.z/10)/200)^(5/8)
  rasterImage2(
    # z = matrix(radar.z,nrow=n.y)[, n.y:1],
    z = matrix(rainfall,nrow=n.y)[, n.y:1],
    zlim = c(0, 120),
    x = grd.x,
    y = grd.y
  )
  Sys.sleep(0.1)
}

##  Read in the data for smaller region  
scan <- readRDS(here::here("data", "radar", "scan_list_34"))
time.id <- seq(1, length(scan), 5)
for (i in time.id){
  radar.z = scan[[i]]$data@data[,1]
  rainfall = (10^(radar.z/10)/200)^(5/8)
  rasterImage2(
    z = matrix(rainfall,nrow=n.y)[, n.y:1][281:380, 301:400],
    #z = matrix(rainfall,nrow=n.y)[, n.y:1][281:380, 291:390],
    zlim = c(0, 120)
  )
  Sys.sleep(0.1)
}

## save the selected dataset 
data.rf <- list()
time.select <- time.id[8:20]
for (i in 1:length(time.select)){
  scan.id = time.select[i]
  radar.z = scan[[scan.id]]$data@data[,1]
  rainfall = (10^(radar.z/10)/200)^(5/8)
  tempt = matrix(rainfall,nrow=n.y)[, n.y:1][281:380, 301:400]
  #tempt = matrix(rainfall,nrow=n.y)[, n.y:1][281:380, 291:390]
  print(range(tempt))
  data.rf[[i]] <- tempt
  rasterImage2(
    z = tempt, 
    zlim = c(0, 120)
  )
  Sys.sleep(0.1)
}
save(data.rf, file = here::here("data", "case", "data.rf.RData"))

data.rd <- list()
time.select <- time.id[8:20]
for (i in 1:length(time.select)){
  scan.id = time.select[i]
  radar.z = scan[[scan.id]]$data@data[,1]
  tempt = matrix(radar.z,nrow=n.y)[, n.y:1][281:380, 301:400]
  print(range(tempt))
  data.rd[[i]] <- tempt
  rasterImage2(
    z = tempt
  )
  Sys.sleep(0.1)
}
save(data.rd, file = here::here("data", "case", "data.rd.RData"))
