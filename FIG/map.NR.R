#!/bin/R
### Map plots from N reduce2grid
library(ggmap)
library(maptools)
library(gpclib)
library(sp)
library(raster)
library(rgdal)
library(dplyr)
library(Cairo)
library(scales)
library(rgeos)

gpclibPermit()


mat <- read.table("symbols_colors.txt", row.names = 1, stringsAsFactors = FALSE)

sites <- rbind(BI = c(41.184326, -71.574127), BP = c(41.2, -73.181154), ER = c(36.807026, 
  -76.290405), F = c(40.9, -73.139791), KC = c(37.3016, -76.4226), NBH = c(41.637174, 
  -70.914284), NYC = c(40.7006, -74.1223), SH = c(40.4, -74.0113))
colnames(sites) <- c("lat", "lon")

# qmap(c(lon= -74.632951,lat=39.433438),zoom=6,source =
# 'google',maptype='satellite')

# mm <- get_map(c(lon= -74.632951,lat=39.433438),zoom=6,source =
# 'google',maptype='satellite')
mm <- ggmap::get_map(c(lon = -74.632951, lat = 39.433438), zoom = 6, source = "stamen")
myloc <- c(-77.25, 36, -69.75, 42)  #left bottom right top
mm <- ggmap::get_map(myloc, zoom = 6, source = "stamen", maptype = "terrain-background")
mm <- ggmap::get_map(myloc, source = "stamen")

box <- as(extent(as.numeric(attr(mm, "bb"))[c(2, 4, 1, 3)] + c(0.001, -0.001, 0.001, 
  -0.001)), "SpatialPolygons")

proj4string(box) <- CRS(summary(shp)[[4]])
shp <- rgdal::readOGR("~/cb_2015_us_state_500k/cb_2017_us_state_500k.shp")
tractSub <- gIntersection(shp, box, byid = TRUE, id = as.character(shp$GEOID))
tractSub <- fortify(tractSub, region = "GEOID")
# plotData <- left_join(tractSub, data, by = 'id')

eastcoast <- ggmap(mm)
eastcoast + geom_polygon(aes(x = long, y = lat, group = group), data = tractSub, 
  colour = "white", fill = "black", alpha = 0.4, size = 0.3) + geom_point(aes(x = lon, 
  y = lat), data = as.data.frame(sites), color = mat[, 5], size = 9) + geom_text(aes(label = c("S1", 
  "T2", "T4", "S2", "S4", "T1", "T3", "S3"), x = lon, y = lat), data = as.data.frame(sites))



shp <- readOGR(dsn = "/Users/noahreid/Downloads/cb_2015_us_state_500k", layer = "cb_2015_us_state_500k")
tractSub <- gIntersection(shp, box, byid = TRUE, id = as.character(shp$GEOID))
par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))
plot(tractSub, col = rgb(0, 0, 0, 0.05))
points(sites[, 2], sites[, 1], pch = 20, cex = 5, col = mat[, 5])
text(sites[, 2], sites[, 1], c("S1", "T2", "T4", "S2", "S4", "T1", "T3", "S3"))
