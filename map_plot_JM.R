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
library(devtools)
gpclibPermit()


mat <- read.table("~/Downloads/symbols_colors.txt",row.names=1,stringsAsFactors=FALSE)

sites <- rbind(
	BI=c(41.184326, -71.574127),
	BP=c(41.2, -73.181154),
	ER=c(36.807026, -76.290405),
	F=c(40.9, -73.139791),
	KC=c(37.3016, -76.4226),
	NBH=c(41.637174, -70.914284),
	NYC=c(40.7006, -74.1223),
	SH=c(40.40, -74.0113)
	)
colnames(sites) <- c("lat","lon")
myloc <- c(-77.25,36,-69.75,42) #left bottom right top
mm <- get_map(myloc,zoom=6,source = "stamen",maptype="terrain-background")
box <- as(extent(as.numeric(attr(mm, 'bb'))[c(2,4,1,3)] + c(.001,-.001,.001,-.001)), "SpatialPolygons")
shp <- readOGR(dsn="~/Dropbox/QTL_Paper/METADATA/cb_2015_us_state_500k",layer="cb_2015_us_state_500k")
shp <- rgdal::readOGR("~/Dropbox/QTL_Paper/METADATA/cb_2015_us_state_500k/cb_2015_us_state_500k.shp")
proj4string(box) <- CRS(summary(shp)[[4]])
tractSub <- gIntersection(shp, box, byid = TRUE, 
                          id = as.character(shp$GEOID))
sites <- sites[c(1,2,3,6,7),]
mat <- mat[rownames(sites),]

png('~/Dropbox/QTL_Paper/FIGURES/Rough Figures/map_only.png')
par(mar=c(0,0,0,0),oma=c(1,1,1,1))
plot(tractSub,col=rgb(0,0,0,.05))
points(sites[,2],sites[,1],pch=c(21),cex=5,col=mat[,5],bg=c('white',mat[,5][2:5]), lwd=5)
text(sites[,2],sites[,1],font=2,c("BLI","BRP","ELR","NBH","NEW"),col=c('black',rep('white',times=4)))
dev.off()


text(sites[,2],sites[,1],font=2,c("BLI","BRP","ELR","NBH","NEW"), cex=0.8)


#qmap(c(lon= -74.632951,lat=39.433438),zoom=6,source = "google",maptype="satellite")

#mm <- get_map(c(lon= -74.632951,lat=39.433438),zoom=6,source = "google",maptype="satellite")
#mm <- get_map(c(lon= -74.632951,lat=39.433438),zoom=6,source = "stamen")


#plotData <- left_join(tractSub, data, by = "id")
#tractSub <- gIntersection(shp, box, byid = TRUE, 
#                        id = as.character(shp$GEOID))
#tractSub <- fortify(tractSub, region = "GEOID")
#eastcoast <- ggmap(mm)
#eastcoast + 
#	geom_polygon(aes(x = long, y = lat, group = group), data = tractSub, colour = 'white', fill = 'black', alpha = .4, size = .3) +
#	geom_point(aes(x = lon, y = lat), data = as.data.frame(sites),color=mat[,5], size=9) +
#	geom_text(aes(label= c("S1","T2","T4","S2","S4","T1","T3","S3"), x = lon, y = lat), data = as.data.frame(sites))





