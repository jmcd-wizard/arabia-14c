############################ Map the spatial intensity of radiocarbon dates #############################################

#load required libraries
library(spatstat)
library(tidyverse)
library(dplyr)
library(raster)
library(sf)
library(maps)

#load the dates and filter
mydates = read.csv("csv/arabia_dates_20.09.24.csv") %>% 
  filter(CRA <=5650 &CRA >= 2750) %>%
  filter(Country=="AE"|Country=="OM")

#turn the dataset into a spatial feature
mydates_sf = st_as_sf(mydates, coords = c("Longitude", "Latitude"))


#provide the EPSG code:
crs_wgs84 = st_crs(4326)
UTM_Zone_37N = st_crs(32637)

#attach the defined crs
st_crs(mydates_sf) = crs_wgs84

#attach the second EPSG code to the geometry of the dates
mydates_sf = st_transform(mydates_sf, UTM_Zone_37N)

#transform it into something usable
mydates_sf=st_as_sfc(mydates_sf)

#load the study area and spatial datasets
window<-st_read("shp/OP.shp")
countries<-st_read("shp/study_area.shp", layer = "study_area")

##tidy spatial wrangling##
##tidy spatial wrangling##
#provide the EPSG code:
crs_wgs84 = st_crs(4326)
#provide the EPSG code:
UTM_Zone_37N = st_crs(32637)
#add UTM_zone_37 to the CRS of the vector files (They are already in wgs84)
window = st_transform(window, UTM_Zone_37N)
countries = st_transform(countries, UTM_Zone_37N)



#Convert the window of analysis into an owin class object
window <- as.owin(window)
#check it has worked!
class(window)

#extract the coordinates from the dates
dates_coords = st_coordinates(mydates_sf)



# Generate a spatial intensity map of dates
sbw <- 50000 # Gaussian bandwdith (1 sd in metres)
cellres <- 1000 # raster cell size (in metres)
mysetppp <- ppp(x=dates_coords[,1],
                y=dates_coords[,2],
                window=as.owin(window))
# you will get the following message "Warning message:data contain duplicated points". Skip it.  
alldens <- density(mysetppp, sigma=sbw, eps=cellres, edge=FALSE)
alldenstc <- plot(alldens, main="", box=FALSE, do.plot=FALSE)


#give maximum long and lat for plotting
lon_min <- 1600000
lon_max <- 2800000
lat_min <- 1950000
lat_max <- 3000000
#long: east/west
#lat: north/south
#top left:(lon_min)1456213/(lat_max)3013154
#bottom right: (lon_max)2816631/(lat_min)1960104




#Plot The Kernal Density Estimate of radiocarbon dates
pdf(file="paper pdf/Figure_1.pdf", width=7.5, height=7)
par(bg = "#4287f5", mar = rep(0, 4))
plot(countries$geometry,col="burlywood3",border="white",xlim = c(lon_min,lon_max),ylim = c(lat_min,lat_max))
plot(alldens, main="", bbox=TRUE, add=TRUE)
plot(countries,col=rgb(128, 128, 128,alpha=0, maxColorValue=255),border="white",xlim = c(lon_min,lon_max),ylim = c(lat_min,lat_max), add=TRUE)
plot(mydates_sf,pch=21, cex=.5,alpha = .01,lwd=.5,bg="grey75", col="grey25",add=TRUE)
legend(1650000, 2125000, legend="Radiocarbon Dates",
       col="grey25",bg="grey75",pch=21, cex=.75, bty = "white")
#density legend
plot(alldenstc,vertical=TRUE, las=2, main="", xlim=c(1700000 ,1800000), ylim=c(2200000,2600000), add=TRUE, cex.axis=0.8, axis=FALSE)
text(x=1750000, y=2625000, labels="High", cex=0.8)
text(x=1750000, y=2175000, labels="Low", cex=0.8)
text(x=1875000, y=2400000, expression(paste(sigma," = 50 km")), cex=.9, col="black", font=1)
#draw the north arrow
xpos <- 2600000
ypos <- 2000000
lines(c(xpos,xpos),c(2090000,2000000),col="black")
polygon(c(xpos,xpos-(10000),xpos,xpos+(10000),xpos),c((ypos+70000),(ypos+70000),(ypos+90000),(ypos+70000),(ypos+70000)), col="black")
text(xpos, ypos+25000, "N", cex=0.9, col="black", font=2)
#draw the scale bar
polygon(c(2450000,2450000,2550000,2550000),c(2050000,2070000,2070000,2050000), col="NA", border="black")
text(2500000, y=2035000, "100 km", cex=0.9, col="black", font=1)
dev.off()
