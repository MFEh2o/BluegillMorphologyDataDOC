# getting catchment information based on lat-long

rm(list=ls())

library(rgdal)
library(rgeos)
library(sf)

# make coordinates into sf object
setwd("~/Documents/Research/People/Students/current/Bishop_Chelsea/lakeWatersheds/")
d=read.csv("Lake_Info_2020.csv",header=TRUE,stringsAsFactors=FALSE)
dsf=st_as_sf(d[,3:4],coords=c("long","lat"),crs=4269)

d$basin=""

# check whether each lake point is in NHD region 04 or 07
#region 04 is the UP and drains into the great lakes
#region 07 is Wisconsin and drains into the Gulf of Mexico
setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusDownloads/")

#check region 04
setwd("NHDPlus04/NHDSnapshot/Hydrography")
#load waterbody shape file
shape<-readOGR(dsn=".",layer="NHDWaterbody")

# makes sf version of shapefile I think...
z=st_as_sf(shape)

# find intersection of lake polygons and coordinate
int04=st_intersects(dsf,z)

setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusDownloads/")

# repeat  process for region 07
setwd("NHDPlus07/NHDSnapshot/Hydrography")
#load waterbody shape file
shape<-readOGR(dsn=".",layer="NHDWaterbody")

# makes sf version of shapefile I think...
z=st_as_sf(shape)

# find intersection of lake polygons and coordinate
int07=st_intersects(dsf,z)

in04=as.data.frame(int04)
in07=as.data.frame(int07)

d$basin[in04$row.id]="04"
d$basin[in07$row.id]="07"

# hummingbird may not be in NHD or the point is missing the polygon
# hummingbird drains into Bay so using its huc for hummingbird
d$basin[5]="04"

d
setwd("~/Documents/Research/People/Students/current/Bishop_Chelsea/lakeWatersheds/")
write.csv(d,"Lake_Info_2020wBasins.csv",row.names=FALSE)
