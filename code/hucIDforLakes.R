# getting catchment information based on lat-long
# Script by Stuart E Jones
# Edited for workflow and file paths by Kaija Gahm, January 2021

# Load packages -----------------------------------------------------------
library(rgdal)
library(rgeos)
library(sf)
library(here)
library(dplyr)

# make coordinates into sf object
setwd("~/Documents/Research/People/Students/current/Bishop_Chelsea/lakeWatersheds/")
d <- read.csv(here("data", "inputs", "lakeInfo.csv"), 
              header = T, stringsAsFactors = F)
dsf <- st_as_sf(d[,c("lat", "long")], 
                coords = c("long", "lat"), crs = 4269)
d$basin <- ""


# Assign lakes to regions -------------------------------------------------
# check whether each lake point is in NHD region 04 or 07
#region 04 is the UP and drains into the great lakes
#region 07 is Wisconsin and drains into the Gulf of Mexico

#check region 04
#load waterbody shape file
shape <- readOGR(dsn = here("data", "inputs", "region4", "NHDSnapshot", "Hydrography"), layer = "NHDWaterbody")

# makes sf version of shapefile I think...
z=st_as_sf(shape)

# find intersection of lake polygons and coordinate
int04=st_intersects(dsf,z)

# repeat  process for region 07
#load waterbody shape file
shape<-readOGR(dsn=here("data", "inputs", "region7", "NHDSnapshot", "Hydrography"),layer="NHDWaterbody")

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
write.csv(d, here("data", "outputs", "Lake_Info_2020wBasins.csv"),row.names=FALSE)
