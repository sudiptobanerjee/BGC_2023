rm(list=ls())

## create an adj matrix from a shapefile in R.
## We load two packages: "maps" and "maptools"
library(maps)
library(maptools)

## Read shapefile and create polygonal list of the map

mn.map.shp = rgdal::readOGR("minnesota.shp")

##Exporting to WinBUGS
sp2WB(mn.map.shp, "mn_bugs.txt")

##Use map library in R
mn.county = map(database="county", region="minnesota", fill=TRUE, plot=FALSE) ## with county boundaries

county.ID <- sapply(strsplit(mn.county$names, ","), function(x) x[2])
mn.poly = map2SpatialPolygons(mn.county, IDs=county.ID)

##Exporting to WinBUGS: How about trying
sp2WB(mn.poly, "mn_bugs_testing.txt")
##Unfortunately, the above does not work as the produced file "mn_bugs_testing.txt" is not in the right format.

##How about converting to a shapefile
library(raster)
shapefile(mn.poly, file="mn_map_r.shp", overwrite=TRUE)

mn.map.r.shp <- rgdal::readOGR("mn_map_r.shp")

##Exporting to WinBUGS
sp2WB(mn.map.r.shp, "mn_bugs_r_shp.txt")

##The output of "mn_bugs_r_shp.txt" will be approximately the same as "mn_bugs.txt" 

