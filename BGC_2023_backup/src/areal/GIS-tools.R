## Author: Sudipto Banerjee, January 07, 2014

rm(list=ls())
## create an adj matrix from a shapefile in R.
## you will need to install two packages: "maptools" and "spdep"

## Using maps in R
library(maps)
library(maptools)
library(spdep)

##Convert map to polygon
mn.county = map(database="county", region="minnesota", fill=TRUE, plot=FALSE) ## with county boundaries
## map("state","minnesota") ## without county boundaries
##mn.SL = map2SpatialLines(mn.map2)
##mn.poly = SpatialLines2PolySet(mn.SL)
##mn.poly = readShapePoly(mn.map2)
county.ID <- sapply(strsplit(mn.county$names, ","), function(x) x[2])
mn.poly = map2SpatialPolygons(mn.county, IDs=county.ID)


##Convert polygon to nb object
mn.nb = poly2nb(mn.poly)
mn.adj.mat = nb2mat(mn.nb, style="B")
## The option style="B" produces the binary adjacency matrix
## Write the 0-1 adjacency matrix 
##W <- mn.adj.mat
##W[(W>0)] <- 1

##Finding neighbors of a county
mn.region.id <- attr(mn.nb, "region.id")
winona.neighbors.index = mn.nb[[match("winona", mn.region.id)]]
winona.neighbors = rownames(mn.adj.mat[winona.neighbors.index,])

##Exporting to WinBUGS
sp2WB(mn.poly, "mn_bugs_testing.txt")

##Creating adjacency list for OpenBUGS
NumCells = length(mn.nb)
num = as.numeric(sapply(mn.nb, length))
mn.adj.list = as.numeric(unlist(mn.nb))
sumNumNeigh = as.numeric(length(unlist(mn.nb)))
mn.bugs <- list(num=num, adj=mn.adj.list, sumNumNeigh=sumNumNeigh)
dput(mn.bugs, "mn_bugs_adjacency.txt")

## Now create adjacency matrices from shapefiles
## Read shapefile and create polygonal list of the map
mn.map.shp = rgdal::readOGR("minnesota.shp")

## Convert the polygonal representation into a neighborhood list 
mn.nb.shp = poly2nb(mn.map.shp)
mn.adj.mat.shp = nb2mat(mn.nb.shp, style="B")

## Write the 0-1 adjacency matrix 
##W2 <- mn.adj.mat.shp
##W2[(W2>0)] <- 1

##Exporting to WinBUGS
sp2WB(mn.map.shp, "mn_bugs_shapefile.txt")


#Now notice that W and W2 are not in the same order!
adj.with.names.no.match = cbind(rownames(mn.adj.mat), as.character(mn.map.shp$NAME))

##o4W2 <-  order(as.character(mn.map$NAME))
ordered.index =  order(as.character(mn.map.shp$NAME))

##W2 <- W2[o4W2, o4W2]
mn.adj.mat.shp = mn.adj.mat.shp[ordered.index, ordered.index]

#now they are in the same order
##adj.with.names.match = cbind(as.character(mn.map.shp$NAME)[as.numeric(rownames(mn.map.shp))+1], rownames(mn.adj.mat))


#W and W2 are very similar: 99.86% 
sum(ifelse(as.numeric(mn.adj.mat)==as.numeric(mn.adj.mat.shp),1,0))/prod(dim(mn.adj.mat))

#The difference is due the polygons for each source are different! 
poly2nb(mn.poly) #It has 458 links
poly2nb(mn.map.shp)  #It has 448 links


