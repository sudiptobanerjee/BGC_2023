## Author: Sudipto Banerjee, January 17, 2014.

rm(list=ls())

library(spdep)
library(maptools)
library(RColorBrewer)

##columbus.poly <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1], proj4string=CRS("+proj=longlat +ellps=clrk66"))
columbus.poly <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1])

##There are two different ways to create a neighbor object for columbus
##Method 1a: Use poly2nb with queen type neighbors
columbus.queen.nb = poly2nb(columbus.poly, queen=TRUE)
##Method 1b: Use poly2nb with rook type neighbors
columbus.rook.nb = poly2nb(columbus.poly, queen =FALSE)
##Method 2: Use the read.gal file to read a GAL file 
columbus.gal.nb <- read.gal(system.file("etc/weights/columbus.gal", package="spdep")[1])

##Compute Moran's I and Geary's C
##We make two distinct listw objects for the rook and gal neighbors
columbus.gal.listw = nb2listw(columbus.gal.nb, style="B", zero.policy=TRUE)
columbus.rook.listw = nb2listw(columbus.rook.nb, style="B", zero.policy=TRUE)
##We now compute Moran's I and Geary's C for each of these neighbors
columbus.gal.moran.out = moran.test(columbus.poly$CRIME, listw=columbus.gal.listw, zero.policy=TRUE)
columbus.rook.moran.out = moran.test(columbus.poly$CRIME, listw=columbus.rook.listw, zero.policy=TRUE)
columbus.gal.geary.out = geary.test(columbus.poly$CRIME, listw=columbus.gal.listw, zero.policy=TRUE)
columbus.rook.geary.out = geary.test(columbus.poly$CRIME, listw=columbus.rook.listw, zero.policy=TRUE)

##SAR model regressing HOVAL+INC
columbus.gal.sar.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="SAR", listw=columbus.gal.listw, zero.policy=TRUE)
columbus.gal.sar.fitted = fitted(columbus.gal.sar.out)
columbus.poly$fitted.gal.sar = columbus.gal.sar.fitted

columbus.rook.sar.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="SAR", listw=columbus.rook.listw, zero.policy=TRUE)
columbus.rook.sar.fitted = fitted(columbus.rook.sar.out)
columbus.poly$fitted.rook.sar = columbus.rook.sar.fitted

##CAR model regressing CRIME on HOVAL + INCOME
columbus.car.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="CAR", listw=columbus.rook.listw, zero.policy=TRUE)
columbus.car.fitted = fitted(columbus.car.out)
columbus.poly$fitted.car = columbus.car.fitted

##Distance based neighbors in spdep
columbus.coords = coordinates(columbus.poly)
columbus.knn = knearneigh(columbus.coords)
columbus.knn2nb = knn2nb(columbus.knn) 
columbus.dist.list = nbdists(columbus.knn2nb, columbus.coords)
columbus.dist.vec = unlist(columbus.dist.list)
columbus.dist.max = max(columbus.dist.vec)
columbus.dnn.nb = dnearneigh(columbus.coords, 0, 0.25*columbus.dist.max)

##Form a listw object using the distance-based nearest neighbors 
columbus.dnn.listw = nb2listw(columbus.dnn.nb, style="B", zero.policy=TRUE)

##SAR model regressing HOUSE_VAL+INCOME using distance-based nearest neighbors
columbus.dnn.sar.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="SAR", listw=columbus.dnn.listw, zero.policy=TRUE)
columbus.dnn.sar.fitted = fitted(columbus.dnn.sar.out)
columbus.poly$fitted.dnn.sar = columbus.dnn.sar.fitted

##CAR model regressing HOUSE_VAL+INCOME using distance-based nearest neighbors
columbus.dnn.car.out = spautolm(CRIME~HOVAL+INC, data=columbus.poly, family="CAR", listw=columbus.dnn.listw, zero.policy=TRUE)
columbus.dnn.car.fitted = fitted(columbus.dnn.car.out)
columbus.poly$fitted.dnn.car = columbus.dnn.car.fitted

##Draw the maps using the spplot (trellis) graphics function
##postscript(file="nc_sids_sar_actual.eps")
##print(spplot(nc.sids, "rates.FT", at = brks, col.regions = rev(brewer.pal(4,"RdBu")), main="a) Actual SIDS rates"))
##dev.off()
##postscript(file="nc_sids_sar_fitted.eps")
##print(spplot(nc.sids, "fitted.sar", at = brks, col.regions = rev(brewer.pal(4,"RdBu")), main="b) Fitted SIDS rates from SAR model"))
##dev.off()

##detach(package:spdep)
##detach(package:maptools)
##detach(package:RColorBrewer)
