library(maps)
library(sf)
library(sp)
library(spdep)

# Get county map for California
ca.county = map("county","california", fill=TRUE, plot=FALSE)
county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.county.sf <- st_as_sf(ca.county)

sf_use_s2(FALSE)
ca.coords = st_coordinates(st_centroid(ca.county.sf))
plot(st_geometry(ca.county.sf), main="County map of California")
text(ca.coords, county.ID, cex=0.5, col = "red")

ca.nb = poly2nb(ca.county.sf)
ca.region.id = attr(ca.nb, "region.id")
ca.region.id = sapply(strsplit(ca.region.id, ","), function(x) x[2])
alameda.neighbors.index = ca.nb[[match("alameda", ca.region.id)]]
ca.adj.mat = nb2mat(ca.nb, style="B")
alameda.neighbors = rownames(ca.adj.mat[alameda.neighbors.index,])
alameda.neighbors = sapply(strsplit(alameda.neighbors, ","), function(x) x[2])


