#' ## Colorado spring temperature data
#' These data were originally part of the `fields` package's `COmonthlyMet` dataset.
## ---- download_data------------------------------------------------------
rm(list=ls())
library(downloader)
load("CO-temp-data.RData")
ls()

#' 
#' Our goal is to create a complete prediction surface of minimum spring temperature with associated estimates of uncertainty. 
#' 
#' We begin by loading the necessary packages.
#' 
## ---- laod_packages, message=FALSE---------------------------------------
library(gstat)
library(raster)
library(leaflet)
library(sp)

#' 
#' Next, set up a `leaflet` basemap to help visualize the data and model output. We'll make heavy use of the pipe operator `%>%` to reduce clutter.
#' 
## ---- leaflet_basemap----------------------------------------------------
blue.red <-  c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")

base.map <- leaflet() %>%
    addProviderTiles("Stamen.Terrain", group="Terrain") %>%
    addProviderTiles("Esri.WorldImagery", group="Satellite") %>%
    addLayersControl(
        baseGroup = c("Terrain", "Satellite"),
        options = layersControlOptions(collapsed = FALSE)
    )

#' 
#' Take a look at station locations and mean minimum spring temperatures across Colorado. This code below produces a clickable dynamic map.
## ---- map_stations-------------------------------------------------------
pal <- colorNumeric(blue.red, domain = temp)

base.map %>%
    addCircleMarkers(lng = coords[,1], lat = coords[,2], col = pal(temp), stroke = FALSE, radius = 5, fillOpacity = 0.9, popup=paste("Mean min temp:",round(temp,1))) %>%
    addLegend("bottomright", pal = pal, values = temp, opacity = 0.9, title = "Temperature C")
    

#' 
#' ### Fit a non-spatial regression
#' Consider the non-spatial regression and a little exploratory data analysis (EDA).
#' 
## ---- fit_lm-------------------------------------------------------------
lm.obj <- lm(temp ~ lon + lat, data=coords)
summary(lm.obj)

#' 
#' We'll reproject the geographic coordinates (i.e., longitude and latitude) to a coordinate system that provides a bit more intuition about distance and reduces spatial distortion. Here we selected Universal Transverse Mercator (UTM) coordinate system. To improve interpretation we convert the UTM distance units from meters to kilometers.
#' 
## ---- reproject_coords---------------------------------------------------
##Promote coords to a sp package SpatialPoints object.
coordinates(coords) <- ~lon+lat
proj4string(coords) <- "+proj=longlat +datum=WGS84"
coords.utm <- spTransform(coords, CRS("+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
coords.utm <- coordinates(coords.utm)/1000

##Changes the sp object coords back to a matrix for later.
coords <- coordinates(coords)

#' 
#' Next let's take a look at the regression model residuals, assess their spatial independence, and start thinking about variogram and covaraince parameters.
## ---- plot_variograms, fig.align="center", fig.width=10, message=FALSE----
d.max <- max(pointDistance(coords.utm, lonlat=FALSE, allpairs=TRUE))
d.max

v.temp <- variog(coords=coords.utm, data=temp, uvec=(seq(0, 0.75*d.max, length=20)))

v.resid <- variog(coords=coords.utm, data=resid(lm.obj), uvec=(seq(0, 0.75*d.max, length=20)))

par(mfrow=c(1,2))
plot(v.temp, xlab="Distance (km)")
plot(v.resid, xlab="Distance (km)")

#' 
#' It is also very helpful to create an interpolated surface of the model residuals to further assess spatial structure and potentially identify missing covariates.
#' 
## ---- map_residuals------------------------------------------------------
resid.surf <- idw(cbind(coords), resid(lm.obj))

proj4string(resid.surf) <- "+proj=longlat +datum=WGS84"

resid.surf <- raster(resid.surf)

pal <- colorNumeric(blue.red, values(resid.surf), na.color = "transparent")

base.map %>%
    addRasterImage(resid.surf, colors = pal, opacity = 0.75, group="Regression residuals") %>%
    addLegend("bottomright", pal = pal, values = values(resid.surf), opacity = 0.75, title = "<center>Regression<br> residuals</center>") %>%
    addLayersControl(
        baseGroup = c("Terrain", "Satellite"),
        overlayGroups = c("Regression residuals"),
        options = layersControlOptions(collapsed = FALSE)
    )

#' 
#' ### Fit some spatial regression models
#' 
## ---- fig.align="center"-------------------------------------------------
#' 
## ---- fit_lm_with_elev_and_variogram, fig.align="center"-----------------
lm.obj <- lm(temp ~ elev + lon + lat, data=data.frame(coords))
summary(lm.obj)

v <- variog(coords=coords.utm, data=resid(lm.obj), uvec=(seq(0, 0.75*d.max, length=20)))
