knitr::opts_chunk$set(comment = NA, tidy = TRUE)

rm(list=ls())
load("../data/CO-temp-data.RData")
ls()

library(spBayes)
library(MBA)
library(geoR)
library(raster)
library(leaflet)
library(sp)

blue.red <-  c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")

base.map <- leaflet(width="100%") %>%
    addProviderTiles("Esri.WorldImagery", group="Satellite") %>%
    addLayersControl(
        baseGroup = c("Satellite"),
        options = layersControlOptions(collapsed = FALSE)
    )

pal <- colorNumeric(blue.red, domain = temp)

base.map %>%
    addCircleMarkers(lng = coords[,1], lat = coords[,2], col = pal(temp), stroke = FALSE, radius = 5, fillOpacity = 0.9, popup=paste("Mean min temp:",round(temp,1))) %>%
    addLegend("bottomright", pal = pal, values = temp, opacity = 0.9, title = "Temperature C")
    

lm.obj <- lm(temp ~ lon + lat, data=coords)
summary(lm.obj)

##Promote coords to a sp package SpatialPoints object.
coordinates(coords) <- ~lon+lat
proj4string(coords) <- "+proj=longlat +datum=WGS84"
coords.utm <- spTransform(coords, CRS("+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
coords.utm <- coordinates(coords.utm)/1000

##Changes the sp object coords back to a matrix for later.
coords <- coordinates(coords)

d.max <- max(iDist(coords.utm))
d.max

v.temp <- variog(coords=coords.utm, data=temp, uvec=(seq(0, 0.75*d.max, length=20)))

v.resid <- variog(coords=coords.utm, data=resid(lm.obj), uvec=(seq(0, 0.75*d.max, length=20)))

par(mfrow=c(1,2))
plot(v.temp, xlab="Distance (km)")
plot(v.resid, xlab="Distance (km)")

resid.surf <- mba.surf(cbind(coords, resid(lm.obj)), no.X=200, no.Y=200, extend=TRUE, sp=TRUE)$xyz.est

proj4string(resid.surf) <- "+proj=longlat +datum=WGS84"

resid.surf <- raster(resid.surf)

pal <- colorNumeric(blue.red, values(resid.surf), na.color = "transparent")

base.map %>%
    addRasterImage(resid.surf, colors = pal, opacity = 0.75, group="Regression residuals") %>%
    addLegend("bottomright", pal = pal, values = values(resid.surf), opacity = 0.75, title = "<center>Regression<br> residuals</center>") %>%
    addLayersControl(
        baseGroup = c("Satellite"),
        overlayGroups = c("Regression residuals"),
        options = layersControlOptions(collapsed = FALSE)
    )

n.samples <- 5000

n <- length(temp)
p <- 3 ##Intercept and spatial locations

phi <- 3/200
alpha <- 1/12

beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

##Assuming an inverse Gamma prior on sigma^2
sigma.sq.prior.shape <- 2.0
sigma.sq.prior.rate <- 1/12 ##1/scale

m.exact  <- bayesGeostatExact(temp~coords, coords=coords.utm,
                          beta.prior.mean=beta.prior.mean,
                          beta.prior.precision=beta.prior.precision,
                          phi=phi, alpha=alpha,
                          sigma.sq.prior.shape=sigma.sq.prior.shape,
                          sigma.sq.prior.rate=sigma.sq.prior.rate, n.samples=n.samples)
names(m.exact)

plot(m.exact$p.samples, density=FALSE)



cov.model <- "exponential"

starting <- list("phi"=3/(0.5*d.max), "sigma.sq"=5, "tau.sq"=5)

tuning <- list("phi"=0.5, "sigma.sq"=0.05, "tau.sq"=0.05)

priors <- list("beta.Flat", "phi.Unif"=c(3/d.max, 3/(0.01*d.max)), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 5))
 
m.1 <- spLM(temp ~ coords, coords=coords.utm, starting=starting,
            tuning=tuning, priors=priors, cov.model=cov.model,
            n.samples=n.samples, n.report=2000)

plot(m.1$p.theta.samples)

burn.in <- floor(0.75*n.samples)

m.1 <- spRecover(m.1, start=burn.in, thin=2)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],3)

round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],3)

quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

##random effects from bayesGeostatExact
w.exact.summary <- apply(m.exact$sp.effects, 1, quants)

w.exact.surf <- mba.surf(cbind(coords, w.exact.summary[1,]), no.X=200, no.Y=200, extend=TRUE, sp=TRUE)$xyz.est
proj4string(w.exact.surf) <- "+proj=longlat +datum=WGS84"
w.exact.surf <- raster(w.exact.surf)

##random effects from spLM
w.m1.summary <- apply(m.1$p.w.recover.samples, 1, quants)

w.m1.surf <- mba.surf(cbind(coords, w.m1.summary[1,]), no.X=200, no.Y=200, extend=TRUE, sp=TRUE)$xyz.est
proj4string(w.m1.surf) <- "+proj=longlat +datum=WGS84"
w.m1.surf <- raster(w.m1.surf)

##make the color ramp
all.values <- c(values(w.exact.surf), values(w.m1.surf), values(resid.surf))
pal <- colorNumeric(blue.red, all.values, na.color = "transparent")

base.map %>%
    addRasterImage(w.exact.surf, colors = pal, opacity = 0.75, group="bayesGeostatExact") %>%
    addRasterImage(w.m1.surf, colors = pal, opacity = 0.75, group="spLM m.1") %>%
    addRasterImage(resid.surf, colors = pal, opacity = 0.75, group="Regression residuals") %>%
    addLegend("bottomright", pal = pal, values = all.values, opacity = 0.75, title = "<center>Random effects &<br>Regression<br>residuals</center>") %>%
    addLayersControl(
        baseGroup = c("Satellite"),
        overlayGroups = c("bayesGeostatExact", "spLM m.1", "Regression residuals"),
        options = layersControlOptions(collapsed = FALSE)
    ) %>% hideGroup(c("spLM m.1", "Regression residuals"))


lm.obj <- lm(temp ~ elev + lon + lat, data=data.frame(coords))
summary(lm.obj)

v <- variog(coords=coords.utm, data=resid(lm.obj), uvec=(seq(0, 0.75*d.max, length=20)))
plot(v, xlab="Distance (km)")

priors$sigma.sq.IG <- c(2, 1)
priors$tau.sq.IG <- c(2, 1)
 
m.2 <- spLM(temp ~ elev + coords, coords=coords.utm, starting=starting,
            tuning=tuning, priors=priors, cov.model=cov.model,
            n.samples=n.samples, n.report=2000)

m.2 <- spRecover(m.2, start=burn.in, thin=2)

round(summary(m.2$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.2$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

w.2.summary <- apply(m.2$p.w.recover.samples, 1, quants)

w.m2.surf <- mba.surf(cbind(coords, w.2.summary[1,]), no.X=200, no.Y=200, extend=TRUE, sp=TRUE)$xyz.est
proj4string(w.m2.surf) <- "+proj=longlat +datum=WGS84"
w.m2.surf <- raster(w.m2.surf)

pal.2 <- colorNumeric(blue.red, values(w.m2.surf), na.color = "transparent")

base.map %>%
    addRasterImage(w.m2.surf, colors = pal.2, opacity = 0.75, group="spLM m.2") %>%
    addRasterImage(w.m1.surf, colors = pal, opacity = 0.75, group="spLM m.1") %>%
    addRasterImage(resid.surf, colors = pal, opacity = 0.75, group="Regression residuals") %>%
    addLegend("bottomleft", pal = pal.2, values = values(w.m2.surf), opacity = 0.75, title = "<center>Random effects (m.2) &<br>Regression<br>residuals</center>") %>%
    addLegend("bottomright", pal = pal, values = all.values, opacity = 0.75, title = "<center>Random effects (m.1) &<br>Regression<br>residuals</center>") %>%
    addLayersControl(
        baseGroup = c("Satellite"),
        overlayGroups = c("spLM m.1", "spLM m.2", "Regression residuals"),
        options = layersControlOptions(collapsed = FALSE)
    ) %>% hideGroup(c("spLM m.1", "Regression residuals"))

pred.X <- cbind(1, pred.elev, pred.coords)

m.2.pred <- spPredict(m.2, pred.covars=pred.X, pred.coords=pred.coords, start=burn.in, thin=20, n.report=10)

y.p.summary <- apply(m.2.pred$p.y.predictive.samples, 1, quants)

y.p.width <- y.p.summary[3,]-y.p.summary[2,]

y.median.surf <- rasterFromXYZ(cbind(pred.coords, y.p.summary[1,]))
proj4string(y.median.surf) <- "+proj=longlat +datum=WGS84"
pal.1 <- colorNumeric(blue.red, c(values(y.median.surf), temp), na.color = "transparent")

y.width.surf <- rasterFromXYZ(cbind(pred.coords, y.p.width)) 
proj4string(y.width.surf) <- "+proj=longlat +datum=WGS84"
pal.2 <- colorNumeric(blue.red, values(y.width.surf), na.color = "transparent")

base.map %>%
    addRasterImage(y.median.surf, colors = pal.1, opacity = 0.9, group="Predicted temperature") %>%
    addRasterImage(y.width.surf, colors = pal.2, opacity = 0.9, group="Prediction 95% CI width") %>%
    addCircleMarkers(lng = coords[,1], lat = coords[,2], col = pal.1(temp), fillOpacity = 0.9, stroke = FALSE,
                     radius = 3, popup=paste("Mean min temp:",round(temp,1)), group="Stations") %>%
      addLegend("bottomleft", pal = pal.1, values = values(y.median.surf), opacity = 0.9, title = "<center>Predicted<br>temperature C</center>") %>%
    addLegend("bottomright", pal = pal.2, values = values(y.width.surf), opacity = 0.9, title = "<center>Prediction<br>95% CI width</center>") %>%
    addLayersControl(
        baseGroup = c("Satellite"),
        overlayGroups = c("Predicted temperature", "Prediction 95% CI width", "Stations"),
        options = layersControlOptions(collapsed = FALSE)
    ) %>% hideGroup(c("Prediction 95% CI width", "Stations"))

lm.ref <- bayesLMRef(lm.obj, 500)

spDiag(lm.ref)$DIC

spDiag(m.2)$DIC
