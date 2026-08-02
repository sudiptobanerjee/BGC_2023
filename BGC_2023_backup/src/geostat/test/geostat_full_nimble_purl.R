knitr::opts_chunk$set(comment = NA, tidy = TRUE)

rm(list=ls())
load("../data/CO-temp-data.RData")
ls()

library(fields)
library(MBA)
library(geoR)
library(sf)
library(leaflet)

blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")

base.map <- leaflet(width = "100%") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addLayersControl(baseGroup = c("Satellite"), options = layersControlOptions(collapsed = FALSE))

pal <- colorNumeric(blue.red, domain = temp)

base.map %>%
    addCircleMarkers(lng = coords[, 1], lat = coords[, 2], col = pal(temp), stroke = FALSE,
        radius = 5, fillOpacity = 0.9, popup = paste("Mean min temp:", round(temp,
            1))) %>%
    addLegend("bottomright", pal = pal, values = temp, opacity = 0.9, title = "Temperature C")

lon2UTM <- function(longitude) {
    (floor((longitude + 180)/6) %% 60) + 1
} ##Function to find the UTM zone of a geographic coordinate in the western hemisphere (only longitude is needed).

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
} ##Function to calculate the mode of the above data to find which zone contains most of the points
utmZone <- Mode(lon2UTM(coords$lon))
utmZone

coordsGeog <- st_as_sf(coords, coords=c("lon","lat"), crs="+proj=longlat +datum=WGS84 +no_defs")
coordsUTM <- st_transform(coordsGeog, crs="+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
coordsUTM <- st_coordinates(coordsUTM)/1000
class(coordsUTM)

elev <- elev/1000
df <- data.frame(temp, elev, coordsUTM)
names(df) <- c("temp", "elev", "xUTM", "yUTM")

library("nimble")

expCov <- nimbleFunction(     
	run = function(dists = double(2), phi=double(0), sigma=double(0)) {
		returnType(double(2))
		sigma2 <- sigma*sigma
		result <- sigma2*exp(-phi*dists)
		return(result)
	}
)
#cExpCov <- compileNimble(expCov) ##Uncomment if you wish to use compiled code in C

geostatCodeFull <- nimbleCode({
        for (i in 1:N) {
                Y[i] ~ dnorm(mu[i], tausq)
                mu[i] <- inprod(X[i,1:p], beta[1:p]) + w[i]
            yFit[i] ~ dnorm(mu[i], tausq)
            yResid[i] <- Y[i] - mu[i]
    }
        
	spCov[1:N, 1:N] <- expCov(Dist[1:N, 1:N], phi, sigmaSp)
	w[1:N] ~ dmnorm(zeros[1:N], cov=spCov[1:N,1:N])
        for (i in 1:p) { beta[i] ~ dflat() }

	phi ~ dunif(aPhi, bPhi)
        sigmasqSp ~ dinvgamma(aSp, bSp)
	sigmaSp <- sqrt(sigmasqSp)	
        sigmasq <- deltasq*sigmasqSp
	deltasq ~ dbeta(a, b)
	tausq <- 1/sigmasq
        sigma <- sqrt(sigmasq)
}
)

nimbleData <- list(Y = df$temp)

modelParameters <- c("beta", "sigmasq", "deltasq", "sigmasqSp", "phi", "w", "yFit", "yResid")

N <- nrow(df)
X <- cbind(rep(1, times=N), df$xUTM, df$yUTM)
p <- ncol(X)
Dist <- rdist(coordsUTM)
aSp <- 2 
bSp <- 10
a <- 1
b <- 9
aPhi <- 3/500
bPhi <- 3/250
zeros <- rep(0, times=N)

nimbleConstants <- list(N=N, X=X, p=p, a=a, b=b, aSp=aSp, bSp=bSp, aPhi=aPhi, bPhi=bPhi, Dist=Dist, zeros=zeros)
nimbleInits <- list(beta=rep(0, times=p), sigmasq=1.0, sigmasqSp=1.0, phi=3/300,  yFit=rep(0, times=N), w=rep(0, times=N))

rModel <- nimbleModel(code=geostatCodeFull, name="coloradoSP", constants=nimbleConstants, data=nimbleData, inits=nimbleInits)
cModel <- compileNimble(rModel)
conf <- configureMCMC(rModel, monitors=modelParameters)
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

n_chains = 1 ## Number of different chains
n_iter = 11000 ## Number of iterations per chain
n_burn = 1000 ## Number of initial iterations ("burn" period) for MCMC to converge. These samples will be discarded and not used in posterior inference. 

samps <- runMCMC(cMCMC, niter = n_iter, nburnin = n_burn, samplesAsCodaMCMC = TRUE)

calculateWAIC(samps, rModel)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_phi <- credibleIntervals[grep("phi", rownames(credibleIntervals)),]
CI_phi

CI_yFit <- credibleIntervals[grep("yFit", rownames(credibleIntervals)),]
##CI_yFit ##Uncomment to print out credible intervals of fitted values
yFitPosteriorMedians <- CI_yFit[, grep("50%", colnames(CI_yFit))]

CI_w <- credibleIntervals[grep("w", rownames(credibleIntervals)),]
wPosteriorMedians <- CI_w[, grep("50%", colnames(CI_w))]

library(sp)
library(raster)

processSurf <- mba.surf(cbind(coords, wPosteriorMedians), no.X = 200, no.Y = 200,
    extend = TRUE, sp = TRUE)$xyz.est

proj4string(processSurf) <- "+proj=longlat +datum=WGS84"

processSurf <- raster(processSurf)

pal <- colorNumeric(blue.red, values(processSurf), na.color = "transparent")

base.map %>%
    addRasterImage(processSurf, colors = pal, opacity = 0.75, group = "Spatial process") %>%
    addLegend("bottomright", pal = pal, values = values(processSurf), opacity = 0.75,
        title = "<center>Spatial<br> process</center>") %>%
    addLayersControl(baseGroup = c("Satellite"), overlayGroups = c("Spatial process"),
        options = layersControlOptions(collapsed = FALSE))
