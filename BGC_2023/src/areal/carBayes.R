rm(list = ls())

library(sp)
library(spdep)
library(maps)
library(maptools)

library(CARBayes)
library(CARBayesdata) 
library(shapefiles) 
library(sp) 
data(lipdata) 
data(lipdbf) 
data(lipshp)

library(RColorBrewer)

library(knitr)

lipdbf$dbf <- lipdbf$dbf[ ,c(2,1)]
data.combined <- combine.data.shapefile(data=lipdata, shp=lipshp, dbf=lipdbf)

W.nb <-poly2nb(data.combined, row.names = rownames(lipdata))
W.mat <- nb2mat(W.nb, style="B")

nb.bound <- poly2nb(data.combined) # shared boundaries
#summary(nb.bound)
coords <- coordinates(data.combined)
#plot(data.combined, border = "gray", main="Scotland")
#plot(nb.bound, coords, pch = 19, cex = 0.6, add = TRUE)
breakpoints <- seq(min(lipdata$observed)-1, max(lipdata$observed)+1, length.out=8)
my.palette <- brewer.pal(n = 7, name = "OrRd")

spplot(data.combined, c("observed", "expected"), main="Scottish Lip Cancer", at=breakpoints, col.regions=my.palette, col="grey")

spplot(data.combined, c("observed", "pcaff"), main="Scottish Lip Cancer", at=breakpoints, col.regions=my.palette, col="black")


glmmodel <- glm(observed~., family="poisson", data=data.combined@data)
#summary(glmmodel)
resid.glmmodel <- residuals(glmmodel)

W.nb <- poly2nb(data.combined, row.names = rownames(data.combined@data))
W.list <- nb2listw(W.nb, style="B")
testglm <- moran.mc(x=resid.glmmodel, listw=W.list, nsim=1000)
W <- nb2mat(W.nb, style="B")

formula <- observed ~ log(expected)+pcaff+latitude+longitude
model.spatial1 <- S.CARbym(formula=formula, data=data.combined@data, family="poisson", W=W, burnin=20000, n.sample=120000, thin=10, verbose=FALSE)
betas1 <- summarise.samples(model.spatial1$samples$beta, quantiles=c(0.5, 0.025, 0.975))
resultsMS1 <- betas1$quantiles
rownames(resultsMS1) <- c("Intercept", "Expected", "Pcaff", "Latitude", "Longitude")
kable(resultsMS1, caption="95% Credible Intervals Model 1")

