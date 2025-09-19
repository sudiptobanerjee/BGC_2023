
## R Programs in BCG Second Edition, pages 44--49.

## Load some libraries we will use
library(spBayes)
library(classInt)
library(RColorBrewer)

## Load a dataset from spBayes. Remove rows with missing values
data(WEF.dat)
WEF.dat <- WEF.dat[!apply(WEF.dat[,c("East_m","North_m", "DBH_cm","Tree_height_m","ELEV_m")], 1, function(x)any(is.na(x))),]
DBH <- WEF.dat$DBH_cm
HT <- WEF.dat$Tree_height_m
coords <- as.matrix(WEF.dat[,c("East_m","North_m")])
plot(coords, pch=1, cex=sqrt(DBH)/10, col="darkgreen", xlab="Easting (m)", ylab="Northing (m)")
leg.vals <- round(quantile(DBH),0)
legend("topleft", pch=1, legend=leg.vals, col="darkgreen",
pt.cex=sqrt(leg.vals)/10, bty="n", title="DBH (cm)")

## Create a color palette for subsequent plots.
col.br <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
col.pal <- col.br(5)

fixed <- classIntervals(DBH, n = 4, style = "fixed", fixedBreaks = c(0, 12.7, 30.48, 60, max(DBH) + 1))
fixed.col <- findColours(fixed, col.pal)
plot(coords, col = fixed.col, pch = 19, cex = 0.5, main = "Forestry tree size classes", xlab = "Easting (m)", ylab = "Northing (m)")
legend("topleft", fill = attr(fixed.col, "palette"), legend = c("sapling", "poletimber", "sawtimber", "large sawtimber"), bty = "n")

## Load the MBA and fields libraries for creating surface interpolation plots
library(MBA)
library(fields) ## For using the image.plot function
x.res <- 100
y.res <- 100
surf <- mba.surf(cbind(coords, DBH), no.X = x.res, no.Y = y.res, h = 5, m = 2, extend = FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab = "Easting (m)", ylab = "Northing (m)", col = col.br(25))
contour(surf, add=T) ## (Optional) Adds contour lines to the plot

## (Optional) Load rgl libraries for drape surface plot
library(rgl)
col <- rbind(0, cbind(matrix(drape.color(surf[[3]],
col = col.br(25)), x.res - 1, y.res - 1), 0))
surface3d(surf[[1]], surf[[2]], surf[[3]], col = col)
axes3d()

title3d(main = "DBH", xlab = "Easting (m)", ylab = "Northing (m)", zlab = "DBH (cm)")
drape.plot(surf[[1]], surf[[2]], surf[[3]], col = col.br(150), theta = 225, phi = 50, border = FALSE, add.legend = FALSE, xlab = "Easting (m)", ylab = "Northing (m)", zlab = "DBH (cm)")
image.plot(zlim = range(surf[[3]], na.rm = TRUE), legend.only = TRUE, horizontal = FALSE)

## Load geoR library for computing variograms
library(geoR)
max.dist <- 0.25 * max(iDist(coords))
bins <- 50

vario.DBH <- variog(coords = coords, data = DBH, uvec = (seq(0, max.dist, length = bins)))
fit.DBH <- variofit(vario.DBH, ini.cov.pars = c(600,200/-log(0.05)), cov.model = "exponential", minimisation.function = "nls", weights = "equal")

## Run an OLS regression and make variograms from the residuals.
lm.DBH <- lm(DBH ~ Species, data = WEF.dat)
summary(lm.DBH)
DBH.resid <- resid(lm.DBH)

vario.DBH.resid <- variog(coords = coords, data = DBH.resid, uvec = (seq(0, max.dist, length = bins)))
fit.DBH.resid <- variofit(vario.DBH.resid, ini.cov.pars = c(300, 200/-log(0.05)), cov.model = "exponential", minimisation.function = "nls", weights = "equal")

par(mfrow = c(1, 2))
plot(vario.DBH, ylim = c(200, 1200), main = "DBH")
lines(fit.DBH)
abline(h = fit.DBH$nugget, col = "blue")
abline(h = fit.DBH$cov.pars[1] + fit.DBH$nugget, col = "green")
abline(v = -log(0.05) * fit.DBH$cov.pars[2], col = "red3")
plot(vario.DBH.resid, ylim = c(200, 500), main = "DBH residuals")
lines(fit.DBH.resid)
abline(h = fit.DBH.resid$nugget, col = "blue")
abline(h = fit.DBH.resid$cov.pars[1] + fit.DBH.resid$nugget, col = "green")
abline(v = -log(0.05) * fit.DBH.resid$cov.pars[2], col = "red3")

### Illustrate variograms using the gstat package
library(gstat)
ELEV_m <- WEF.dat$ELEV_m
sp.dat <- as.data.frame(cbind(DBH, HT, coords, ELEV_m))
coordinates(sp.dat) <- c("East_m", "North_m")
vario.DBH <- variogram(DBH ~ ELEV_m, data = sp.dat,utoff = max.dist, width = 5, alpha = (0:3) * 45)
fit.DBH <- fit.variogram(vario.DBH, vgm(1000, "Exp", 200/-log(0.05), 600))
print(plot(vario.DBH, fit.DBH))
