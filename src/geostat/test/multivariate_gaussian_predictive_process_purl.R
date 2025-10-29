rm(list=ls())

knitr::opts_chunk$set(comment = NA, tidy = TRUE)

library(spBayes)
library(fields)
library(cluster)
library(maptools)

data(BEF.dat)

BEF.dat <- BEF.dat[BEF.dat$ALLBIO02_KGH>0,]
bio <- BEF.dat$ALLBIO02_KGH*0.001

bio <- BEF.dat$ALLBIO02_KGH*0.001;
log.bio <- as.matrix(log(0.001*BEF.dat[, c("BOLE02_KGH","BRANCH02_KGH","FOLIAGE02_KGH")]))
colnames(log.bio) <- c("log.bole.mt", "log.branch.mt","log.foliage.mt")
coords <- as.matrix(BEF.dat[,c("XUTM","YUTM")])
plot(coords, pch=19, cex=0.5, xlab="Easting (m)", ylab="Northing (m)")

print(cov(log.bio))

m <- 50
km.knots <- kmeans(coords, m)$centers
cl.knots <- clara(coords, m)$medoids
cd.knots <- cover.design(coords, nd=m)$design

plot(coords, pch=19, cex=0.5, xlab="Easting (m)", ylab="Northing (m)")
points(km.knots, pch=5, cex=1, col="blue")
points(cl.knots, pch=6, cex=1, col="green")
points(cd.knots, pch=7, cex=1, col="red")
legend("bottomright", cex=1, pch=c(19,5,6,7), bty="n", col=c("black","blue","green","red"), legend=c("observations","kmeans","clara", "cover.design"))

q <- 3
n.samples <- 2000
n.ltr <- q*(q+1)/2
model <- list(log.bio[,"log.bole.mt"] ~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3, log.bio[,"log.branch.mt"]~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3, log.bio[,"log.foliage.mt"] ~ ELEV + SLOPE + SUM_02_TC1 + SUM_02_TC2 + SUM_02_TC3)

bef.spMvLM <- spMvLM(model, coords=coords, knots=km.knots, data=BEF.dat, starting=list("phi"=rep(3/500,q), "A"=rep(0.05,n.ltr), "Psi"=rep(0.05,q)), tuning=list("phi"=rep(0.3,q), "A"=rep(0.0001,n.ltr), "Psi"=rep(0.001,q)), priors=list("phi.Unif"=list(rep(3/2500,q), rep(3/100,q)), modified.pp=TRUE, "K.IW"=list(q+1, diag(0.1,q)), "Psi.IG"=list(rep(2,q), rep(0.05,q))), cov.model="exponential", n.samples=n.samples, verbose=TRUE, n.report=100)

burn.in <- floor(0.75*n.samples)
bef.spMvLM <- spRecover(bef.spMvLM, start=burn.in)

library(coda)
round(summary(mcmc(cbind(bef.spMvLM$p.beta.recover.samples, bef.spMvLM$p.theta.recover.samples)))$quantiles[,c(1,3,5)],3)

fitted <- spPredict(bef.spMvLM, start=burn.in, thin=10, pred.covars=bef.spMvLM$X, pred.coords=bef.spMvLM$coords)
y.hat <- rowMeans(fitted$p.y.predictive.samples)
bole <- y.hat[seq(1,length(y.hat),q)]
branch <- y.hat[seq(2,length(y.hat),q)]
foliage <- y.hat[seq(3,length(y.hat),q)]

library(MBA)
res <- 100
par(mfrow=c(2,2))
surf <- mba.surf(cbind(coords,bole), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Bole fitted values")
points(coords)
surf <- mba.surf(cbind(coords,branch), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Branch fitted values")
points(coords)
surf <- mba.surf(cbind(coords,foliage), no.X=res, no.Y=res, extend=FALSE)$xyz.est
image.plot(surf, main="Foliage fitted values")
points(coords)
