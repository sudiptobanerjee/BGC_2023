rm(list = ls())

library(rjags)
library(coda)

source("../src/data_prep_arealModeling.R")

inits = list(
	     list(beta=rep(0,times=p), tausq=1.0, tausqSp=1.0, w=rep(0, times=N), yFit = rep(0, times=N)),
	     list(beta=rep(-100, times=p), tausq=1.0, tausqSp=1.0, w=rep(0, times=N), yFit = rep(0, times=N)),
	     list(beta=rep(100,times=p), tausq=1.0, tausqSp=1.0, w=rep(0, times=N), yFit = rep(0, times=N))
)

nSamp = 10000 ## Number of posterior samples to be used for inference
nChains = 3 ## Number of different chains
nAdapt = 500 ## The initial number of runs to tune parameters
nBurn = 1000 ## Number of initial runs ("burn" period) for MCMC to converge. 

model.parameters = c("beta", "tausq", "w", "tausqSp", "sigma", "sigmaSp", "yFit")

m1 <- jags.model("../src/model_jags_arealModeling.txt", data=jagsData, inits = inits, n.chains=nChains, n.adapt = nAdapt, quiet=TRUE)

update(m1, n.iter=nBurn)

m1.out <- coda.samples(model=m1, variable.names=model.parameters, n.iter=nSamp)

sm <- summary(m1.out)

samps <- do.call(rbind, m1.out)
dim(samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))
credibleIntervals

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_sigmaSp <- credibleIntervals[grep("sigmaSp", rownames(credibleIntervals)),]
CI_sigmaSp

CI_w <- credibleIntervals[grep("w", rownames(credibleIntervals)),]
CI_w

CI_yFit <- credibleIntervals[grep("yFit", rownames(credibleIntervals)),]
CI_yFit

library(maps)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(classInt)

caCounty = map("county","california", fill=TRUE, plot=FALSE)
countyID <- sapply(strsplit(caCounty$names, ","), function(x) x[2])
caCountySF <- st_as_sf(caCounty)

yFitPosteriorMedians <- CI_yFit[,grep("50%", colnames(CI_yFit))]
caCountySF$ModelFitted <- yFitPosteriorMedians
caCountySF$ObsData <- Y
brksModelFitted = quantile(caCountySF$ModelFitted, c(0, 0.2, 0.4, 0.6, 0.8, 1))
brksObsData = quantile(caCountySF$ObsData, c(0, 0.2, 0.4, 0.6, 0.8, 1))

colorPalette = rev(brewer.pal(5,"RdBu"))

par(mfrow=c(1,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(caCountySF["ObsData"], pal = colorPalette, breaks = brksObsData, main = "Age adjusted lung cancer rates")

plot(caCountySF["ModelFitted"], pal = colorPalette, breaks = brksModelFitted, main = "Model fitted lung cancer rates")

cor(Y, yFitPosteriorMedians) ## correlation between observed and model fitted values 
plot(Y, yFitPosteriorMedians, xlab="Observed data", ylab="Model fitted data")
abline(0,1) ## adds the 45 degree (y=x) line to identify departures of fitted values from observations.

yFitSamps <- samps[,grep("yFit", colnames(samps))]

G <- sum((Y - colMeans(yFitSamps))^2)
G

P <- sum(apply(yFitSamps, 2, function (x) {var(x)} ))
P

D <- G+P
D

inits = list(
	     list(beta=rep(0,times=p), tausq=1.0, yFit = rep(0, times=N)),
	     list(beta=rep(-100, times=p), tausq=1.0, yFit = rep(0, times=N)),
	     list(beta=rep(100,times=p), tausq=1.0, yFit = rep(0, times=N))
)

model.parameters = c("beta", "tausq", "sigma", "yFit")

m2 <- jags.model("../src/model_jags_simpleLM.txt", data=jagsDataSimpleLM, inits = inits, n.chains=nChains, n.adapt = nAdapt, quiet=TRUE)

update(m2, n.iter=nBurn)

m2.out <- coda.samples(model=m2, variable.names=model.parameters, n.iter=nSamp)

samps <- do.call(rbind, m2.out)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))
CI_yFit <- credibleIntervals[grep("yFit", rownames(credibleIntervals)),]
yFitPosteriorMedians <- CI_yFit[,grep("50%", colnames(CI_yFit))]
cor(Y, yFitPosteriorMedians) ## correlation between observed and model fitted values 
plot(Y, yFitPosteriorMedians, xlab="Observed data", ylab="Model fitted data")
abline(0,1) ## adds the 45 degree (y=x) line to identify departures of fitted values from observations.


yFitSamps <- samps[,grep("yFit", colnames(samps))]

G <- sum((Y - colMeans(yFitSamps))^2)
G

P <- sum(apply(yFitSamps, 2, function (x) {var(x)} ))
P

D <- G+P
D
