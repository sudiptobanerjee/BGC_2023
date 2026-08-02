rm(list = ls())

source("../src/data_prep_multivariate_areal.R")

library(nimble)
library(coda)

spCode <- nimbleCode({

	for(i in 1:N) {
		Y1[i] ~ dnorm(mu1[i], tau1)               
		Y2[i] ~ dnorm(mu2[i], tau2)
		mu1[i] <- inprod(X[i,1:p], beta1[1:p]) + w1[i]
		mu2[i] <- inprod(X[i,1:p], beta2[1:p]) + w2[i]
		
		mu[i] <- eta_0*w2[i] + eta_1*inprod(A[i,1:N], w2[1:N]) #To be used in CAR model mean for w1 | w2
	
		y1Fit[i] ~ dnorm(mu1[i], tau1)
		y2Fit[i] ~ dnorm(mu2[i], tau2)
	}

	w2[1:N] ~ dcar_proper(zeros[1:N], adj=adj[1:L], num=num[1:N], tau=tau1Sp, gamma=gamma1)
 	w1[1:N] ~ dcar_proper(mu[1:N], adj=adj[1:L], num=num[1:N], tau=tau2Sp, gamma=gamma2) 
	
	tau1Sp ~ dgamma(0.01, 0.01)
        sigma1Sp <- 1/sqrt(tau1Sp)
	tau2Sp ~ dgamma(0.01, 0.01)
	sigma2Sp <- 1/sqrt(tau2Sp)

        for (i in 1:p) { 
		beta1[i] ~ dflat()
		beta2[i] ~ dflat()
	}
        tau1 ~ dgamma(0.001, 0.001)
        sigma1 <- 1/sqrt(tau1)
	tau2 ~ dgamma(0.001, 0.001)
        sigma2 <- 1/sqrt(tau2)

	eta_0 ~ dnorm(0, 0.01)                      
	eta_1 ~ dnorm(0, 0.01)

	gamma1 ~ dunif(0, 0.999)       # avoid improper CAR to maintain identifiability     
	gamma2 ~ dunif(0, 0.999)
})

nimbleInits <- list(beta1=rep(0,times=p), beta2=rep(0,times=p), tau1=1, tau2=1, eta_0=0.3, eta_1=0.3, gamma1=0.5, gamma2=0.5, w1=rep(0, times=N), w2=rep(0, times=N), tau1Sp=1, tau2Sp=1, y1Fit=rep(0, times=N), y2Fit=rep(0, times=N))

model_parameters = c("beta1", "beta2", "sigma1", "sigma2", "sigma1Sp", "sigma2Sp", "eta_0", "eta_1", "gamma1", "gamma2", "w1", "w2", "y1Fit", "y2Fit")

adjInfo <- as.carAdjacency(A)
L <- length(adjInfo$adj)
nimbleConstants <- c(nimbleConstants, adjInfo, L=L)

n_iter = 11000 ## Number of iterations
n_chains = 1 ## Number of different chains
n_burn = 1000 ## Number of initial runs ("burn" period) for MCMC to converge. 

rModel <- nimbleModel(code = spCode, constants=nimbleConstants, data=nimbleData)
mcmc.out <- nimbleMCMC(model = rModel, nchains=n_chains, inits=nimbleInits, niter = n_iter, nburnin=n_burn, monitor=model_parameters)

samps <- as.matrix(mcmc.out)
dim(samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta1 <- credibleIntervals[grep("beta1", rownames(credibleIntervals)),]
CI_beta1

CI_beta2 <- credibleIntervals[grep("beta2", rownames(credibleIntervals)),]
CI_beta2

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_eta <- credibleIntervals[grep("eta_", rownames(credibleIntervals)),]
CI_eta

CI_w1 <- credibleIntervals[grep("w1", rownames(credibleIntervals)),]
CI_w2 <- credibleIntervals[grep("w2", rownames(credibleIntervals)),]

CI_y1Fit <- credibleIntervals[grep("y1Fit", rownames(credibleIntervals)),]
CI_y2Fit <- credibleIntervals[grep("y2Fit", rownames(credibleIntervals)),]

library(maps)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(classInt)

caCounty = map("county","california", fill=TRUE, plot=FALSE)
countyID <- sapply(strsplit(caCounty$names, ","), function(x) x[2])
caCountySF <- st_as_sf(caCounty)

y1FitPosteriorMedians <- CI_y1Fit[,grep("50%", colnames(CI_y1Fit))]
caCountySF$ModelFitted <- y1FitPosteriorMedians
caCountySF$ObsData <- Y1
brksModelFitted = quantile(caCountySF$ModelFitted, c(0, 0.2, 0.4, 0.6, 0.8, 1))
brksObsData = quantile(caCountySF$ObsData, c(0, 0.2, 0.4, 0.6, 0.8, 1))

colorPalette = rev(brewer.pal(5,"RdBu"))

par(mfrow=c(1,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)
plot(caCountySF["ObsData"], pal = colorPalette, breaks = brksObsData, main = "Age adjusted lung cancer rates")

plot(caCountySF["ModelFitted"], pal = colorPalette, breaks = brksModelFitted, main = "Model fitted lung cancer rates")

y2FitPosteriorMedians <- CI_y2Fit[,grep("50%", colnames(CI_y2Fit))]
caCountySF$ModelFitted <- y2FitPosteriorMedians
caCountySF$ObsData <- Y2
brksModelFitted = quantile(caCountySF$ModelFitted, c(0, 0.2, 0.4, 0.6, 0.8, 1))
brksObsData = quantile(caCountySF$ObsData, c(0, 0.2, 0.4, 0.6, 0.8, 1))

colorPalette = rev(brewer.pal(5,"RdBu"))

par(mfrow=c(1,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)
plot(caCountySF["ObsData"], pal = colorPalette, breaks = brksObsData, main = "Age adjusted larynx cancer rates")

plot(caCountySF["ModelFitted"], pal = colorPalette, breaks = brksModelFitted, main = "Model fitted larynx cancer rates")

w1PosteriorMedians <- CI_w1[,grep("50%", colnames(CI_w1))]
w2PosteriorMedians <- CI_w2[,grep("50%", colnames(CI_w2))]

caCountySF$w1PosteriorMedians <- w1PosteriorMedians
caCountySF$w2PosteriorMedians <- w2PosteriorMedians

brks1 = quantile(w1PosteriorMedians, c(0, 0.2, 0.4, 0.6, 0.8, 1))
brks2 = quantile(w2PosteriorMedians, c(0, 0.2, 0.4, 0.6, 0.8, 1))

colorPalette = rev(brewer.pal(5,"RdBu"))

par(mfrow=c(1,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)
plot(caCountySF["w1PosteriorMedians"], pal = colorPalette, breaks = brks1, main = "Posterior medians of residual lung cancer effects")

plot(caCountySF["w2PosteriorMedians"], pal = colorPalette, breaks = brks2, main = "Posterior medians of residual larynx cancer effects")
