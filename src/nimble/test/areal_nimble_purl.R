rm(list = ls())

source("../src/data_prep_areal.R")

library(nimble)
library(coda)

attributes(nimbleData)
attributes(nimbleConstants)
attributes(nimbleInits)

n_iter = 11000 ## Number of iterations
n_chains = 3 ## Number of different chains
n_burn = 1000 ## Number of initial runs ("burn" period) for MCMC to converge. 

model_parameters = c("beta", "tausq", "sigma", "w", "sigma_sp", "yFit")

adjInfo <- as.carAdjacency(A)
L <- length(adjInfo$adj)
nimbleConstants <- c(nimbleConstants, adjInfo, L=L)

spCode <- nimbleCode({
        for (i in 1:N) {
                Y[i] ~ dnorm(mu[i], tausq)
                mu[i] <- inprod(X[i,1:p], beta[1:p]) + w[i]
        	yFit[i] ~ dnorm(mu[i], tausq)
	}

	w[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tausq_sp, zero_mean = 1)
	tausq_sp ~ dgamma(0.01, 0.01)
	sigma_sp <- 1/sqrt(tausq_sp)

	for (i in 1:p) { beta[i] ~ dflat()}
	tausq ~ dgamma(0.001, 0.001)
	sigma <- 1/sqrt(tausq)
}
)   

rModel <- nimbleModel(code = spCode, constants=nimbleConstants, data=nimbleData)
mcmc.out <- nimbleMCMC(model = rModel, nchains=n_chains, inits=nimbleInits, niter = n_iter, nburnin=n_burn, monitor=model_parameters)
##mcmc.out <- nimbleMCMC(code = spCode, constants=nimbleConstants, data=nimbleData, inits=nimbleInits, nchains=n_chains, niter = n_iter, nburnin=n_burn, monitor=model_parameters) ##The above two steps can be combined into one step


attributes(mcmc.out)

samps <- rbind(mcmc.out$chain1, mcmc.out$chain2, mcmc.out$chain3)
dim(samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_sigma_sp <- credibleIntervals[grep("sigma_sp", rownames(credibleIntervals)),]
CI_sigma_sp

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
