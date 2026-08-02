rm(list = ls())

source("../src/data_prep.R")

library(nimble)
library(coda)
set.seed(1234) ##Optional

attributes(nimbleData)
attributes(nimbleConstants)
attributes(nimbleInits)

n_iter = 11000 ## Number of iterations
n_chains = length(nimbleInits) ## Number of different chains
n_burn = 1000 ## Number of initial runs ("burn" period) for MCMC to converge. 

model_parameters = c("beta", "tausq", "sigma", "yFit")

lmCode <- nimbleCode({
        for (i in 1:N) {
                Y[i] ~ dnorm(mu[i], tausq)
                mu[i] <- inprod(X[i,1:p], beta[1:p])
        	yFit[i] ~ dnorm(mu[i], tausq)
	}

                for (i in 1:p) { beta[i] ~ dflat()}
                tausq ~ dgamma(0.001, 0.001)
                sigma <- 1/sqrt(tausq)

}
)   

mcmc.out <- nimbleMCMC(code = lmCode, constants=nimbleConstants, data=nimbleData, inits = nimbleInits, nchains=n_chains, niter = n_iter, nburnin=n_burn, monitor=model_parameters)

attributes(mcmc.out)

samps <- rbind(mcmc.out$chain1, mcmc.out$chain2, mcmc.out$chain3)
dim(samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_yFit <- credibleIntervals[grep("yFit", rownames(credibleIntervals)),]
CI_yFit

	betaSamps <- samps[,1:p]
	sigmaSamps <- samps[,(p+1)]
	yFitSamps <- samps[,(p+3):ncol(samps)]

	muRep <- matrix(tcrossprod(betaSamps, X), ncol=1)
	yRep <- rnorm(nrow(betaSamps)*N, muRep, rep(sigmaSamps, times=N))
	yRep <- matrix(yRep, ncol=58)
	CI_yRep <- t(apply(yRep, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))


plot(CI_yRep[,grep("50%", colnames(CI_yRep))], CI_yFit[,grep("50%", colnames(CI_yFit))], xlab="Posterior medians using yRep", ylab="Posterior medians using yFit")
abline(0,1)

par(mfrow=c(1,2))
plot(CI_yRep[,grep("2.5%", colnames(CI_yRep))], CI_yFit[,grep("2.5%", colnames(CI_yFit))], xlab="Posterior 2.5th quantiles using yRep", ylab="Posterior 2.5th quantiles using yFit")
abline(0,1)
plot(CI_yRep[,grep("97.5%", colnames(CI_yRep))], CI_yFit[,grep("97.5%", colnames(CI_yFit))], xlab="Posterior 97.5th quantiles using yRep", ylab="Posterior 97.5th quantiles using yFit")
abline(0,1)

par(mfrow=c(1,3))
qqplot(yRep[,1], yFitSamps[,1])
abline(0,1)
qqplot(yRep[,37], yFitSamps[,37])
abline(0,1)
qqplot(yRep[,58], yFitSamps[,58])
abline(0,1)

G <- sum((Y - colMeans(yFitSamps))^2)
G

P <- sum(apply(yFitSamps, 2, function (x) {var(x)} ))
P

D <- G+P
D

rModel <- nimbleModel(code=lmCode, name="lm_simple", constants=nimbleConstants, data=nimbleData)
cModel <- compileNimble(rModel)
conf <- configureMCMC(rModel, monitors=model_parameters)
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)

samps <- runMCMC(cMCMC, inits=nimbleInits, nchains=length(nimbleInits), niter=n_iter, nburnin = n_burn, samplesAsCodaMCMC = TRUE)
samps <- do.call(rbind, samps)

credibleIntervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credibleIntervals[grep("beta", rownames(credibleIntervals)),]
CI_beta

CI_sigma <- credibleIntervals[grep("sigma", rownames(credibleIntervals)),]
CI_sigma

CI_yFit <- credibleIntervals[grep("yFit", rownames(credibleIntervals)),]
## CI_yFit ##Uncomment to see posterior credible intervals of all fitted values.

yFitSamps <- samps[,grep("yFit", colnames(samps))]

G <- sum((Y - colMeans(yFitSamps))^2)
G

P <- sum(apply(yFitSamps, 2, function (x) {var(x)} ))
P

D <- G + P
D
