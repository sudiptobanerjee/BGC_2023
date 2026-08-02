rm(list = ls())

library(rjags)
library(coda)

source("../src/data_prep.R")

inits = list(
	     list(beta=c(0,0,0,0,0,0,0), tausq=1.0, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
	     list(beta=c(-100,-100,-100,-100,-100,-100,-100), tausq=10.0, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
	     list(beta=c(100,100,100,100,100,100,100), tausq=0.10, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
	     )

n_samp = 10000 ## Number of posterior samples to be used for inference
n_chains = 3 ## Number of different chains
n_adapt = 500 ## The initial number of runs to tune parameters
n_burn = 1000 ## Number of initial runs ("burn" period) for MCMC to converge. 

model.parameters = c("beta", "tausq", "sigma", "Y.tilde")

m1 <- jags.model("../src/model_jags.txt", data=jags.data, inits = inits, n.chains=n_chains, n.adapt = n_adapt)

update(m1, n.iter=n_burn)

m1.out <- coda.samples(model=m1, variable.names=model.parameters, n.iter=n_samp)

sm <- summary(m1.out)

samps <- do.call(rbind, m1.out)
dim(samps)

credible.intervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))
credible.intervals

CI_beta <- credible.intervals[grep("beta", rownames(credible.intervals)),]
CI_beta
