rm(list = ls())

library(rjags)
library(coda)

source("../src/data_prep.R")

inits = list(
	     list(beta=c(0,0,0,0,0,0,0), tausq=1.0, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
	     list(beta=c(-100,-100,-100,-100,-100,-100,-100), tausq=10.0, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
	     list(beta=c(100,100,100,100,100,100,100), tausq=0.10, Y.tilde = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
	     )

n_samp = 10000
n_chains = 3
n_adapt = 500
n_burn = 1000

model.parameters = c("beta", "tausq", "sigma", "Y.tilde")

#—————————————————–
# compile 
#—————————————————–
m1 <- jags.model("../src/model_jags.txt", data=jags.data, inits = inits, n.chains=n_chains, n.adapt = n_adapt)

#—————————————————–
#update and burn 
#—————————————————–
update(m1, n.iter=n_burn)

#—————————————————–
# posterior sampling
#—————————————————–
m1.out <- coda.samples(model=m1, variable.names=model.parameters, n.iter=n_samp)

#—————————————————–
# descriptive statistics of posterior densities
sm <- summary(m1.out)
#—————————————————–

#——————————————————
# extracting full posterior samples n
#—————————————————–
samps <- do.call(rbind, m1.out)

#——————————————————
# medians and 95% credible intervals
#—————————————————–
credible.intervals <- t(apply(samps, 2, function(x){quantile(x, c(0.50, 0.025, 0.975))}))

CI_beta <- credible.intervals[grep("beta", rownames(credible.intervals)),]

