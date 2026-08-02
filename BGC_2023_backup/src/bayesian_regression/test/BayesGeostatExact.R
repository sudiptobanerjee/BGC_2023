set.seed(1234)

library(spBayes)

data(FORMGMT.dat)

n = nrow(FORMGMT.dat)
p = 5 ##an intercept an four covariates

nsample <- 500

phi <- 0.0012

coords <- cbind(FORMGMT.dat$Longi, FORMGMT.dat$Lat)
coords <- coords*(pi/180)*6378

beta.prior.mean <- rep(0, times=p)
beta.prior.precision <- matrix(0, nrow=p, ncol=p)

alpha <- 1/1.5

signal.prior.shape <- 2.0
signal.prior.rate <- 10.0

simple.lm <- lm(Y ~ X1+X2+X3+X4, data = FORMGMT.dat)

source("../src/BayesianRegression.R")

sp.lm.conjugate.marginalized <- bayes.geostat.simple.marginalized(simple.lm, nsample, beta.prior.mean, beta.prior.precision, coords, phi, alpha, signal.prior.shape, signal.prior.rate)

beta.posterior.samples = as.matrix(sp.lm.conjugate.marginalized[[1]])
signal.posterior.samples = as.vector(sp.lm.conjugate.marginalized[[2]])

w.recover <- sp.effects.recover(simple.lm, coords, beta.posterior.samples, signal.posterior.samples, alpha)

