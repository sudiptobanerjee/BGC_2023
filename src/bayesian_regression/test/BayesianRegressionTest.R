set.seed(1234)

source("../src/BayesianRegression.R")

dataObject <- read.table("../data/LinearModelExample.txt", header=T)

dataObject.lm <- lm(Y ~ X1+X2+X3+X4+X5+X6, data = dataObject)

lm.frame = model.frame(dataObject.lm)
lm.frame.cov = vcov(dataObject.lm)

n = nrow(lm.frame)
p = ncol(lm.frame.cov)

nsample <- 500


##Here we obtain posterior samples using improper priors
out.bayes.lm.ref <- bayes.lm.ref(dataObject.lm, nsample)


##Below we demonstrate the conjugate function in the special case with improper priors. The results are the same as for the above, up to MC error. 

beta.prior.mean <- rep(0, times=p)

beta.prior.precision <- matrix(0, nrow=p, ncol=p)

prior.shape <- -p/2

prior.rate <- 0

out.bayes.lm.conjugate <- bayes.lm.conjugate(dataObject.lm, nsample, beta.prior.mean, beta.prior.precision, prior.shape, prior.rate)

##Below we demonstrate exact posterior sampling for a simple geostatistical model by casting it within a conjugate Bayesian linear model 

phi <- 2500
alpha <- 1.5
coords <- cbind(dataObject$Longi, dataObject$Lat)
coords <- coords*(pi/180)*6378

nugget.prior.shape <- 2.0
nugget.prior.rate <- 10.0

#sp.lm.conjugate <- bayes.geostat.exp.simple(dataObject.lm, nsample, beta.prior.mean, beta.prior.precision, coords, phi, alpha, nugget.prior.shape, nugget.prior.rate)

alpha <- 1/1.5
signal.prior.shape <- 2.0
signal.prior.rate <- 10.0

sp.lm.conjugate.marginalized <- bayes.geostat.simple.marginalized(dataObject.lm, nsample, beta.prior.mean, beta.prior.precision, coords, phi, alpha, signal.prior.shape, signal.prior.rate) 

