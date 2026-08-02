rm(list = ls())

library(tidyr)

rate_lung <- read.csv("../data/rate_lung_cancer_with_smoking.csv")


## Prepare variables for nimble
## N denotes number of observations
N <- nrow(rate_lung)

## Y denotes Age_Adjusted_rate
Y <- rate_lung$Age_Adjusted_Rate

## X denotes design matrix with intercept and Smoking_Rate
X <- cbind(rep(1, times=N), rate_lung$Smoking_Rate)

## p denotes number of covariates
p <- ncol(X)


## Create adjacency matrix from map of California
library(maps)
library(sf)
library(spdep)

## Get county map for California from R “maps” library and create an “sf” object. Turn off the spherical geometry option as the “maps” library uses planar projections.
ca.county = map("county", "california", fill=TRUE, plot=FALSE)
ca.county.sf <- st_as_sf(ca.county)
sf_use_s2(FALSE)


## We first create an “nb” object from the sf object. 
caNb = poly2nb(ca.county.sf)

## Create the binary adjacency matrix.
A = nb2mat(caNb, style="B")

## Create nimbleConstants and nimbleData
nimbleConstants <- list(N=N, p=p, X=X)
nMissing <- floor(0.25*N)
missingIndex <- sample(1:N, nMissing)
Y[missingIndex] <- NA
nimbleData <- list(Y=Y)

#Initialize model parameters. We will create a variable to store posterior predictive samples yFit (or replicated data) which will sample from p(Y_{rep}|Y) and also initialize spatial random effects to zero

nimbleInits <- list(beta=rep(0, times=p), tausq=1.0, yFit=rep(0, times=N), w=rep(0, times=N), tausq_sp=1.0)

