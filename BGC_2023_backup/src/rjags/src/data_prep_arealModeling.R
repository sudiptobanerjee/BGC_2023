rm(list = ls())

library(tidyr)

rate_5y <- read.csv("../data/age_adjusted.csv")
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]

rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]

##We will need to extract polished county names from the data base
##First inspect the State_county variable.
##Get rid of "CA:" from the names
rate_lung$State_county = sapply(strsplit(rate_lung$State_county, ":"), function(x) x[2])
##Get rid of "County" inside the names
rate_lung$State_county = sub("County", "", rate_lung$State_county)
##Get rid of "Registry" inside the names (only "Los Angeles")
rate_lung$State_county = sub("Registry", "", rate_lung$State_county)
##Get rid of "the last 7 characters" from the names (they represent FIPS codes)
rate_lung$State_county = substr(rate_lung$State_county, 1, nchar(rate_lung$State_county)-7)
##Get rid of unnecessary white spaces
rate_lung$State_county = trimws(rate_lung$State_county)

names(rate_lung)[names(rate_lung) == "State_county"] <- "County"

##Remove row numbers from the data frame:
row.names(rate_lung) <- NULL


smoking <- read.csv("../data/smoking.csv")
names(smoking) <- c("County", "Smoking_Rate")

##We will merge the "smoking" data frame with the rate_lung data frame by column "County". It will be important to check for hidden whitespaces in this column for both data frames.

smoking$County <- trimws(smoking$County, whitespace="[\t\r\n\\h\\v]")

##Also make sure that Smoking_Rate is numeric
smoking$Smoking_Rate <- as.numeric(trimws(sub("%", "", smoking$Smoking_Rate), whitespace="[\t\r\n\\h\\v]"))

##Now we will merge the data into a single data frame 
rate_lung <- merge(rate_lung, smoking, by="County")

##If you wish to write this cleaned up dataframe as a new csv file:
write.csv(rate_lung, file="../data/rate_lung_cancer_with_smoking.csv", row.names = FALSE)

## Prepare variables for JAGS
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

## Create diagonal matrix with number of neighbors along the diagonal
D = diag(rowSums(A))

## Create proper precision matrix
rho <- 0.99 ## Spatial autocorrelation (smoothness) parameter; rho <- 1 max smoothing
Q <- D - rho*A

## Create a vector of zeros to be used for the prior mean
zeros <- rep(0, times=N)

## Create data for JAGS
jagsData <- list(N=N, p=p, Y=Y, X=X, Q=Q, zeros=zeros)

## Also create data for JAGS to run a linear regression without spatial effects
jagsDataSimpleLM <- list(N=N, p=p, Y=Y, X=X)


