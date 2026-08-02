rm(list = ls())

rate_5y <- read.csv("../data/age_adjusted.csv")
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]

rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]
rate_lung$State_county = sapply(strsplit(rate_lung$State_county, ":"), function(x) x[2])
rate_lung$State_county = sub("County", "", rate_lung$State_county)
rate_lung$State_county = sub("Registry", "", rate_lung$State_county)
rate_lung$State_county = substr(rate_lung$State_county, 1, nchar(rate_lung$State_county)-7)
rate_lung$State_county = trimws(rate_lung$State_county)
names(rate_lung)[names(rate_lung) == "State_county"] <- "County"
row.names(rate_lung) <- NULL

rate_larynx = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Larynx",]
rate_larynx = rate_larynx[order(readr::parse_number(as.character(rate_larynx$State_county))),]
rate_larynx$State_county = sapply(strsplit(rate_larynx$State_county, ":"), function(x) x[2])
rate_larynx$State_county = sub("County", "", rate_larynx$State_county)
rate_larynx$State_county = sub("Registry", "", rate_larynx$State_county)
rate_larynx$State_county = substr(rate_larynx$State_county, 1, nchar(rate_larynx$State_county)-7)
rate_larynx$State_county = trimws(rate_larynx$State_county)
names(rate_larynx)[names(rate_larynx) == "State_county"] <- "County"
row.names(rate_larynx) <- NULL

smoking <- read.csv("../data/smoking.csv")

names(smoking) <- c("County", "Smoking_Rate")

smoking$County <- trimws(smoking$County, whitespace="[\t\r\n\\h\\v]")

smoking$Smoking_Rate <- as.numeric(trimws(sub("%", "", smoking$Smoking_Rate), whitespace="[\t\r\n\\h\\v]"))

rate_lung_larynx <- merge(rate_lung, rate_larynx, by="County")
rate_lung_larynx <- merge(rate_lung_larynx, smoking, by="County")

N <- nrow(rate_lung_larynx)

Y1 <- rate_lung_larynx$Age_Adjusted_Rate.x
Y2 <- rate_lung_larynx$Age_Adjusted_Rate.y

X <- rate_lung_larynx$Smoking_Rate
X <- cbind(rep(1, times=N), X)

p <- ncol(X)

zeros <- rep(0, times=N)

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


nimbleConstants <- list(N=N, p=p, X=X, A=A, zeros=zeros) 
nimbleData <- list(Y1=Y1, Y2=Y2)

