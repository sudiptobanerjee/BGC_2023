library(maps)
library(spdep)
library(maptools)
library(classInt)
library(RColorBrewer)
library(CARBayes)

# Get county map for California
ca.county = map("county","california", fill=TRUE, plot=FALSE)
county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords = coordinates(ca.poly)

plot(ca.poly)
text(ca.coords, county.ID, cex=0.5, col = "red")

# Plot data on map
# Get data for lung cancer
rate_5y <- read.csv("age_adjusted.csv")
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]
rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]

# Import lung cancer data in map data
ca.poly$rate_lung = rate_lung$Age_Adjusted_Rate

##Convert polygon to nb object
ca.nb = poly2nb(ca.poly)
ca.adj.mat = nb2mat(ca.nb, style="B")

model.spatial <- S.CARleroux(Age_Adjusted_Rate~1, data=rate_lung, family="gaussian", W=ca.adj.mat, rho=1, burnin=2000, n.sample=10000, thin=1, verbose=FALSE)

