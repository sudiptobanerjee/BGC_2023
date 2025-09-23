rm(list = ls())

library(tidyr)

rate_5y <- read.csv("data/age_adjusted.csv")
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


smoking <- read.csv("data/smoking.csv")
names(smoking) <- c("County", "Smoking_Rate")

##We will merge the "smoking" data frame with the rate_lung data frame by column "County". It will be important to check for hidden whitespaces in this column for both data frames.

smoking$County <- trimws(smoking$County, whitespace="[\t\r\n\\h\\v]")

##Also make sure that Smoking_Rate is numeric
smoking$Smoking_Rate <- as.numeric(trimws(sub("%", "", smoking$Smoking_Rate), whitespace="[\t\r\n\\h\\v]"))

##Now we will merge the data into a single data frame 
rate_lung <- merge(rate_lung, smoking, by="County")

##Run a simple linear regression of Age_Adjusted_Rate on Smoking_Rate
out.lm <- lm(Age_Adjusted_Rate ~ Smoking_Rate, data = rate_lung)
summary(out.lm)

##Obtain the residuals from the linear regression model
resid.lm <- residuals(out.lm)

##Plot the residuals on a map
library(maps)
library(sf)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(classInt)

##Get county map for California from R “maps” library and create an “sf” object.
ca.county = map("county","california", fill=TRUE, plot=FALSE)
ca.county.sf <- st_as_sf(ca.county)

ca.county.sf$resid.lm <- resid.lm

brks_resid.lm = quantile(ca.county.sf$resid.lm, c(0, 0.2, 0.4, 0.6, 0.8, 1))

color.palette = rev(brewer.pal(5,"RdBu"))
class(color.palette)
plot(ca.county.sf["resid.lm"], pal = color.palette, breaks = brks_resid.lm, main = "Lung cancer rates after accounting for age and smoking")


