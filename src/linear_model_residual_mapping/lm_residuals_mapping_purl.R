rm(list = ls())

library(tidyr)

rate_5y <- read.csv("data/age_adjusted.csv")
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]
rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]

head(rate_lung)
rate_lung$State_county

rate_lung$State_county = sapply(strsplit(rate_lung$State_county, ":"), function(x) x[2])

rate_lung$State_county = sub("County", "", rate_lung$State_county)

rate_lung$State_county = sub("Registry", "", rate_lung$State_county)

rate_lung$State_county = substr(rate_lung$State_county, 1, nchar(rate_lung$State_county)-7)

rate_lung$State_county = trimws(rate_lung$State_county)

names(rate_lung)[names(rate_lung) == "State_county"] <- "County"

row.names(rate_lung) <- NULL

smoking <- read.csv("data/smoking.csv")
head(smoking)

names(smoking) <- c("County", "Smoking_Rate")

utf8::utf8_print(smoking$County, utf8 = FALSE)
smoking$County <- trimws(smoking$County, whitespace="[\t\r\n\\h\\v]")

smoking$Smoking_Rate <- as.numeric(trimws(sub("%", "", smoking$Smoking_Rate), whitespace="[\t\r\n\\h\\v]"))

rate_lung <- merge(rate_lung, smoking, by="County")

out.lm <- lm(Age_Adjusted_Rate ~ Smoking_Rate, data = rate_lung)
summary(out.lm)

resid.lm <- residuals(out.lm)

library(maps)
library(sf)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(classInt)

ca.county = map("county","california", fill=TRUE, plot=FALSE)
ca.county.sf <- st_as_sf(ca.county)

ca.county.sf$resid.lm <- resid.lm

brks_resid.lm = quantile(ca.county.sf$resid.lm, c(0, 0.2, 0.4, 0.6, 0.8, 1))

color.palette = rev(brewer.pal(5,"RdBu"))
class(color.palette)

plot(ca.county.sf["resid.lm"], pal = color.palette, breaks = brks_resid.lm, main = "Lung cancer rates after accounting for age and smoking")

knitr::purl(input = "lm_residuals_mapping.Rmd", output = "lm_residuals_mapping_purl.R",documentation = 0)
