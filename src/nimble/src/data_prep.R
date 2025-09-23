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

smoking <- read.csv("../data/smoking.csv")

names(smoking) <- c("County", "Smoking_Rate")

smoking$County <- trimws(smoking$County, whitespace="[\t\r\n\\h\\v]")

smoking$Smoking_Rate <- as.numeric(trimws(sub("%", "", smoking$Smoking_Rate), whitespace="[\t\r\n\\h\\v]"))

rate_lung <- merge(rate_lung, smoking, by="County")

N <- nrow(rate_lung)

Y <- rate_lung$Age_Adjusted_Rate

X <- rate_lung$Smoking_Rate

X <- cbind(rep(1, times=N), X)

p <- ncol(X)

nimbleConstants <- list(N=N, p=p, X=X) 
nimbleData <- list(Y=Y)

#We will create a variable to store posterior predictive samples yFit (or replicated data) which will sample from p(Y_{rep}|Y)

nimbleInits = list(
             list(beta=rep(0, times=p), tausq=1.0, yFit=rep(0, times=N)),
             list(beta=rep(-100, times=p), tausq=10.0, yFit=rep(0, times=N)),
             list(beta=rep(100, times=p), tausq=0.10, yFit=rep(0, times=N))
             )


