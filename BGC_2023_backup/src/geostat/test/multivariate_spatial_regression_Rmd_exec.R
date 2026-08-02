rm(list=ls())

library(rmarkdown)

render("multivariate_spatial_regression.Rmd")

knitr::purl(input = "multivariate_spatial_regression.Rmd", output = "multivariate_spatial_regression_purl.R", documentation = 0)

