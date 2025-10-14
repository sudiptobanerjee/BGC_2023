rm(list=ls())

library(rmarkdown)

render("univariate_spatial_regression.Rmd")

knitr::purl(input = "univariate_spatial_regression.Rmd", output = "univariate_spatial_regression__purl.R", documentation = 0)

