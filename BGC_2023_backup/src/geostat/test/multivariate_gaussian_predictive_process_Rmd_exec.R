rm(list=ls())

library(rmarkdown)

render("multivariate_gaussian_predictive_process.Rmd")

knitr::purl(input = "multivariate_gaussian_predictive_process.Rmd", output = "multivariate_gaussian_predictive_process_purl.R", documentation = 0)

