rm(list=ls())

library(rmarkdown)

render("gaussian_predictive_process.Rmd")

knitr::purl(input = "gaussian_predictive_process.Rmd", output = "gaussian_predictive_process_purl.R", documentation = 0)

