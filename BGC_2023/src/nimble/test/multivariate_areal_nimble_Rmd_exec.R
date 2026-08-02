rm(list=ls())

library(rmarkdown)

render("multivariate_areal_nimble.Rmd")

knitr::purl(input = "multivariate_areal_nimble.Rmd", output = "multivariate_areal_nimble_purl.R", documentation = 0)

