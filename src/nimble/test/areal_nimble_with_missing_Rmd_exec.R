rm(list=ls())

library(rmarkdown)

render("areal_nimble_with_missing.Rmd")

knitr::purl(input = "areal_nimble_with_missing.Rmd", output = "areal_nimble_with_missing_purl.R", documentation = 0)

