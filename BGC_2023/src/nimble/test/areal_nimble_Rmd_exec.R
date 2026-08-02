rm(list=ls())

library(rmarkdown)

render("areal_nimble.Rmd")

knitr::purl(input = "areal_nimble.Rmd", output = "areal_nimble_purl.R", documentation = 0)

