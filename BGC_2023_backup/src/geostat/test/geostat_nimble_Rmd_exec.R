rm(list=ls())

library(rmarkdown)

render("geostat_nimble.Rmd")

knitr::purl(input = "geostat_nimble.Rmd", output = "geostat_nimble_purl.R", documentation = 0)

