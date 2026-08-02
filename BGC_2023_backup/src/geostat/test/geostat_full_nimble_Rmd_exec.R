rm(list=ls())

library(rmarkdown)

render("geostat_full_nimble.Rmd")

knitr::purl(input = "geostat_full_nimble.Rmd", output = "geostat_full_nimble_purl.R", documentation = 0)

