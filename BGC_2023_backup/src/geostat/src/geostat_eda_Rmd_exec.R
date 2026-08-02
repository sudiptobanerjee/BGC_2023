rm(list=ls())

library(rmarkdown)

render("geostat_eda.Rmd")

knitr::purl(input = "geostat_eda.Rmd", output = "geostat_eda_purl.R", documentation = 0)

