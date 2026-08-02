rm(list=ls())

library(rmarkdown)

render("CO-temp.Rmd")

knitr::purl(input = "CO-temp.Rmd", output = "CO-temp_purl.R", documentation = 0)

