rm(list=ls())

library(rmarkdown)

render("lm_nimble.Rmd")

knitr::purl(input = "lm_nimble.Rmd", output = "lm_nimble_purl.R", documentation = 0)

