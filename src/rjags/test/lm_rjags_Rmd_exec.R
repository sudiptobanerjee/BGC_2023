library(rmarkdown)

render("lm_rjags.Rmd")

knitr::purl(input = "lm_rjags.Rmd", output = "lm_rjags_purl.R",documentation = 0)

