library(rmarkdown)

render("arealModeling_rjags.Rmd")

knitr::purl(input = "arealModeling_rjags.Rmd", output = "arealModeling_rjags_purl.R",documentation = 0)

