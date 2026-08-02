#!/bin/bash
# ==========================================
# File: run/renderRaw.sh
# ==========================================

# Print status
echo "Rendering src/theoryAndRaw.Rmd into HTML format..."

# Use Rscript to call rmarkdown::render
# We point to the Rmd file in the src directory and explicitly set the output directory to '.' (the run directory)
Rscript -e "rmarkdown::render('theoryAndRaw.Rmd', output_dir = '.')"

Rscript -e "rmarkdown::render('theoryPredictions.Rmd', output_dir = '.')"

Rscript -e "rmarkdown::render('modelComparisons.Rmd', output_dir = '.')"

Rscript -e "rmarkdown::render('theoryAndRawWithQR.Rmd', output_dir = '.')"

#Use purl command in knitr package to create R programs from Rmd files

Rscript -e "knitr::purl(input = 'theoryAndRaw.Rmd', output = 'theoryAndRaw_purl.R', documentation = 0)"

Rscript -e "knitr::purl(input = 'theoryPredictions.Rmd', output = 'theoryPredictions_purl.R', documentation = 0)"

Rscript -e "knitr::purl(input = 'modelComparisons.Rmd', output = 'modelComparisons_purl.R', documentation = 0)"

Rscript -e "knitr::purl(input = 'theoryAndRawWithQR.Rmd', output = 'theoryAndRawWithQR_purl.R', documentation = 0)"


# Notify when finished
echo "Rendering and purl complete. The HTML and R files are all located in the run directory."

