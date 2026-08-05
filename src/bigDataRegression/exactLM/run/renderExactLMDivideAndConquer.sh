#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

RMD_FILE="exactLMDivideAndConquer.Rmd"
R_OUTPUT="exactLMDivideAndConquer.R"

echo "=================================================="
echo "Processing: $RMD_FILE"
echo "=================================================="

# 1. Extract pure R code from the Rmd file
echo "[1/2] Extracting pure R code via knitr::purl()..."
Rscript -e "knitr::purl(input = '$RMD_FILE', output = '$R_OUTPUT', documentation = 1)"

# 2. Render Rmd file into an HTML report
echo "[2/2] Rendering HTML document via rmarkdown::render()..."
Rscript -e "rmarkdown::render(input = '$RMD_FILE', output_format = 'html_document')"

echo "=================================================="
echo "Successfully generated:"
echo " - HTML Report: ${RMD_FILE%.Rmd}.html"
echo " - Pure R Code: $R_OUTPUT"
echo "=================================================="

