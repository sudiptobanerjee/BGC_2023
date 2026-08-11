#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# List your Rmd files here
RMD_FILES=(
  "exactMvLM.Rmd"
  "exactMvLMFederated.Rmd"
)

for RMD_FILE in "${RMD_FILES[@]}"; do
  # Extract base name without .Rmd extension
  BASE_NAME="${RMD_FILE%.Rmd}"
  R_OUTPUT="${BASE_NAME}.R"
  HTML_OUTPUT="${BASE_NAME}.html"

  echo "=================================================="
  echo "Processing: $RMD_FILE"
  echo "=================================================="

  # 1. Extract pure R code from the Rmd file
  echo "[1/2] Extracting pure R code via knitr::purl()..."
  Rscript -e "knitr::purl(input = '$RMD_FILE', output = '$R_OUTPUT', documentation = 1)"

  # 2. Render Rmd file into an HTML report
  echo "[2/2] Rendering HTML document via rmarkdown::render()..."
  Rscript -e "rmarkdown::render(input = '$RMD_FILE', output_format = 'html_document', output_file = '$HTML_OUTPUT')"

  echo "=================================================="
  echo "Successfully generated:"
  echo " - HTML Report: $HTML_OUTPUT"
  echo " - Pure R Code: $R_OUTPUT"
  echo "=================================================="
  echo ""
done
