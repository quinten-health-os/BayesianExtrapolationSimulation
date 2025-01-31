#!/bin/bash

# Optionally set the 'envs' environment variable if provided
# Example: ./run_main.sh prod
if [ -n "$1" ]; then
  export envs="$1"
fi

# Run the R script
Rscript ./inst/scripts/main.R
