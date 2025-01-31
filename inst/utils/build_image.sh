#! /bin/bash
echo "Building the singularity container ./singularity/container.sif from ./inst/singularity/recipe.txt"
sudo singularity build --force ./inst/singularity/container.sif ./inst/singularity/recipe.txt
