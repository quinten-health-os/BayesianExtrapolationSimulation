#!/bin/bash
#SBATCH --account=cad15186
#SBATCH --job-name=singularity_run_main
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=23:59:59
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.fauvel@quinten-health.com    # Email to send notifications

# Path to your Singularity container (either local or a container registry)
CONTAINER_PATH="/opt/software/containers/images/users/cad15186/container.sif"

# Directory where you want to bind host directories to the container
# Adjust paths for your input/output directories if necessary
INPUT_DIR="/lus/home/CT3/cad15186/pgodbillot/EMA_ROC02_Bayesian_simulation"
OUTPUT_DIR="/lus/home/CT3/cad15186/pgodbillot/EMA_ROC02_Bayesian_simulation"

# Run the Singularity container 
salloc --account=cad15186 --job-name="interactive" --constraint=GENOA --nodes=1 --time=23:59:59 --exclusive

singularity run --writable-tmpfs -B/lus/home/CT3/cad15186/pgodbillot/EMA_ROC02_Bayesian_simulation -B/sys -B/dev -B/proc -B/lus /opt/software/containers/images/users/cad15186/container.sif fast_cases_config

singularity shell --writable-tmpfs -B/lus/home/CT3/cad15186/pgodbillot/EMA_ROC02_Bayesian_simulation -B/sys -B/dev -B/proc -B/lus /opt/software/containers/images/users/cad15186/container.sif 
