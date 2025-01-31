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

# Run the Singularity container with bind paths for I/O
singularity exec --bind $INPUT_DIR --bind $OUTPUT_DIR \
    $CONTAINER_PATH \
    /bin/bash -c "echo 'Running Singularity job'; singularity run ./container.sif"
