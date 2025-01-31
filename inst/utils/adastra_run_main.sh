#!/bin/bash
#SBATCH --account=cad15186
#SBATCH --job-name=run_main
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=23:59:59

module purge

# A CrayPE environment version
module load cpe/23.12
# An architecture
module load craype-x86-genoa
# A compiler to target the architecture
module load PrgEnv-cray
# R environment
module load cray-R

module list

# Optionally set the 'envs' environment variable if provided
# Example: ./adastra_run_main.sh prod
if [ -n "$1" ]; then
  export envs="$1"
fi

srun --ntasks-per-node=1 --cpus-per-task=96 --threads-per-core=1 -- Rscript /lus/home/CT3/cad15186/pgodbillot/EMA_ROC02_Bayesian_simulation/inst/scripts/main.R
