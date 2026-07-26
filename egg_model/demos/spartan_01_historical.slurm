#!/bin/bash
#SBATCH --job-name=pa_forecast_historical
#SBATCH --account=REPLACE_ME
#SBATCH --partition=REPLACE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Solves and caches the historical (SILO) leg only -- must complete before
# spartan_02_forecast_array.slurm, since every forecast member needs its
# final soil state as init=. Submit that with --dependency=afterok:<this job id>.

set -euo pipefail

module load Julia   # check `module avail julia` for the exact module name/version

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100   # headless plotting -- no display on a compute node
export RASTERDATASOURCES_PATH=/data/projects/REPLACE_ME/spatial_data

PROJECT_DIR=/home/REPLACE_ME/BiophysicalEcologyEnv
SCRIPT=/home/REPLACE_ME/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_microclimate.jl

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT 0
