#!/bin/bash
#SBATCH --job-name=pa_forecast_eggmodel
#SBATCH --account=REPLACE_ME
#SBATCH --partition=REPLACE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G                       # egg-model runs are light; microclimate data is loaded from cache, not re-solved
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Loads every cached microclimate leg/member (written by spartan_01/02) and
# runs the egg model over all of them, producing the stats CSV and maps.
# Submit with --dependency=afterok:<spartan_02_forecast_array job id> -- for
# a job array, afterok waits for every task in the array to finish.

set -euo pipefail

module load Julia

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/projects/REPLACE_ME/spatial_data

PROJECT_DIR=/home/REPLACE_ME/BiophysicalEcologyEnv
SCRIPT=/home/REPLACE_ME/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_eggmodel.jl

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT
