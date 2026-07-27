#!/bin/bash
#SBATCH --job-name=pa_forecast_egghistorical
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Egg model over the historical (SILO) leg only, cached to disk -- needed by
# every forecast member's egg-model job (spartan_04_eggmodel_array.slurm) for
# forecast_point_indices and each surviving point's final_state. Only needs
# the historical microclimate leg (spartan_01), not the forecast members, so
# submit this in parallel with spartan_02's forecast-member array (both with
# --dependency=afterok:<spartan_01 job id>) rather than after it.

set -euo pipefail

module load Julia/1.12.6-linux-x86_64

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/scratch/projects/punim0593/rasters
export LOCUST_FORECAST_OUTPUT_DIR=/data/scratch/projects/punim0593/locust_forecast
export JULIA_PKG_PRECOMPILE_AUTO=0   # cache was already warmed by spartan_01

PROJECT_DIR=/home/mrke/git/BiophysicalEcologyEnv
SCRIPT=/home/mrke/git/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_eggmodel_historical.jl

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT
