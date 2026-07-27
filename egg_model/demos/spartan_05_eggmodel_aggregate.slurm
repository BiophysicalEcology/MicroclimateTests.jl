#!/bin/bash
#SBATCH --job-name=pa_forecast_eggaggregate
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Loads the historical egg-model cache and every member's cached egg-model
# outcomes (written by spartan_03/spartan_04) and produces the combined
# stats CSV + maps. No egg-model solving here -- CSV write + a few heatmaps,
# doesn't benefit from more cores, so this deliberately stays small rather
# than matching the 32 used by the compute-heavy stages.
# Submit with --dependency=afterok:<spartan_04 job id> -- for a job array,
# afterok waits for every task in the array to finish.

set -euo pipefail

module load Julia/1.12.6-linux-x86_64

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/scratch/projects/punim0593/rasters
export LOCUST_FORECAST_OUTPUT_DIR=/data/scratch/projects/punim0593/locust_forecast
export JULIA_PKG_PRECOMPILE_AUTO=0

PROJECT_DIR=/home/mrke/git/BiophysicalEcologyEnv
SCRIPT=/home/mrke/git/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_eggmodel_aggregate.jl

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT
