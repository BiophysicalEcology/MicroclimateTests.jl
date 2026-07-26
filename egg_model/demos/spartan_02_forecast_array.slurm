#!/bin/bash
#SBATCH --job-name=pa_forecast_members
#SBATCH --account=REPLACE_ME
#SBATCH --partition=REPLACE_ME
#SBATCH --array=1-25                    # match n_ensembles in points_australia_forecast_setup.jl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# One array task per ACCESS-S2 member -- these overlap across the cluster
# instead of queueing behind each other on one node, which is the actual
# point of the split (see the earlier discussion on core allocation). Submit
# with --dependency=afterok:<spartan_01_historical job id>, since every
# member needs the historical leg's cached final soil state.

set -euo pipefail

module load Julia

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/projects/REPLACE_ME/spatial_data

PROJECT_DIR=/home/REPLACE_ME/BiophysicalEcologyEnv
SCRIPT=/home/REPLACE_ME/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_microclimate.jl

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT $SLURM_ARRAY_TASK_ID
