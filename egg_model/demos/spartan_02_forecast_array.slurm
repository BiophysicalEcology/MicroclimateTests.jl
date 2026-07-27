#!/bin/bash
#SBATCH --job-name=pa_forecast_members
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --array=0-24                     # default for n_ensembles=25, members_per_task=1 -- spartan_submit_all.sh overrides this via `sbatch --array=...` when N_ENSEMBLES/MEMBERS_PER_TASK_FORECAST differ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Each array task solves a *chunk* of members_per_task ACCESS-S2 members in
# one julia process -- chunks still overlap across the cluster instead of
# queueing behind each other on one node (the point of the split), but this
# way the historical leg is only loaded once per chunk instead of once per
# member, and there are fewer array tasks competing for queue slots. Submit
# with --dependency=afterok:<spartan_01_historical job id>, since every
# member needs the historical leg's cached final soil state.

set -euo pipefail

module load Julia/1.12.6-linux-x86_64

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/scratch/projects/punim0593/rasters
export LOCUST_FORECAST_OUTPUT_DIR=/data/scratch/projects/punim0593/locust_forecast
export JULIA_PKG_PRECOMPILE_AUTO=0   # cache was already warmed by spartan_01 -- without this, tasks starting at once can race to recompile

PROJECT_DIR=/home/mrke/git/BiophysicalEcologyEnv
SCRIPT=/home/mrke/git/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_microclimate.jl

members_per_task=${MEMBERS_PER_TASK_FORECAST:-1}   # set by spartan_submit_all.sh; falls back to 1 (one member/task) when run standalone
n_ensembles=${LOCUST_N_ENSEMBLES:-25}
start=$((SLURM_ARRAY_TASK_ID * members_per_task + 1))
end=$((start + members_per_task - 1))
# n_ensembles doesn't always divide evenly by members_per_task -- without
# this the last task's range overshoots n_ensembles (e.g. 99 members,
# members_per_task=10 -> last task computes end=100) and fails outright.
(( end > n_ensembles )) && end=$n_ensembles

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT $start $end
