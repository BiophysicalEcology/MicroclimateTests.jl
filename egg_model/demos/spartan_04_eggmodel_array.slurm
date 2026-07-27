#!/bin/bash
#SBATCH --job-name=pa_forecast_eggmembers
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --array=0-4                     # default for n_ensembles=25, members_per_task=5 -- spartan_submit_all.sh overrides this via `sbatch --array=...` when N_ENSEMBLES/MEMBERS_PER_TASK_EGG differ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# Each array task solves a *chunk* of members_per_task ACCESS-S2 members --
# same members_per_task pattern as spartan_02_forecast_array.slurm. Fewer
# array tasks means fewer separate jobs competing for a priority slot on a
# busy shared partition (25 one-member-each tasks were queuing slowly under
# `(Priority)`, not a resource shortage); cpus-per-task is also shrunk from
# 32 to 8 since the egg model is cheap per point (microclimate is loaded
# from a single-depth NetCDF cache, not re-solved) -- smaller resource asks
# backfill faster. Each member's outcome is cache-skipped if already on
# disk, so re-chunking after a partial run is safe -- already-done members
# are skipped, not resolved.
# Submit with --dependency=afterok:<spartan_02 job id>:<spartan_03 job id>,
# since every member needs its own cached forecast microclimate (spartan_02)
# and the historical egg-model cache (spartan_03).

set -euo pipefail

module load Julia/1.12.6-linux-x86_64

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100
export RASTERDATASOURCES_PATH=/data/scratch/projects/punim0593/rasters
export LOCUST_FORECAST_OUTPUT_DIR=/data/scratch/projects/punim0593/locust_forecast
export JULIA_PKG_PRECOMPILE_AUTO=0   # cache was already warmed by spartan_01

PROJECT_DIR=/home/mrke/git/BiophysicalEcologyEnv
SCRIPT=/home/mrke/git/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_eggmodel_members.jl

members_per_task=${MEMBERS_PER_TASK_EGG:-5}   # set by spartan_submit_all.sh; falls back to 5 when run standalone
n_ensembles=${LOCUST_N_ENSEMBLES:-25}
start=$((SLURM_ARRAY_TASK_ID * members_per_task + 1))
end=$((start + members_per_task - 1))
# n_ensembles doesn't always divide evenly by members_per_task -- without
# this the last task's range overshoots n_ensembles (e.g. 99 members,
# members_per_task=5 -> last task computes end=100) and fails outright.
(( end > n_ensembles )) && end=$n_ensembles

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT $start $end
