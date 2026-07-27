#!/bin/bash
#SBATCH --job-name=pa_forecast_eggmembers
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --array=1-25                    # one array task per ACCESS-S2 member -- update to match n_ensembles
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# One array task per ACCESS-S2 member -- spread across as many nodes as
# Slurm has free, rather than chunking several members onto one node like
# spartan_02_forecast_array.slurm does for the (much heavier) microclimate
# solve. The egg model itself is cheap per member (microclimate is loaded
# from a single-depth NetCDF cache, not re-solved), so per-member overhead
# here is dominated by node/job startup, not compute -- worth paying that to
# get full member-level concurrency across nodes.
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

julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID
