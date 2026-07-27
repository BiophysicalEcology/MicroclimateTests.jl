#!/bin/bash
#SBATCH --job-name=pa_forecast_historical
#SBATCH --account=punim0593
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Solves and caches the historical (SILO) leg only -- must complete before
# spartan_02_forecast_array.slurm, since every forecast member needs its
# final soil state as init=. Submit that with --dependency=afterok:<this job id>.

set -euo pipefail

module load Julia/1.12.6-linux-x86_64   # check `module avail julia` if this version disappears

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GKSwstype=100   # headless plotting -- no display on a compute node
export RASTERDATASOURCES_PATH=/data/scratch/projects/punim0593/rasters
export LOCUST_FORECAST_OUTPUT_DIR=/data/scratch/projects/punim0593/locust_forecast
export JULIA_PKG_PRECOMPILE_AUTO=0   # this job precompiles explicitly below; the array job (spartan_02) inherits the warm cache and must not race to recompile it

PROJECT_DIR=/home/mrke/git/BiophysicalEcologyEnv
SCRIPT=/home/mrke/git/MicroclimateTests.jl/egg_model/demos/points_australia_forecast_microclimate.jl

julia --project=$PROJECT_DIR -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=$PROJECT_DIR --threads=$SLURM_CPUS_PER_TASK $SCRIPT 0
