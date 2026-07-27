#!/bin/bash
# Submits the full points_australia_forecast pipeline as a dependency DAG:
#   historical microclimate (job1)
#     |-> forecast-member microclimate array (job2) --\
#     \-> egg-model historical (job3)  ----------------+-> egg-model member array (job4) -> aggregate: stats + maps (job5)
# job2 and job3 both only need job1, so they run in parallel rather than
# chained one after the other.
#
# Grid extent/spacing and ensemble count are overridable from the environment
# (fall back to points_australia_forecast_setup.jl's own defaults if unset)
# -- e.g.:
#   N_ENSEMBLES=50 LOCUST_GRID_SPACING_DEG=0.2 ./spartan_submit_all.sh
#
# Usage: ./spartan_submit_all.sh  (run from this directory)

set -euo pipefail

N_ENSEMBLES=${N_ENSEMBLES:-99}                    # how many of the up-to-99 ACCESS-S2 members to run
MEMBERS_PER_TASK_FORECAST=${MEMBERS_PER_TASK_FORECAST:-1}   # microclimate array chunk size (spartan_02)
MEMBERS_PER_TASK_EGG=${MEMBERS_PER_TASK_EGG:-5}             # egg-model array chunk size (spartan_04)

(( N_ENSEMBLES >= 1 && N_ENSEMBLES <= 99 )) || { echo "N_ENSEMBLES=$N_ENSEMBLES out of bounds 1:99 -- ACCESS-S2 only has 99 members" >&2; exit 1; }

# LOCUST_N_ENSEMBLES (and LOCUST_LON_MIN/MAX, LOCUST_LAT_MIN/MAX,
# LOCUST_GRID_SPACING_DEG, if the caller set them) reach the Julia jobs via
# sbatch's default environment passthrough -- no --export needed, just
# export them here before submitting.
export LOCUST_N_ENSEMBLES=$N_ENSEMBLES

# array upper bound = ceil(N_ENSEMBLES / members_per_task) - 1
n_forecast_tasks=$(( (N_ENSEMBLES + MEMBERS_PER_TASK_FORECAST - 1) / MEMBERS_PER_TASK_FORECAST ))
n_egg_tasks=$(( (N_ENSEMBLES + MEMBERS_PER_TASK_EGG - 1) / MEMBERS_PER_TASK_EGG ))

job1=$(sbatch --parsable spartan_01_historical.slurm)
echo "Historical microclimate job: $job1"

job2=$(MEMBERS_PER_TASK_FORECAST=$MEMBERS_PER_TASK_FORECAST sbatch --parsable --dependency=afterok:$job1 \
    --array=0-$((n_forecast_tasks - 1)) spartan_02_forecast_array.slurm)
echo "Forecast-member microclimate array job: $job2 ($n_forecast_tasks tasks, $MEMBERS_PER_TASK_FORECAST members/task)"

job3=$(sbatch --parsable --dependency=afterok:$job1 spartan_03_eggmodel_historical.slurm)
echo "Egg-model historical job: $job3"

job4=$(MEMBERS_PER_TASK_EGG=$MEMBERS_PER_TASK_EGG sbatch --parsable --dependency=afterok:$job2:$job3 \
    --array=0-$((n_egg_tasks - 1)) spartan_04_eggmodel_array.slurm)
echo "Egg-model member array job: $job4 ($n_egg_tasks tasks, $MEMBERS_PER_TASK_EGG members/task)"

job5=$(sbatch --parsable --dependency=afterok:$job4 spartan_05_eggmodel_aggregate.slurm)
echo "Aggregate (stats + maps) job: $job5"

echo "Submitted: $job1 -> {$job2, $job3} -> $job4 -> $job5"
