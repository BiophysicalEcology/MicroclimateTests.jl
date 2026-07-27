#!/bin/bash
# Submits the full points_australia_forecast pipeline as a dependency DAG:
#   historical microclimate (job1)
#     |-> forecast-member microclimate array (job2) --\
#     \-> egg-model historical (job3)  ----------------+-> egg-model member array (job4) -> aggregate: stats + maps (job5)
# job2 and job3 both only need job1, so they run in parallel rather than
# chained one after the other.
# Usage: ./spartan_submit_all.sh  (run from this directory)

set -euo pipefail

job1=$(sbatch --parsable spartan_01_historical.slurm)
echo "Historical microclimate job: $job1"

job2=$(sbatch --parsable --dependency=afterok:$job1 spartan_02_forecast_array.slurm)
echo "Forecast-member microclimate array job: $job2"

job3=$(sbatch --parsable --dependency=afterok:$job1 spartan_03_eggmodel_historical.slurm)
echo "Egg-model historical job: $job3"

job4=$(sbatch --parsable --dependency=afterok:$job2:$job3 spartan_04_eggmodel_array.slurm)
echo "Egg-model member array job: $job4"

job5=$(sbatch --parsable --dependency=afterok:$job4 spartan_05_eggmodel_aggregate.slurm)
echo "Aggregate (stats + maps) job: $job5"

echo "Submitted: $job1 -> {$job2, $job3} -> $job4 -> $job5"
