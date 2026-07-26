#!/bin/bash
# Submits the full points_australia_forecast pipeline as a dependency chain:
# historical -> forecast member array -> egg model + maps.
# Usage: ./spartan_submit_all.sh  (run from this directory)

set -euo pipefail

job1=$(sbatch --parsable spartan_01_historical.slurm)
echo "Historical job: $job1"

job2=$(sbatch --parsable --dependency=afterok:$job1 spartan_02_forecast_array.slurm)
echo "Forecast member array job: $job2"

job3=$(sbatch --parsable --dependency=afterok:$job2 spartan_03_eggmodel.slurm)
echo "Egg-model + maps job: $job3"

echo "Submitted chain: $job1 -> $job2 -> $job3"
