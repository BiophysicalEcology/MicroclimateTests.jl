#!/bin/bash
# Re-runs just the egg-model stage (historical -> member array -> aggregate)
# against already-cached microclimate .nc files, skipping spartan_01/02
# entirely -- for when only the egg model itself changed (e.g. a parameter
# or bugfix) and the microclimate solve is still valid and doesn't need
# re-running.
#
# IMPORTANT: delete the *stale* egg-result caches first (egg_dir's
# *_egg_n*.jls and *_eggout_n*.jls files) -- these are separate from the
# microclimate cache and won't be recomputed on their own if left in place,
# since solve_batched-style caching means "file exists" short-circuits
# straight past any code changes.
#
# Same env-var overrides as spartan_submit_all.sh -- e.g.:
#   N_ENSEMBLES=50 ./spartan_resimulate_eggs.sh
#
# Usage: ./spartan_resimulate_eggs.sh  (run from this directory)

set -euo pipefail

N_ENSEMBLES=${N_ENSEMBLES:-99}
MEMBERS_PER_TASK_EGG=${MEMBERS_PER_TASK_EGG:-5}

(( N_ENSEMBLES >= 1 && N_ENSEMBLES <= 99 )) || { echo "N_ENSEMBLES=$N_ENSEMBLES out of bounds 1:99 -- ACCESS-S2 only has 99 members" >&2; exit 1; }

export LOCUST_N_ENSEMBLES=$N_ENSEMBLES

n_egg_tasks=$(( (N_ENSEMBLES + MEMBERS_PER_TASK_EGG - 1) / MEMBERS_PER_TASK_EGG ))

job3=$(sbatch --parsable spartan_03_eggmodel_historical.slurm)
echo "Egg-model historical job: $job3"

job4=$(MEMBERS_PER_TASK_EGG=$MEMBERS_PER_TASK_EGG sbatch --parsable --dependency=afterok:$job3 \
    --array=0-$((n_egg_tasks - 1)) spartan_04_eggmodel_array.slurm)
echo "Egg-model member array job: $job4 ($n_egg_tasks tasks, $MEMBERS_PER_TASK_EGG members/task)"

job5=$(sbatch --parsable --dependency=afterok:$job4 spartan_05_eggmodel_aggregate.slurm)
echo "Aggregate (stats + maps) job: $job5"

echo "Submitted: $job3 -> $job4 -> $job5"
