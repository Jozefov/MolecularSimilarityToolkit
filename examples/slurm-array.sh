#!/usr/bin/env bash
#SBATCH --job-name=mol-sim
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=02:00:00
# Set --array=0-(task_count-1) at submission, after inspecting plan.json.
set -euo pipefail
: "${MST_PLAN:?Export MST_PLAN as the absolute shared plan.json path}"
: "${SLURM_ARRAY_TASK_ID:?Run as a Slurm array task}"
mol-sim task "$MST_PLAN" "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
