#!/usr/bin/env bash
#PBS -N mol-sim
#PBS -l select=1:ncpus=4:mem=6gb
#PBS -l walltime=02:00:00
# Set -J 0-(task_count-1) at submission, after inspecting plan.json.
set -euo pipefail
: "${MST_PLAN:?Export MST_PLAN as the absolute shared plan.json path}"
: "${PBS_ARRAY_INDEX:?Run as a PBS Pro array task}"
# Activate the same environment used by the coordinator before this command.
mol-sim task "$MST_PLAN" "$PBS_ARRAY_INDEX" --workers 4
