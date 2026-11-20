#!/bin/bash
#SBATCH --job-name=boundary_array
#SBATCH --output=logs/boundary_array_%A_%a.out
#SBATCH --error=logs/boundary_array_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=0-3
#SBATCH --time=01:00:00
#SBATCH --mem=50G

set -euo pipefail

ARRAY_FIXED_PARAMETER="Delta"
ARRAY_VALUES=(0.01 0.02 0.05 0.1)

if [[ $# -gt 1 || -z ${SLURM_ARRAY_TASK_ID:-} || -z ${SLURM_ARRAY_JOB_ID:-} ]]; then
  echo "Usage: sbatch runscripts/boundary_array.sbatch [base_config_path]" >&2
  exit 1
fi
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
BOUNDARY_DIR="$(cd "${BOUNDARY_SEARCH_DIR:-$SUBMIT_DIR}" && pwd)"
BASE_CONFIG="${1:-$BOUNDARY_DIR/config/default_config.jl}"
if [[ "$BASE_CONFIG" != /* ]]; then
  BASE_CONFIG="$SUBMIT_DIR/$BASE_CONFIG"
fi
if [[ ! -f "$BASE_CONFIG" ]]; then
  echo "Configuration not found: $BASE_CONFIG" >&2
  exit 1
fi
if [[ ! $SLURM_ARRAY_TASK_ID =~ ^[0-9]+$ ]] || (( 10#$SLURM_ARRAY_TASK_ID >= ${#ARRAY_VALUES[@]} )); then
  echo "Array task index is outside ARRAY_VALUES: $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi
export BOUNDARY_ARRAY_BASE_CONFIG="$(realpath "$BASE_CONFIG")"
export BOUNDARY_ARRAY_PARAMETER="$ARRAY_FIXED_PARAMETER"
export BOUNDARY_ARRAY_VALUE="${ARRAY_VALUES[10#$SLURM_ARRAY_TASK_ID]}"
export BOUNDARY_SEARCH_DIR="$BOUNDARY_DIR"
echo "[boundary_array.sbatch] Task $SLURM_ARRAY_TASK_ID: $BOUNDARY_ARRAY_PARAMETER=$BOUNDARY_ARRAY_VALUE"
exec bash "$BOUNDARY_DIR/runscripts/boundary_master.sbatch" "$BOUNDARY_DIR/config/array_task_config.jl"