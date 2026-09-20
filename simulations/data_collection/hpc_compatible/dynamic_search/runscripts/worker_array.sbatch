#!/bin/bash
#SBATCH --job-name=dysearch_worker
#SBATCH --output=logs/worker_%A_%a.out
#SBATCH --error=logs/worker_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# Defaults below are intentionally conservative and can be overridden by sbatch CLI args.
# master.jl can pass calibrated --time/--mem and cluster routing args.
#SBATCH --time=02:00:00
#SBATCH --mem=4G

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: sbatch worker_array.sbatch <config_path> <iter_dir>"
  exit 1
fi

CONFIG_PATH="$1"
ITER_DIR="$2"
TASK_ID="${SLURM_ARRAY_TASK_ID}"

CHUNK_FILE=$(printf "%s/points/points_chunk_%04d.csv" "$ITER_DIR" "$TASK_ID")
RESULTS_DIR="$ITER_DIR/results"
LOGS_DIR="$ITER_DIR/logs"
OUT_FILE=$(printf "%s/worker_results_%04d.csv" "$RESULTS_DIR" "$TASK_ID")

# Derive DYN_DIR from CONFIG_PATH (more robust than BASH_SOURCE on compute nodes)
# CONFIG_PATH is DYN_DIR/config/first_test_config.jl, so go up 2 levels
DYN_DIR="$(dirname "$(dirname "$CONFIG_PATH")")"
ROOT_LOGS_DIR="$DYN_DIR/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$ROOT_LOGS_DIR"

if [ ! -f "$CHUNK_FILE" ]; then
  echo "Chunk file not found: $CHUNK_FILE"
  exit 2
fi

# Temporarily disable set -u for sourcing /etc/profile (cluster-prompt.sh uses PS1)
set +u
source /etc/profile
set -u
module load languages/julia/1.10.3

julia --project="$DYN_DIR" "$DYN_DIR/worker.jl" "$CONFIG_PATH" "$CHUNK_FILE" "$OUT_FILE"
