#!/bin/bash
#SBATCH --job-name=dyslice
#SBATCH --output=logs/slice_master_%j.out
#SBATCH --error=logs/slice_master_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00
#SBATCH --mem=64G

set -euo pipefail

if [[ "$#" -lt 1 ]]; then
  echo "Usage: sbatch slice_master.sbatch <config_path>"
  exit 1
fi

CONFIG_PATH="$1"
if [[ ! "$CONFIG_PATH" = /* ]]; then
  CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"
fi

SLICE_DIR="$(dirname "$(dirname "$CONFIG_PATH")")"
mkdir -p "$SLICE_DIR/logs"

set +u
source /etc/profile
set -u
module load languages/julia/1.10.3

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

julia --project="$SLICE_DIR" "$SLICE_DIR/setup_env.jl"
echo "[slice_master.sbatch] Config: $CONFIG_PATH"
echo "[slice_master.sbatch] JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
julia --project="$SLICE_DIR" "$SLICE_DIR/master_slice.jl" "$CONFIG_PATH"