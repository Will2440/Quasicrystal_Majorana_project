#!/bin/bash
#SBATCH --job-name=dysearch_single
#SBATCH --output=logs/single_master_%j.out
#SBATCH --error=logs/single_master_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mem=64G

set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: sbatch single_master.sbatch <config_path>"
  exit 1
fi

CONFIG_PATH="$1"
if [[ ! "$CONFIG_PATH" = /* ]]; then
  CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"
fi

# Derive DYN_DIR from CONFIG_PATH (config is in DYN_DIR/config/)
DYN_DIR="$(dirname "$(dirname "$CONFIG_PATH")")"
LOGS_DIR="$DYN_DIR/logs"
mkdir -p "$LOGS_DIR"

# Temporarily disable set -u while sourcing /etc/profile.
set +u
source /etc/profile
set -u
module load languages/julia/1.10.3

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

julia --project="$DYN_DIR" "$DYN_DIR/setup_env.jl"

echo "[single_master.sbatch] Using config file: $CONFIG_PATH"
echo "[single_master.sbatch] JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

julia --project="$DYN_DIR" "$DYN_DIR/master_single.jl" "$CONFIG_PATH"
