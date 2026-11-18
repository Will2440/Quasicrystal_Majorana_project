#!/bin/bash
#SBATCH --job-name=boundary_search
#SBATCH --output=boundary_search_%j.out
#SBATCH --error=boundary_search_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --mem=64G

set -euo pipefail

if [[ $# -gt 1 ]]; then
  echo "Usage: sbatch runscripts/boundary_master.sbatch [config_path]" >&2
  exit 1
fi

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
BOUNDARY_DIR="${BOUNDARY_SEARCH_DIR:-$SUBMIT_DIR}"
BOUNDARY_DIR="$(cd "$BOUNDARY_DIR" && pwd)"
CONFIG_PATH="${1:-$BOUNDARY_DIR/config/default_config.jl}"
if [[ "$CONFIG_PATH" != /* ]]; then
  CONFIG_PATH="$SUBMIT_DIR/$CONFIG_PATH"
fi

if [[ ! -f "$BOUNDARY_DIR/Project.toml" || ! -f "$BOUNDARY_DIR/master.jl" ]]; then
  echo "Expected Project.toml and master.jl in $BOUNDARY_DIR; submit from boundary_search or set BOUNDARY_SEARCH_DIR." >&2
  exit 1
fi
if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Configuration not found: $CONFIG_PATH" >&2
  exit 1
fi

set +u
source /etc/profile
set -u
module load languages/julia/1.10.3

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-24}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

julia --project="$BOUNDARY_DIR" "$BOUNDARY_DIR/setup_env.jl"
echo "[boundary_master.sbatch] Config: $CONFIG_PATH"
echo "[boundary_master.sbatch] JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
julia --project="$BOUNDARY_DIR" "$BOUNDARY_DIR/master.jl" "$CONFIG_PATH"