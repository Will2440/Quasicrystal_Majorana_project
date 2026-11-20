#!/bin/bash
#SBATCH --job-name=boundary_kappa_window
#SBATCH --output=logs/boundary_kappa_window_%j.out
#SBATCH --error=logs/boundary_kappa_window_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --mem=8G
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
set +u
source /etc/profile
set -u
module load languages/julia/1.10.3
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
julia --startup-file=no --project=. optimisation_notes/kappa_map.jl fine