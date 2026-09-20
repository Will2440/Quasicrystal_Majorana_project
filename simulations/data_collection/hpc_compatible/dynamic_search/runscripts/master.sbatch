#!/bin/bash
#SBATCH --job-name=dysearch_master
#SBATCH --output=logs/master_%j.out
#SBATCH --error=logs/master_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# Master script that orchestrates entire adaptive search and submits worker jobs.
# 
# USAGE:
#   sbatch --partition=<part> --account=<acct> --time=<time> --mem=<mem> \
#     runscripts/master.sbatch <config_path>
#
# The master will:
# 1. Run calibration on the same partition/account to estimate worker resources.
# 2. Submit worker array jobs with matching partition/account + calibrated time/mem.
# 3. Manage the full iterative adaptive search until completion.
#
# All worker jobs will inherit the partition/account from the master's Slurm allocation.
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: sbatch master.sbatch <config_path>"
  exit 1
fi

CONFIG_PATH="$1"
# Convert relative path to absolute
if [[ ! "$CONFIG_PATH" = /* ]]; then
  CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"
fi

# Derive DYN_DIR from CONFIG_PATH (config file is in DYN_DIR/config/)
# This is more robust than BASH_SOURCE on compute nodes after Slurm staging
DYN_DIR="$(dirname "$(dirname "$CONFIG_PATH")")"
RUNS_DIR="$DYN_DIR/runs"
GLOBAL_STOP_FILE="$DYN_DIR/STOP_REQUESTED"
ACTIVE_JOBS_FILE="$RUNS_DIR/active_worker_jobs.txt"
LOGS_DIR="$DYN_DIR/logs"

# Ensure logs directory exists for Slurm output redirection
mkdir -p "$LOGS_DIR"

# Temporarily disable set -u for sourcing /etc/profile (cluster-prompt.sh uses PS1)
set +u
source /etc/profile
set -u
module load languages/julia/1.10.3

# Capture Slurm-allocated partition and account for worker submissions.
# These are set by Slurm if the user supplied --partition and --account to sbatch.
MJOB_PARTITION="${SLURM_JOB_PARTITION:-}"
MJOB_ACCOUNT="${SLURM_JOB_ACCOUNT:-}"

julia --project="$DYN_DIR" "$DYN_DIR/setup_env.jl"

cleanup() {
  echo "[master.sbatch] Termination signal received. Requesting graceful stop."
  mkdir -p "$RUNS_DIR"
  touch "$GLOBAL_STOP_FILE"
  if [ -f "$ACTIVE_JOBS_FILE" ]; then
    while IFS= read -r jid; do
      if [ -n "$jid" ]; then
        echo "[master.sbatch] Cancelling active worker job: $jid"
        scancel "$jid" || true
      fi
    done < "$ACTIVE_JOBS_FILE"
  fi
  if [ -n "${MASTER_PID:-}" ]; then
    kill -TERM "$MASTER_PID" 2>/dev/null || true
    wait "$MASTER_PID" || true
  fi
  exit 0
}

trap cleanup TERM INT

rm -f "$GLOBAL_STOP_FILE"

echo "[master.sbatch] Master job allocated to partition=$MJOB_PARTITION, account=$MJOB_ACCOUNT"
echo "[master.sbatch] Using config file: $CONFIG_PATH"
echo "[master.sbatch] Passing to master.jl for worker submissions..."

julia --project="$DYN_DIR" "$DYN_DIR/master.jl" "$CONFIG_PATH" "$MJOB_PARTITION" "$MJOB_ACCOUNT" &
MASTER_PID=$!
wait "$MASTER_PID"
