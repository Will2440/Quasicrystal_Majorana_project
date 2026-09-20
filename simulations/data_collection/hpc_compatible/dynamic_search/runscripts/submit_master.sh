#!/bin/bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: ./runscripts/submit_master.sh <config_path> [partition] [account]"
  exit 1
fi

CONFIG_PATH="$1"
PARTITION="${2:-}"
ACCOUNT="${3:-}"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
MASTER_SCRIPT="$SCRIPT_DIR/master.sbatch"

CMD=(sbatch)
if [ -n "$PARTITION" ]; then
  CMD+=(--partition "$PARTITION")
fi
if [ -n "$ACCOUNT" ]; then
  CMD+=(--account "$ACCOUNT")
fi
CMD+=("$MASTER_SCRIPT" "$CONFIG_PATH")

echo "Submitting: ${CMD[*]}"
"${CMD[@]}"
