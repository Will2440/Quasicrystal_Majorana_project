# Single `sbatch` Command Integration Flow

## User Command
```bash
sbatch --partition=compute --account=myaccount --time=48:00:00 --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

## What Happens (Step by Step)

### 1. Slurm Accepts the Command
- Slurm allocates a master job on `compute` partition with `myaccount`.
- Sets environment variables: `SLURM_JOB_PARTITION=compute`, `SLURM_JOB_ACCOUNT=myaccount`.

### 2. master.sbatch Runs
- Sources module environment and loads Julia 1.10.3.
- Captures the Slurm environment variables:
  ```bash
  MJOB_PARTITION="${SLURM_JOB_PARTITION:-}"  # = "compute"
  MJOB_ACCOUNT="${SLURM_JOB_ACCOUNT:-}"      # = "myaccount"
  ```
- Passes them as CLI arguments to master.jl:
  ```bash
  julia master.jl config/default_config.jl compute myaccount
  ```

### 3. master.jl Receives and Applies Overrides
- Parses CLI args: `cli_partition="compute"`, `cli_account="myaccount"`
- Overrides config values:
  ```julia
  cfg[:slurm_partition] = "compute"
  cfg[:slurm_account] = "myaccount"
  ```
- Logs: `[timestamp] CLI override: slurm_partition=compute`
- Logs: `[timestamp] CLI override: slurm_account=myaccount`

### 4. Auto Calibration Stage (if enabled)
- Builds a tiny test chunk (default 3 points).
- Submits ONE calibration worker with:
  ```bash
  sbatch --array=1-1%1 --partition=test --account="" --time=00:30:00 --mem=4096M worker_array.sbatch
  ```
  (Note: uses `:calibration_partition=test` by default to save time.)
- Waits for calibration to complete.
- Reads Slurm accounting to get measured runtime and peak memory.
- Computes estimated per-chunk time and memory with safety factors.

### 5. Main Adaptive Loop
- Iteration 1:
  - Generates candidate points.
  - Submits worker array:
    ```bash
    sbatch --array=1-N%M --partition=compute --account=myaccount --time=<calibrated> --mem=<calibrated> worker_array.sbatch
    ```
  - Waits for workers to complete.
  - Processes results and refines mesh.
  
- Iterations 2, 3, ...: Same pattern until stopping condition.

### 6. Output and Checkpoints
- All data saved in: `runs/dynamic_mb_search_<timestamp>/`
- Includes: `point_results_compact.csv`, `search_history.csv`, `run_summary.jld2`
- Calibration summary at: `calibration/calibration_summary.jld2`

## Key Points for BluePebble Compatibility

✓ **Single `sbatch` command** (no multi-step submission)
✓ **Partition/account flow-through** (auto-detected from Slurm allocation)
✓ **Worker jobs use same partition/account** as master
✓ **Auto resource estimation** (no manual time/mem tuning for workers)
✓ **Graceful stop support** (cancels active workers on master termination)
✓ **All output in single run directory** (easy to find and archive)

## Slurm Directives Timeline

| Stage | Partition | Account | Time | Memory |
|-------|-----------|---------|------|--------|
| Master | user's arg | user's arg | user's arg | user's arg |
| Calibration | `test` (or override) | empty or override | `00:30:00` | `4096M` |
| Workers (iter 1-N) | compute (from master) | myaccount (from master) | ~calibrated | ~calibrated |

All worker job IDs are tracked in `runs/active_worker_jobs.txt` for safe cancellation if master is terminated.
