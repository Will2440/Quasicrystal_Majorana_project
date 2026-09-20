# BluePebble Usage Guide

This document shows how to run the dynamic adaptive search on BluePebble using standard Slurm `sbatch` commands.

## Single `sbatch` Command Workflow

Submit the entire adaptive search pipeline with a single command:

```bash
sbatch \
  --partition=compute \
  --account=your_account_name \
  --time=48:00:00 \
  --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

### What Happens

1. **Master job** is allocated on `compute` partition with `your_account_name`.
2. **Calibration stage** (if enabled):
   - Runs a tiny test worker on the same partition/account.
   - Measures runtime and memory from actual worker code.
   - Estimates per-chunk resources (time + memory) for full run.
3. **Main loop**:
   - Each iteration generates candidate points.
   - Master submits worker array jobs with:
     - Same partition/account as master.
     - Calibrated `--time` and `--mem` per worker task.
   - Master waits for workers and processes results.
   - Refines mesh and repeats.
4. **Graceful completion** or early stop preserves all data.

### Customizing for Your Account/Partition

Replace placeholders:

```bash
sbatch \
  --partition=PARTITION_NAME \
  --account=ACCOUNT_NAME \
  --time=HH:MM:SS \
  --mem=XG \
  runscripts/master.sbatch config/default_config.jl
```

Standard BluePebble recommendations:
- `--partition`: typically `compute`, `gpu`, or `test` depending on your needs.
- `--account`: your project account (use `sacctmgr show assoc user=$USER | grep account` to list).
- `--time`: total walltime for the master job (should cover all iterations + calibration).
- `--mem`: memory for the master process itself (not worker memory, which is auto-calibrated).

## Disabling or Tuning Calibration

### Disable auto-calibration (use fallback values)

Edit `config/default_config.jl`:

```julia
:auto_calibrate_worker_resources => false,
:fallback_worker_time => "02:00:00",
:fallback_worker_mem_mb => 4096,
```

### Use a different partition for calibration

By default, calibration runs on the `test` partition to save time. To use your main partition:

```julia
:calibration_partition => "",  # empty string = use :slurm_partition
```

### Adjust calibration safety factors

If workers are timing out or running out of memory, increase safety factors:

```julia
:calibration_time_safety_factor => 2.5,   # more conservative
:calibration_mem_safety_factor => 2.0,
```

## Monitor Progress

While the master job is running, check:

```bash
# Current job status
squeue -u $USER -j <job_id>

# Master and worker Slurm logs
ls runs/*/master_*.out runs/*/iter_*/worker_*.out

# Real-time output from master
tail -f runs/*/master_*.out
```

## Early Stop / Graceful Cancel

While the master job is running, gracefully stop it:

```bash
touch dynamic_search/STOP_REQUESTED
```

Or cancel the master job (which auto-cancels active workers):

```bash
scancel <master_job_id>
```

All collected data is saved and can be recovered from `runs/<run_name>_<timestamp>/`.

## Troubleshooting

### "Could not parse sbatch output"
- Ensure Slurm is available and `sbatch` works: `sbatch --version`
- Check master Slurm output for the full error: `cat runs/*/master_*.err`

### Workers timing out or running out of memory
- Increase safety factors in config
- Re-run with explicit `--time` and `--mem` to override estimates
- Check calibration output: `cat runs/<run_name>/calibration/calibration_summary.jld2`

### Calibration fails or is skipped
- Check if `:auto_calibrate_worker_resources => true` in config
- Look for error message in master output about calibration
- Falls back to `:fallback_worker_time` and `:fallback_worker_mem_mb`

### Different partition for calibration vs workers
- Calibration defaults to `test` partition (fast, short-time).
- Worker jobs use the same partition as the master job.
- Adjust `:calibration_partition` if your test partition is unavailable.
