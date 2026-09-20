# BluePebble Slurm Compatibility Verification

## ✓ Single `sbatch` Command Works Correctly

The implementation now fully supports the standard Slurm workflow you're familiar with:

```bash
sbatch --partition=compute --account=myaccount --time=48:00:00 --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

### Step-by-Step Verification

**Step 1: Slurm Directive Passthrough**
- ✓ `master.sbatch` is a valid Slurm script with SBATCH directives
- ✓ Accepts `--partition`, `--account`, `--time`, `--mem` on the command line
- ✓ These override the #SBATCH defaults in the script

**Step 2: Environment Variable Capture**
- ✓ `master.sbatch` captures `SLURM_JOB_PARTITION` and `SLURM_JOB_ACCOUNT` (set by Slurm)
- ✓ Logs them for transparency: `[master.sbatch] Master job allocated to partition=compute, account=myaccount`
- ✓ Passes them to Julia as CLI arguments

**Step 3: Config Override in master.jl**
- ✓ `master.jl` receives partition/account as ARGS[2] and ARGS[3]
- ✓ Overrides config: `cfg[:slurm_partition] = cli_partition`
- ✓ Logs the override: `[timestamp] CLI override: slurm_partition=compute`
- ✓ All subsequent code uses these values

**Step 4: Worker Submission**
- ✓ When submitting calibration job, `build_common_sbatch_args(cfg)` includes partition/account
- ✓ When submitting worker arrays, same routing args are used
- ✓ Fallback values in config are only used if CLI args are empty

**Step 5: Output Format**
- ✓ All worker jobs receive: `--partition=compute --account=myaccount --time=<est> --mem=<est>`
- ✓ Format is standard Slurm sbatch compatible
- ✓ No special parsing or Slurm vendor-specific syntax

## Example Command for BluePebble

Replace with your actual account:

```bash
# First time: check available partitions
sinfo --Format=partitionname,statelong --all

# Then submit:
sbatch \
  --partition=compute \
  --account=uob-example-account \
  --time=48:00:00 \
  --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

## Files Modified for Slurm Compatibility

| File | Change |
|------|--------|
| `runscripts/master.sbatch` | Capture SLURM env vars, pass to master.jl |
| `master.jl` | Accept partition/account CLI args, override config |
| `README.md` | Document single-sbatch-command workflow |
| `BLUEPEBBLE_USAGE.md` | New: BluePebble-specific guide |
| `SLURM_FLOW.md` | New: Visual flow diagram |

## What Gets Passed Through

```
User sbatch command
  ├─ --partition=compute
  ├─ --account=myaccount
  └─ --time=48:00:00
       ↓
  [Slurm allocates job with these values]
       ↓
  master.sbatch detects SLURM_JOB_PARTITION, SLURM_JOB_ACCOUNT
       ↓
  Passes to: julia master.jl config.jl compute myaccount
       ↓
  master.jl overrides config:
    - cfg[:slurm_partition] = "compute"
    - cfg[:slurm_account] = "myaccount"
       ↓
  All subsequent job submissions (calibration + workers):
    sbatch --partition=compute --account=myaccount --time=<est> --mem=<est> ...
```

## Backward Compatibility

- ✓ If you omit `--partition` and `--account` from sbatch, the script still runs
- ✓ Falls back to config file values (if set) or empty strings (no routing)
- ✓ `prepare_only` mode still works for manual submission
- ✓ Existing configs continue to work unchanged

## Tested Scenarios

| Scenario | Behavior |
|----------|----------|
| `sbatch master.sbatch config.jl` (no routing args) | Uses config values; falls back to fallback defaults |
| `sbatch --partition=test master.sbatch config.jl` | Partition=test, account=empty or from config |
| `sbatch --partition=compute --account=proj master.sbatch config.jl` | ✓ Fully routed as expected |
| Master job cancelled while workers running | ✓ Trap catches TERM, cancels active worker jobs |
| Config has empty partition/account | ✓ No Slurm routing args injected (OK) |

## Key Safety Features Intact

1. **Graceful cancellation**: `scancel <master_job_id>` auto-cancels workers
2. **Early stop**: `touch dynamic_search/STOP_REQUESTED` exits cleanly
3. **Data persistence**: All results saved in `runs/<run>/` before exit
4. **Worker tracking**: Active job IDs in `runs/active_worker_jobs.txt`
5. **Auto calibration**: Runs on `test` partition by default (fast)
6. **Timeout handling**: Workers get calibrated time estimates + safety margins

## Recommended Usage Pattern for BluePebble

```bash
#!/bin/bash
# save as: submit_search.sh

PARTITION="compute"
ACCOUNT="your-account-name"
TIME="48:00:00"
MEM="8G"
CONFIG="config/default_config.jl"

sbatch \
  --partition=$PARTITION \
  --account=$ACCOUNT \
  --time=$TIME \
  --mem=$MEM \
  runscripts/master.sbatch $CONFIG

echo "Master job submitted. Check squeue for job status."
```

Usage:
```bash
chmod +x submit_search.sh
./submit_search.sh
```
