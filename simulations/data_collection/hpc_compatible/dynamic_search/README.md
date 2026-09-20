# Dynamic Adaptive Search (HPC)

This directory is a standalone adaptive search pipeline for the `(mu, rho)` phase map.
It leaves `serial_batch` untouched and reuses the same spectral-localiser physics kernels.

## Design Goals

- Keep spectral-localiser computations unchanged from the validated implementation.
- Avoid uniform high-resolution gridding of the full domain.
- Refine only near mixed/near-critical regions.
- Track termination with physically meaningful metrics, including trivial-gap count on the largest-rho cut.

## Quick Start (HPC Production)

For a full autonomous adaptive search on BluePebble or similar Slurm clusters:

```bash
sbatch --partition=compute --account=<your_account> --time=48:00:00 --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

The master will:
1. Auto-calibrate worker resources using the same partition/account.
2. Run the full iterative adaptive search, submitting worker arrays for each iteration.
3. All workers inherit the partition/account and calibrated time/memory from the master.
4. Outputs and checkpoints saved in `runs/<run_name>_<timestamp>/`.

No manual worker job submissions needed—everything is orchestrated end-to-end.

## Directory Layout

- `core/physics.jl`: unchanged localiser compute kernels copied from `serial_batch/solvers.jl`.
- `core/sequences.jl`: load/select sequence by slope/phason from the usual BSON sequence files.
- `core/bandwidth.jl`: `W_gamma ≈ 4 tbar_gamma` with `tbar_gamma = (1-gamma)t1 + gamma t2`, and `|mu_c| ≈ W_gamma/2`.
- `core/adaptive.jl`: adaptive mesh scheduling and refinement logic.
- `core/io_utils.jl`: chunk IO and config loading helpers.
- `worker.jl`: evaluates one chunk of candidate points.
- `master.jl`: iterative coordinator for adaptive search.
- `config/default_config.jl`: run-time settings.
- `runscripts/worker_array.sbatch`: Slurm array script for workers.

## Execution Model

1. `master.jl` builds coarse cells in `(rho, u)` where `u = mu / mu_c_est`.
2. It schedules corner+center samples only for unresolved cells.
3. Workers evaluate spectral localiser invariants/gaps at scheduled points.
4. Master refines only cells with mixed phase labels or low localiser gaps.
5. Iteration terminates by one or more criteria:
   - no new refinements,
   - max points budget,
   - stable trivial-gap count at `rho_max` over recent iterations,
   - optional expected-gap cap.

## Sequence Selection

Sequence extraction uses the same BSON source style as `serial_batch`:

- `:sequence_bson_path`
- `:target_slope`, `:slope_tolerance`
- `:target_phason`, `:phason_tolerance`
- `:sequence_selection_mode` (`:closest` or `:first`)

## Sweep Modes

- `:vary_t2` (default): `t1` fixed, `t2 = rho * t1`
- `:vary_t1`: `t2` fixed, `t1 = t2 / rho`

Both use the same adaptive framework.

## Quick Start

Edit `config/default_config.jl`, then:

```bash
cd simulations/data_collection/hpc_compatible/dynamic_search
julia --project=. setup_env.jl
julia --project=. master.jl config/default_config.jl
```

Default mode is `:prepare_only` so it creates iteration chunks and prints the `sbatch` command.

### Slurm run for one prepared iteration

```bash
sbatch --array=1-N runscripts/worker_array.sbatch config/default_config.jl <iter_dir>
```

Where `N` is the number of generated chunks and `<iter_dir>` is printed by the master.

## Fully Iterative HPC Mode

Set in config:

- `:execution_mode => :submit_and_wait`

Then `master.jl` will submit each iteration as an array job and wait for completion before refining.

### Single `sbatch` Command (Recommended)

Run the entire adaptive search in a single Slurm submission:

```bash
sbatch --partition=<partition> --account=<account> --time=12:00:00 --mem=8G \
  runscripts/master.sbatch config/default_config.jl
```

Where:
- `--partition` and `--account` are passed to the master job and flow through to all worker submissions.
- Master internally runs calibration on the same partition/account to estimate per-chunk resources.
- All worker arrays are submitted with matching partition/account and calibrated time/mem.
- No separate manual `sbatch` commands needed.

### For convenience:

```bash
./runscripts/submit_master.sh config/default_config.jl <partition> <account>
```

This wrapper constructs the full command for you.

## Auto Calibration

In `:submit_and_wait` mode, master automatically estimates worker resources on startup:

1. **Calibration job** runs a tiny test chunk (`:calibration_points`) to measure:
   - Elapsed runtime per point
   - Peak memory (MaxRSS)

2. **Scaling**: these metrics are scaled to your full `:points_per_job` using safety factors and buffers.

3. **Result**: all worker jobs in the run use `--time` and `--mem` computed from calibration.

4. **Routing**: calibration runs on `:calibration_partition` (default `"test"`, can be overridden by Slurm allocation).

### Key Config Knobs

- `:auto_calibrate_worker_resources` - enable/disable calibration
- `:calibration_points` - how many points to evaluate in the calibration chunk (smaller = faster calibration)
- `:calibration_partition`, `:calibration_account` - routing for calibration job (if empty, uses main `:slurm_partition`/`:slurm_account`)
- `:calibration_time_seed`, `:calibration_mem_seed_mb` - resource requests for the calibration job itself
- `:calibration_time_safety_factor`, `:calibration_mem_safety_factor` - scale-up factors for robustness (default 1.8x time, 1.6x memory)
- `:calibration_time_buffer_seconds`, `:calibration_mem_buffer_mb` - fixed buffers added to estimates
- `:min_worker_time_seconds`, `:min_worker_mem_mb` - floor values to avoid underestimating
- `:fallback_worker_time`, `:fallback_worker_mem_mb` - used if calibration fails or is disabled
- `:slurm_partition`, `:slurm_account` - routing for actual worker jobs (CLI args override config)

## Safety Features

- Hard worker concurrency cap via config: `:hard_max_workers => 199`.
- Adjustable cap below the hard limit: `:max_workers_per_array`.
- Master submits Slurm arrays as `--array=1-N%M`, where `M = min(N, hard_max_workers, max_workers_per_array)`.
- Optional strict completeness check: `:require_all_worker_results`.
- Worker resources are no longer fixed: master injects calibrated/fallback `--time` and `--mem` into worker submissions.

## Graceful Termination

- Cooperative early stop: create either of these files while the run is active:
  - `dynamic_search/STOP_REQUESTED`
  - `<run_root>/STOP_REQUESTED`
- Master detects the stop request, cancels active worker array jobs, writes checkpoints, and exits cleanly.
- `runscripts/master.sbatch` now traps `TERM/INT` and attempts to cancel active workers from `runs/active_worker_jobs.txt`.

## Progress and ETA Reporting

- Master Slurm output now reports:
  - iteration sizing (`cells`, `new_points`, `chunks`, `max_parallel_workers`),
  - worker queue state (`done/running/pending`),
  - elapsed time and ETA (when runtime estimate is available),
  - stop/cancel/checkpoint events.
- Worker Slurm output now reports per-chunk progress and ETA every `:worker_log_every` points.
- Master learns chunk runtime from Slurm accounting (`sacct`) and updates ETA estimates online.

## Outputs

Inside run folder:

- `search_history.csv`: per-iteration summary
- `point_results_compact.csv`: all evaluated point statistics
- `run_checkpoint.jld2`: checkpoint snapshot (history + point tables + gap history)
- `run_summary.jld2`: final run metadata and stop reason
- `iter_XXX/points/*.csv`: scheduled points/chunks
- `iter_XXX/results/worker_results_*.csv`: worker outputs
- `iter_XXX/cells_resolved.csv`: classified cells
- `iter_XXX/cells_next.csv`: next-iteration mesh

## Notes

- This pipeline is intentionally independent of `serial_batch` orchestration.
- It reuses localiser compute kernels but changes scheduler/search logic.
- If you want tighter/looser refinement, tune:
  - `:gap_refine_tol`, `:rho_cell_tol`, `:u_cell_tol`, `:max_depth`.
- For gap counting on the largest-rho cut, tune:
  - `:rho_max_probe_points`, `:min_gap_run_points`, `:gap_count_stable_iters`.
- For HPC safety and throughput, tune:
  - `:points_per_job`, `:max_workers_per_array`, `:poll_seconds`.
