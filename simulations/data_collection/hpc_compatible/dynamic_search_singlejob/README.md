# Dynamic Search (Single-Job Variant)

This directory contains a clean single-job implementation of the adaptive spectral-localizer search.

## Why this variant exists

The legacy array-based implementation remains in:
- `../dynamic_search`

This single-job variant avoids repeated Slurm queue waits between iterations by keeping all adaptive evaluation inside one long-lived allocation.

## Design goals

- Safety: keep checkpoint and stop-file behavior.
- Simplicity: remove array-job orchestration and per-iteration `sbatch` complexity.
- Efficiency: evaluate candidate points in-process using Julia parallel workers.

## Main files

- `master_single.jl`: adaptive orchestrator + local parallel evaluation.
- `config/default_single_config.jl`: default N=500 configuration with `rho_range=(1.0, 5.0)`.
- `runscripts/single_master.sbatch`: one-job Slurm launcher for up to 64 CPU cores.

## Initial runtime assumptions

- Prior observed throughput guidance for N=500 was approximately 44 points/s.
- Config starts with `estimated_points_per_second = 44.0` for rough ETA logs.
- Actual achieved points/s is measured per iteration and written to outputs.

## Usage

From this directory:

```bash
sbatch --partition=compute --account=phys030424 runscripts/single_master.sbatch config/default_single_config.jl
```

## Outputs

- Slurm logs: `logs/single_master_<jobid>.out/.err`
- Run outputs: `runs/dynamic_mb_search_singlejob_<timestamp>/`

Each run directory includes:
- `iter_*/points/points_all.csv`
- `iter_*/results/worker_results_local.csv`
- `iter_*/cells_input.csv`, `cells_next.csv`, `cells_resolved.csv`
- `point_results_compact.csv`, `search_history.csv`, `run_checkpoint.jld2`, `run_summary.jld2`

## Stop behavior

Stop file support remains:
- local stop: `<run_root>/STOP_REQUESTED`
- global stop: `STOP_REQUESTED` in this project root
