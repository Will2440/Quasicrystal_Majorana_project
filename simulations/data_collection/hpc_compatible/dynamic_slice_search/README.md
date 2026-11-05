# Dynamic Slice Search

This is a standalone, slice-oriented spectral-localizer search. It scans rho from a starting value, discovers phase transitions along each full mu slice, refines observed transition brackets, and uses conservative transition tracking to add samples near previously observed boundaries.

The original adaptive-cell implementations remain unchanged in `../dynamic_search` and `../dynamic_search_singlejob`.

## Search behavior

Each rho slice receives a full coarse baseline scan over the configured normalized chemical-potential interval. Predicted transition bands add local samples, low-localizer-gap intervals receive bounded midpoint refinement, and every observed phase-change bracket is bisected to `transition_mu_tol`.

Transition motion is only extrapolated after at least three observed slice positions. The extrapolation is damped and capped by `momentum_max_shift_fraction`; every slice still receives its full baseline scan. Predictions therefore focus additional work but never remove the discovery scan.

After the configured rho grid, additional midpoint slices are evaluated around detected gap-count changes until their brackets are narrower than `rho_onset_tol` or the extra-slice budget is exhausted.

The baseline scan cannot guarantee detection of a phase interval narrower than its spacing. Treat `mu_baseline_points` as a declared discovery resolution and test convergence by increasing it. The localizer-gap refinement is a supplemental signal, not a proof that unsampled intervals contain no transition.

## Configuration

Edit `config/default_slice_config.jl` for the default run or pass another config file to the launcher.

- `rho_start`, `rho_max`, `n_rho_slices`: inclusive sweep range and number of initial slices.
- `rho_spacing_mode`: `:linear` or `:log` spacing in rho.
- `u_range`: normalized mu range; physical mu bounds are recalculated at each rho.
- `mu_baseline_points`: full-range discovery samples on every slice.
- `transition_mu_tol`: stopping width for each observed transition bracket in physical mu.
- `gap_refine_tol`, `gap_refine_mu_tol`, `gap_refine_max_depth`: bounded refinement for adjacent samples with small localizer gaps.
- `momentum_min_history`, `momentum_damping`, `momentum_max_shift_fraction`: conservative prediction controls.
- `prediction_band_half_width_fraction`, `mu_refine_points_per_prediction`: extra sample bands around tracked edges.
- `track_match_max_fraction`, `track_max_missed_slices`: edge matching and temporary disappearance behavior.
- `rho_onset_tol`, `rho_onset_max_extra_slices`: gap-count change bracket width and total midpoint-probe budget.
- `max_points_per_slice`, `max_total_points`: point budgets; execution stops cleanly at a budget and saves completed work.
- `n_eval_workers`, `blas_threads`: Julia CPU parallelism controls.

All hopping, sequence-selection, localizer, and disorder-free physical parameters are also in the config. Transition tracks are evaluated in physical mu; stored rows include both `mu` and `u`.

## Run

From this directory on the cluster:

```bash
sbatch --partition=compute --account=phys030424 runscripts/slice_master.sbatch config/default_slice_config.jl
```

The job uses Julia 1.10.3, instantiates this standalone project, and runs all slices inside one allocation. Adjust Slurm time, memory, account, and CPU requests in `runscripts/slice_master.sbatch` for the target cluster and planned baseline resolution.

## Outputs

Each run is written under `runs/dynamic_mb_slice_search_<timestamp>/`:

- `slice_summary.csv`: gap count, transition count, sample counts, and timing for each rho slice.
- `transition_boundaries.csv`: transition positions and mu brackets, with track IDs for the initial sweep.
- `onset_brackets.csv`: sampled rho intervals where the gap count changes and the bracket meets the configured tolerance.
- `slices/slice_NNNN/point_results.csv`: all localizer results for that slice.
- `slices/slice_NNNN/predictions.csv`: prediction centers and sampling bands used on the slice.
- `run_checkpoint.jld2`, `run_summary.jld2`: cumulative checkpoint and final run metadata.
- `STOP_REQUESTED`: create this file inside the run directory to stop before the next slice.

Run the synthetic algorithm tests with:

```bash
julia --project=. test/runtests.jl
```