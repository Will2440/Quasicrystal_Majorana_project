# Dynamic Slice Search

This is a standalone, slice-oriented spectral-localizer search. It scans rho from a starting value, discovers phase transitions along each full mu slice, refines observed transition brackets, and uses conservative transition tracking to add samples near previously observed boundaries.

The original adaptive-cell implementations remain unchanged in `../dynamic_search` and `../dynamic_search_singlejob`.

## Search behavior

Each rho slice uses two sampling scales. A sparse full-range discovery grid remains active on every slice. `mu_baseline_points` is its count at `rho_start`; physical-span scaling preserves its starting physical-mu spacing, and `discovery_rho_growth_factor` adds a modest extra density toward `rho_max`. Predicted transition bands receive more samples, increasing from `mu_refine_points_per_prediction` to its configured maximum. Their padding is driven by observed edge motion and curvature, with configured minimum and maximum fractions of the mu span. Newly detected transitions become tracks and receive the same focused sampling on subsequent slices. Low-localizer-gap intervals still receive bounded midpoint refinement, and every observed phase-change bracket is bisected to `transition_mu_tol`.

Transition motion is only extrapolated after at least three observed slice positions. The extrapolation is damped and capped by `momentum_max_shift_fraction`; every slice still receives its full baseline scan. Predictions therefore focus additional work but never remove the discovery scan.

Transition IDs are assigned by an order-preserving sequence alignment in mu, using predicted positions and phase-transition orientation. The matcher can leave old edges unmatched or start new tracks, and its allowed match distance is based on a small baseline tolerance plus predicted movement, capped by `track_match_max_tolerance_fraction`. This prevents a newly appearing interior edge from taking the ID of an established outer edge. Midpoint probes for onset bracketing are evaluated out of sweep order and are deliberately saved with track ID zero rather than used to update motion histories.

After the configured rho grid, additional midpoint slices are evaluated around detected gap-count changes until their brackets are narrower than `rho_onset_tol` or the extra-slice budget is exhausted.

The sparse discovery scan cannot guarantee detection of a phase interval narrower than its spacing. Treat the per-slice `baseline_mu_spacing` in `slice_summary.csv` as the discovery-resolution record and test convergence by increasing `mu_baseline_points` or `discovery_rho_growth_factor`. The localizer-gap refinement is a supplemental signal, not proof that unsampled intervals contain no transition.

## Configuration

Edit `config/default_slice_config.jl` for the default run or pass another config file to the launcher.

- `rho_start`, `rho_max`, `n_rho_slices`: inclusive sweep range and number of initial slices.
- `rho_spacing_mode`: `:linear` or `:log` spacing in rho.
- `u_range`: normalized mu range; physical mu bounds are recalculated at each rho.
- `mu_baseline_points`: sparse full-range discovery count at `rho_start`, scaled with physical mu span.
- `discovery_rho_growth_factor`, `discovery_rho_growth_power`: gentle extra discovery-grid growth from `rho_start` to `rho_max`.
- `transition_mu_tol`: stopping width for each observed transition bracket in physical mu.
- `gap_refine_tol`, `gap_refine_mu_tol`, `gap_refine_max_depth`: bounded refinement for adjacent samples with small localizer gaps.
- `momentum_min_history`, `momentum_damping`, `momentum_max_shift_fraction`: conservative prediction controls.
- `prediction_band_min_fraction`, `prediction_band_max_fraction`, `prediction_band_motion_multiplier`, `prediction_band_curvature_multiplier`: data-informed padding around each predicted edge, bounded by the configured mu-span fractions. The default 3% minimum preserves the prior low-rho band width.
- `mu_refine_points_per_prediction`, `mu_refine_points_per_prediction_max`, `boundary_resolution_growth_power`: local boundary-band sample density from the start to the high-rho end.
- `track_match_base_tolerance_fraction`, `track_match_step_scale`, `track_match_max_tolerance_fraction`: prediction-centered association tolerance, scaled by the predicted edge movement and capped relative to the physical mu span.
- `track_unmatched_penalty`: sequence-alignment cost for leaving an old edge unmatched or treating a detected transition as new.
- `track_max_missed_slices`: how many sequential slices an unmatched track remains eligible for matching.
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

- `slice_summary.csv`: gap count, transition count, baseline and total sample counts, actual baseline mu spacing, number of predicted edges, local points per edge, and timing for each rho slice.
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