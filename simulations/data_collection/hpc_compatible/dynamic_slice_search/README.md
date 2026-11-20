# Dynamic Slice Search

This is a standalone, slice-oriented spectral-localizer search. It sweeps either rho at fixed Delta or Delta at fixed hoppings, discovers phase transitions along each full mu slice, refines observed transition brackets, and uses conservative transition tracking to add samples near previously observed boundaries.

The original adaptive-cell implementations remain unchanged in `../dynamic_search` and `../dynamic_search_singlejob`.

## Search behavior

Each slice uses two sampling scales. A sparse full-range discovery grid remains active on every slice. `mu_baseline_points` is its count at the sweep start; physical-span scaling preserves its starting physical-mu spacing, and `discovery_rho_growth_factor` adds density toward the sweep endpoint, including toward decreasing Delta. Predicted transition bands receive more samples, increasing from `mu_refine_points_per_prediction` to its configured maximum. Their padding is driven by observed edge motion and curvature, with configured minimum and maximum fractions of the mu span. Newly detected transitions become tracks and receive the same focused sampling on subsequent slices. Low-localizer-gap intervals still receive bounded midpoint refinement, and every observed phase-change bracket is bisected to `transition_mu_tol`.

Transition motion is only extrapolated after at least three observed slice positions. The extrapolation is damped and capped by `momentum_max_shift_fraction`; every slice still receives its full baseline scan. Predictions therefore focus additional work but never remove the discovery scan.

Transition IDs are assigned by an order-preserving sequence alignment in mu, using predicted positions and phase-transition orientation. The matcher can leave old edges unmatched or start new tracks, and its allowed match distance is based on a small baseline tolerance plus predicted movement, capped by `track_match_max_tolerance_fraction`. This prevents a newly appearing interior edge from taking the ID of an established outer edge. Midpoint probes for onset bracketing are evaluated out of sweep order and are deliberately saved with track ID zero rather than used to update motion histories.

After the configured sweep grid, additional midpoint slices are evaluated around detected gap-count changes until their brackets meet the corresponding onset tolerance or the extra-slice budget is exhausted. Tracking and onset refinement use the selected sweep parameter and work in either Delta direction. Changes that occur and disappear between slices, or leave the endpoint gap counts equal, are not necessarily detected.

The sparse discovery scan cannot guarantee detection of a phase interval narrower than its spacing. Treat the per-slice `baseline_mu_spacing` in `slice_summary.csv` as the discovery-resolution record and test convergence by increasing `mu_baseline_points` or `discovery_rho_growth_factor`. The localizer-gap refinement is a supplemental signal, not proof that unsampled intervals contain no transition.

## Configuration

Edit `config/default_slice_config.jl` for the default run or pass another config file to the launcher.

- `sweep_parameter`: `:rho` or `:Delta`. Older configs without this key default to `:rho`; the provided default config selects `:rho`.
- For a rho sweep, `Delta` is fixed. `sweep_mode = :vary_t2` fixes `t1_fixed` and sets `t2 = rho * t1`; `:vary_t1` fixes `t2_fixed` and sets `t1 = t2 / rho`.
- `rho_start`, `rho_max`, `n_rho_slices`: inclusive sweep range and number of initial slices.
- `rho_spacing_mode`: `:linear` or `:log` spacing in rho. `:log` interpolates in log(rho) and exponentially warps the sample positions. It clusters toward the smaller-rho endpoint when the range's geometric midpoint is at or below 1, and toward the larger-rho endpoint when it is above 1. Thus ranges below 1 cluster toward low rho; ranges above 1 cluster toward high rho. For a range crossing 1, the geometric midpoint selects the biased endpoint.
- `rho_log_concentration`: positive concentration parameter for `:log` spacing (default 4.0); larger values cluster more strongly toward the selected endpoint. This setting has no effect with `:linear` spacing.
- For a Delta sweep, set both `t1_fixed` and `t2_fixed` to the desired hopping values. Both remain fixed, and `rho = t2_fixed / t1_fixed`; `sweep_mode` and the fixed `Delta` value are ignored.
- `Delta_start`, `Delta_max`, `n_Delta_slices`, `Delta_spacing_mode`: inclusive Delta sweep in either direction. Despite its historical name, `Delta_max` is the final endpoint, not necessarily the maximum. `:linear` uses equal Delta steps and permits zero. `:log` uses geometric spacing (equal steps in log(Delta)) and requires both endpoints to be positive. `:exponential` permits zero and concentrates samples toward the smaller Delta endpoint, whether the sweep is increasing or decreasing.
- `Delta_exponential_concentration`: positive exponent controlling `:exponential` spacing (default 4.0; accepted range `(0, 700]`). For endpoints 1 and 0, the grid is `expm1(c * x) / expm1(c)` with `x` evenly spaced from 1 down to 0. Larger values concentrate more samples near the smaller endpoint. This is a sampling choice, not a guarantee of resolving an exponential gap density.
- `u_range`: normalized mu range. In rho mode physical bounds are recalculated at each rho; in Delta mode they are calculated once for the selected sequence slope and fixed hoppings and reused on every slice.
- `kappa_mode`: `:user_defined` evaluates every entry in `kappas` at every `x0`; `:delta_over_L` uses one value `Delta / N`, where `N` is the number of chain sites (unit lattice spacing). Older configs without this setting use `:user_defined`. In rho sweeps, Delta is fixed, so `:delta_over_L` gives the same kappa on every slice. Onset probes use their own Delta.
- `kappas`: a nonempty collection of finite, nonnegative values for `:user_defined`; for example, `collect(range(1e-5, 1e-4; length=5))`. Every kappa and position contributes to the same phase vote and summary statistics, not independent kappa-resolved runs. The list is ignored in `:delta_over_L` mode.
- `mu_baseline_points`: full-range discovery samples on every slice.
- `transition_mu_tol`: stopping width for each observed transition bracket in physical mu.
- `gap_refine_tol`, `gap_refine_mu_tol`, `gap_refine_max_depth`: bounded refinement for adjacent samples with small localizer gaps.
- `momentum_min_history`, `momentum_damping`, `momentum_max_shift_fraction`: conservative prediction controls.
- `prediction_band_min_fraction`, `prediction_band_max_fraction`, `prediction_band_motion_multiplier`, `prediction_band_curvature_multiplier`: data-informed padding around each predicted edge, bounded by the configured mu-span fractions. The default 3% minimum preserves the prior low-rho band width.
- `mu_refine_points_per_prediction`, `mu_refine_points_per_prediction_max`, `boundary_resolution_growth_power`: local boundary-band sample density from the start to the sweep endpoint.
- `track_match_base_tolerance_fraction`, `track_match_step_scale`, `track_match_max_tolerance_fraction`: prediction-centered association tolerance, scaled by the predicted edge movement and capped relative to the physical mu span.
- `track_unmatched_penalty`: sequence-alignment cost for leaving an old edge unmatched or treating a detected transition as new.
- `track_max_missed_slices`: how many sequential slices an unmatched track remains eligible for matching.
- `rho_onset_tol`, `rho_onset_max_extra_slices`: gap-count change bracket width and total midpoint-probe budget.
- `Delta_onset_tol`, `Delta_onset_max_extra_slices`: corresponding controls for a Delta sweep, in Delta units.
- `max_points_per_slice`, `max_total_points`: point budgets; execution stops cleanly at a budget and saves completed work.
- `n_eval_workers`, `blas_threads`: Julia CPU parallelism controls.

All hopping, sequence-selection, localizer, and disorder-free physical parameters are also in the config. Transition tracks are evaluated in physical mu; stored rows include both `mu` and `u`.

The provided default config sweeps rho from 0.001 to 1 with `t1` fixed at 1, `Delta` fixed at 0.05, and linear rho spacing. For example, to use the logarithmic bias toward the low-rho endpoint, set:

```julia
:sweep_parameter => :rho,
:sweep_mode => :vary_t2,
:t1_fixed => 1.0,
:rho_start => 1e-3,
:rho_max => 1.0,
:n_rho_slices => 101,
:rho_spacing_mode => :log,
:rho_log_concentration => 4.0,
```

A Delta-sweep configuration can instead fix the hoppings and use exponential sampling toward smaller Delta:

```julia
:sweep_parameter => :Delta,
:t1_fixed => 1.0,
:t2_fixed => 1.5,
:Delta_start => 1.0,
:Delta_max => 0.001,
:n_Delta_slices => 21,
:Delta_spacing_mode => :exponential,
:Delta_exponential_concentration => 4.0,
:kappa_mode => :delta_over_L,
```

At `Delta = 0`, `:delta_over_L` gives `kappa = 0`, so the localizer's position-dependent term vanishes and it may become singular. Treat the zero endpoint with care; choose a positive endpoint or a user-defined positive kappa if a nonzero position term is required. There is no automatic kappa floor.

The fixed estimate is `mu_c = 2 * ((1 - gamma) * t1 + gamma * t2)`, where `gamma` is the slope of the selected sequence. This is a bandwidth-based estimate, not a numerically measured critical point or a guarantee that every Delta-dependent transition lies inside the chosen `u_range`. Check the physical bounds and expand `u_range` or run convergence checks when needed.

Near Delta zero, gaps may become too narrow in mu or too closely spaced in Delta for a finite scan. Increase `n_Delta_slices`, adjust the exponential concentration, and increase the mu discovery density to test convergence; neither the onset refinement nor the extra boundary samples guarantee discovery. More slices and tracks increase computational cost, and the point budgets may stop a run before it reaches zero. Inspect `stop_reason`, completed slices, and `baseline_mu_spacing` in the output.

## Run

From this directory on the cluster:

```bash
sbatch --partition=compute --account=phys030424 runscripts/slice_master.sbatch config/default_slice_config.jl
```

The job uses Julia 1.10.3, instantiates this standalone project, and runs all slices inside one allocation. Adjust Slurm time, memory, account, and CPU requests in `runscripts/slice_master.sbatch` for the target cluster and planned baseline resolution.

## Outputs

Each run is written under `runs/dynamic_mb_slice_search_<timestamp>/`:

- `slice_summary.csv`: gap count, transition count, sample counts, timing, physical mu bounds, mu-critical estimate, and kappa mode and values for each slice.
- `transition_boundaries.csv`: transition positions and mu brackets, with track IDs for the initial sweep.
- `onset_brackets.csv`: sampled intervals where the gap count changes and the bracket meets the configured tolerance; columns are `rho_lo`/`rho_hi` or `Delta_lo`/`Delta_hi`.
- `slices/slice_NNNN/point_results.csv`: all localizer results for that slice.
- `slices/slice_NNNN/predictions.csv`: prediction centers and sampling bands used on the slice.
- `run_checkpoint.jld2`, `run_summary.jld2`: cumulative checkpoint and final run metadata.
- `STOP_REQUESTED`: create this file inside the run directory to stop before the next slice.

The summary, transition, and point tables contain `sweep_parameter`, `sweep_value`, `rho`, and `Delta`. `rho` always means the physical hopping ratio, even in a Delta sweep; use `Delta` or `sweep_value` for its horizontal axis. The JLD2 metadata also records the sweep parameter and values. The scripts in `simulations/data_processing/boundary_processing/slice_boundary_processing` select the plotting axis from the sweep parameter and support legacy rho-only runs.

The summary and point tables record `kappa_mode` and `kappa_values`. Multiple kappa values are stored in one semicolon-separated field in the same order as the config; numerical results are aggregated across them and `x0` positions.

## Merge rho ranges

For rho sweeps only, combine consecutive runs while giving the second run priority throughout its observed rho interval:

```bash
julia merge_slice_runs.jl <first_run_directory> <priority_run_directory> [output_directory]
```

If no output directory is given, a timestamped `dynamic_mb_slice_search_merged_...` directory is created beside the source runs. First-run slices inside the priority run's inclusive `[min(rho), max(rho)]` range are discarded; all priority-run slices and first-run slices outside that interval are retained. The merger sorts and reindexes the summary rows, transition rows, and slice folders together. It writes the summary, transitions, and slice data required by the post-processing scripts, plus `merge_info.md`; source plots and restart checkpoints are not copied because they would not describe the merged run.

Run the synthetic algorithm tests with:

```bash
julia --project=. test/runtests.jl
julia --project=. test/integration.jl
```

The integration test uses a synthetic phase classifier to exercise complete runs in both modes and includes a small real-localizer smoke test.