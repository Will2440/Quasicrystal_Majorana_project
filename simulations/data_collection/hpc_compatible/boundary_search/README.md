# Boundary search

This standalone, experimental program searches phase boundaries in slices of a rho or Delta sweep. It copies the localizer and sequence-selection routines and leaves `../dynamic_slice_search` unchanged.

## Coordinates and sweeps

The search coordinate is `u = mu / mu_critical`, where `mu_critical = 2 * ((1 - gamma) * t1 + gamma * t2)` is a bandwidth estimate, not a measured critical point. Tracking, prediction, padding and transition tolerances are in u; the localizer evaluates physical `mu = u * mu_critical`. `u_range` sets the search bounds and may extend above 1. Boundaries outside it cannot be detected.

The intended gap-appearance directions are increasing rho above 1, decreasing rho below 1 and decreasing Delta. Trivial regions and their boundaries may move in u; a fixed u is not assumed to remain trivial. The current predictor uses damped, capped linear extrapolation from recent motion, not a quadratic-well fit.

Set `sweep_parameter` to `:rho` or `:Delta`. Despite their historical names, `rho_max` and `Delta_max` are the final endpoints, not necessarily numerical maxima. Both endpoints must be finite and positive:

- For rho above 1: `1 <= rho_start <= rho_max`; for rho below 1: `0 < rho_max <= rho_start <= 1`. A sweep cannot cross 1.
- For Delta: `0 < Delta_max <= Delta_start`.
- Equal endpoints are allowed only with one slice. With one slice and unequal endpoints, the grid contains the start value only.
- Rho spacing is `:linear` or `:log`; the latter retains the legacy endpoint bias, controlled by `rho_log_concentration`. Delta spacing is `:linear`, `:log`, or `:exponential`; `Delta_exponential_concentration` controls the last mode. Concentrations must be finite and in `(0, 700]` when used.

## Configuration

The default [configuration](config/default_config.jl) returns a `Dict{Symbol,Any}` and resolves paths relative to this directory. The supplied sequence path points to the shared `serial_batch` data; ensure that file exists on the execution host. Custom configurations should likewise return a dictionary and resolve their own paths explicitly.

- `run_name`, `output_root`, `global_stop_file`, `run_root_pointer_file`: run naming, output and stop locations.
- `sequence_bson_path`, `target_slope`, `slope_tolerance`, `target_phason`, `phason_tolerance`, `sequence_selection_mode`: sequence input and selection.
- `N`, `sweep_parameter`, `sweep_mode`, `t1_fixed`, `t2_fixed`, `Delta`: system size and physics. For rho, `:vary_t2` fixes `t1` and `:vary_t1` fixes `t2`; for Delta, both hoppings are fixed.
- `rho_start`, `rho_max`, `n_rho_slices`, `rho_spacing_mode`, `rho_log_concentration`, `Delta_start`, `Delta_max`, `n_Delta_slices`, `Delta_spacing_mode`, `Delta_exponential_concentration`, `u_range`: sweep and search grids.
- `x0_alpha`, `x0_count`, `kappa_mode`, `kappas`, `n_lowest_evals`, `kk_tol`, `kk_maxiter`, `topological_threshold`, `topological_fraction`: localiser sampling, solver and phase classification. Kappa mode is `:user_defined` or `:delta_over_L`.
- `u_baseline_points`, `u_refine_points_per_prediction`, `transition_u_tol`, `transition_max_refine_iters`: baseline sampling target, samples per predicted boundary band (0 disables them; otherwise at least 3), and transition refinement in u. Audits, band points and bisection points are additional, so the baseline count is not a hard budget.
- `rho_onset_tol`, `rho_onset_max_extra_slices`, `Delta_onset_tol`, `Delta_onset_max_extra_slices`: onset resolution and global extra-slice budgets. Budget exhaustion stops frontier advancement with pending brackets recorded.
- `max_points_per_slice`, `max_total_points`: evaluation budgets.
- `confidence_min_history`, `confidence_max_residual`, `boundary_padding_min`, `boundary_padding_max`, `audit_points_per_interval`, `skip_trivial_interiors`: confidence, boundary padding and coverage controls. Residual and padding are absolute u distances. The default does not skip trivial interiors; enabling it requires at least one audit point per excluded interval.
- `momentum_enabled`, `momentum_min_history`, `momentum_damping`, `momentum_max_shift_fraction`, `track_match_base_tolerance_fraction`, `track_match_max_tolerance_fraction`, `track_match_step_scale`, `track_unmatched_penalty`, `track_max_missed_slices`: tracking prediction and matching controls.
- `n_eval_workers`, `progress_log_every_points`, `blas_threads`: execution and logging controls.

## Search method

The first slice scans the whole u range. Every later nominal grid slice uses predictions based on previously incorporated observations. If the count of phase-0 runs changes, midpoint probes immediately refine the interval between frontiers to the onset tolerance or global probe budget. Probes are evaluated out of order but incorporated into tracking in sweep order before the next frontier. Track IDs are assigned to probe boundaries as well as frontier boundaries.

The tracker matches ordered transitions by orientation and predicted u, allowing unmatched edges and new IDs. Confidence requires consecutive converged, unambiguous matches with residuals below `confidence_max_residual`; misses, ambiguity or large residuals reset it. Boundaries are bisected from actual phase-label changes. Predictions get two-sided sampling bands. A single `x0` is at the chain midpoint; multiple positions span `N/2 +/- x0_alpha*N`.

With `skip_trivial_interiors = true`, the padded interior between two confident predicted boundaries of an observed phase-0 interval may be excluded. The candidate interval is recomputed from moving predictions and clipped to the previously observed interval, not kept as a permanent mask. Each exclusion receives audit samples, and both bounding phase transitions must be found in their predicted bands. A phase-1 audit or failed boundary confirmation triggers a full-baseline rescan. First slices and onset probes always use full baseline scans. Point rows are actual evaluations, not inferred labels.

Skipping is experimental: an audit can miss a narrow topological region, and a confident extrapolation can still be wrong. Even a full baseline is finite and can miss narrow features. Equal endpoint gap counts can hide appearance, disappearance or merger events between slices. A resolved onset bracket is an observed gap-count change, not proof of the exact physical critical point. Compare against full-scan runs at increasing u and sweep resolutions before relying on skipping.

## Outputs

Each run gets a unique directory. `slice_summary.csv` is in sweep order and distinguishes frontier and onset-probe slices; execution IDs give their computation order. `transition_boundaries.csv` records track IDs and measured u and physical-mu brackets. `onset_brackets.csv` contains resolved and pending brackets. Each `slices/slice_NNNN/` contains actual `point_results.csv`, `predictions.csv` and proposed `excluded_intervals.csv`, including fallback reasons. JLD2 checkpoint and summary files contain config, tracks, counts, pending brackets and stop reason. Completed work and partial slice point files are saved on stop or budget exhaustion; automatic restart is not implemented.

## Launch and tests

Submit from this directory with `sbatch --partition=compute --account=phys030424 runscripts/boundary_master.sbatch config/default_config.jl`, adjusting account, partition, time, memory and CPUs for the system. The launcher loads Julia 1.10.3 and instantiates this project. Relative config paths are resolved from the submission directory; for submission from elsewhere, set `BOUNDARY_SEARCH_DIR` to this directory's absolute path. Slurm logs are written in the submission directory.

Run the lightweight grid, coverage, tracking and integration tests, including synthetic sweeps in all three directions and a tiny real-localizer smoke test:

```bash
module load languages/julia/1.10.3
julia --project=. test/sweep.jl
julia --project=. test/coverage.jl
julia --project=. test/tracking.jl
JULIA_NUM_THREADS=2 julia --project=. test/integration.jl
```