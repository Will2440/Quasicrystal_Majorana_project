# Returns Dict{Symbol,Any}

project_root = @__DIR__
slice_dir = normpath(joinpath(project_root, ".."))
hpc_dir = normpath(joinpath(slice_dir, ".."))

Dict{Symbol,Any}(
    :run_name => "dynamic_mb_slice_search",
    :output_root => joinpath(slice_dir, "runs"),

    :sequence_bson_path => joinpath(hpc_dir, "serial_batch", "batch_params", "sequences", "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1.bson"),
    :target_slope => 2 / (1 + sqrt(5)),
    :slope_tolerance => 0.001,
    :target_phason => 0.0,
    :phason_tolerance => 1e-4,
    :sequence_selection_mode => :closest,

    :N => 500,
    :Delta => 0.05,
    :sweep_mode => :vary_t2,
    :t1_fixed => 1.0,
    :t2_fixed => 1.0,
    :rho_start => 1.0,
    :rho_max => 3.0,
    :n_rho_slices => 41,
    :rho_spacing_mode => :linear,
    :u_range => (0.0, 1.05),

    :x0_alpha => 0.125 / 2,
    :x0_count => 5,
    :kappas => [1e-4],
    :n_lowest_evals => 4,
    :kk_tol => 1e-8,
    :kk_maxiter => 300,
    :topological_threshold => 0.5,
    :topological_fraction => 0.5,

    # Full-range baseline at every rho is the discovery safeguard. Features
    # narrower than this sampling interval are not guaranteed to be detected.
    :mu_baseline_points => 201,
    :mu_refine_points_per_prediction => 9,
    :transition_mu_tol => 1e-4,
    :transition_max_refine_iters => 40,
    :gap_refine_tol => 1e-3,
    :gap_refine_mu_tol => 5e-4,
    :gap_refine_max_depth => 2,

    # Momentum is off until at least three positions are observed. Predictions
    # are damped and capped; baseline discovery remains active on every slice.
    :momentum_enabled => true,
    :momentum_min_history => 3,
    :momentum_damping => 0.2,
    :momentum_max_shift_fraction => 0.02,
    :prediction_band_half_width_fraction => 0.03,
    :track_match_max_fraction => 0.15,
    :track_max_missed_slices => 1,

    :rho_onset_tol => 0.01,
    :rho_onset_max_extra_slices => 40,
    :max_points_per_slice => 20_000,
    :max_total_points => 500_000,
    :n_eval_workers => 32,
    :progress_log_every_points => 25,
    :blas_threads => 1,

    :global_stop_file => joinpath(slice_dir, "STOP_REQUESTED"),
    :run_root_pointer_file => joinpath(slice_dir, "runs", "latest_run_root.txt"),
)