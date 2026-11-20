# Returns Dict{Symbol,Any}

project_root = @__DIR__
slice_dir = normpath(joinpath(project_root, ".."))
hpc_dir = normpath(joinpath(slice_dir, ".."))

Dict{Symbol,Any}(
    :run_name => "dynamic_mb_slice_search",
    :output_root => joinpath(slice_dir, "runs"),

    :sequence_bson_path => joinpath(hpc_dir, "serial_batch", "batch_params", "sequences", "literal", "slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1).bson"),
    :target_slope => 2 / (1 + sqrt(5)), #1 / (1 + sqrt(2)),
    :slope_tolerance => 0.001,
    :target_phason => 2 / (1 + sqrt(5)), #1 / (1 + sqrt(2)),
    :phason_tolerance => 1e-4,
    :sequence_selection_mode => :closest,

    :N => 500,
    :sweep_parameter => :rho,           # :rho or :Delta
    :Delta => 0.05,                     # fixed value of Delta when :sweep_parameter=:rho
    :sweep_mode => :vary_t2,            # :vary_t2 or :vary_t1
    :t1_fixed => 1.0,                   # fixed value of t1 when :sweep_parameter=:rho & :sweep_mode=:vary_t2
    :t2_fixed =>2.0,                    # fixed value of t2 when :sweep_parameter=:rho & :sweep_mode=:vary_t1
    :rho_start => 1.0,
    :rho_max => 2.0,
    :n_rho_slices => 11,
    :rho_spacing_mode => :linear,       # spacing mode for rho slices: :linear or :log
    :rho_log_concentration => 4.0,       # concentration for :log rho spacing
    :Delta_start => 1.0,                    # starting value of Delta when :sweep_parameter=:Delta
    :Delta_max => 0.001,                    # end value of Delta when :sweep_parameter=:Delta
    :n_Delta_slices => 21,
    :Delta_spacing_mode => :exponential,    # spacing mode for Delta slices: :linear or :exponential
    :Delta_exponential_concentration => 4.0,
    :u_range => (0.0, 1.05),             # range of u (normalsied mu/mu_c) values to consider

    :x0_alpha => 0.125 / 4,
    :x0_count => 1,
    :kappa_mode => :user_defined,   # :user_defined, :delta_over_L
    :kappas => [1e-4],
    :n_lowest_evals => 4,
    :kk_tol => 1e-8,
    :kk_maxiter => 300,
    :topological_threshold => 0.5,
    :topological_fraction => 0.5,

    # Sparse full-range discovery scan: preserve physical-mu spacing and add
    # only modest extra density toward the sweep endpoint. Narrower features may be missed.
    :mu_baseline_points => 501,
    :discovery_rho_growth_factor => 2.0,
    :discovery_rho_growth_power => 1.0,

    # Dense local sampling follows tracked boundaries. Point count increases
    # toward the sweep endpoint; preserve the old 3% minimum band and expand it from observed
    # motion/curvature, up to a configured cap.
    :mu_refine_points_per_prediction => 75,
    :mu_refine_points_per_prediction_max => 201,
    :boundary_resolution_growth_power => 1.5,
    :prediction_band_min_fraction => 0.03,
    :prediction_band_max_fraction => 0.08,
    :prediction_band_motion_multiplier => 2.5,
    :prediction_band_curvature_multiplier => 1.5,
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
    :track_match_base_tolerance_fraction => 0.01,
    :track_match_max_tolerance_fraction => 0.05,
    :track_match_step_scale => 1.5,
    :track_unmatched_penalty => 0.75,
    :track_max_missed_slices => 1,

    :rho_onset_tol => 0.01,
    :rho_onset_max_extra_slices => 40,
    :Delta_onset_tol => 0.0005,
    :Delta_onset_max_extra_slices => 200,
    :max_points_per_slice => 100_000,
    :max_total_points => 10_000_000,
    :n_eval_workers => 32,
    :progress_log_every_points => 100,
    :blas_threads => 1,

    :global_stop_file => joinpath(slice_dir, "STOP_REQUESTED"),
    :run_root_pointer_file => joinpath(slice_dir, "runs", "latest_run_root.txt"),
)