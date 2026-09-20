# Returns Dict{Symbol,Any}

project_root = @__DIR__
dyn_dir = normpath(joinpath(project_root, ".."))
hpc_dir = normpath(joinpath(dyn_dir, ".."))

Dict{Symbol,Any}(
    # -----------------------------------------------------------------
    # Run identity and output locations
    # -----------------------------------------------------------------
    :run_name => "dynamic_mb_search_singlejob",
    :output_root => joinpath(dyn_dir, "runs"),

    # -----------------------------------------------------------------
    # Sequence source and selection
    # -----------------------------------------------------------------
    :sequence_bson_path => joinpath(hpc_dir, "serial_batch", "batch_params", "sequences", "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1.bson"),
    :target_slope => 2 / (1 + sqrt(5)),
    :slope_tolerance => 0.001,
    :target_phason => 0.0,
    :phason_tolerance => 1e-4,
    :sequence_selection_mode => :closest,

    # -----------------------------------------------------------------
    # Physical parameters
    # -----------------------------------------------------------------
    :N => 500,
    :Delta => 0.05,

    :sweep_mode => :vary_t2,
    :t1_fixed => 1.0,
    :t2_fixed => 1.0,
    :rho_range => (1.0, 3.0),
    :u_range => (0.0, 1.05),

    # -----------------------------------------------------------------
    # Spectral localiser sampler settings
    # -----------------------------------------------------------------
    :x0_alpha => 0.125 / 2,
    :x0_count => 5,
    :kappas => [1e-4],
    :n_lowest_evals => 4,
    :kk_tol => 1e-8,
    :kk_maxiter => 300,

    :topological_threshold => 0.5,
    :topological_fraction => 0.5,

    # -----------------------------------------------------------------
    # Adaptive search controls
    # -----------------------------------------------------------------
    # rho range is wider than legacy [1,2], so increase bins to keep similar rho resolution.
    :init_rho_bins => 32,
    :init_u_bins => 24,
    :max_depth => 7,
    :rho_cell_tol => 5e-4,
    :u_cell_tol => 5e-4,
    :gap_refine_tol => 5e-4,
    :point_round_digits => 12,

    :max_iters => 12,
    :max_points => 250_000,
    :gap_count_stable_iters => 3,
    :gap_count_tol => 0,
    :max_expected_gaps => 0,
    :max_expected_gaps_tol => 0,

    :rho_max_probe_points => 2000,
    :min_gap_run_points => 3,

    # -----------------------------------------------------------------
    # Single-job parallel execution controls
    # -----------------------------------------------------------------
    :n_eval_workers => 56,
    :progress_log_every_points => 25,
    :require_all_worker_results => true,

    # Conservative initial estimate from prior N=500 throughput discussions.
    # The run records actual achieved points/s each iteration.
    :estimated_points_per_second => 44.0,

    # Avoid BLAS oversubscription when many Julia workers are active.
    :blas_threads => 1,

    # Graceful stop/cancel controls
    :global_stop_file => joinpath(dyn_dir, "STOP_REQUESTED"),
    :run_root_pointer_file => joinpath(dyn_dir, "runs", "latest_run_root.txt"),
)
