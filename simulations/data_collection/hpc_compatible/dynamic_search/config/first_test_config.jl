# Returns Dict{Symbol,Any}

project_root = @__DIR__
dyn_dir = normpath(joinpath(project_root, ".."))
hpc_dir = normpath(joinpath(dyn_dir, ".."))

Dict{Symbol,Any}(
    # -----------------------------------------------------------------
    # Run identity and output locations
    # -----------------------------------------------------------------
    :run_name => "dynamic_mb_search",
    :output_root => joinpath(dyn_dir, "runs"),

    # -----------------------------------------------------------------
    # Sequence source and selection
    # -----------------------------------------------------------------
    :sequence_bson_path => joinpath(hpc_dir, "serial_batch", "batch_params", "sequences", "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1.bson"),
    :target_slope => 2 / (1 + sqrt(5)),
    :slope_tolerance => 0.001,
    :target_phason => 0.0,
    :phason_tolerance => 1e-4,
    :sequence_selection_mode => :closest,   # :closest or :first

    # -----------------------------------------------------------------
    # Physical parameters
    # -----------------------------------------------------------------
    :N => 500,
    :Delta => 0.05,

    # rho = t2/t1 exploration mode
    :sweep_mode => :vary_t2,                # :vary_t2 or :vary_t1
    :t1_fixed => 1.0,
    :t2_fixed => 1.0,
    :rho_range => (1.0, 2.0),

    # Normalized chemical potential coordinate u = mu/mu_c_est
    :u_range => (0.0, 1.05),

    # -----------------------------------------------------------------
    # Spectral localiser sampler settings
    # -----------------------------------------------------------------
    :x0_alpha => 0.125 / 2,
    :x0_count => 5,
    :kappas => [1e-4], #collect(exp.(range(log(1e-6), log(1e-1), length=10))),
    :n_lowest_evals => 4,
    :kk_tol => 1e-8,
    :kk_maxiter => 300,

    # Classification robustness
    :topological_threshold => 0.5,
    :topological_fraction => 0.5,

    # -----------------------------------------------------------------
    # Adaptive search controls
    # -----------------------------------------------------------------
    :init_rho_bins => 20,
    :init_u_bins => 24,
    :max_depth => 7,
    :rho_cell_tol => 5e-4,
    :u_cell_tol => 5e-4,
    :gap_refine_tol => 5e-4,
    :point_round_digits => 12,

    # Stop conditions
    :max_iters => 12,
    :max_points => 250_000,
    :gap_count_stable_iters => 3,
    :gap_count_tol => 0,
    :max_expected_gaps => 0,              # 0 disables this criterion
    :max_expected_gaps_tol => 0,

    # rho_max cut diagnostics
    :rho_max_probe_points => 2000,
    :min_gap_run_points => 3,

    # -----------------------------------------------------------------
    # Batch execution controls
    # -----------------------------------------------------------------
    :execution_mode => :submit_and_wait,      # :prepare_only, :submit_and_wait
    :points_per_job => 2000,
    :sbatch_script => joinpath(dyn_dir, "runscripts", "worker_array.sbatch"),
    :sbatch_extra_args => "",
    :poll_seconds => 10,

    # Slurm routing for worker jobs submitted by master.jl
    :slurm_partition => "compute",               # e.g. "compute"
    :slurm_account => "phys030424",                 # e.g. "my_account"
    :slurm_qos => "",
    :slurm_constraint => "",
    :slurm_reservation => "",

    # Startup auto-calibration of worker resources on a tiny chunk
    :auto_calibrate_worker_resources => true,
    :calibration_require_success => false,
    :calibration_points => 3,
    :calibration_poll_seconds => 15,
    :calibration_partition => "compute",     # set "" to use :slurm_partition
    :calibration_account => "phys030424",           # set if your test partition requires account
    :calibration_time_seed => "00:30:00",
    :calibration_mem_seed_mb => 4096,
    :calibration_time_safety_factor => 1.8,
    :calibration_time_buffer_seconds => 120,
    :calibration_mem_safety_factor => 1.6,
    :calibration_mem_buffer_mb => 512,
    :min_worker_time_seconds => 300,
    :min_worker_mem_mb => 2048,
    :fallback_worker_time => "00:45:00",
    :fallback_worker_mem_mb => 4096,

    # Absolute and user-tunable worker concurrency controls.
    # Effective Slurm array spec is: --array=1-N%min(hard_max_workers, max_workers_per_array, N)
    :hard_max_workers => 199,
    :max_workers_per_array => 100,

    # If true, require a full set of worker results for each iteration unless a stop is requested.
    :require_all_worker_results => true,

    # Runtime estimate seed (seconds/chunk). Master updates this online from sacct p90 values.
    :initial_chunk_eta_seconds => 0.0,

    # Graceful stop/cancel controls
    :global_stop_file => joinpath(dyn_dir, "STOP_REQUESTED"),
    :run_root_pointer_file => joinpath(dyn_dir, "runs", "latest_run_root.txt"),
    :active_jobs_file => joinpath(dyn_dir, "runs", "active_worker_jobs.txt"),

    # Worker logging verbosity
    :worker_log_every => 10,
)
