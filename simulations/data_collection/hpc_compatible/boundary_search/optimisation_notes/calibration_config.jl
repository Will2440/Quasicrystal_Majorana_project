root = normpath(joinpath(@__DIR__, ".."))
cfg = include(joinpath(root, "config", "default_config.jl"))
merge!(cfg, Dict{Symbol,Any}(
    :run_name => "calibration_rho3",
    :output_root => joinpath(@__DIR__, "runs"),
    :run_root_pointer_file => joinpath(@__DIR__, "runs", "calibration_latest.txt"),
    :global_stop_file => joinpath(@__DIR__, "STOP_REQUESTED"),
    :N => 500,
    :rho_start => 3.0,
    :rho_max => 3.0,
    :n_rho_slices => 1,
    :Delta => 0.05,
    :u_baseline_points => 201,
    :skip_trivial_interiors => false,
    :n_eval_workers => 2,
    :blas_threads => 1,
    :max_points_per_slice => 1000,
    :max_total_points => 1000,
))
cfg