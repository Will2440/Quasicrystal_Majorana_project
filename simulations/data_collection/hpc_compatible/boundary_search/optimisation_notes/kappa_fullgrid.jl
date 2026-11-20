root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))

cases = isempty(ARGS) ? (
    ("rho5_k316e5", 5.0, 0.05, 10.0^(-4.5)),
    ("rho5_k562e5", 5.0, 0.05, 10.0^(-4.25)),
    ("rho35_delta002_heuristic", 3.5, 0.02, 0.02 / 500),
    ("rho35_delta002_k1e4", 3.5, 0.02, 1e-4),
) : ARGS == ["fine"] ? (
    ("fine_rho5_k316e5", 5.0, 0.05, 10.0^(-4.5)),
    ("fine_rho35_delta002_heuristic", 3.5, 0.02, 0.02 / 500),
) : error("Expected no arguments or fine")
point_count = isempty(ARGS) ? 3201 : 6401
for (label, rho, delta, kappa) in cases
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :run_name => "kappa_fullgrid_" * label,
        :run_root_pointer_file => joinpath(@__DIR__, "runs", "kappa_fullgrid_" * label * "_latest.txt"),
        :rho_start => rho, :rho_max => rho, :Delta => delta, :kappas => [kappa],
        :u_baseline_points => point_count, :u_refine_points_per_prediction => 0,
        :transition_max_refine_iters => 0, :n_eval_workers => 4,
        :progress_log_every_points => 1000, :max_points_per_slice => point_count + 100,
        :max_total_points => point_count + 100,
    ))
    result = BoundarySearch.run_search(cfg; config_path=@__FILE__)
    result.stop_reason == "completed" || error("Incomplete scan: $label")
    println("RESULT ", label, " kappa=", kappa, " gaps=", only(result.results).gap_count,
        " points=", result.total_points, " path=", result.run_root)
    flush(stdout)
end