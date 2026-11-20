root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))

for (label, rho, delta, count, kappa) in (
    ("rho5_grid1601_k3e5", 5.0, 0.05, 1601, 3e-5),
    ("rho5_grid1601_k1e5", 5.0, 0.05, 1601, 1e-5),
    ("rho5_grid3201_k1e4", 5.0, 0.05, 3201, 1e-4),
    ("rho35_delta002", 3.5, 0.02, 1601, 1e-4),
    ("rho35_delta010", 3.5, 0.10, 1601, 1e-4),
)
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :run_name => label,
        :run_root_pointer_file => joinpath(@__DIR__, "runs", label * "_latest.txt"),
        :rho_start => rho,
        :rho_max => rho,
        :Delta => delta,
        :u_baseline_points => count,
        :kappas => [kappa],
        :max_points_per_slice => 8000,
        :max_total_points => 8000,
    ))
    result = BoundarySearch.run_search(cfg; config_path=@__FILE__)
    println("RESULT ", label, " status=", result.stop_reason, " points=", result.total_points,
        " gaps=", isempty(result.results) ? "missing" : only(result.results).gap_count,
        " output=", result.run_root)
    flush(stdout)
    result.stop_reason == "completed" || error("Incomplete comparison: $label")
end