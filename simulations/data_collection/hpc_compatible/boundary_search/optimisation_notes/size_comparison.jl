root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))

for size in (300, 400)
    label = "rho5_N$(size)_grid1601_k1e4"
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :run_name => label,
        :run_root_pointer_file => joinpath(@__DIR__, "runs", label * "_latest.txt"),
        :N => size,
        :rho_start => 5.0,
        :rho_max => 5.0,
        :u_baseline_points => 1601,
        :kappas => [1e-4],
        :max_points_per_slice => 5000,
        :max_total_points => 5000,
    ))
    result = BoundarySearch.run_search(cfg; config_path=@__FILE__)
    println("RESULT ", label, " status=", result.stop_reason, " points=", result.total_points,
        " gaps=", isempty(result.results) ? "missing" : only(result.results).gap_count,
        " output=", result.run_root)
    flush(stdout)
    result.stop_reason == "completed" || error("Incomplete comparison: $label")
end