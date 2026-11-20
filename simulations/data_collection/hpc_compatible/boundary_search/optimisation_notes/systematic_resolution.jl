using CSV, DataFrames
root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))
include(joinpath(@__DIR__, "resolution_analysis.jl"))
using .ResolutionAnalysis

cases = isempty(ARGS) ? (("rho2_delta005", 2.0, 0.05),
    ("rho35_delta005", 3.5, 0.05), ("rho5_delta005", 5.0, 0.05),
    ("rho35_delta002", 3.5, 0.02)) : ARGS == ["fine"] ?
    (("fine_rho5_delta005", 5.0, 0.05), ("fine_rho35_delta002", 3.5, 0.02)) :
    ARGS == ["rho-range"] ? (("rho25_delta005", 2.5, 0.05), ("rho3_delta005", 3.0, 0.05),
        ("rho4_delta005", 4.0, 0.05), ("rho45_delta005", 4.5, 0.05)) : error("Expected no arguments, fine or rho-range")
fine_count = ARGS == ["fine"] ? 12801 : 6401
for (label, rho, delta) in cases
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :run_name => "systematic_" * label,
        :run_root_pointer_file => joinpath(@__DIR__, "runs", "systematic_" * label * "_latest.txt"),
        :rho_start => rho, :rho_max => rho, :Delta => delta,
        :u_baseline_points => fine_count, :u_refine_points_per_prediction => 0,
        :transition_max_refine_iters => 0, :n_eval_workers => 4,
        :progress_log_every_points => 1000, :max_points_per_slice => fine_count + 100,
        :max_total_points => fine_count + 100,
    ))
    result = BoundarySearch.run_search(cfg; config_path=@__FILE__)
    result.stop_reason == "completed" || error("Incomplete scan: $label")
    points = CSV.read(joinpath(result.run_root, "slices", "slice_0001", "point_results.csv"), DataFrame)
    sort!(points, :u)
    nrow(points) == fine_count || error("Unexpected fine grid size")
    labels = Int.(points.phase_label)
    step = (last(points.mu) - first(points.mu)) / (fine_count - 1)
    rows = NamedTuple[]
    for stride in (16, 8, 4, 2, 1)
        for offset in (stride == 1 ? (0,) : (0, stride ÷ 2))
            metrics = grid_metrics(labels, stride, offset, step)
            push!(rows, merge((rho=rho, Delta=delta, kappa=1e-4, stride=stride,
                offset=offset, fine_mu_spacing=step), metrics))
        end
    end
    CSV.write(joinpath(result.run_root, "resolution_comparison.csv"), DataFrame(rows))
    println("RESULT ", label, " path=", result.run_root, " metrics=", rows)
    flush(stdout)
end