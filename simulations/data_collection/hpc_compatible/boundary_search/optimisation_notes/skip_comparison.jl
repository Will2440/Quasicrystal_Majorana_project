root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))

for skipping in (false, true)
    label = skipping ? "skip_on" : "skip_off"
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :run_name => label,
        :run_root_pointer_file => joinpath(@__DIR__, "runs", label * "_latest.txt"),
        :rho_start => 2.0,
        :rho_max => 2.04,
        :n_rho_slices => 5,
        :rho_spacing_mode => :linear,
        :rho_onset_tol => 1.0,
        :u_baseline_points => 401,
        :u_refine_points_per_prediction => 5,
        :skip_trivial_interiors => skipping,
        :confidence_min_history => 1,
        :confidence_max_residual => 0.003,
        :boundary_padding_min => 0.005,
        :boundary_padding_max => 0.01,
        :track_match_base_tolerance_fraction => 0.003,
        :track_match_max_tolerance_fraction => 0.01,
        :max_points_per_slice => 3000,
        :max_total_points => 15000,
    ))
    result = BoundarySearch.run_search(cfg; config_path=@__FILE__)
    println("RESULT ", label, " status=", result.stop_reason, " points=", result.total_points,
        " gaps=", [row.gap_count for row in result.results],
        " exclusions=", [row.n_exclusions for row in result.results],
        " fallbacks=", [row.fallback for row in result.results], " output=", result.run_root)
    flush(stdout)
    result.stop_reason == "completed" || error("Incomplete comparison: $label")
end