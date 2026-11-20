using CSV, DataFrames
include("resolution_analysis.jl")
using .ResolutionAnalysis

isempty(ARGS) && error("Usage: survey_summary.jl run_directory [...]")
for directory in ARGS
    table = CSV.read(joinpath(directory, "resolution_comparison.csv"), DataFrame)
    rows = collect(eachrow(table))
    reference = only(filter(row -> row.stride == 1 && row.offset == 0, rows))
    passing = coarsest_passing_stride(rows)
    points = CSV.read(joinpath(directory, "slices", "slice_0001", "point_results.csv"), DataFrame)
    points = points[points.stage .== "planned", :]
    sort!(points, :u)
    interior = filter(run -> run.first > 1 && run.last < nrow(points), phase_runs(points.phase_label))
    shortest(phase) = minimum([run.last - run.first + 1 for run in interior if run.phase == phase]; init=typemax(Int))
    println("rho=", reference.rho, " Delta=", reference.Delta, " kappa=", reference.kappa,
        " finest_points=", reference.n_points, " fine_mu_step=", reference.fine_mu_spacing,
        " reference_runs=", reference.phase0_runs,
        " coarsest_passing_mu_step=", isnothing(passing) ? "unresolved" : passing * reference.fine_mu_spacing,
        " coarsest_passing_points=", isnothing(passing) ? "unresolved" : (reference.n_points-1) ÷ passing + 1,
        " min_interior_phase0_samples=", shortest(0), " min_interior_phase1_samples=", shortest(1))
end