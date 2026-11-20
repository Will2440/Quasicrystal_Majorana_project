using CSV, DataFrames, JLD2
include("kappa_analysis.jl")
include("resolution_analysis.jl")
using .KappaAnalysis, .ResolutionAnalysis

length(ARGS) == 2 || error("Usage: compare_aligned.jl reference_run comparison_run")
function planned_points(directory)
    data = CSV.read(joinpath(directory, "slices", "slice_0001", "point_results.csv"), DataFrame)
    sort!(data[data.stage .== "planned", :], :u)
end
reference, comparison = planned_points.(ARGS)
validate_comparable_configs(JLD2.load(joinpath(ARGS[1], "run_summary.jld2"), "config"),
    JLD2.load(joinpath(ARGS[2], "run_summary.jld2"), "config"))
for column in (:rho, :Delta, :t1, :t2, :mu)
    all(isapprox.(reference[!, column], comparison[!, column]; rtol=1e-12, atol=1e-14)) ||
        error("Physical coordinates differ in $column")
end
changes = aligned_changes(reference.u, reference.phase_label, comparison.u, comparison.phase_label)
metrics = agreement_counts(reference.phase_label, comparison.phase_label)
ref_runs = phase_runs(reference.phase_label)
cmp_runs = phase_runs(comparison.phase_label)
println("aligned_points=", nrow(reference), " changed=", length(changes),
    " phase0_to_phase1=", metrics.trivial_changed, " phase1_to_phase0=", metrics.topological_lost,
    " phase0_runs=", count(run -> run.phase == 0, ref_runs), " -> ", count(run -> run.phase == 0, cmp_runs))
for index in changes
    println("u=", reference.u[index], " mu=", reference.mu[index], " label=", reference.phase_label[index],
        " -> ", comparison.phase_label[index], " localizer_gaps=", reference.gap_min[index], ",", comparison.gap_min[index])
end