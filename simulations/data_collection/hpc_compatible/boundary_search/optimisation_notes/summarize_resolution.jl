using CSV, DataFrames
include("resolution_analysis.jl")
using .ResolutionAnalysis

isempty(ARGS) && error("Usage: summarize_resolution.jl run_directory [...]")
for directory in ARGS
    points = CSV.read(joinpath(directory, "slices", "slice_0001", "point_results.csv"), DataFrame)
    points = points[points.stage .== "planned", :]
    sort!(points, :u)
    runs = phase_runs(points.phase_label)
    rows = NamedTuple[]
    for (number, run) in enumerate(runs)
        left_edge = run.first == 1
        right_edge = run.last == nrow(points)
        push!(rows, (run=number, phase=run.phase, first_u=points.u[run.first],
            last_u=points.u[run.last], first_mu=points.mu[run.first], last_mu=points.mu[run.last],
            n_samples=run.last-run.first+1, sampled_span=points.mu[run.last]-points.mu[run.first],
            left_censored=left_edge, right_censored=right_edge,
            outer_span=left_edge || right_edge ? missing : points.mu[run.last+1]-points.mu[run.first-1]))
    end
    CSV.write(joinpath(directory, "phase_run_widths.csv"), DataFrame(rows))
    interior = filter(row -> !row.left_censored && !row.right_censored, rows)
    println(basename(directory), " rho=", points.rho[1], " Delta=", points.Delta[1],
        " fine_mu_step=", (last(points.mu)-first(points.mu))/(nrow(points)-1),
        " phase0_runs=", count(row -> row.phase == 0, rows),
        " minimum_interior_phase0_samples=", minimum([row.n_samples for row in interior if row.phase == 0]; init=typemax(Int)),
        " minimum_interior_phase1_samples=", minimum([row.n_samples for row in interior if row.phase == 1]; init=typemax(Int)),
        " interior_runs_with_at_most_two_samples=", count(row -> row.n_samples <= 2, interior))
end