#!/usr/bin/env julia

include(joinpath(@__DIR__, "slice_data.jl"))
using .SliceBoundaryData
using CSV
using DataFrames

function main()
    length(ARGS) >= 1 || error("Usage: julia export_slice_boundaries.jl <run_directory> [output.csv]")
    run_dir = abspath(ARGS[1])
    transitions = load_transitions(run_dir)
    summary = load_slice_summary(run_dir)

    if !isempty(transitions) && !isempty(summary)
        select!(summary, [:slice_index, :rho, :gap_count, :n_points])
        transitions = leftjoin(transitions, summary; on=[:slice_index, :rho])
    end

    out_path = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(run_dir, "boundary_data.csv")
    mkpath(dirname(out_path))
    CSV.write(out_path, transitions)
    println("Exported $(nrow(transitions)) transition boundaries to $out_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end