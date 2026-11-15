#!/usr/bin/env julia

using CSV
using DataFrames

function usage()
    println("Usage: julia merge_slice_runs.jl <first_run> <priority_run> [output_directory]")
end

function require_columns(df::AbstractDataFrame, required, path::AbstractString)
    all(name -> name in propertynames(df), required) || error("Unexpected columns in $path")
end

function slice_directory(run_dir::AbstractString, slice_index::Integer)
    return joinpath(run_dir, "slices", string("slice_", lpad(slice_index, 4, '0')))
end

function load_run(run_dir::AbstractString)
    isdir(run_dir) || error("Run directory not found: $run_dir")
    summary_path = joinpath(run_dir, "slice_summary.csv")
    transitions_path = joinpath(run_dir, "transition_boundaries.csv")
    isfile(summary_path) || error("Slice summary not found: $summary_path")
    isfile(transitions_path) || error("Transition table not found: $transitions_path")

    summary = CSV.read(summary_path, DataFrame)
    transitions = CSV.read(transitions_path, DataFrame)
    require_columns(summary, (:slice_index, :rho, :gap_count, :n_points), summary_path)
    require_columns(transitions, (:slice_index, :rho, :track_id, :mu), transitions_path)
    isempty(summary) && error("Slice summary is empty: $summary_path")
    if :sweep_parameter in propertynames(summary) && any(value -> !ismissing(value) && value != "rho", summary.sweep_parameter)
        error("This merger only supports rho sweeps, not Delta sweeps: $summary_path")
    end
    all(isfinite, Float64.(summary.rho)) || error("Non-finite rho values in $summary_path")
    length(unique(summary.rho)) == nrow(summary) || error("Duplicate rho values in $summary_path")
    length(unique(summary.slice_index)) == nrow(summary) || error("Duplicate slice indices in $summary_path")

    summary_slice_indices = Set(Int.(summary.slice_index))
    all(index -> Int(index) in summary_slice_indices, transitions.slice_index) ||
        error("Transition rows reference slices absent from $summary_path")
    for slice_index in summary.slice_index
        path = slice_directory(run_dir, Int(slice_index))
        isdir(path) || error("Slice directory not found: $path")
        isfile(joinpath(path, "point_results.csv")) || error("point_results.csv not found in $path")
    end
    return (path=String(run_dir), summary=summary, transitions=transitions)
end

function run_tag(run_dir::AbstractString)
    matched = match(r"_(\d{8}_\d{6})$", basename(run_dir))
    return matched === nothing ? basename(run_dir) : matched.captures[1]
end

function is_within(path::AbstractString, parent::AbstractString)
    relative = relpath(path, parent)
    return relative == "." || !(relative == ".." || startswith(relative, "../"))
end

function write_merge_info(path, first_run, priority_run, priority_range, first_kept, first_dropped, priority_kept)
    open(path, "w") do io
        println(io, "# Merged dynamic slice-search run")
        println(io)
        println(io, "- First source: `$(first_run.path)`")
        println(io, "- Priority source: `$(priority_run.path)`")
        println(io, "- Priority rho interval: [$(priority_range[1]), $(priority_range[2])] ")
        println(io, "- Retained from first source: $first_kept slices")
        println(io, "- Dropped from first source inside priority interval: $first_dropped slices")
        println(io, "- Retained from priority source: $priority_kept slices")
        println(io)
        println(io, "The slice summary, transition table, and `slices/slice_NNNN` folders are merged and reindexed together in increasing rho order. Nonzero transition track IDs from the priority source are offset to prevent unrelated run-local IDs from being connected by plotting tools. Source plots and restart checkpoints are not copied because they do not represent the merged data.")
    end
end

function merge_runs(first_path::AbstractString, priority_path::AbstractString, output_path::AbstractString)
    first_run = load_run(abspath(first_path))
    priority_run = load_run(abspath(priority_path))
    output = abspath(output_path)
    output == first_run.path && error("Output directory cannot be the first source run")
    output == priority_run.path && error("Output directory cannot be the priority source run")
    is_within(output, first_run.path) && error("Output directory cannot be inside the first source run")
    is_within(output, priority_run.path) && error("Output directory cannot be inside the priority source run")
    ispath(output) && error("Output path already exists: $output")

    priority_range = extrema(Float64.(priority_run.summary.rho))
    first_keep = (Float64.(first_run.summary.rho) .< priority_range[1]) .|
        (Float64.(first_run.summary.rho) .> priority_range[2])
    first_selected = first_run.summary[first_keep, :]
    first_dropped = nrow(first_run.summary) - nrow(first_selected)

    combined_summary = vcat(first_selected, priority_run.summary; cols=:union)
    source_ids = vcat(fill(1, nrow(first_selected)), fill(2, nrow(priority_run.summary)))
    sort_order = sortperm(Float64.(combined_summary.rho))
    combined_summary = combined_summary[sort_order, :]
    source_ids = source_ids[sort_order]
    original_indices = Int.(combined_summary.slice_index)

    index_maps = (Dict{Int,Int}(), Dict{Int,Int}())
    for (new_index, (source_id, original_index)) in enumerate(zip(source_ids, original_indices))
        index_maps[source_id][original_index] = new_index
    end
    combined_summary.slice_index = collect(1:nrow(combined_summary))

    first_transitions = filter(row -> haskey(index_maps[1], Int(row.slice_index)), first_run.transitions)
    priority_transitions = copy(priority_run.transitions)
    first_transitions.slice_index = [index_maps[1][Int(index)] for index in first_transitions.slice_index]
    priority_transitions.slice_index = [index_maps[2][Int(index)] for index in priority_transitions.slice_index]

    retained_track_ids = filter(!=(0), Int.(first_transitions.track_id))
    track_offset = isempty(retained_track_ids) ? 0 : maximum(retained_track_ids)
    priority_transitions.track_id = [Int(id) == 0 ? 0 : Int(id) + track_offset for id in priority_transitions.track_id]
    combined_transitions = vcat(first_transitions, priority_transitions; cols=:union)
    sort!(combined_transitions, [:rho, :mu])

    mkpath(dirname(output))
    mktempdir(dirname(output); prefix=".slice_merge_") do staging
        mkpath(joinpath(staging, "slices"))
        for (new_index, (source_id, original_index)) in enumerate(zip(source_ids, original_indices))
            source_run = source_id == 1 ? first_run : priority_run
            source_slice = slice_directory(source_run.path, original_index)
            target_slice = slice_directory(staging, new_index)
            cp(source_slice, target_slice)
        end
        CSV.write(joinpath(staging, "slice_summary.csv"), combined_summary)
        CSV.write(joinpath(staging, "transition_boundaries.csv"), combined_transitions)
        write_merge_info(
            joinpath(staging, "merge_info.md"),
            first_run,
            priority_run,
            priority_range,
            nrow(first_selected),
            first_dropped,
            nrow(priority_run.summary),
        )
        mv(staging, output)
    end

    println("Created merged run: $output")
    println("Retained $(nrow(first_selected)) first-run slices and $(nrow(priority_run.summary)) priority-run slices; dropped $first_dropped first-run slices in the priority rho interval.")
    println("Combined slices: $(nrow(combined_summary)); transitions: $(nrow(combined_transitions))")
    return output
end

function main()
    length(ARGS) in (2, 3) || (usage(); error("Two source run directories and an optional output directory are required"))
    first_path, priority_path = ARGS[1], ARGS[2]
    output_path = if length(ARGS) == 3
        ARGS[3]
    else
        parent = dirname(abspath(first_path))
        name = "dynamic_mb_slice_search_merged_$(run_tag(first_path))_$(run_tag(priority_path))"
        joinpath(parent, name)
    end
    merge_runs(first_path, priority_path, output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
