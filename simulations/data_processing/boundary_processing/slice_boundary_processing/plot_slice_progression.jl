#!/usr/bin/env julia

using Plots
include(joinpath(@__DIR__, "slice_data.jl"))
using .SliceBoundaryData

function usage()
    println("Usage: julia plot_slice_progression.jl <run_directory> [--save output_dir] [--interval N] [--max-slice N]")
end

function make_frame(samples, transitions, summary, slice_index::Int)
    current_summary = filter(row -> row.slice_index <= slice_index, summary)
    current_transitions = filter(row -> row.slice_index <= slice_index, transitions)
    p = plot(
        xlabel="ρ = t₂/t₁",
        ylabel="μ",
        title="Slice search progression: through slice $slice_index",
        legend=:topright,
        size=(950, 650),
        dpi=110,
    )

    visible_samples = filter(row -> row.slice_index <= slice_index, samples)
    for (phase, color, label) in ((0, :blue, "Trivial"), (1, :red, "Topological"))
        subset = filter(row -> row.phase_label == phase, visible_samples)
        scatter!(p, subset.rho, subset.mu; label=label, color=color, alpha=0.35, markersize=2)
    end

    track_colors = [:orange, :green, :purple, :brown, :cyan, :magenta, :gray, :black]
    for track_id in sort(unique(Int.(current_transitions.track_id)))
        subset = sort(filter(row -> row.track_id == track_id, current_transitions), :rho)
        color = track_id == 0 ? :black : track_colors[mod1(track_id, length(track_colors))]
        label = track_id == 0 ? "Onset probes" : "Boundary track $track_id"
        if nrow(subset) > 1 && track_id != 0
            plot!(p, subset.rho, subset.mu; color=color, linewidth=1.5, label=label)
        else
            scatter!(p, subset.rho, subset.mu; color=color, markersize=3, label=label)
        end
    end

    if !isempty(current_summary) && !isempty(visible_samples)
        latest = argmax(current_summary.slice_index)
        annotate!(p,
            minimum(current_summary.rho), maximum(visible_samples.mu),
            text("latest gap count: $(current_summary[latest, :gap_count])", 9, :left))
    end
    return p
end

function main()
    isempty(ARGS) && (usage(); error("Run directory required"))
    run_dir = abspath(ARGS[1])
    save_dir = nothing
    interval = 1
    max_slice = typemax(Int)
    i = 2
    while i <= length(ARGS)
        if ARGS[i] == "--save"
            i += 1
            i <= length(ARGS) || error("--save requires a directory")
            save_dir = abspath(ARGS[i])
        elseif ARGS[i] == "--interval"
            i += 1
            i <= length(ARGS) || error("--interval requires an integer")
            interval = parse(Int, ARGS[i])
            interval >= 1 || error("--interval must be positive")
        elseif ARGS[i] == "--max-slice"
            i += 1
            i <= length(ARGS) || error("--max-slice requires an integer")
            max_slice = parse(Int, ARGS[i])
        else
            error("Unknown option: $(ARGS[i])")
        end
        i += 1
    end

    slices = filter(pair -> pair[1] <= max_slice, available_slices(run_dir))
    isempty(slices) && error("No slice directories found in $run_dir")
    samples = load_samples(run_dir)
    transitions = load_transitions(run_dir)
    summary = load_slice_summary(run_dir)
    save_dir === nothing || mkpath(save_dir)

    for (position, (slice_index, _)) in enumerate(slices)
        if mod(position - 1, interval) != 0 && position != length(slices)
            continue
        end
        p = make_frame(samples, transitions, summary, slice_index)
        if save_dir === nothing
            display(p)
        else
            path = joinpath(save_dir, "slice_progress_$(lpad(slice_index, 4, '0')).png")
            savefig(p, path)
            println("Saved $path")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end