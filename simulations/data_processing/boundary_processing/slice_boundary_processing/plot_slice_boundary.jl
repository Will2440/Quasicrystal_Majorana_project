#!/usr/bin/env julia

using Plots
include(joinpath(@__DIR__, "slice_data.jl"))
using .SliceBoundaryData

function usage()
    println("Usage: julia plot_slice_boundary.jl <run_directory> [--save output.png] [--points]")
end

function main()
    isempty(ARGS) && (usage(); error("Run directory required"))
    run_dir = abspath(ARGS[1])
    save_path = nothing
    show_points = false
    i = 2
    while i <= length(ARGS)
        if ARGS[i] == "--save"
            i += 1
            i <= length(ARGS) || error("--save requires a path")
            save_path = abspath(ARGS[i])
        elseif ARGS[i] == "--points"
            show_points = true
        else
            error("Unknown option: $(ARGS[i])")
        end
        i += 1
    end

    transitions = load_transitions(run_dir)
    summary = load_slice_summary(run_dir)
    p = plot(
        xlabel="ρ = t₂/t₁",
        ylabel="μ",
        title="Slice-search phase boundaries",
        legend=:topright,
        size=(1000, 700),
        dpi=150,
    )

    if show_points
        samples = load_samples(run_dir)
        for (phase, color, label) in ((0, :blue, "Trivial samples"), (1, :red, "Topological samples"))
            subset = filter(row -> row.phase_label == phase, samples)
            scatter!(p, subset.rho, subset.mu; label=label, color=color, alpha=0.25, markersize=2)
        end
    end

    track_ids = sort(unique(Int.(transitions.track_id)))
    colors = [:black, :orange, :green, :purple, :brown, :cyan, :magenta, :gray]
    for (idx, track_id) in enumerate(track_ids)
        subset = sort(filter(row -> row.track_id == track_id, transitions), :rho)
        label = track_id == 0 ? "Onset-refinement transitions" : "Transition track $track_id"
        color = track_id == 0 ? :black : colors[mod1(idx, length(colors))]
        yerr = (subset.mu .- subset.mu_lo, subset.mu_hi .- subset.mu)
        scatter!(p, subset.rho, subset.mu; yerror=yerr, label=label, color=color, markersize=3)
        if track_id != 0 && nrow(subset) > 1
            plot!(p, subset.rho, subset.mu; color=color, alpha=0.65, label="")
        end
    end

    gap_plot = plot(
        summary.rho,
        summary.gap_count;
        xlabel="ρ = t₂/t₁",
        ylabel="Trivial intervals",
        label="Gap count",
        color=:darkgreen,
        linewidth=2,
        marker=:diamond,
        markersize=3,
        legend=:topright,
        size=(1000, 300),
        dpi=150,
    )
    p = plot(p, gap_plot; layout=(2, 1), link=:x, size=(1000, 900))

    if save_path !== nothing
        mkpath(dirname(save_path))
        savefig(p, save_path)
        println("Saved plot to $save_path")
    else
        display(p)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end