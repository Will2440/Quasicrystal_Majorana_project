#!/usr/bin/env julia

using CSV
using DataFrames
using Plots
include(joinpath(@__DIR__, "slice_data.jl"))
using .SliceBoundaryData

function usage()
    println("Usage: julia plot_gap_plateaus.jl <run_directory> [max_gap_index]")
end

function trivial_intervals(transitions::AbstractDataFrame, slice_index::Integer)
    rows = sort(filter(row -> row.slice_index == slice_index, transitions), :mu)
    intervals = NamedTuple[]
    for index in 1:(nrow(rows) - 1)
        if rows.phase_hi[index] == 0 && rows.phase_lo[index + 1] == 0
            push!(intervals, (
                mu_lo=Float64(rows.mu[index]),
                mu_hi=Float64(rows.mu[index + 1]),
                mu_center=(Float64(rows.mu[index]) + Float64(rows.mu[index + 1])) / 2,
            ))
        end
    end
    return intervals
end

function unmatched_interval_indices(previous, current)
    old_count, current_count = length(previous), length(current)
    current_count >= old_count || error("Trivial interval count decreased unexpectedly during an increasing gap-count step")
    old_count == 0 && return collect(1:current_count)

    costs = fill(Inf, old_count + 1, current_count + 1)
    matched = falses(old_count + 1, current_count + 1)
    costs[1, :] .= 0.0
    for old_index in 1:old_count, current_index in 1:current_count
        skip_cost = costs[old_index + 1, current_index]
        match_cost = costs[old_index, current_index] +
            abs(previous[old_index].mu_center - current[current_index].mu_center)
        if match_cost <= skip_cost
            costs[old_index + 1, current_index + 1] = match_cost
            matched[old_index + 1, current_index + 1] = true
        else
            costs[old_index + 1, current_index + 1] = skip_cost
        end
    end

    matched_current = Int[]
    old_index, current_index = old_count, current_count
    while old_index > 0
        if matched[old_index + 1, current_index + 1]
            push!(matched_current, current_index)
            old_index -= 1
            current_index -= 1
        else
            current_index -= 1
        end
        current_index >= 0 || error("Could not match trivial intervals across rho slices")
    end
    return setdiff(1:current_count, matched_current)
end

function mirrored_transitions(transitions::AbstractDataFrame)
    original = copy(transitions)
    original[!, :is_mirrored] = falses(nrow(original))
    mirrored = copy(transitions)
    mirrored[!, :mu] = -Float64.(transitions.mu)
    mirrored[!, :mu_lo] = -Float64.(transitions.mu_hi)
    mirrored[!, :mu_hi] = -Float64.(transitions.mu_lo)
    mirrored[!, :phase_lo] = Int.(transitions.phase_hi)
    mirrored[!, :phase_hi] = Int.(transitions.phase_lo)
    mirrored[!, :is_mirrored] = trues(nrow(mirrored))
    result = vcat(original, mirrored; cols=:setequal)
    sort!(result, [:slice_index, :mu])
    return result
end

function extract_gap_plateaus(summary::AbstractDataFrame, transitions::AbstractDataFrame, max_gap_index::Int)
    max_gap_index >= 1 || error("Maximum gap index must be positive")
    ordered = sort(summary, :rho)
    counts = trivial_gap_counts(ordered)
    stop_index = findfirst(index -> index > 1 && counts[index] < counts[index - 1], eachindex(counts))
    last_index = stop_index === nothing ? length(counts) : stop_index - 1
    last_index > 0 || error("No increasing gap-count data found")
    counts = counts[1:last_index]
    ordered = ordered[1:last_index, :]
    maximum(counts) >= max_gap_index || error("Adjusted trivial gap count reaches only $(maximum(counts)); cannot extract $max_gap_index gaps")

    intervals_by_slice = [trivial_intervals(transitions, index) for index in ordered.slice_index]
    for index in eachindex(counts)
        length(intervals_by_slice[index]) == counts[index] || error(
            "Slice $(ordered.slice_index[index]) at rho $(ordered.rho[index]) has $(length(intervals_by_slice[index])) bounded trivial intervals but adjusted count $(counts[index])")
    end

    births = NamedTuple[]
    previous_count = 0
    for slice_position in eachindex(counts)
        count = counts[slice_position]
        count <= previous_count && continue
        new_indices = unmatched_interval_indices(
            intervals_by_slice[slice_position - 1 + (slice_position == 1 ? 1 : 0)],
            intervals_by_slice[slice_position])
        if slice_position == 1
            new_indices = collect(1:length(intervals_by_slice[slice_position]))
        end
        count_delta = count - previous_count
        length(new_indices) == count_delta || error(
            "Gap count rose by $count_delta at rho $(ordered.rho[slice_position]), but interval matching found $(length(new_indices)) new gaps")

        plateau_end_position = slice_position
        while plateau_end_position < last_index && counts[plateau_end_position + 1] == count
            plateau_end_position += 1
        end
        rho_start = Float64(ordered.rho[slice_position])
        rho_end = Float64(ordered.rho[plateau_end_position])
        for (offset, interval_index) in enumerate(new_indices)
            push!(births, (
                gap_index=previous_count + offset,
                mu_center=intervals_by_slice[slice_position][interval_index].mu_center,
                rho_onset=rho_start,
                rho_plateau_end=rho_end,
                rho_plateau_width=rho_end - rho_start,
                count_plateau=count,
            ))
        end
        previous_count = count
    end
    length(births) >= max_gap_index || error("Only $(length(births)) gap births could be identified")
    filter!(birth -> birth.gap_index <= max_gap_index, births)
    return births
end

function make_plateau_table(births)
    gap_count = length(births)
    rows = NamedTuple[]
    for birth in births
        for (sign, mu_center) in (("-", -birth.mu_center), ("+", birth.mu_center))
            push!(rows, (
                gap_index=birth.gap_index,
                symmetry=sign,
                mu_center=mu_center,
                rho_onset=birth.rho_onset,
                rho_plateau_end=birth.rho_plateau_end,
                rho_plateau_width=birth.rho_plateau_width,
                count_plateau=birth.count_plateau,
            ))
        end
    end
    sort!(rows; by=row -> row.mu_center)
    table = DataFrame(rows)
    table.normalised_gap_rank = gap_count == 1 ? fill(0.5, nrow(table)) :
        collect(range(0.0, 1.0; length=nrow(table)))
    return table
end

function normalise_boundary_mu(transitions::AbstractDataFrame)
    scales = Dict(
        Float64(group.rho[1]) => maximum(abs.(Float64.(group.mu)))
        for group in groupby(transitions, :rho)
    )
    all(>(0), values(scales)) || error("Cannot normalize boundary mu values: an outermost boundary is zero")
    normalized = copy(transitions)
    row_scales = [scales[Float64(rho)] for rho in normalized.rho]
    normalized[!, :mu_normalisation_scale] = row_scales
    normalized[!, :mu_normalised] = Float64.(transitions.mu) ./ row_scales
    normalized[!, :mu_lo_normalised] = Float64.(transitions.mu_lo) ./ row_scales
    normalized[!, :mu_hi_normalised] = Float64.(transitions.mu_hi) ./ row_scales
    return normalized, scales
end

function save_symmetric_boundary_plot(
    transitions::AbstractDataFrame,
    output_path::AbstractString;
    inverse_rho::Bool=false,
    log_rho::Bool=false,
    normalised_mu::Bool=false,
)
    inverse_rho && log_rho && error("Choose inverse-rho or log-rho scaling, not both")
    x_values = inverse_rho ? 1.0 ./ Float64.(transitions.rho) : Float64.(transitions.rho)
    x_label = inverse_rho ? "1/ρ = t₁/t₂" : "ρ = t₂/t₁"
    title = normalised_mu ? "Symmetric phase boundaries (normalized mu)" :
        inverse_rho ? "Symmetric phase boundaries versus inverse rho" :
        log_rho ? "Symmetric phase boundaries versus rho (log scale)" : "Symmetric phase boundaries"
    mu_column = normalised_mu ? :mu_normalised : :mu
    p = plot(
        xlabel=x_label,
        ylabel=normalised_mu ? "Normalised μ (per-rho outer-boundary scale)" : "μ",
        title=title,
        xscale=log_rho ? :log10 : :identity,
        legend=:topright,
        size=(1100, 700),
        dpi=150,
    )
    for (is_mirrored, color, label) in ((false, :navy, "+μ boundaries"), (true, :darkorange, "Mirrored -μ boundaries"))
        subset = filter(row -> row.is_mirrored == is_mirrored, transitions)
        x_subset = inverse_rho ? 1.0 ./ Float64.(subset.rho) : Float64.(subset.rho)
        scatter!(p, x_subset, subset[!, mu_column]; label=label, color=color, alpha=0.65, markersize=2)
        for track_id in filter(!=(0), sort(unique(Int.(subset.track_id))))
            track = sort(filter(row -> row.track_id == track_id, subset), :rho)
            track_x = inverse_rho ? 1.0 ./ Float64.(track.rho) : Float64.(track.rho)
            if inverse_rho
                order = sortperm(track_x)
                track_x = track_x[order]
                track_mu = Float64.(track[!, mu_column])[order]
            else
                track_mu = Float64.(track[!, mu_column])
            end
            nrow(track) > 1 && plot!(p, track_x, track_mu; color=color, alpha=0.4, linewidth=1, label="")
        end
    end
    savefig(p, output_path)
end

function save_gap_plateau_plot(
    table::AbstractDataFrame,
    max_gap_index::Int,
    output_path::AbstractString;
    normalised_mu::Bool=false,
)
    mu_column = normalised_mu ? :mu_center_normalised : :mu_center
    rank_column = normalised_mu ? :normalised_mu_gap_rank : :normalised_gap_rank
    x_label = normalised_mu ? "Normalised μ centre (per-rho boundary scale)" : "μ centre of trivial gap"
    plot_title = normalised_mu ? "First $max_gap_index gaps on normalized μ scale" :
        "First $max_gap_index trivial gaps and their rho plateau widths"
    p = plot(
        xlabel=x_label,
        ylabel=normalised_mu ? "Normalised gap rank (ordered by normalized μ)" : "Normalised gap rank (ordered by μ)",
        title=plot_title,
        legend=false,
        size=(1100, 750),
        dpi=150,
        ylims=(-0.04, 1.04),
        yticks=0.0:0.1:1.0,
    )
    colors = palette(:default, max_gap_index)
    for row in eachrow(table)
        color = colors[row.gap_index]
        half_width = row.rho_plateau_width / 2
        center = row[mu_column]
        plot!(p, [center - half_width, center + half_width],
            [row[rank_column], row[rank_column]]; color=color, linewidth=3)
        scatter!(p, [center], [row[rank_column]]; color=color, markersize=3, label="")
        annotate!(p, center, row[rank_column] + 0.018,
            text("$(row.gap_index)$(row.symmetry)", 7, color, :center))
    end
    savefig(p, output_path)
end

function main()
    length(ARGS) in (1, 2) || (usage(); error("Run directory and optional maximum gap index required"))
    run_dir = abspath(ARGS[1])
    max_gap_index = length(ARGS) == 2 ? something(tryparse(Int, ARGS[2]), 0) : 20
    max_gap_index > 0 || error("Maximum gap index must be a positive integer")
    summary = load_slice_summary(run_dir)
    transitions = load_transitions(run_dir)
    births = extract_gap_plateaus(summary, transitions, max_gap_index)
    plateaus = make_plateau_table(births)

    symmetric = mirrored_transitions(transitions)
    normalized_symmetric, boundary_mu_scales = normalise_boundary_mu(symmetric)
    plateau_scales = [boundary_mu_scales[Float64(rho)] for rho in plateaus.rho_onset]
    plateaus[!, :mu_center_normalised] = plateaus.mu_center ./ plateau_scales
    plateaus[!, :mu_normalisation_scale] = plateau_scales
    normalized_order = sortperm(plateaus.mu_center_normalised)
    plateaus.normalised_mu_gap_rank = similar(plateaus.normalised_gap_rank)
    plateaus.normalised_mu_gap_rank[normalized_order] = collect(range(0.0, 1.0; length=nrow(plateaus)))
    symmetric_path = joinpath(run_dir, "symmetric_transition_boundaries.csv")
    normalized_symmetric_path = joinpath(run_dir, "symmetric_transition_boundaries_normalised.csv")
    plateau_path = joinpath(run_dir, "gap_plateau_data.csv")
    boundary_plot_path = joinpath(run_dir, "symmetric_boundary_overview.png")
    normalized_boundary_plot_path = joinpath(run_dir, "symmetric_boundary_normalised_overview.png")
        inverse_boundary_plot_path = joinpath(run_dir, "symmetric_boundary_inverse_rho.png")
        log_boundary_plot_path = joinpath(run_dir, "symmetric_boundary_log_rho.png")
    plateau_plot_path = joinpath(run_dir, "gap_plateau_overview.png")
    normalized_plateau_plot_path = joinpath(run_dir, "gap_plateau_normalised_mu_overview.png")
    CSV.write(symmetric_path, symmetric)
    CSV.write(normalized_symmetric_path, normalized_symmetric)
    CSV.write(plateau_path, plateaus)
    save_symmetric_boundary_plot(symmetric, boundary_plot_path)
    save_symmetric_boundary_plot(symmetric, inverse_boundary_plot_path; inverse_rho=true)
    save_symmetric_boundary_plot(symmetric, log_boundary_plot_path; log_rho=true)
    save_symmetric_boundary_plot(normalized_symmetric, normalized_boundary_plot_path; normalised_mu=true)
    save_gap_plateau_plot(plateaus, max_gap_index, plateau_plot_path)
    save_gap_plateau_plot(plateaus, max_gap_index, normalized_plateau_plot_path; normalised_mu=true)

    println("Saved mirrored transitions: $symmetric_path ($(nrow(symmetric)) rows)")
    println("Saved per-rho normalized mirrored transitions: $normalized_symmetric_path ($(length(boundary_mu_scales)) rho scales)")
    println("Saved gap plateau data: $plateau_path ($(nrow(plateaus)) rows)")
    println("Saved mirrored boundary plot: $boundary_plot_path")
    println("Saved inverse-rho mirrored boundary plot: $inverse_boundary_plot_path")
    println("Saved log-rho mirrored boundary plot: $log_boundary_plot_path")
    println("Saved normalized mirrored boundary plot: $normalized_boundary_plot_path")
    println("Saved gap plateau plot: $plateau_plot_path")
    println("Saved normalized-mu gap plateau plot: $normalized_plateau_plot_path")
    println("Extracted $(length(births)) gap births through rho $(maximum(plateaus.rho_onset)); analysis stopped before the first gap-count decrease")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
