#!/usr/bin/env julia

using DataFrames
using LinearAlgebra
using Printf
using Plots
include(joinpath(@__DIR__, "slice_data.jl"))
using .SliceBoundaryData

function usage()
    println("Usage: julia fit_slice_gap_count.jl <run_directory> <max_gap_count> <lin|exp|poly> [poly_order] [--plot-excluded true|false] [normalise=true] [half_norm=true]")
end

function parse_bool_option(value::AbstractString, option::AbstractString)
    parsed = lowercase(value)
    parsed in ("true", "false") || error("$option requires true or false")
    return parsed == "true"
end

format_number(value::Real) = @sprintf("%.4g", value)

function format_polynomial_equation(coefficients::Vector{Float64}, variable::String)
    terms = String[]
    for (index, coefficient) in enumerate(coefficients)
        coefficient == 0 && continue
        degree = index - 1
        monomial = degree == 0 ? "" : degree == 1 ? variable : "$variable^$degree"
        term = format_number(abs(coefficient)) * monomial
        prefix = isempty(terms) ? (coefficient < 0 ? "-" : "") : coefficient < 0 ? " - " : " + "
        push!(terms, prefix * term)
    end
    return "y = " * (isempty(terms) ? "0" : join(terms))
end

function fit_polynomial(x::Vector{Float64}, y::Vector{Float64}, order::Int)
    order >= 1 || error("Polynomial order must be at least 1")
    length(unique(x)) >= order + 1 || error("Polynomial order $order requires at least $(order + 1) distinct rho values")
    center = (minimum(x) + maximum(x)) / 2
    scale = (maximum(x) - minimum(x)) / 2
    scale > 0 || error("Fit data must contain at least two distinct rho values")
    scaled_x = (x .- center) ./ scale
    design = hcat([scaled_x .^ degree for degree in 0:order]...)
    coefficients = design \ y
    model = x_values -> hcat([((x_values .- center) ./ scale) .^ degree for degree in 0:order]...) * coefficients
    return model, coefficients, center, scale
end

function fit_exponential(x::Vector{Float64}, y::Vector{Float64})
    length(y) >= 3 || error("Exponential fit requires at least 3 included points")
    minimum(x) < maximum(x) || error("Exponential fit requires at least two distinct rho values")
    x_min = minimum(x)
    x_range = maximum(x) - x_min
    scaled_x = (x .- x_min) ./ x_range

    function fit_at_exponent(exponent)
        exponential_values = exp.(exponent .* scaled_x)
        design = hcat(exponential_values, ones(length(y)))
        coefficients = design \ y
        residuals = design * coefficients - y
        return sum(abs2, residuals), coefficients
    end

    grid = collect(range(-30.0, 30.0; length=601))
    grid_errors = first.(fit_at_exponent.(grid))
    best_index = argmin(grid_errors)
    left = grid[max(1, best_index - 1)]
    right = grid[min(length(grid), best_index + 1)]
    ratio = (sqrt(5) - 1) / 2
    first_point = right - ratio * (right - left)
    second_point = left + ratio * (right - left)
    first_error = first(fit_at_exponent(first_point))
    second_error = first(fit_at_exponent(second_point))
    for _ in 1:100
        if first_error < second_error
            right = second_point
            second_point = first_point
            second_error = first_error
            first_point = right - ratio * (right - left)
            first_error = first(fit_at_exponent(first_point))
        else
            left = first_point
            first_point = second_point
            first_error = second_error
            second_point = left + ratio * (right - left)
            second_error = first(fit_at_exponent(second_point))
        end
    end
    exponent = (left + right) / 2
    _, coefficients = fit_at_exponent(exponent)
    model = x_values -> coefficients[1] .* exp.(exponent .* ((x_values .- x_min) ./ x_range)) .+ coefficients[2]
    return model, coefficients, exponent, x_min, x_range
end

function fit_curve(x::Vector{Float64}, y::Vector{Float64}, fit_type::String, poly_order::Union{Nothing,Int})
    if fit_type == "lin"
        length(y) >= 2 || error("Linear fit requires at least 2 included points")
        model, coefficients, center, scale = fit_polynomial(x, y, 1)
        intercept = coefficients[1] - coefficients[2] * center / scale
        slope = coefficients[2] / scale
        equation = format_polynomial_equation([intercept, slope], "ρ")
        return model, "Linear fit", equation
    elseif fit_type == "exp"
        model, coefficients, exponent, x_min, x_range = fit_exponential(x, y)
        offset = coefficients[2]
        offset_text = offset < 0 ? " - $(format_number(abs(offset)))" : " + $(format_number(offset))"
        equation = "y = $(format_number(coefficients[1])) exp($(format_number(exponent))z)$offset_text\nz = (ρ - $(format_number(x_min))) / $(format_number(x_range))"
        return model, "Exponential fit", equation
    elseif fit_type == "poly"
        poly_order === nothing && error("Polynomial fit requires an order, for example: poly 2")
        length(y) >= poly_order + 1 || error("Polynomial order $poly_order requires at least $(poly_order + 1) included points")
        model, coefficients, center, scale = fit_polynomial(x, y, poly_order)
        equation = format_polynomial_equation(coefficients, "z") *
            "\nz = (ρ - $(format_number(center))) / $(format_number(scale))"
        return model, "Polynomial fit (order $poly_order)", equation
    end
    error("Unknown fit type '$fit_type'; choose lin, exp, or poly")
end

function main()
    length(ARGS) >= 3 || (usage(); error("Run directory, maximum gap count, and fit type are required"))
    run_dir = abspath(ARGS[1])
    max_gap_count = tryparse(Int, ARGS[2])
    max_gap_count === nothing && error("Maximum gap count must be an integer: $(ARGS[2])")
    max_gap_count >= 0 || error("Maximum gap count must be non-negative")
    fit_type = ARGS[3]
    poly_order = nothing
    plot_excluded = true
    normalise = false
    half_norm = false
    i = 4
    if fit_type == "poly"
        i <= length(ARGS) || error("Polynomial fit requires an order, for example: poly 2")
        poly_order = tryparse(Int, ARGS[i])
        poly_order === nothing && error("Polynomial order must be an integer: $(ARGS[i])")
        i += 1
    end
    while i <= length(ARGS)
        argument = ARGS[i]
        if argument == "--plot-excluded"
            i += 1
            i <= length(ARGS) || error("--plot-excluded requires true or false")
            plot_excluded = parse_bool_option(ARGS[i], "--plot-excluded")
        elseif argument in ("--normalise", "--normalize")
            i += 1
            i <= length(ARGS) || error("$argument requires true or false")
            normalise = parse_bool_option(ARGS[i], argument)
        elseif startswith(argument, "--normalise=") || startswith(argument, "--normalize=")
            option, value = split(argument, "="; limit=2)
            normalise = parse_bool_option(value, option)
        elseif startswith(argument, "normalise=")
            _, value = split(argument, "="; limit=2)
            normalise = parse_bool_option(value, "normalise")
        elseif argument in ("--half-norm", "--half_norm")
            i += 1
            i <= length(ARGS) || error("$argument requires true or false")
            half_norm = parse_bool_option(ARGS[i], argument)
        elseif startswith(argument, "--half-norm=") || startswith(argument, "--half_norm=")
            option, value = split(argument, "="; limit=2)
            half_norm = parse_bool_option(value, option)
        elseif startswith(argument, "half_norm=")
            _, value = split(argument, "="; limit=2)
            half_norm = parse_bool_option(value, "half_norm")
        else
            error("Unknown option: $argument")
        end
        i += 1
    end
    normalise && half_norm && error("Choose either normalise=true or half_norm=true, not both")

    summary = load_slice_summary(run_dir)
    gap_counts = trivial_gap_counts(summary)
    included = findall(<=(max_gap_count), gap_counts)
    isempty(included) && error("No data points have trivial gap count <= $max_gap_count")
    step_edges = leading_step_edge_indices(gap_counts)
    fit_indices = filter(index -> gap_counts[index] <= max_gap_count, step_edges)
    isempty(fit_indices) && error("No step-edge points have trivial gap count <= $max_gap_count")
    x_fit = Float64.(summary.rho[fit_indices])
    normalisation_scale = normalise ? max(maximum(gap_counts), 1) :
        half_norm ? 2 * max(maximum(gap_counts), 1) : 1
    plotted_counts = Float64.(gap_counts) ./ normalisation_scale
    y_fit = plotted_counts[fit_indices]
    model, fit_label, equation = fit_curve(x_fit, y_fit, fit_type, poly_order)
    x_curve = collect(range(minimum(x_fit), maximum(x_fit); length=400))
    y_curve = model(x_curve)
    plotted_gap_counts = plot_excluded ? plotted_counts : plotted_counts[included]
    y_min, y_max = extrema(vcat(Float64.(plotted_gap_counts), y_curve))
    y_padding = max(0.5, 0.1 * (y_max - y_min))
    y_ticks = normalise ? (0.0:0.1:1.0) :
        half_norm ? (0.0:0.05:0.5) : (minimum(plotted_gap_counts):maximum(plotted_gap_counts))
    y_limits = normalise ? (0.0, 1.0) :
        half_norm ? (0.0, 0.5) : (y_min - y_padding, y_max + y_padding)
    y_label = normalise ? "Normalised trivial gap count" :
        half_norm ? "Half-normalised trivial gap count" : "Trivial gap count"

    plot_data = plot(
        xlabel="ρ = t₂/t₁",
        ylabel=y_label,
        title="Trivial gap count fit (included count ≤ $max_gap_count)",
        legend=:outerright,
        legendfontsize=8,
        size=(1000, 600),
        dpi=150,
        yticks=y_ticks,
        ylims=y_limits,
    )
    excluded = setdiff(1:length(gap_counts), included)
    if plot_excluded && !isempty(excluded)
        scatter!(plot_data, summary.rho[excluded], plotted_counts[excluded];
            label="Excluded points", color=:gray, alpha=0.55, markersize=4)
    end
    other_included = setdiff(included, fit_indices)
    if !isempty(other_included)
        scatter!(plot_data, summary.rho[other_included], plotted_counts[other_included];
            label="Other included points", color=:darkgreen, markersize=4)
    end
    scatter!(plot_data, x_fit, y_fit; label="Step-edge fit points", color=:orange, markersize=6, markerstrokecolor=:black)
    plot!(plot_data, x_curve, y_curve; label="$fit_label\n$equation", color=:black, linewidth=2)

    order_tag = fit_type == "poly" ? "_order$(poly_order)" : ""
    normalise_tag = normalise ? "_normalised" : half_norm ? "_half_norm" : ""
    output_path = joinpath(run_dir, "trivial_gap_fit_max$(max_gap_count)_$(fit_type)$(order_tag)$(normalise_tag).png")
    savefig(plot_data, output_path)
    println("Saved fit plot to $output_path")
    println("Fit includes $(length(fit_indices)) step-edge points from $(length(included)) points below the cutoff")
    normalise && println("Normalised plotted and fitted gap counts by $normalisation_scale")
    half_norm && println("Half-normalised plotted and fitted gap counts by $normalisation_scale")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end