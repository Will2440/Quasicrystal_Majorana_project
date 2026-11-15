using Colors
using BSON
using DataFrames
using JLD2
using Plots
using Printf
using Statistics

const MP_FLOOR_QUANTILE = 0.05
const DEFAULT_SHADE_ALPHA = 1.0
const TARGET_RHO = 1.5
const Y_AXIS_PADDING_BELOW = 0.05
const Y_AXIS_PADDING_ABOVE = 0.02
const DATA_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "MP_verying_Delta"))
const DATASET_NAME = "hof_style_slopes_N500_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p01-0p1-10)_x0(218p75-281p25-5)_kappa(0p0001-0p0001-1)"
const RAW_DATA_DIR = normpath(joinpath(
    DATA_ROOT, DATASET_NAME
))
const OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..", "results", "MP_verying_Delta", DATASET_NAME,
    "mp_vs_normalized_mu_delta_rho"
))

function load_results(folder::AbstractString)
    files = sort(filter(name -> endswith(name, ".jld2") || endswith(name, ".bson"), readdir(folder)))
    isempty(files) && error("No JLD2 or BSON files found in $folder")

    frames = DataFrame[]
    for filename in files
        path = joinpath(folder, filename)
        data = endswith(filename, ".bson") ? BSON.load(path) : JLD2.load(path)
        key = endswith(filename, ".bson") ? :results_df : "results_df"
        haskey(data, key) || error("Missing results_df in $filename")
        push!(frames, DataFrame(data[key]))
    end

    df = vcat(frames...; cols=:union)
    required = ("N", "t_n", "mu", "Delta", "sequence_id", "mp")
    missing_columns = filter(column -> !(column in names(df)), required)
    isempty(missing_columns) || error("Missing required columns: $(join(missing_columns, ", "))")
    return df
end

analytic_mu_c(slope::Real, t1::Real, t2::Real) =
    2.0 * ((1.0 - slope) * t1 + slope * t2)

function gradient_colors(count::Int, start_color, end_color)
    count == 1 && return [start_color]
    return [RGB(
        (1 - fraction) * red(start_color) + fraction * red(end_color),
        (1 - fraction) * green(start_color) + fraction * green(end_color),
        (1 - fraction) * blue(start_color) + fraction * blue(end_color),
    ) for fraction in range(0.0, 1.0; length=count)]
end

format_value(value::Real) = @sprintf("%.3f", value)

function contiguous_true_ranges(mask)
    ranges = UnitRange{Int}[]
    range_start = nothing
    for index in eachindex(mask)
        if mask[index] && isnothing(range_start)
            range_start = index
        elseif !mask[index] && !isnothing(range_start)
            push!(ranges, range_start:(index - 1))
            range_start = nothing
        end
    end
    !isnothing(range_start) && push!(ranges, range_start:lastindex(mask))
    return ranges
end

function make_overlay(curve_group, fixed_parameter::Symbol, fixed_value::Float64,
                      varying_parameter::Symbol,
                      colors, shade_alpha::Float64, output_path::AbstractString)
    sort!(curve_group; by=curve -> getproperty(curve, varying_parameter))
    normalized_curves = [merge(curve, (
        y_normalized=(curve.log_values .- curve.floor_value) ./ -curve.floor_value,
    )) for curve in curve_group]
    p = plot(
        xlabel="",
        ylabel="",
        title="",
        xlims=(0.0, 1.0 + Y_AXIS_PADDING_BELOW),
        ylims=(-Y_AXIS_PADDING_BELOW, 1.0 + Y_AXIS_PADDING_ABOVE),
        yticks=(collect(range(0.0, 1.0; length=6)), string.(collect(range(-1.0, 0.0; length=6)))),
        label="",
        size=(1050, 680),
        dpi=170,
        framestyle=:box,
        grid=true,
        gridalpha=0.22,
        legend=false,
    )

    draw_order = fixed_parameter == :Delta ? reverse(eachindex(normalized_curves)) : eachindex(normalized_curves)
    for index in draw_order
        curve = normalized_curves[index]
        color = colors[index]
        finite_curve = isfinite.(curve.y_normalized)
        for segment in contiguous_true_ranges(finite_curve)
            length(segment) < 2 && continue
            plot!(
                p,
                curve.mu_normalized[segment],
                curve.y_normalized[segment];
                fillrange=-Y_AXIS_PADDING_BELOW,
                fillcolor=color,
                fillalpha=shade_alpha,
                linealpha=0.0,
                linewidth=0,
                label="",
            )
        end
        plot!(
            p,
            curve.mu_normalized,
            curve.y_normalized;
            color=color,
            linewidth=1.8,
            label="",
        )
    end

    hline!(p, [0.0]; color=:gray35, linestyle=:dot, linewidth=0.9, label="")
    hline!(p, [1.0]; color=:black, linestyle=:dash, linewidth=1.0, label="")
    savefig(p, output_path)
    println("Saved $(output_path)")
end

function main()
    length(ARGS) <= 1 || error("Usage: julia plot_mp_vs_normalized_mu_delta_rho.jl [shade_alpha]")
    shade_alpha = isempty(ARGS) ? DEFAULT_SHADE_ALPHA : parse(Float64, only(ARGS))
    isfinite(shade_alpha) && 0.0 <= shade_alpha <= 1.0 ||
        error("shade_alpha must be a finite value between 0 and 1")

    isdir(RAW_DATA_DIR) || error("Raw data directory not found: $RAW_DATA_DIR")
    df = load_results(RAW_DATA_DIR)
    isempty(df) && error("No rows found in $RAW_DATA_DIR")
    n_values = unique(Int.(df.N))
    length(n_values) == 1 || error("Expected fixed N, found $(sort(n_values))")
    n_value = only(n_values)
    df.t1 = getindex.(df.t_n, 1)
    df.t2 = getindex.(df.t_n, 2)
    df.rho = df.t2 ./ df.t1
    df = filter(:rho => rho -> isapprox(rho, TARGET_RHO; atol=1e-10, rtol=0.0), df)
    isempty(df) && error("No rows found for rho=$TARGET_RHO")
    slope = only(unique(Float64.(getindex.(df.sequence_id, 1))))

    grouped = combine(
        groupby(df, [:Delta, :rho, :mu]),
        :mp => mean => :mp,
        :t1 => first => :t1,
        :t2 => first => :t2,
    )

    curves = NamedTuple[]
    for group in groupby(grouped, [:Delta, :rho])
        sort!(group, :mu)
        t1 = Float64(first(group.t1))
        t2 = Float64(first(group.t2))
        mu_c = analytic_mu_c(slope, t1, t2)
        isfinite(mu_c) && mu_c > 0 || error("Computed invalid mu_c=$mu_c for rho=$(first(group.rho))")
        mu_normalized = Float64.(group.mu) ./ mu_c
        log_values = [
            isfinite(mp) && 1.0 + mp > 0.0 ? log10(1.0 + mp) : NaN for mp in group.mp
        ]
        visible = (0.0 .<= mu_normalized .<= 1.0) .& isfinite.(log_values)
        floor_samples = log_values[visible]
        isempty(floor_samples) && error(
            "No finite floor samples for Delta=$(first(group.Delta)), rho=$(first(group.rho))"
        )
        floor_value = quantile(floor_samples, MP_FLOOR_QUANTILE)
        floor_value < 0.0 || error(
            "Estimated floor must be below zero for Delta=$(first(group.Delta)), rho=$(first(group.rho)); got $floor_value"
        )
        push!(curves, (
            Delta=Float64(first(group.Delta)),
            rho=Float64(first(group.rho)),
            mu_c=mu_c,
            mu_normalized=mu_normalized,
            floor_value=floor_value,
            log_values=log_values,
        ))
    end

    mkpath(OUTPUT_DIR)
    println("Loaded $(nrow(df)) rows; N=$n_value, $(length(curves)) Delta/rho curves")
    println("Plotting fixed rho=$TARGET_RHO over Delta values: $(join(format_value.(sort(unique(getproperty.(curves, :Delta)))), ", "))")
    println("Using slope=$slope and mu_c(rho)=2*((1-slope)*t1+slope*t2)")
    println("Per-curve floor estimates use quantile=$MP_FLOOR_QUANTILE: $(join(["Delta=$(format_value(curve.Delta)), rho=$(format_value(curve.rho)): $(format_value(curve.floor_value))" for curve in curves], "; "))")
    println("Each curve is normalized from its own estimated floor to log10(MP + 1)=0")
    println("Saving plots to $OUTPUT_DIR")

    deltas = sort(unique(getproperty.(curves, :Delta)))
    delta_colors = gradient_colors(length(deltas), RGB(0.61, 0.04, 0.08), RGB(1.0, 0.87, 0.83))
    sorted_curves = sort(curves; by=curve -> curve.Delta)
    path = joinpath(OUTPUT_DIR, "fixed_rho_$(replace(format_value(TARGET_RHO), "." => "p"))_overlay.svg")
    make_overlay(sorted_curves, :rho, TARGET_RHO, :Delta, delta_colors, shade_alpha, path)
end

main()