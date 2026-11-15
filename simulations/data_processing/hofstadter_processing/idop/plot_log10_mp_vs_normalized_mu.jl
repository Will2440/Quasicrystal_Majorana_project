using Colors
using BSON
using DataFrames
using Plots
using Printf
using Statistics

const MP_ZERO_TOL = 1e-3
const MU_C_REFERENCE_N = 1000
const MAX_PLOT_N = 800
const MP_FLOOR_QUANTILE = 0.05
const Y_AXIS_PADDING_FRACTION = 0.05
const X_AXIS_MAX = 1.05
const NORMALIZED_Y_LIMITS = (-0.05, 1.02)
const DATASET_NAME = "hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1-0.0_nbins1000_npb1_N(100-1000-10)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"
const RAW_DATA_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..", "raw_data", "Ns100-1000_for_mp_tol_plots", DATASET_NAME
))
const OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..", "results", "Ns100-1000_for_mp_tol_plots",
    DATASET_NAME, "log10_mp_vs_normalized_mu"
))

function load_results(folder::AbstractString)
    files = sort(filter(name -> endswith(name, ".bson"), readdir(folder)))
    isempty(files) && error("No BSON files found in $folder")

    frames = DataFrame[]
    for filename in files
        data = BSON.load(joinpath(folder, filename))
        haskey(data, :results_df) || error("Missing :results_df in $filename")
        push!(frames, DataFrame(data[:results_df]))
    end

    df = vcat(frames...; cols=:union)
    required = ("N", "t_n", "mu", "sequence_id", "mp")
    missing_columns = filter(column -> !(column in names(df)), required)
    isempty(missing_columns) || error("Missing required columns: $(join(missing_columns, ", "))")
    return df
end

function single_parameter_value(values, name::AbstractString; atol::Float64=1e-10)
    numeric_values = Float64.(values)
    reference = first(numeric_values)
    all(value -> isapprox(value, reference; atol=atol, rtol=0.0), numeric_values) ||
        error("Expected one $name value, found $(unique(numeric_values))")
    return reference
end

function find_sustained_zero_mu(mu_values, mp_values; tol::Float64=MP_ZERO_TOL)
    order = sortperm(mu_values)
    mus = Float64.(mu_values[order])
    mps = Float64.(mp_values[order])

    last_outside = findlast(mp -> !isfinite(mp) || abs(mp) > tol, mps)
    candidate = isnothing(last_outside) ? 1 : last_outside + 1
    candidate <= length(mus) || error("MP does not settle within tolerance $tol before the end of the mu range")
    all(mp -> isfinite(mp) && abs(mp) <= tol, @view mps[candidate:end]) ||
        error("MP is not continuously within tolerance $tol above the candidate mu")

    return mus[candidate]
end

function main()
    isdir(RAW_DATA_DIR) || error("Raw data directory not found: $RAW_DATA_DIR")
    df = load_results(RAW_DATA_DIR)
    isempty(df) && error("No rows found in $RAW_DATA_DIR")

    slope = single_parameter_value(getindex.(df.sequence_id, 1), "slope")
    t1 = single_parameter_value(getindex.(df.t_n, 1), "t1")
    t2 = single_parameter_value(getindex.(df.t_n, 2), "t2")
    rho = t2 / t1

    mkpath(OUTPUT_DIR)
    println("Loaded $(nrow(df)) rows from $(RAW_DATA_DIR)")
    println("slope=$slope, rho=$rho, MP zero tolerance=$MP_ZERO_TOL")
    println("Saving plots to $(OUTPUT_DIR)")

    reference_rows = filter(:N => ==(MU_C_REFERENCE_N), df)
    isempty(reference_rows) && error("Reference size N=$MU_C_REFERENCE_N not found in dataset")
    reference_curve = combine(groupby(reference_rows, :mu), :mp => mean => :mp)
    sort!(reference_curve, :mu)
    mu_c = find_sustained_zero_mu(reference_curve.mu, reference_curve.mp)
    isfinite(mu_c) && mu_c > 0 || error("Computed invalid mu_c=$mu_c from N=$MU_C_REFERENCE_N")
    println("Using mu_c=$mu_c from N=$MU_C_REFERENCE_N (MP tolerance=$MP_ZERO_TOL) for all sizes")

    overlay = plot(
        xlabel="",
        ylabel="",
        title="",
        xlims=(0.0, X_AXIS_MAX),
        label="",
        size=(1100, 680),
        dpi=160,
        framestyle=:box,
        grid=true,
        gridalpha=0.25,
        legend=false,
    )
    range_overlay = plot(
        xlabel="",
        ylabel="",
        title="",
        xlims=(0.0, X_AXIS_MAX),
        ylims=NORMALIZED_Y_LIMITS,
        yticks=(collect(range(0.0, 1.0; length=6)), string.(collect(range(-1.0, 0.0; length=6)))),
        label="",
        size=(1100, 680),
        dpi=160,
        framestyle=:box,
        grid=true,
        gridalpha=0.25,
        legend=false,
    )
    range_y_min = Inf
    range_y_max = -Inf

    n_values = sort(filter(n_value -> n_value <= MAX_PLOT_N, unique(df.N)))
    isempty(n_values) && error("No system sizes found at or below N=$MAX_PLOT_N")
    grayscale = [RGB(shade, shade, shade) for shade in range(0.82, 0.12; length=length(n_values))]
    for (n_index, n_value) in enumerate(n_values)
        grouped_rows = filter(:N => ==(n_value), df)
        plot_rows = combine(groupby(grouped_rows, :mu), :mp => mean => :mp)
        sort!(plot_rows, :mu)

        mu_normalized = plot_rows.mu ./ mu_c
        log10_mp_plus_one = [
            isfinite(mp) && 1.0 + mp > 0.0 ? log10(1.0 + mp) : NaN for mp in plot_rows.mp
        ]
        visible_mask = (0.0 .<= mu_normalized .<= 1.0) .& isfinite.(log10_mp_plus_one)
        visible_log_values = log10_mp_plus_one[visible_mask]
        isempty(visible_log_values) && error("No finite plotted values for N=$n_value in mu/mu_c ∈ [0, 1]")
        floor_value = quantile(visible_log_values, MP_FLOOR_QUANTILE)
        floor_normalized_values = visible_log_values .- floor_value
        floor_value < 0.0 || error("Estimated floor must be below zero for N=$n_value, got $floor_value")
        range_normalized_values = floor_normalized_values ./ -floor_value
        range_y_min = min(range_y_min, minimum(range_normalized_values))
        range_y_max = max(range_y_max, maximum(range_normalized_values))
        plot!(
            overlay,
            mu_normalized[visible_mask],
            floor_normalized_values;
            label="N=$n_value",
            linewidth=1.6,
            color=grayscale[n_index],
        )
        plot!(
            range_overlay,
            mu_normalized[visible_mask],
            range_normalized_values;
            label="N=$n_value",
            linewidth=1.6,
            color=grayscale[n_index],
        )
        println("N=$n_value: estimated log10(MP + 1) floor (q=$(MP_FLOOR_QUANTILE)) = $floor_value")

        y_lower_limit = n_value == 400 ? -7.0 :
                        n_value == 600 ? -10.5 :
                        n_value == 800 ? -14.0 : nothing
        finite_y_values = filter(isfinite, log10_mp_plus_one)
        data_y_min = minimum(finite_y_values)
        data_y_max = maximum(finite_y_values)
        y_range_min = isnothing(y_lower_limit) ? data_y_min : y_lower_limit
        y_range = max(data_y_max - y_range_min, eps(Float64))
        y_padding = Y_AXIS_PADDING_FRACTION * y_range
        plot_y_limits = isnothing(y_lower_limit) ?
            (data_y_min - y_padding, data_y_max + y_padding) :
            (y_range_min, data_y_max + y_padding)
        y_tick_limits = plot_y_limits
        y_ticks = floor(Int, y_tick_limits[1]):ceil(Int, y_tick_limits[2])
        p = plot(
            mu_normalized,
            log10_mp_plus_one;
            xlabel="",
            ylabel="",
            yticks=y_ticks,
            ylims=plot_y_limits,
            xlims=(0.0, 1.0),
            title="",
            label="",
            legend=false,
            linewidth=1.5,
            color=:steelblue,
            size=(900, 560),
            dpi=160,
            framestyle=:box,
            grid=true,
            gridalpha=0.25,
        )
        # vline!(p, [1.0]; color=:firebrick, linestyle=:dash, linewidth=1.2, label="mu_c")

        output_path = joinpath(OUTPUT_DIR, "log10_mp_vs_mu_normalized_N$(n_value).svg")
        savefig(p, output_path)
        println("Saved $(output_path)")
    end

    hline!(overlay, [0.0]; color=:black, linestyle=:dash, linewidth=1.0, label="Estimated floor")
    overlay_path = joinpath(OUTPUT_DIR, "log10_mp_floor_normalized_overlay.svg")
    savefig(overlay, overlay_path)
    println("Saved combined floor-normalized overlay to $(overlay_path)")

    hline!(range_overlay, [0.0]; color=:black, linestyle=:dash, linewidth=1.0, label="Estimated floor")
    hline!(range_overlay, [1.0]; color=:gray, linestyle=:dot, linewidth=1.0, label="log10(MP + 1) = 0")
    range_overlay_path = joinpath(OUTPUT_DIR, "log10_mp_range_normalized_overlay.svg")
    savefig(range_overlay, range_overlay_path)
    println("Saved combined range-normalized overlay to $(range_overlay_path)")
end

main()