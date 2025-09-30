
include("bson_unpacker.jl")

using DataFrames
using CSV

# """
#     filter_diagonal_band(df, outfile; lower_slope=0.5, upper_slope=2.0)

# Filter rows of a DataFrame `df` (with columns `:mu_t` and `:rho`)
# to keep only points within a diagonal band around the line `rho ≈ mu_t`.

# Arguments:
# - `outfile`: string, path to the CSV file to write
# - `lower_slope`: minimum ratio rho/mu_t to keep
# - `upper_slope`: maximum ratio rho/mu_t to keep

# Returns: the filtered DataFrame.
# Also writes a CSV with just the `:mu_t` and `:rho` values.
# """
# function filter_diagonal_band(df, outfile; lower_slope=0.5, upper_slope=2.0)
#     df_filtered = filter(row -> begin
#         mt, r = row.mu_t, row.rho
#         if mt == 0  # avoid division by zero
#             return false
#         end
#         ratio = r / mt
#         lower_slope ≤ ratio ≤ upper_slope
#     end, df)

#     # keep only the two columns we need
#     df_out = select(df_filtered, [:mu_t, :rho])

#     # write to CSV
#     CSV.write(outfile, df_out)

#     return df_filtered
# end

using DataFrames, CSV, Statistics, Plots

"""
    auto_band_filter(df, outfile; threshold=-0.05, quantiles=(0.05, 0.95))

Automatically find a diagonal band in (mu_t, rho) space where `:mp`
is below a given threshold (i.e. interesting region).

Steps:
- Keep only rows where mp < threshold
- Compute ratio rho/mu_t
- Use given quantiles of ratio distribution to define band
- Save mu_t, rho pairs to CSV

Returns: (filtered DataFrame, (q_low, q_high))
"""
function auto_band_filter(df, outfile; threshold=-0.05, quantiles=(0.05, 0.95))
    # Step 1: keep interesting rows
    df_interest = filter(:mp => x -> x < threshold, df)

    # Step 2: compute ratio distribution (avoid mu_t=0)
    ratios = [row.rho / row.mu_t for row in eachrow(df_interest) if row.mu_t != 0]

    if isempty(ratios)
        @warn "No interesting points found under threshold"
        return DataFrame(), (NaN, NaN)
    end

    # Step 3: find band bounds
    q_low, q_high = quantile(ratios, quantiles)

    # Step 4: filter whole dataframe into band
    df_band = filter(row -> begin
        mt, r = row.mu_t, row.rho
        mt == 0 && return false
        ratio = r / mt
        q_low ≤ ratio ≤ q_high
    end, df)

    # Step 5: save parameter pairs
    CSV.write(outfile, select(df_band, [:mu_t, :rho]))

    return df_band, (q_low, q_high)
end

"""
    plot_band(df, band; threshold=-0.05)

Plot mu_t vs rho heatmap colored by :mp,
with overlayed lines indicating the slope band.
"""
function plot_band(df, band; threshold=-0.05)
    q_low, q_high = band

    x = df.mu_t
    y = df.rho
    z = df.mp

    # Create the heatmap
    plt = Plots.heatmap(
        x, 
        y, 
        z, 
        xlabel=L"\mu/t_A", #string(x_variable), 
        ylabel=L"\rho", #string(y_variable), 
        color=:viridis,
        labelfontsize = 35,
        tickfontsize = 24,
        left_margin=10mm,
        bottom_margin=5mm,
        right_margin=10mm,
        top_margin=5mm,
        # clims=(minimum(value_matrix), maximum(value_matrix)),
        colorbar=true,
        grid=false,
        size=(1600, 1200)
    )

    # # overlay band lines (rho = slope * mu_t)
    # xline = range(minimum(x), stop=maximum(x), length=200)
    # plot!(plt, xline, q_low .* xline, color=:red, lw=2, label="band lower")
    # plot!(plt, xline, q_high .* xline, color=:red, lw=2, label="band upper")

    Plots.savefig(plt, "mu_rho_band_plot.png")
end



filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/mu_vs_rho_mp_heatmaps/GQC_N(200-200-1)_t1(1.0-1.0-51__t2(0.0-10.0-51)_mu(0.0-10.0-51)_Delta(0.1-0.1-1)/"
mp_tol = 0.1
df = unpack_bason_standard(filepath; mp_tol=mp_tol)

# filtered_df = filter_diagonal_band(df, "filtered_mu_rho.csv"; lower_slope=0.5, upper_slope=2.0)

# Filter and save
df_band, band = auto_band_filter(df, "refine_params.csv"; threshold=-0.1)

# Visualize
plot_band(df, band; threshold=-0.1)
