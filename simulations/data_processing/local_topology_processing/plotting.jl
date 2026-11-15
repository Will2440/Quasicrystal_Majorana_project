module Plotting

using Plots
using DataFrames
using Printf

function plt_spectral_function(df::DataFrame, omega_range::Vector{Float64}, savepath::String; fixed_values...)
    """
    Plots the spectral function as a function of energy for each site.

    Parameters:
    df           -- DataFrame containing the data (must include :spectral_function column).
    omega_range  -- Vector of energy values for the spectral function.
    savepath     -- String path to save the plot.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    if hasproperty(filtered_df, :mp)
        mp_value = filtered_df[!, :mp][1]
        title_str = @sprintf("Spectral Function A(ω) for μ=%.3f, MP=%.2f", fixed_values[:mu], mp_value)
    else
        title_str = @sprintf("Spectral Function A(ω) for μ=%.3f", fixed_values[:mu])
    end

    spectral_function = filtered_df[!, :spectral_function][1]  # Assuming one matching row
    site_indices = 1:size(spectral_function, 1)

    heatmap(omega_range, site_indices, spectral_function,
            xlabel="ω",
            ylabel="Site Index",
            title=title_str,
            color=:viridis,
            clims = (0, maximum(spectral_function)),
            colorbar=true
    )

    Plots.savefig(savepath)
end

function plt_local_mp_scales_middle_site(df::DataFrame, savepath::String; fixed_values...)
    """
    For each scale, plots the local MP for all local chains of that scale.
    x axis is the middle site of the local chain, y axis is the scale, colour is the local MP value.
    """
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    # if isempty(filtered_df)
    #     error("No results match the specified fixed values.")
    # end

    isempty(filtered_df) && error("No results match fixed_values=$(fixed_values)")
    nrow(filtered_df) == 1 || error("Expected exactly 1 row, got $(nrow(filtered_df)); fixed_values=$(fixed_values)")

    local_mp_data = filtered_df[!, :local_mp][1]

    scales = [data.scale for data in local_mp_data]
    middle_sites = [(data.start_site + data.end_site) ÷ 2 for data in local_mp_data]
    overall_mps = [data.overall_MP for data in local_mp_data]

    clims = (-1.0, 0.0)

    # println("scales: ", scales)
    # println("middle_sites: ", middle_sites)
    # println("overall_mps: ", overall_mps)
    # println("clims: ", clims)
    # println("lengths: ", (length(scales), length(middle_sites), length(overall_mps)))

    plt = Plots.scatter(
        middle_sites,
        scales,
        zcolor=overall_mps,
        color=:viridis,
        markershape=:square,
        markersize=4,
        markerstrokewidth=0,
        xlabel="Middle Site",
        ylabel="Scale (Subchain Length)",
        title="Local MP across Scales\n$(join(["$k=$v" for (k,v) in fixed_values], ", "))",
        colorbar_title="Overall MP",
        size=(800, 600),
        grid=false,
        clims=clims,
        framestyle=:box
    )

    Plots.savefig(plt, savepath)
end

# function plt_spectral_and_mp_scales(
#     df::DataFrame,
#     omega_range::Vector{Float64},
#     savepath::String;
#     fixed_values...
# )
#     # Filter for the desired row
#     filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
#     isempty(filtered_df) && error("No results match fixed_values=$(fixed_values)")
#     nrow(filtered_df) == 1 || error("Expected exactly 1 row, got $(nrow(filtered_df)); fixed_values=$(fixed_values)")

#     row = filtered_df[1, :]
#     # Spectral function: site index on x, omega on y
#     spectral_function = row.spectral_function
#     site_indices = 1:size(spectral_function, 1)
#     ω = omega_range

#     # Local MP data
#     local_mp_data = row.mp_local_scales
#     scales = [data.scale for data in local_mp_data]
#     middle_sites = [(data.start_site + data.end_site) ÷ 2 for data in local_mp_data]
#     overall_mps = [data.overall_MP for data in local_mp_data]

#     # Gather parameters for title/annotation
#     param_str = join(["$k=$(getproperty(row, k))" for k in propertynames(row) if k in (:N, :mu, :Delta, :t_n, :sequence_name, :mp)], ", ")
#     global_mp = hasproperty(row, :mp) ? row.mp : NaN

#     # Plot
#     plt = plot(layout = (1, 2), size=(1200, 500))

#     # Spectral function subplot (site index on x, omega on y)
#     heatmap!(
#         plt[1],
#         site_indices,
#         ω,
#         spectral_function',
#         xlabel = "Site Index",
#         ylabel = "ω",
#         title = "Spectral Function\n$param_str",
#         color = :viridis,
#         clims = (0, maximum(spectral_function)),
#         colorbar = true,
#         flipy = false
#     )

#     # Local MP scales subplot
#     scatter!(
#         plt[2],
#         middle_sites,
#         scales,
#         zcolor = overall_mps,
#         color = :viridis,
#         markershape = :square,
#         markersize = 4,
#         markerstrokewidth = 0,
#         xlabel = "Middle Site",
#         ylabel = "Scale (Subchain Length)",
#         title = "Local MP across Scales",
#         colorbar_title = "Overall MP",
#         grid = false,
#         clims = (-1.0, 1.0),
#         framestyle = :box
#     )

#     # Annotate global MP on the MP subplot
#     annotate!(plt[2], (maximum(middle_sites), maximum(scales)), text("Global MP = $(round(global_mp, digits=3))", :right, 12, :black))

#     savefig(plt, savepath)
# end

function plt_spectral_and_mp_scales(
    df::DataFrame,
    omega_range::Vector{Float64},
    savepath::String;
    fixed_values...
)
    # Filter for the desired row
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    isempty(filtered_df) && error("No results match fixed_values=$(fixed_values)")
    nrow(filtered_df) == 1 || error("Expected exactly 1 row, got $(nrow(filtered_df)); fixed_values=$(fixed_values)")

    row = filtered_df[1, :]
    # Spectral function: site index on x, omega on y
    spectral_function = row.spectral_function
    site_indices = 1:size(spectral_function, 1)
    ω = omega_range

    # Local MP data
    local_mp_data = row.mp_local_scales
    scales = [data.scale for data in local_mp_data]
    middle_sites = [(data.start_site + data.end_site) ÷ 2 for data in local_mp_data]
    overall_mps = [data.overall_MP for data in local_mp_data]

    # Gather parameters for title/annotation
    param_str = join(["$k=$(getproperty(row, k))" for k in propertynames(row) if k in (:N, :mu, :Delta, :t_n, :sequence_name)], ", ")
    global_mp = hasproperty(row, :mp) ? row.mp : NaN

    # Plot (Added margins and a global plot title)
    plt = plot(
        layout = (1, 2), 
        size = (1200, 550),
        plot_title = "Parameters: $param_str",
        plot_titlevspan = 0.1,
        margin = 10Plots.mm  # Gives enough breathing room for labels and colorbars
    )

    # Spectral function subplot (site index on x, omega on y)
    heatmap!(
        plt[1],
        site_indices,
        ω,
        spectral_function',
        xlabel = "Site Index",
        ylabel = "ω",
        title = "Spectral Function",
        color = :viridis,
        clims = (0, maximum(spectral_function)),
        colorbar = true,
        flipy = false
    )

    # Local MP scales subplot
    scatter!(
        plt[2],
        middle_sites,
        scales,
        zcolor = overall_mps,
        color = :viridis,
        markershape = :square,
        markersize = 4,
        markerstrokewidth = 0,
        xlabel = "Middle Site",
        ylabel = "Scale (Subchain Length)",
        title = "Local MP across Scales\n(Global MP = $(round(global_mp, digits=3)))", # Moved annotation here
        colorbar_title = "Overall MP",
        grid = false,
        # clims = (-1.0, 1.0),
        framestyle = :box
    )

    savefig(plt, savepath)
end

end # module Plotting