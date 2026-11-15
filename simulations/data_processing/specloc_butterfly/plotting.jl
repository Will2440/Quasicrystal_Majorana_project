"""
    file name:   plotting.jl
    created:     2026-09-22

    overview:
        Simplified plotting module for spectral localizer butterfly plots.
        Creates 2D heatmaps of specloc invariant on mu vs phi (slope) axes.
"""

module SpecLocPlotting

using Plots
using Colors
using ColorSchemes
using DataFrames
using Plots.PlotMeasures

export plot_specloc_butterfly, plot_specloc_winding


"""
    plot_specloc_butterfly(df::AbstractDataFrame, output_path::String; 
                           title_suffix="", plt_xlims=nothing, plt_ylims=nothing,
                           value_col::Symbol=:specloc_invariant,
                           mu_col::Symbol=:mu)

Plot spectral localizer invariant vs mu and phi (slope).
Creates a 2D heatmap with mu on x-axis, phi on y-axis, colored by value_col.

Args:
    value_col: Which column to use for colorbar values (default: :specloc_invariant)
    mu_col: Which column to use for x-axis (default: :mu, use :mu_specloc_invariant_norm for normalized)
"""
function plot_specloc_butterfly(
    df::AbstractDataFrame, 
    output_path::String; 
    title_suffix="", 
    plt_xlims=nothing, 
    plt_ylims=nothing,
    value_col::Symbol=:specloc_invariant,
    mu_col::Symbol=:mu,
    nonwhite_colour::AbstractString="#072f6b4d",
    nxn::Union{Nothing, Int}=nothing,
    linewidth::Real=1.0,
    marker::Symbol=:square
)
    
    # Get unique values and sort
    mus = sort(unique(df[!, mu_col]))
    phis = sort(unique(df.phi))
    
    # Create matrix for heatmap
    specloc_matrix = fill(NaN, length(phis), length(mus))
    
    # Pre-create lookup dicts (faster than repeated unique calls)
    mu_to_idx = Dict(mu => i for (i, mu) in enumerate(mus))
    phi_to_idx = Dict(phi => i for (i, phi) in enumerate(phis))
    
    # Check which column to use (do this once, not per row)
    use_invariant = "specloc_invariant" in names(df)
    
    # Extract column vectors for efficient access (avoid repeated column lookups)
    mu_col_data = df[!, mu_col]
    phi_col = df.phi
    
    # Use the specified value_col if it exists, otherwise fallback
    if string(value_col) in names(df)
        value_col_actual = value_col
        col_used = string(value_col)
    elseif use_invariant
        value_col_actual = :specloc_invariant
        col_used = "specloc_invariant"
    else
        value_col_actual = :specloc_gap
        col_used = "specloc_gap"
    end
    
    value_col_data = df[!, value_col_actual]
    colour_hex = startswith(nonwhite_colour, "#") ? nonwhite_colour[2:end] : nonwhite_colour
    occursin(r"^(?:[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$", colour_hex) ||
        throw(ArgumentError("nonwhite_colour must be a 6- or 8-digit hexadecimal RGB(A) value"))
    nonwhite_color = parse(Colorant, "#" * colour_hex)
    specloc_colormap = cgrad([:white, nonwhite_color])
    plot_size = isnothing(nxn) ? (1000, 800) : (nxn, nxn)

    if mu_col != :mu
        x_points = Float64[]
        y_points = Float64[]
        color_points = Float64[]
        for i in eachindex(mu_col_data, phi_col, value_col_data)
            x = mu_col_data[i]
            y = phi_col[i]
            value = value_col_data[i]
            if !ismissing(x) && !ismissing(y) && !ismissing(value) && isfinite(x) && isfinite(y) && isfinite(value)
                push!(x_points, Float64(x))
                push!(y_points, Float64(y))
                push!(color_points, Float64(value))
            end
        end

        p = scatter(
            x_points, y_points;
            marker_z=color_points,
            color=specloc_colormap,
            clims=(0.0, 1.0),
            colorbar=true,
            colorbar_title="Specloc Invariant",
            marker=marker,
            markersize=linewidth,
            markerstrokewidth=0,
            legend=false,
            xlabel="μ",
            ylabel="γ (slope)",
            title="Spectral Localizer Invariant $(title_suffix)",
            grid=false,
            size=plot_size,
            margin=10mm
        )
        if !isnothing(plt_xlims)
            xlims!(p, plt_xlims)
        end
        if !isnothing(plt_ylims)
            ylims!(p, plt_ylims)
        end
        savefig(p, output_path)
        println("Saved plot to: $output_path")
        return p
    end
    
    # Fill matrix using vectorized column access
    n_rows = nrow(df)
    @inbounds for i in 1:n_rows
        mu_idx = mu_to_idx[mu_col_data[i]]
        phi_idx = phi_to_idx[phi_col[i]]
        value = value_col_data[i]
        if !ismissing(value)
            specloc_matrix[phi_idx, mu_idx] = value
        end
    end
    
    # Calculate clim once (not twice)
    valid_mask = .!isnan.(specloc_matrix)
    if any(valid_mask)
        min_val = minimum(specloc_matrix[valid_mask])
        max_val = maximum(specloc_matrix[valid_mask])
        clim_val = (min_val, max_val)
    else
        clim_val = (-1.0, 1.0)  # fallback if all NaN
    end
    
    # Create plot with pre-computed clim
    p = heatmap(
        mus, phis, specloc_matrix,
        xlabel="μ",
        ylabel="γ (slope)",
        title="Spectral Localizer Invariant $(title_suffix)",
        color=specloc_colormap,
        clim=clim_val,
        colorbar_title="Specloc Invariant",
        size=plot_size,
        margin=10mm
    )
    
    # Set custom limits if provided
    if !isnothing(plt_xlims)
        xlims!(p, plt_xlims)
    end
    if !isnothing(plt_ylims)
        ylims!(p, plt_ylims)
    end
    
    # Save plot
    savefig(p, output_path)
    println("Saved plot to: $output_path")
    
    return p
end


"""
    plot_specloc_winding(df::DataFrame, output_path::String; 
                         title_suffix="", plt_xlims=nothing, plt_ylims=nothing)

Plot spectral localizer winding number vs mu and phi (slope).
Creates a 2D heatmap with mu on x-axis, phi on y-axis, colored by winding number.

NOTE: This is a utility function for potential future use. Currently, specloc_invariant
is the primary quantity. If specloc_winding data becomes available, this can be used.
"""
function plot_specloc_winding(
    df::AbstractDataFrame, 
    output_path::String; 
    title_suffix="", 
    plt_xlims=nothing, 
    plt_ylims=nothing
)
    
    # Get unique values and sort
    mus = sort(unique(df.mu))
    phis = sort(unique(df.phi))
    
    # Create matrix for heatmap
    winding_matrix = fill(NaN, length(phis), length(mus))
    
    # Pre-create lookup dicts
    mu_to_idx = Dict(mu => i for (i, mu) in enumerate(mus))
    phi_to_idx = Dict(phi => i for (i, phi) in enumerate(phis))
    
    # Extract column vectors for efficient access
    mu_col = df.mu
    phi_col = df.phi
    winding_col = df.specloc_winding
    
    # Fill matrix using vectorized column access
    n_rows = nrow(df)
    @inbounds for i in 1:n_rows
        mu_idx = mu_to_idx[mu_col[i]]
        phi_idx = phi_to_idx[phi_col[i]]
        value = winding_col[i]
        if !ismissing(value)
            winding_matrix[phi_idx, mu_idx] = Float64(value)
        end
    end
    
    # Create plot with discrete colormap for integer winding numbers
    p = heatmap(
        mus, phis, winding_matrix,
        xlabel="μ",
        ylabel="γ (slope)",
        title="Spectral Localizer Winding Number $(title_suffix)",
        color=:Set3,
        colorbar_title="Winding #",
        size=(1000, 800),
        margin=10mm,
        dpi=150
    )
    
    # Set custom limits if provided
    if !isnothing(plt_xlims)
        xlims!(p, plt_xlims)
    end
    if !isnothing(plt_ylims)
        ylims!(p, plt_ylims)
    end
    
    # Save plot
    savefig(p, output_path)
    println("Saved plot to: $output_path")
    
    return p
end


"""
    plot_discrete_specloc_projections(df::DataFrame, 
                                      value_col::Symbol,
                                      fixed_param::Symbol,
                                      output_path::String;
                                      kwargs...)

Create a discrete phase projection plot similar to IDOP phase projections.
"""
function plot_discrete_specloc_projections(
    df::AbstractDataFrame,
    value_col::Symbol,
    fixed_param::Symbol,
    output_path::String;
    plt_xlims=(-1.0, 1.0),
    plt_ylims=(0.0, 1.0),
    N=nothing,
    Delta=nothing,
    t_n=nothing,
    phason=nothing
)
    
    # Get unique fixed param values
    if fixed_param ∉ names(df)
        @warn "Parameter $fixed_param not in DataFrame"
        return nothing
    end
    
    unique_fixed = unique(df[!, fixed_param])
    
    # Create subplot
    p = scatter(
        df[!, :mu], df[!, value_col],
        xlabel="μ",
        ylabel=string(value_col),
        title="Specloc $(value_col) Projection",
        markersize=1,
        alpha=0.6,
        legend=false,
        size=(1000, 600),
        margin=10mm,
        dpi=150
    )
    
    if !isnothing(plt_xlims)
        xlims!(p, plt_xlims)
    end
    if !isnothing(plt_ylims)
        ylims!(p, plt_ylims)
    end
    
    savefig(p, output_path)
    println("Saved projection plot to: $output_path")
    
    return p
end


end  # end module SpecLocPlotting
