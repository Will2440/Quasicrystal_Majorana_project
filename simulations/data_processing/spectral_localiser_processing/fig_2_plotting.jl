using DataFrames
using JLD2
using Glob
using Plots
using ProgressMeter
using Statistics


data_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/spec_loc_final/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(200-200-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-301)_Delta(0.05-0.05-1)_kappa(0.005-0.015-11)_x0s(100.0-100.0-1)"

output_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/spec_loc_final/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(200-200-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-301)_Delta(0.05-0.05-1)_kappa(0.005-0.015-11)_x0s(100.0-100.0-1)"
isdir(output_path) || mkpath(output_path)
println("Output path: $output_path")

function process_jld2_files_fluid(folder_path::String)
    file_paths = glob("*.jld2", folder_path)
    
    # Pre-allocate a vector to hold all the dataframes
    dfs = DataFrame[]
    failed_files = String[]

    @showprogress for file_path in file_paths
        jld_data = nothing
        try
            jld_data = load(file_path)
        catch err
            push!(failed_files, file_path)
            @warn "Failed to load JLD2 file (skipping)" file=file_path error=err
            continue
        end

        if haskey(jld_data, "results_df")
            push!(dfs, DataFrame(jld_data["results_df"]))
        else
            @warn "JLD2 file missing results_df key (skipping)" file=file_path
        end
    end

    println("Number of files read: $(length(dfs)) / $(length(file_paths))")
    if !isempty(failed_files)
        @warn "Some JLD2 files could not be processed" failed_count=length(failed_files) failed_files=failed_files
    end
    
    # Combine everything dynamically
    combined_dataframe = isempty(dfs) ? DataFrame() : vcat(dfs..., cols=:union)

    # Print the resulting column names dynamically
    println("\n=== DataFrame Loading Complete ===")
    println("DataFrame column names:")
    for col in names(combined_dataframe)
        println("  - :$col")
    end
    println("==================================\n")

    return combined_dataframe
end

df = process_jld2_files_fluid(data_path)

function plt_spectrum_and_sig(
    df::DataFrame,
    filename::String;
    x_var::Symbol=:mu,
    log_spectrum::Bool=false,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    no_of_evals::Union{Int, Nothing}=nothing,
    fixed_values...
)

    # approx-aware filtering on row-level context columns
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    filtered_rows = isempty(fixed_values) ? df :
                    filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_rows)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
        end
        error("No data points match the specified fixed values.")
    end

    # Sort by x_var
    sorted_rows = sort(filtered_rows, x_var)
    
    # Extract x-axis data
    x_data = sorted_rows[!, x_var]
    
    # Extract spectrum data and invariant
    spectrum_data = sorted_rows.specloc_spectrum
    invariant_col = :specloc_winding #:specloc_gap #
    invariant_data = sorted_rows[!, invariant_col]
    
    # Transform spectrum to log scale if requested
    if log_spectrum
        spectrum_data = [s isa Vector ? sign.(s) .* log.(abs.(s) .+ 1e-15) : s for s in spectrum_data]
    end

    if no_of_evals !== nothing
        spectrum_data = [s isa Vector ? s[1:min(no_of_evals, length(s))] : s for s in spectrum_data]
    end
    
    # Create plot with two y-axes
    y_label = log_spectrum ? "log(|Spectrum|) (signed)" : "Spectral Localiser Spectrum"
    plt = plot(
        xlabel=string(x_var),
        ylabel=y_label,
        title="Spectral Localiser Spectrum and $(invariant_col) vs $(x_var)"
    )
    
    # Plot spectrum eigenvalues on primary y-axis
    for (i, spectrum) in enumerate(spectrum_data)
        if spectrum isa Vector
            scatter!(plt, fill(x_data[i], length(spectrum)), spectrum, 
                  label=i == 1 ? "Spectrum" : "", 
                #   linewidth=1,
                  markerstrokewidth=0,
                  marker=:circle,
                  markersize=1,
                  color=:blue,
                  alpha=0.6,
                  ylims = (-0.5, 0.5)
                )
        end
    end
    
    # Create secondary y-axis and plot invariant
    ax2 = twinx(plt)
    plot!(ax2, x_data, invariant_data,
        ylabel=string(invariant_col),
        label=string(invariant_col),
        # markerstrokewidth=0,
        # marker=:circle,
        # markersize=1,
        linewidth=2,
        color=:red,
        legend=:topright)
    
    savefig(plt, filename)

end

function plt_spectrum_and_sig_shaded(
    df::DataFrame,
    filename::String;
    x_var::Symbol=:mu,
    log_spectrum::Bool=false,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    fixed_values...
)

    # approx-aware filtering on row-level context columns
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    filtered_rows = isempty(fixed_values) ? df :
                    filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_rows)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
        end
        error("No data points match the specified fixed values.")
    end

    # Sort by x_var
    sorted_rows = sort(filtered_rows, x_var)
    
    # Extract x-axis data
    x_data = sorted_rows[!, x_var]
    
    # Extract spectrum data and invariant
    spectrum_data = sorted_rows.specloc_spectrum
    invariant_col = :specloc_winding
    invariant_data = sorted_rows[!, invariant_col]
    
    # Transform spectrum to log scale if requested
    if log_spectrum
        spectrum_data = [s isa Vector ? sign.(s) .* log.(abs.(s) .+ 1e-15) : s for s in spectrum_data]
    end
    
    # Determine y-axis range from all spectrum data
    all_spectrum_vals = vcat([s for s in spectrum_data if s isa Vector]...)
    y_min, y_max = minimum(all_spectrum_vals), maximum(all_spectrum_vals)
    y_range = range(y_min, y_max, length=100)
    
    # Create shading matrix: 1 where invariant==1, 0 otherwise
    shading_matrix = fill(0.0, length(y_range), length(x_data))
    for (j, inv_val) in enumerate(invariant_data)
        if isapprox(inv_val, 1.0; atol=atol, rtol=rtol)
            shading_matrix[:, j] .= 1.0
        end
    end
    
    # Create plot with heatmap background showing invariant==1 regions
    y_label = log_spectrum ? "log(|Energy|) (signed)" : "Energy"
    plt = heatmap(
        x_data, y_range, shading_matrix,
        xlabel="", #string(x_var),
        ylabel="", #y_label,
        title="", #"Spectral Localiser Spectrum and $(invariant_col) vs $(x_var)",
        colorbar=false,
        color=:Blues,
        ylims=(-0.5, 0.5),
        xlims=(0.0, 3.0),
        axis=false,
        legend=false,
        size=(1000,200)
    )
    
    # Sort eigenvalues at each x_var value and prepare for line plotting
    sorted_spectra = [sort(s) for s in spectrum_data]
    
    # Get maximum number of eigenvalues
    max_evals = maximum(length(s) for s in sorted_spectra if s isa Vector)
    
    # Create line for each eigenvalue branch (index)
    for eval_idx in 1:max_evals
        x_branch = Float64[]
        y_branch = Float64[]
        
        for (i, spectrum) in enumerate(sorted_spectra)
            if spectrum isa Vector && eval_idx <= length(spectrum)
                push!(x_branch, x_data[i])
                push!(y_branch, spectrum[eval_idx])
            end
        end
        
        # Plot line for this eigenvalue branch
        if !isempty(x_branch)
            plot!(plt, x_branch, y_branch,
                  label=eval_idx == 1 ? "Eigenvalue Bands" : "",
                  linewidth=2.0,
                  color=:black,
                  alpha=0.7)
        end
    end
    
    # Draw dashed line at y=0
    hline!(plt, [0.0], linestyle=:dash, color=:red, label="y=0")

    savefig(plt, filename)

end

function plt_log_gap(
    df::DataFrame,
    filename::String;
    x_var::Symbol=:mu,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    fixed_values...
)

    # approx-aware filtering on row-level context columns
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    filtered_rows = isempty(fixed_values) ? df :
                    filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_rows)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
        end
        error("No data points match the specified fixed values.")
    end

    # Sort by x_var
    sorted_rows = sort(filtered_rows, x_var)
    
    # Extract x-axis data
    x_data = sorted_rows[!, x_var]
    
    # Extract spectrum data and invariant
    spectrum_data = sorted_rows.specloc_spectrum
    gap_col = :specloc_gap
    gap_data = sorted_rows[!, gap_col]
    log_gap_data = log.(abs.(gap_data) .+ 1e-15)

    # create plot of log gap vs x_var
    plt = plot(
        x_data, log_gap_data,
        xlabel=string(x_var),
        ylabel="", #"log(|Gap|)",
        title="", #"Log Gap vs $(x_var)",
        linewidth=2.0,
        color=:blue,
        legend=false,
        size=(1000,200),
        axis=false,
        xlims=(0.0,3.0),
        ylims=(-8.0, 0.0)
    )
    savefig(plt, filename)

end

## for trying to average over x0s, but not useful
# function plt_spectrum_and_sig_x0_averaged(
#     df::DataFrame,
#     filename::String;
#     x_var::Symbol=:mu,
#     type::Symbol=:gap,
#     log_spectrum::Bool=false,
#     show_std::Bool=true,
#     atol::Float64=1e-8,
#     rtol::Float64=1e-5,
#     verbose::Bool=true,
#     fixed_values...
# )

#     # approx-aware filtering on row-level context columns (excluding x0)
#     _kv_match(row, key, val) = begin
#         rv = getproperty(row, key)
#         (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
#     end
#     _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

#     filtered_rows = isempty(fixed_values) ? df :
#                     filter(row -> _row_matches(row, fixed_values), df)

#     if isempty(filtered_rows)
#         if verbose
#             @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
#         end
#         error("No data points match the specified fixed values.")
#     end

#     # Group by x_var to collect all x0 samples at each x_var value
#     grouped = groupby(filtered_rows, x_var)
    
#     # Extract unique x_var values
#     x_values = unique(filtered_rows[!, x_var])
#     sort!(x_values)
    
#     # Collect all spectrum points for density plot
#     all_x_points = Float64[]
#     all_y_points = Float64[]
    
#     # Collect statistics for invariant
#     x_data = Float64[]
#     invariant_mean = Float64[]
#     invariant_std = Float64[]
    
#     invariant_col = :specloc_winding
    
#     for x_val in x_values
#         # Get all rows with this x_var value
#         matching_rows = filter(row -> isapprox(getproperty(row, x_var), x_val; atol=atol, rtol=rtol), filtered_rows)
        
#         if isempty(matching_rows)
#             continue
#         end
        
#         # Collect all spectrum values across x0s at this x_var value
#         for row in eachrow(matching_rows)
#             spectrum = row.specloc_spectrum
#             if spectrum isa Vector
#                 # Transform to log scale if requested
#                 if log_spectrum
#                     spectrum = sign.(spectrum) .* log.(abs.(spectrum) .+ 1e-15)
#                 end
#                 # Add all spectrum points
#                 append!(all_x_points, fill(x_val, length(spectrum)))
#                 append!(all_y_points, spectrum)
#             end
#         end
        
#         # Compute invariant statistics
#         inv_values = [row[invariant_col] for row in eachrow(matching_rows)]
#         push!(x_data, x_val)
#         push!(invariant_mean, mean(inv_values))
#         push!(invariant_std, std(inv_values))
#     end
    
#     # Determine y-axis range
#     y_min, y_max = minimum(all_y_points), maximum(all_y_points)
    
#     # Create 2D histogram for density
#     x_bins = length(x_values)
#     y_bins = 100
#     y_edges = range(y_min, y_max, length=y_bins+1)
#     x_edges = vcat(x_values[1] - (x_values[2]-x_values[1])/2, 
#                    [(x_values[i] + x_values[i+1])/2 for i in 1:length(x_values)-1],
#                    x_values[end] + (x_values[end]-x_values[end-1])/2)
    
#     # Build density matrix
#     density_matrix = zeros(y_bins, x_bins)
#     for (x_pt, y_pt) in zip(all_x_points, all_y_points)
#         x_idx = searchsortedfirst(x_values, x_pt)
#         y_idx = min(searchsortedfirst(y_edges, y_pt), y_bins)
#         if x_idx >= 1 && x_idx <= x_bins && y_idx >= 1 && y_idx <= y_bins
#             density_matrix[y_idx, x_idx] += 1
#         end
#     end
    
#     # Normalize density matrix by column (each x_var value)
#     for i in 1:x_bins
#         col_sum = sum(density_matrix[:, i])
#         if col_sum > 0
#             density_matrix[:, i] ./= col_sum
#         end
#     end
    
#     # Create plot with density heatmap
#     y_label = log_spectrum ? "log(|Energy|) (signed)" : "Energy"
#     plt = heatmap(
#         x_values,
#         collect(range(y_min, y_max, length=y_bins)),
#         density_matrix,
#         xlabel=string(x_var),
#         ylabel=y_label,
#         title="Spectral Localiser (x0-averaged) and $(invariant_col) vs $(x_var)",
#         colorbar_title="Density",
#         color=:viridis,
#         alpha=0.8,
#         ylims=(-0.5, 0.5)
#     )
    
#     # Create secondary y-axis for invariant
#     ax2 = twinx(plt)
    
#     # Plot mean invariant
#     plot!(ax2, x_data, invariant_mean,
#           ylabel=string(invariant_col),
#           label="Mean $(invariant_col)",
#           linewidth=2,
#           color=:red,
#           legend=:topright)
    
#     # Add shaded region for standard deviation if requested
#     if show_std && !all(iszero, invariant_std)
#         plot!(ax2, x_data, invariant_mean,
#               ribbon=invariant_std,
#               fillalpha=0.3,
#               fillcolor=:red,
#               label="",
#               linewidth=0)
#     end
#     println("made plot: $(filename)")
    
#     savefig(plt, filename)

# end

kappas = unique(df.kappa)
t_ns = unique(df.t_n)

x_var = :mu
# kappas=[0.01]

# for t_n in t_ns
#     for kappa in kappas
#         plt_spectrum_and_sig(
#             df, 
#             joinpath(output_path, "spectrum_and_sig_vs_$(x_var)_kappa_$(kappa)_t_n_$(t_n).png"), 
#             x_var=x_var,
#             log_spectrum=false,
#             kappa=kappa,
#             t_n=t_n
#         )
#     end
# end

for t_n in t_ns
    for kappa in kappas
        plt_spectrum_and_sig_shaded(
            df, 
            joinpath(output_path, "spectrum_and_sig_shaded_vs_$(x_var)_kappa_$(kappa)_t_n_$(t_n).svg"), 
            x_var=x_var,
            log_spectrum=false,
            kappa=kappa,
            t_n=t_n
        )

        plt_log_gap(
            df, 
            joinpath(output_path, "bare_log_gap_vs_$(x_var)_kappa_$(kappa)_t_n_$(t_n).svg"), 
            x_var=x_var,
            kappa=kappa,
            t_n=t_n
        )
    end
end
