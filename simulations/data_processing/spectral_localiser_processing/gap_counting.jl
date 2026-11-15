using DataFrames
using CSV
using JLD2
using Glob
using Plots
using ProgressMeter
using Statistics
using Printf

data_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/spec_loc_qc-sc_comp/hof_style_slopes_N500_phason_0.0-1-0.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p0-2p5-501)_mu(0p0-4p5-45001)_Delta(0p05-0p05-1)"

output_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/spec_loc_qc-sc_comp/hof_style_slopes_N500_phason_0.0-1-0.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p0-2p5-501)_mu(0p0-4p5-45001)_Delta(0p05-0p05-1)"
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

function add_rho_to_df!(df::DataFrame)
    df.rho = map(x -> x[2] / x[1], df.t_n)
end

add_rho_to_df!(df)

kappas = unique(df.kappa)
Ns = unique(df.N)
x0s = unique(df.x0)
t_ns = unique(df.t_n)
rhos = unique(df.rho)
mus = unique(df.mu)
Deltas = unique(df.Delta)
sequence_names = unique(df.sequence_name)
sequence_ids = unique(df.sequence_id)
t_1s = unique(map(x -> x[1], df.t_n))
t_2s = unique(map(x -> x[2], df.t_n))

println("Unique values: min-length-max")
println("  kappas: $(minimum(kappas)) - $(length(kappas)) - $(maximum(kappas))")
println("  Ns: $(minimum(Ns)) - $(length(Ns)) - $(maximum(Ns))")
println("  x0s: $(minimum(x0s)) - $(length(x0s)) - $(maximum(x0s))")
println(" t_1s: $(minimum(t_1s)) - $(length(t_1s)) - $(maximum(t_1s))")
println(" t_2s: $(minimum(t_2s)) - $(length(t_2s)) - $(maximum(t_2s))")
println("  rhos: $(minimum(rhos)) - $(length(rhos)) - $(maximum(rhos))")
println("  mus: $(minimum(mus)) - $(length(mus)) - $(maximum(mus))")
println("  Deltas: $(minimum(Deltas)) - $(length(Deltas)) - $(maximum(Deltas))")
println("  sequence_names: $sequence_names")
println("  sequence_ids: $sequence_ids")


# function group_and_concatenate_data(
#     df::DataFrame;
#     group_cols::Vector{Symbol}=[:N, :t_n, :Delta, :sequence_name, :sequence_id, :kappa, :x0],
#     verbose::Bool=true
# )
#     """
#     Group dataframe by fixed parameters and concatenate mu-dependent columns into vectors.
    
#     For each unique combination of group parameters, all mu values and corresponding
#     specloc_gap and specloc_winding values are concatenated into Vector{Float64},
#     maintaining their original order relationship (index i of mus corresponds to 
#     index i of specloc_gaps and specloc_windings).
    
#     Additionally, counts the number of contiguous regions where specloc_gap == 0.0.
    
#     Args:
#         df: DataFrame with columns :mu, :specloc_gap, :specloc_winding, and grouping columns
#         group_cols: Columns to group by
#         verbose: Print information messages
    
#     Returns:
#         DataFrame with one row per group, containing:
#         - Group parameters (from group_cols)
#         - :mus (Vector{Float64}) - sorted mu values for this group
#         - :specloc_gaps (Vector{Float64}) - corresponding specloc_gap values
#         - :specloc_windings (Vector{Float64}) - corresponding specloc_winding values
#         - :gap_count (Int64) - number of contiguous regions where specloc_gap == 0.0
#     """
    
#     # Filter to only include group_cols that exist in the dataframe
#     # Note: names(df) returns Strings, but group_cols are Symbols, so convert for comparison
#     df_col_names_symbols = Symbol.(names(df))
#     group_cols_filtered = Symbol[col for col in group_cols if col in df_col_names_symbols]
    
#     if verbose
#         @info "Requested grouping columns: $group_cols"
#         @info "Available columns in dataframe: $(names(df))"
#         @info "Grouping by columns: $group_cols_filtered"
#         if isempty(group_cols_filtered)
#             @warn "WARNING: No grouping columns found! Groupby will create only 1 group with all data combined."
#         end
#     end
    
#     # Group by the specified columns
#     grouped = groupby(df, group_cols_filtered)
    
#     if verbose
#         @info "Groupby created $(length(grouped)) groups"
#     end
    
#     # Warning if we got only 1 group (usually indicates empty group_cols_filtered)
#     if length(grouped) == 1 && !isempty(group_cols_filtered)
#         if verbose
#             @warn "Only 1 group created despite specifying group columns. This may indicate a data issue."
#         end
#     elseif length(grouped) == 1 && isempty(group_cols_filtered)
#         if verbose
#             @warn "CRITICAL: No grouping columns specified! All data will be combined into 1 row."
#         end
#     end
    
#     # Pre-allocate results
#     results = Dict[]
    
#     for group in grouped
#         # Sort by mu within this group to maintain order
#         sort!(group, :mu)
        
#         # Concatenate mu, specloc_gap, and specloc_winding into vectors
#         mus = Vector{Float64}(group.mu)
#         specloc_gaps = Vector{Float64}(group.specloc_gap)
#         specloc_windings = Vector{Float64}(group.specloc_winding)
        
#         # Count contiguous regions where specloc_gap == 0.0
#         gap_count = 0
#         in_gap = false
        
#         for gap_val in specloc_gaps
#             is_gap = gap_val == 0.0
            
#             # If we just entered a gap region, increment counter
#             if is_gap && !in_gap
#                 gap_count += 1
#                 in_gap = true
#             elseif !is_gap && in_gap
#                 in_gap = false
#             end
#         end
        
#         # Create a row for this group with its parameters and concatenated data
#         row_dict = Dict(col => group[1, col] for col in group_cols_filtered)
#         row_dict[:mus] = mus
#         row_dict[:specloc_gaps] = specloc_gaps
#         row_dict[:specloc_windings] = specloc_windings
#         row_dict[:gap_count] = gap_count
        
#         # Debug on first row
#         if verbose && length(results) == 0
#             @info "First row_dict keys before push: $(collect(keys(row_dict)))"
#             @info "row_dict grouping column values:"
#             for col in group_cols_filtered
#                 @info "  :$col => $(row_dict[col])"
#             end
#         end
        
#         push!(results, row_dict)
#     end
    
#     # Convert to dataframe by unpacking the results vector of dicts properly
#     if isempty(results)
#         result_df = DataFrame()
#     else
#         # Debug: inspect first result dict to see what keys we have
#         first_result = results[1]
#         col_names = collect(keys(first_result))
        
#         if verbose
#             @info "First result dict keys: $col_names"
#             @info "First result dict values types: $([(k, typeof(first_result[k])) for k in col_names])"
#             @info "Number of results to convert: $(length(results))"
#         end
        
#         # Build a clean dictionary of column arrays by iterating through each column name
#         data_dict = Dict()
#         for col_name in col_names
#             # Extract this column from all result dicts
#             col_data = [results[i][col_name] for i in 1:length(results)]
#             data_dict[col_name] = col_data
            
#             if verbose
#                 first_val = col_data[1]
#                 if isa(first_val, Vector)
#                     @info "  Column $col_name: Vector{$(eltype(first_val))} ($(length(results)) rows with vectors of lengths [$(length(first_val))])"
#                 else
#                     @info "  Column $col_name: $(typeof(first_val))"
#                 end
#             end
#         end
        
#         # Create DataFrame from the clean dictionary
#         result_df = DataFrame(data_dict)
        
#         # Verify column count
#         if verbose
#             @info "After DataFrame construction:"
#             @info "  Expected $(length(col_names)) columns, got $(length(names(result_df)))"
#             @info "  DataFrame columns: $(names(result_df))"
#         end
        
#         if length(names(result_df)) < length(col_names)
#             missing_cols = setdiff(col_names, names(result_df))
#             @warn "Some columns were lost: $missing_cols"
#         end
#     end
    
#     if verbose
#         @info "Grouped into $(nrow(result_df)) unique parameter combinations"
#         @info "Result DataFrame columns: $(names(result_df))"
#         if !isempty(result_df)
#             @info "First few rows details:"
#             for col in names(result_df)
#                 val = result_df[1, col]
#                 if isa(val, Vector)
#                     @info "  :$col => Vector{$(eltype(val))} of length $(length(val))"
#                 else
#                     @info "  :$col => $val ($(typeof(val)))"
#                 end
#             end
#         end
#     end
    
#     return result_df
# end








"""
    count_contiguous_gaps(gaps::Vector{Float64}; tol::Float64=1e-5)

Count contiguous regions where values fall below `tol`.
"""
function count_contiguous_gaps(gaps::Vector{Float64}; tol::Float64=1e-5)
    isempty(gaps) && return 0
    is_gap = gaps .<= tol
    
    # Region starts at index 1 if true, plus any false -> true transition
    count = is_gap[1] ? 1 : 0
    for i in 2:length(is_gap)
        if is_gap[i] && !is_gap[i-1]
            count += 1
        end
    end
    return count
end

function group_and_concatenate_data(
    df::DataFrame;
    group_cols::Vector{Symbol}=[:N, :t_n, :Delta, :sequence_name, :sequence_id, :kappa, :x0],
    gap_tol::Float64=1e-5,
    smoothing_halfwidth::Int=25,
    min_gap_width::Int=5,
    verbose::Bool=true
)
    df_props = propertynames(df)
    group_cols_filtered = Symbol[col for col in group_cols if col in df_props]
    
    if verbose
        @info "Grouping by columns: $group_cols_filtered"
    end
    
    grouped = groupby(df, group_cols_filtered)
    results = Dict{Symbol, Any}[]
    
    for group in grouped
        # Materialize into DataFrame and sort explicitly by mu
        group_df = sort(DataFrame(group), :mu)
        
        mus = Vector{Float64}(group_df.mu)
        specloc_gaps = Vector{Float64}(group_df.specloc_gap)
        specloc_windings = Vector{Float64}(group_df.specloc_winding)
        
        # NEW: Find critical mu first (where trivial phase becomes permanent)
        critical_mu = find_critical_mu_for_group(mus, specloc_windings; smoothing_halfwidth=smoothing_halfwidth)
        
        # NEW: Use robust counting with smoothing and minimum gap width
        gap_count, _ = count_robust_gaps(
            mus, 
            specloc_windings; 
            critical_mu=critical_mu,
            smoothing_halfwidth=smoothing_halfwidth,
            min_gap_width=min_gap_width
        )
        
        row_dict = Dict{Symbol, Any}(col => group_df[1, col] for col in group_cols_filtered)
        row_dict[:mus] = mus
        row_dict[:specloc_gaps] = specloc_gaps
        row_dict[:specloc_windings] = specloc_windings
        row_dict[:critical_mu] = critical_mu
        row_dict[:gap_count] = gap_count
        
        push!(results, row_dict)
    end
    
    if isempty(results)
        return DataFrame()
    end
    
    # Construct output DataFrame
    col_names = collect(keys(results[1]))
    data_dict = Dict(col => [r[col] for r in results] for col in col_names)
    result_df = DataFrame(data_dict)
    
    if verbose
        @info "Grouped into $(nrow(result_df)) parameter combinations"
        @info "Gap count range observed: $(extrema(result_df.gap_count))"
        @info "Critical mu range observed: $(extrema(result_df.critical_mu))"
    end
    
    return result_df
end









"""
    find_critical_mu_for_group(
        mus::Vector{Float64},
        specloc_windings::Vector{Float64};
        min_sustained_length::Int=100,
        smoothing_halfwidth::Int=25
    )

Find the critical mu above which trivial phase (winding=0) never ends.
Works backwards from high mu to find where trivial region starts.
"""
function find_critical_mu_for_group(
    mus::Vector{Float64},
    specloc_windings::Vector{Float64};
    min_sustained_length::Int=100,
    smoothing_halfwidth::Int=25
)
    if length(mus) < min_sustained_length
        return Inf
    end
    
    # Simple moving average smoothing
    is_trivial_raw = abs.(specloc_windings) .< 0.01  # or use your gap threshold
    is_trivial_smooth = similar(is_trivial_raw)
    
    for i in eachindex(is_trivial_raw)
        lo = max(1, i - smoothing_halfwidth)
        hi = min(length(is_trivial_raw), i + smoothing_halfwidth)
        is_trivial_smooth[i] = mean(is_trivial_raw[lo:hi]) > 0.5
    end
    
    # Scan backwards: find first position where everything after it is trivial
    for i in length(is_trivial_smooth):-1:min_sustained_length
        if all(is_trivial_smooth[i:end])
            return mus[i]
        end
    end
    
    return Inf
end

"""
    count_robust_gaps(
        mus::Vector{Float64},
        specloc_windings::Vector{Float64};
        critical_mu::Float64=Inf,
        smoothing_halfwidth::Int=25,
        min_gap_width::Int=5
    )

Count distinct gap regions below critical_mu with smoothing for robustness.
"""
function count_robust_gaps(
    mus::Vector{Float64},
    specloc_windings::Vector{Float64};
    critical_mu::Float64=Inf,
    smoothing_halfwidth::Int=25,
    min_gap_width::Int=5
)
    # Filter to regions below critical mu
    mask = mus .<= critical_mu
    mus_filtered = mus[mask]
    windings_filtered = specloc_windings[mask]
    
    if isempty(mus_filtered)
        return 0, Inf
    end
    
    # Smooth to eliminate noise
    is_trivial = similar(windings_filtered, Bool)
    for i in eachindex(windings_filtered)
        lo = max(1, i - smoothing_halfwidth)
        hi = min(length(windings_filtered), i + smoothing_halfwidth)
        mean_winding = mean(abs.(windings_filtered[lo:hi]))
        is_trivial[i] = mean_winding < 0.01
    end
    
    # Count distinct gap regions (ignoring very short noise)
    gap_count = 0
    in_gap = false
    gap_width = 0
    
    for is_gap in is_trivial
        if is_gap
            if !in_gap
                in_gap = true
                gap_width = 1
            else
                gap_width += 1
            end
        else
            if in_gap && gap_width >= min_gap_width
                gap_count += 1
            end
            in_gap = false
            gap_width = 0
        end
    end
    
    # Don't forget the last gap
    if in_gap && gap_width >= min_gap_width
        gap_count += 1
    end
    
    return gap_count, critical_mu
end







println("Started grouping and concatenating data...")
@time grouped_df = group_and_concatenate_data(
    df;
    group_cols=[:N, :t_n, :Delta, :sequence_name, :sequence_id, :kappa, :x0],
    smoothing_halfwidth=25,
    min_gap_width=5
)

# Debug: print the structure of grouped_df
println("\n=== Grouped DataFrame Structure ===")
println("Columns: $(names(grouped_df))")
println("Size: $(size(grouped_df))")
println("\nFirst row:")
for col in names(grouped_df)
    val = grouped_df[1, col]
    if isa(val, Vector)
        println("  :$col => Vector of length $(length(val))")
    else
        println("  :$col => $val")
    end
end
println("==================================\n")

# Print statistics on the new robust gap counting
println("\n=== Robust Gap Counting Statistics ===")
println("Gap count range: $(extrema(grouped_df.gap_count))")
println("Mean gap count: $(mean(grouped_df.gap_count))")
println("Critical mu range: $(extrema(grouped_df.critical_mu))")
println("Fraction with Inf critical mu: $(sum(isinf.(grouped_df.critical_mu)) / nrow(grouped_df))")
println("=========================================\n")



function plt_gap_count_vs_rho(
    df::DataFrame,
    filename::String;
    log_rho::Bool=false,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    fixed_values...
)
    """
    Plot gap count vs rho from the grouped dataframe with concatenated vectors.
    Uses robust gap counting that accounts for critical mu and smoothing.
    
    Args:
        df: DataFrame with columns :t_n and :gap_count (from group_and_concatenate_data)
        filename: Path to save the plot
        log_rho: Whether to log-scale the rho axis
        atol/rtol: Tolerance for approximate float matching in fixed_values filter
        verbose: Print information messages
        fixed_values: Keyword arguments to filter rows (e.g., Delta=0.05, N=200)
    """
    
    if verbose
        @info "Available columns in input dataframe: $(names(df))"
    end
    
    # Convert column names to symbols for reliable checking (names() returns Strings)
    df_col_names_symbols = Symbol.(names(df))
    
    # Check for required columns
    if !(:t_n in df_col_names_symbols)
        error("Column :t_n is required but not found in dataframe. Available columns: $(names(df))")
    end
    
    if !(:gap_count in df_col_names_symbols)
        error("Column :gap_count is required but not found in dataframe. Available columns: $(names(df))")
    end
    
    # Approx-aware filtering on row-level context columns
    _kv_match(row, key, val) = begin
        rv = row[key]  
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    # Filter rows by fixed_values if provided
    filtered_df = isempty(fixed_values) ? df :
                  filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_df)
        if verbose
            @info "No rows match fixed_values: $(Dict(fixed_values))"
        end
        error("No data matches the specified fixed values.")
    end
    
    if verbose
        @info "After filtering: $(nrow(filtered_df)) rows match criteria"
    end

    # Calculate rho from t_n (t_n should be a tuple/pair)
    rho_vals = map(filtered_df.t_n) do t_n_val
        if isa(t_n_val, Tuple) || isa(t_n_val, Vector)
            # t_n is stored as (t1, t2) or similar
            t_n_val[2] / t_n_val[1]
        else
            error("Expected :t_n to be a tuple or vector, got $(typeof(t_n_val))")
        end
    end
    
    gap_counts = filtered_df.gap_count
    
    # Create dataframe for sorting
    plot_df = DataFrame(rho=rho_vals, gap_count=gap_counts)
    sort!(plot_df, :rho)
    
    if verbose
        @info "Prepared $(nrow(plot_df)) points for plotting"
        @info "Rho range: $(extrema(plot_df.rho))"
        @info "Gap count range: $(extrema(plot_df.gap_count))"
    end
    
    # Apply log scale if requested
    x_vals = log_rho ? log10.(plot_df.rho) : plot_df.rho
    y_vals = plot_df.gap_count

    # Create line plot with improved styling
    plt = plot(
        x_vals,
        y_vals,
        xlabel=log_rho ? "log₁₀(ρ)" : "ρ",
        ylabel="Gap Count",
        title="Gap Count vs ρ (Robust Detection with Critical μ)",
        marker=:circle,
        markersize=4,
        linewidth=2,
        legend=false,
        size=(900, 600),
        ylims=(0, maximum(y_vals) + 1),
        xlims=(minimum(x_vals) - 0.05, maximum(x_vals) + 0.05)
    )

    savefig(plt, filename)
    if verbose
        @info "Plot saved to $filename"
    end
end

plt_gap_count_vs_rho(
    grouped_df,
    joinpath(output_path, "gap_count_vs_rho.png");
    log_rho=false
)

# Extract and save gap count ranges
function extract_gap_count_ranges(
    df::DataFrame;
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    fixed_values...
)
    """
    Extract rho ranges for each gap count region.
    
    For each contiguous region where gap_count is constant, returns the
    rho_min, rho_max, and width of that region.
    
    Args:
        df: Grouped DataFrame with :t_n and :gap_count columns
        atol/rtol: Tolerance for float matching in filtering
        verbose: Print information
        fixed_values: Keyword args to filter rows (e.g., Delta=0.05, N=200)
    
    Returns:
        DataFrame with columns: gap_count, rho_min, rho_max, rho_width
    """
    
    # Approx-aware filtering
    _kv_match(row, key, val) = begin
        rv = row[key]  
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    # Filter rows by fixed_values if provided
    filtered_df = isempty(fixed_values) ? df :
                  filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_df)
        if verbose
            @info "No rows match fixed_values: $(Dict(fixed_values))"
        end
        return DataFrame()
    end

    # Calculate rho from t_n
    rho_vals = map(filtered_df.t_n) do t_n_val
        if isa(t_n_val, Tuple) || isa(t_n_val, Vector)
            t_n_val[2] / t_n_val[1]
        else
            error("Expected :t_n to be a tuple or vector, got $(typeof(t_n_val))")
        end
    end
    
    gap_counts = filtered_df.gap_count
    
    # Create sorted dataframe
    data_df = DataFrame(rho=rho_vals, gap_count=gap_counts)
    sort!(data_df, :rho)
    
    if verbose
        @info "Extracted $(nrow(data_df)) rho points for gap range analysis"
    end
    
    # Find boundaries where gap_count changes
    ranges = Dict{Int, Tuple{Float64, Float64}}()
    
    if nrow(data_df) == 0
        return DataFrame()
    end
    
    current_gap_count = data_df[1, :gap_count]
    rho_start = data_df[1, :rho]
    
    for i in 2:nrow(data_df)
        new_gap_count = data_df[i, :gap_count]
        
        if new_gap_count != current_gap_count
            # Save the previous region
            rho_end = data_df[i-1, :rho]
            if !haskey(ranges, current_gap_count)
                ranges[current_gap_count] = (rho_start, rho_end)
            else
                # If the same gap_count appears again, extend the range (shouldn't happen with clean data)
                old_start, old_end = ranges[current_gap_count]
                ranges[current_gap_count] = (min(old_start, rho_start), max(old_end, rho_end))
            end
            
            # Start new region
            current_gap_count = new_gap_count
            rho_start = data_df[i, :rho]
        end
    end
    
    # Don't forget the last region
    rho_end = data_df[end, :rho]
    if !haskey(ranges, current_gap_count)
        ranges[current_gap_count] = (rho_start, rho_end)
    else
        old_start, old_end = ranges[current_gap_count]
        ranges[current_gap_count] = (min(old_start, rho_start), max(old_end, rho_end))
    end
    
    # Convert to sorted dataframe
    gap_counts_list = sort(collect(keys(ranges)))
    rho_mins = [ranges[gc][1] for gc in gap_counts_list]
    rho_maxs = [ranges[gc][2] for gc in gap_counts_list]
    rho_widths = rho_maxs .- rho_mins
    
    result_df = DataFrame(
        gap_count=gap_counts_list,
        rho_min=rho_mins,
        rho_max=rho_maxs,
        rho_width=rho_widths
    )
    
    if verbose
        @info "Found $(nrow(result_df)) distinct gap count regions:"
        for row in eachrow(result_df)
            @info "  Gap count $(row.gap_count): ρ ∈ [$(round(row.rho_min; digits=6)), $(round(row.rho_max; digits=6))], width = $(round(row.rho_width; digits=6))"
        end
    end
    
    return result_df
end

# Extract gap ranges and save to CSV
gap_ranges_df = extract_gap_count_ranges(grouped_df; verbose=true)

if !isempty(gap_ranges_df)
    gap_ranges_csv_path = joinpath(output_path, "gap_count_ranges.csv")
    CSV.write(gap_ranges_csv_path, gap_ranges_df)
    println("Gap count ranges saved to: $gap_ranges_csv_path")
else
    @warn "No gap count ranges could be extracted"
end

# Apply gap labeling theorem: N = p + φq where φ = 1/golden_ratio
function apply_gap_labeling_theorem(
    gap_ranges_df::DataFrame;
    phi::Float64=(sqrt(5) - 1) / 2,  # 1/golden_ratio ≈ 0.618...
    p_range::Vector{Int}=collect(-5:5),
    q_max::Int=5,
    verbose::Bool=true
)
    """
    Apply quasicrystal gap labeling theorem to gap count ranges.
    
    For each gap region with gap_count value N, find (p, q) integers such that
    N ≈ p + φq minimizes the error, where φ is the quasicrystal generating number.
    
    Args:
        gap_ranges_df: DataFrame with columns gap_count, rho_min, rho_max, rho_width
        phi: Quasicrystal generating number (default: 1/golden_ratio)
        p_range: Range of p values to search (default: -5 to 5)
        q_max: Maximum absolute q value (default: 5)
        verbose: Print information
    
    Returns:
        DataFrame with additional columns: p, q, label_error
    """
    
    if verbose
        @info "Applying gap labeling theorem with φ = $(round(phi; digits=6))"
        @info "Searching p ∈ [$( p_range[1]), $(p_range[end])], |q| ≤ $q_max"
    end
    
    result_df = copy(gap_ranges_df)
    p_labels = Int[]
    q_labels = Int[]
    label_errors = Float64[]
    
    for row in eachrow(result_df)
        gap_count = row.gap_count
        
        best_err = Inf
        best_p = 0
        best_q = 0
        
        # Search over (p, q) pairs for best match
        for p_try in p_range
            # For given p, calculate which q gives closest match
            q_candidate = round(Int, (gap_count - p_try) / phi)
            
            # Enforce q_max constraint
            if abs(q_candidate) > q_max
                continue
            end
            
            # Compute label value
            label_val = p_try + q_candidate * phi
            err = abs(gap_count - label_val)
            
            if err < best_err
                best_err = err
                best_p = p_try
                best_q = q_candidate
            end
        end
        
        push!(p_labels, best_p)
        push!(q_labels, best_q)
        push!(label_errors, best_err)
        
        if verbose && nrow(result_df) <= 20
            @info "Gap count $(gap_count): (p, q) = ($(best_p), $(best_q)), " *
                  "p + φq = $(round(best_p + best_q * phi; digits=4)), error = $(round(best_err; digits=6))"
        end
    end
    
    result_df[!, :p] = p_labels
    result_df[!, :q] = q_labels
    result_df[!, :label_error] = label_errors
    
    # Add computed label value for reference
    result_df[!, :label_value] = result_df.p .+ result_df.q .* phi
    
    if verbose
        @info "Gap labeling complete. Summary:"
        @info "  Mean labeling error: $(round(mean(label_errors); digits=6))"
        @info "  Max labeling error: $(round(maximum(label_errors); digits=6))"
        @info "  Unique (p, q) labels: $(length(unique(zip(p_labels, q_labels))))"
    end
    
    return result_df
end

# Apply gap labeling
labeled_gaps_df = apply_gap_labeling_theorem(gap_ranges_df; verbose=true)

# Save labeled gaps to CSV
if !isempty(labeled_gaps_df)
    labeled_gaps_csv_path = joinpath(output_path, "gap_count_ranges_labeled.csv")
    CSV.write(labeled_gaps_csv_path, labeled_gaps_df)
    println("Labeled gap ranges saved to: $labeled_gaps_csv_path")
    println("  Columns in CSV: $(join(names(labeled_gaps_df))), \", \"))")
    
    # Print summary table
    println("\n=== Gap Labeling Summary ===")
    println("gap_count | (p, q) | ρ_min    | ρ_max    | ρ_width  | p+φq     | error")
    println("-" ^ 75)
    for row in eachrow(labeled_gaps_df)
        @printf("%9d | (%2d,%2d) | %8.4f | %8.4f | %8.4f | %8.4f | %.2e\n",
            row.gap_count, row.p, row.q, row.rho_min, row.rho_max, row.rho_width, row.label_value, row.label_error)
    end
    println("=" ^ 75)
else
    @warn "No labeled gaps could be created"
end

# Optional: Visualize smoothing effect on a single rho slice for validation
function visualize_gap_labels_on_plot(
    grouped_df::DataFrame,
    labeled_gaps_df::DataFrame,
    output_path::String;
    log_rho::Bool=false,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    fixed_values...
)
    """
    Create an enhanced plot showing gap count vs rho with (p, q) labels annotated.
    """
    
    # Approx-aware filtering
    _kv_match(row, key, val) = begin
        rv = row[key]  
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    # Filter grouped_df by fixed_values if provided
    filtered_df = isempty(fixed_values) ? grouped_df :
                  filter(row -> _row_matches(row, fixed_values), grouped_df)

    if isempty(filtered_df)
        if verbose
            @info "No rows match fixed_values for gap label plot"
        end
        return
    end

    # Calculate rho from t_n
    rho_vals = map(filtered_df.t_n) do t_n_val
        if isa(t_n_val, Tuple) || isa(t_n_val, Vector)
            t_n_val[2] / t_n_val[1]
        else
            error("Expected :t_n to be a tuple or vector")
        end
    end
    
    gap_counts = filtered_df.gap_count
    
    # Create sorted dataframe
    plot_df = DataFrame(rho=rho_vals, gap_count=gap_counts)
    sort!(plot_df, :rho)
    
    # Create main plot
    x_vals = log_rho ? log10.(plot_df.rho) : plot_df.rho
    y_vals = plot_df.gap_count

    plt = plot(
        x_vals,
        y_vals,
        xlabel=log_rho ? "log₁₀(ρ)" : "ρ",
        ylabel="Gap Count",
        title="Gap Count vs ρ with (p,q) Labels",
        marker=:circle,
        markersize=4,
        linewidth=2,
        legend=false,
        size=(1000, 600),
        ylims=(0, maximum(y_vals) + 1),
        xlims=(minimum(x_vals) - 0.05, maximum(x_vals) + 0.05)
    )

    # Add (p, q) label annotations for each gap region
    for row in eachrow(labeled_gaps_df)
        rho_center = (row.rho_min + row.rho_max) / 2
        gap_count_val = row.gap_count
        p_val = row.p
        q_val = row.q
        
        x_pos = log_rho ? log10(rho_center) : rho_center
        y_pos = gap_count_val + 0.3
        
        # Annotate with (p, q) label
        annotate!(plt, x_pos, y_pos, 
            text("($(p_val),$(q_val))", 10, :center, color=:red))
    end

    savefig(plt, joinpath(output_path, "gap_count_vs_rho_labeled.png"))
    if verbose
        @info "Labeled gap count plot saved to: $(joinpath(output_path, "gap_count_vs_rho_labeled.png"))"
    end
end

# Create labeled visualization
try
    visualize_gap_labels_on_plot(grouped_df, labeled_gaps_df, output_path; verbose=true)
catch e
    @warn "Could not generate labeled gap plot: $e"
end

# Optional: Visualize smoothing effect on a single rho slice for validation
function visualize_smoothing_example(
    df::DataFrame,
    output_path::String;
    row_idx::Int=1,
    smoothing_halfwidth::Int=25,
    verbose::Bool=true
)
    """
    Plot raw winding number vs smoothed trivial detection for one parameter set.
    Useful for validating smoothing parameters.
    """
    if row_idx > nrow(df)
        @warn "Row index $row_idx exceeds dataframe size $(nrow(df))"
        return
    end
    
    row = df[row_idx, :]
    mus = row.mus
    windings = row.specloc_windings
    critical_mu = row.critical_mu
    
    # Compute smoothed is_trivial for this row
    is_trivial_raw = abs.(windings) .< 0.01
    is_trivial_smooth = similar(is_trivial_raw)
    
    for i in eachindex(is_trivial_raw)
        lo = max(1, i - smoothing_halfwidth)
        hi = min(length(is_trivial_raw), i + smoothing_halfwidth)
        is_trivial_smooth[i] = mean(is_trivial_raw[lo:hi]) > 0.5
    end
    
    # Create figure with subplots
    p1 = plot(mus, abs.(windings), 
        xlabel="μ", ylabel="|Winding Number|",
        title="Raw Spectral Localizer Winding Number",
        marker=:dot, markersize=2, legend=false, linewidth=1)
    hline!(p1, [0.01], line=:dash, color=:red, label="Threshold", legend=true)
    
    # Mark critical mu
    if !isinf(critical_mu)
        vline!(p1, [critical_mu], line=:dash, color=:orange, label="Critical μ", legend=true)
    end
    
    p2 = plot(mus, Float64.(is_trivial_raw),
        xlabel="μ", ylabel="Is Trivial (Raw)",
        title="Raw Trivial Detection (threshold < 0.01)",
        marker=:dot, markersize=2, legend=false, linewidth=1,
        ylim=(-0.1, 1.1))
    
    p3 = plot(mus, Float64.(is_trivial_smooth),
        xlabel="μ", ylabel="Is Trivial (Smoothed)",
        title="Smoothed Trivial Detection (moving avg)",
        marker=:dot, markersize=2, legend=false, linewidth=1,
        ylim=(-0.1, 1.1))
    
    # Mark critical mu
    if !isinf(critical_mu)
        vline!(p2, [critical_mu], line=:dash, color=:orange, alpha=0.5)
        vline!(p3, [critical_mu], line=:dash, color=:orange, alpha=0.5)
    end
    
    plt_combined = plot(p1, p2, p3, layout=(3,1), size=(1000, 800), dpi=100)
    
    filename = joinpath(output_path, "smoothing_validation_example.png")
    savefig(plt_combined, filename)
    
    if verbose
        @info "Smoothing validation plot saved to $filename"
        @info "Row $row_idx: N=$(row.N), t_n=$(row.t_n), Δ=$(row.Delta)"
        @info "  Critical μ: $(critical_mu)"
        @info "  Gap count: $(row.gap_count)"
    end
end

# Generate smoothing validation plot
try
    visualize_smoothing_example(grouped_df, output_path; 
        row_idx=div(nrow(grouped_df), 2),  # Use middle row as example
        smoothing_halfwidth=25,
        verbose=true
    )
catch e
    @warn "Could not generate smoothing validation plot: $e"
end
