using DataFrames
using BSON
using Glob
using Plots
using ProgressMeter
using Statistics

data_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/N_range_data/hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1-0.0_nbins1000_npb1_N(100-1000-10)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"

output_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/spec_loc_testing/N_range_data/hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1-0.0_nbins1000_npb1_N(100-1000-10)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"
isdir(output_path) || mkpath(output_path)

# function process_bson_files_fluid(folder_path::String)
#     file_paths = glob("*.bson", folder_path)
    
#     # Pre-allocate a vector to hold all the dataframes
#     dfs = DataFrame[]
#     failed_files = String[]

#     @showprogress for file_path in file_paths
#         bson_data = nothing
#         try
#             bson_data = load(file_path)
#         catch err
#             push!(failed_files, file_path)
#             @warn "Failed to load bson file (skipping)" file=file_path error=err
#             continue
#         end

#         if haskey(bson_data, "results_df")
#             push!(dfs, DataFrame(bson_data["results_df"]))
#         else
#             @warn "bson file missing results_df key (skipping)" file=file_path
#         end
#     end

#     println("Number of files read: $(length(dfs)) / $(length(file_paths))")
#     if !isempty(failed_files)
#         @warn "Some bson files could not be processed" failed_count=length(failed_files) failed_files=failed_files
#     end
    
#     # Combine everything dynamically
#     combined_dataframe = isempty(dfs) ? DataFrame() : vcat(dfs..., cols=:union)

#     # Print the resulting column names dynamically
#     println("\n=== DataFrame Loading Complete ===")
#     println("DataFrame column names:")
#     for col in names(combined_dataframe)
#         println("  - :$col")
#     end
#     println("==================================\n")

#     return combined_dataframe
# end

function process_bson_files(folder_path::String)
    combined_dataframe = DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        sequence_id = Any[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Vector{BigFloat}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Matrix{BigFloat}, Missing}[]
    )
    
    file_paths = glob("*.bson", folder_path)
    files_read = 0
    failed_files = String[]

    @showprogress for file_path in file_paths
        bson_data = nothing
        try
            bson_data = BSON.load(file_path)
        catch err
            push!(failed_files, file_path)
            @warn "Failed to load BSON file (skipping)" file=file_path error=err
            continue
        end

        if haskey(bson_data, :results_df)
            results_df = DataFrame(bson_data[:results_df])
            append!(combined_dataframe, results_df)
            files_read += 1
        else
            @warn "BSON file missing :results_df key (skipping)" file=file_path
        end
    end

    println("number of files read: $files_read / $(length(file_paths))")
    if !isempty(failed_files)
        @warn "Some BSON files could not be processed" failed_count=length(failed_files) failed_files=failed_files
    end
    
    return combined_dataframe
end

df = process_bson_files(data_path)

# println("DataFrame column names & content indication:")
# for col in names(df)
#     println("  - :$col => $(eltype(df[!, col])) | First 3 entries: $(first(df[!, col], 3))")
# end

function plt_mp_vs_mu(
    df::DataFrame,
    filename::String;
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
        error("No rows match the specified fixed values.")
    end

    # Sort by mu for plotting
    sort!(filtered_rows, :mu)
    
    # Create the plot
    p = plot(
        filtered_rows.mu,
        filtered_rows.mp,
        xlabel="μ",
        ylabel="MP",
        title="Majorana Polarization vs Chemical Potential",
        legend=false,
        linewidth=2,
        # marker=:circle,
        grid=false
    )
    
    # Save the plot
    savefig(p, filename)
    
    println("Plot saved to: $filename")
    
    return p
end

# Ns = unique(df.N)
# for N in Ns
#     plt_mp_vs_mu(df, joinpath(output_path, "mp_vs_mu_$N.png"), N=N)
# end

function plt_mp_vs_mu_multiline_L(
    df::DataFrame,
    filename::String;
    L_col::Symbol=:N,
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
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k, v) in kvs)

    filtered_rows = isempty(fixed_values) ? df :
                    filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_rows)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
        end
        error("No rows match the specified fixed values.")
    end

    L_values = unique(filtered_rows[!, L_col])
    sort!(L_values)

    if isempty(L_values)
        error("No values found for $(L_col) after filtering.")
    end

    nL = length(L_values)
    tvals = nL == 1 ? [0.5] : collect(range(0.0, 1.0; length=nL))
    line_colors = [get(cgrad(:Blues), t) for t in tvals]

    p = plot(
        xlabel="μ",
        ylabel="MP",
        title="Majorana Polarization vs Chemical Potential (by $(String(L_col)))",
        legend=:outerright,
        linewidth=2,
        grid=false
    )

    for (i, L_val) in enumerate(L_values)
        L_rows = filter(row -> _kv_match(row, L_col, L_val), filtered_rows)
        isempty(L_rows) && continue
        sort!(L_rows, :mu)

        plot!(
            p,
            L_rows.mu,
            L_rows.mp,
            label="$(String(L_col))=$(L_val)",
            color=line_colors[i]
        )
    end

    savefig(p, filename)

    if verbose
        println("Plot saved to: $filename")
    end

    return p
end

plt_mp_vs_mu_multiline_L(df, joinpath(output_path, "mp_vs_mu_all_L.png"))

