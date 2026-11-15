using DataFrames
using JLD2
using Glob
using Plots
using ProgressMeter
using Statistics


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


function plt_specloc_mu_vs_kappa(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap, # choose :gap or :winding
    log_kappa::Bool=false,
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
        error("No gaps match the specified fixed values.")
    end

    # Extract unique slopes (first element of sequence_id) and mu values, sort them
    kappas = unique(filtered_rows.kappa)
    sort!(kappas)
    mu_values = unique(filtered_rows.mu)
    sort!(mu_values)

    # Choose the column to plot
    z_col = type == :gap ? :specloc_gap : :specloc_winding
    
    # Create matrix for heatmap indexed by (slope, mu)
    z_matrix = fill(NaN, length(kappas), length(mu_values))

    for (i, kappa) in enumerate(kappas)
        for (j, mu) in enumerate(mu_values)
            # Find row matching this (kappa, mu) pair
            matching_rows = filter(row -> 
                isapprox(row.kappa, kappa; atol=atol, rtol=rtol) && 
                isapprox(row.mu, mu; atol=atol, rtol=rtol), 
                filtered_rows)
            if !isempty(matching_rows)
                val = matching_rows[1, z_col]
                z_matrix[i, j] = type == :gap ? log(val) : val
            end
        end
    end

    # Determine title and colorbar label based on type
    title_str = type == :gap ? "Spectral Localiser Gap vs Mu" : "Spectral Localiser Winding Number vs Mu"
    colorbar_title_str = type == :gap ? "log(Gap)" : "Winding Number"

    if log_kappa
        kappas = log10.(kappas)
    end

    plt = heatmap(
        mu_values,
        kappas,
        z_matrix,
        xlabel="mu",
        ylabel=log_kappa ? "log(Kappa)" : "Kappa",
        title=title_str,
        colorbar_title=colorbar_title_str
    )

    # display(plt)
    savefig(plt, filename)

end

function plt_specloc_mu_vs_rho(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap,
    log_rho::Bool=false,
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    mu_norm::Bool=true,
    fixed_values...
)
    # Approx-aware filtering on row-level context columns
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
        error("No data matches the specified fixed values.")
    end

    # Extract and sort unique values for matrix dimensions
    mu_values = sort(unique(filtered_rows.mu))
    rho_values = sort(unique(filtered_rows.rho))

    # Build O(1) index lookup maps
    rho_to_i = Dict(val => i for (i, val) in enumerate(rho_values))
    mu_to_j  = Dict(val => j for (j, val) in enumerate(mu_values))

    z_col = type == :gap ? :specloc_gap : :specloc_winding
    z_matrix = fill(NaN, length(rho_values), length(mu_values))

    # Single-pass population (O(N_rows) complexity)
    for row in eachrow(filtered_rows)
        i = rho_to_i[row.rho]
        j = mu_to_j[row.mu]
        val = getproperty(row, z_col)
        # Prevent DomainError on log(0.0)
        z_matrix[i, j] = (type == :gap) ? (val > 0.0 ? log(val) : NaN) : val
    end

    title_str = type == :gap ? "Spectral Localiser Gap vs Mu" : "Spectral Localiser Winding Number vs Mu"
    colorbar_title_str = type == :gap ? "log(Gap)" : "Winding Number"

    y_axis_vals = log_rho ? log10.(rho_values) : rho_values

    plt = heatmap(
        mu_values,
        y_axis_vals,
        z_matrix,
        xlabel="mu",
        ylabel=log_rho ? "log₁₀(Rho)" : "Rho",
        title=title_str,
        colorbar_title=colorbar_title_str
    )

    savefig(plt, filename)
    if verbose
        @info "Plot saved to $filename"
    end
end



plt_types = [:gap, :winding]
mu_norm = true
for type in plt_types
    plt_specloc_mu_vs_rho(
        df,
        joinpath(output_path, "specloc_$(type)_rho_vs_mu$(mu_norm ? "_norm" : "").png");
        type=type,
        log_rho=false,
        atol=1e-8,
        rtol=1e-5,
        verbose=true,
        mu_norm=mu_norm
    )
end

