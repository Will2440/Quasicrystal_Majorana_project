using DataFrames
using JLD2
using Glob
using Plots
using ProgressMeter
using Statistics


data_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/test_new_specloc/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p236004390779364-4238)_Delta(0p05-0p05-1)_x0(218p75-281p25-5)_kappa(1p0000000000000004e-6-0p10000000000000002-10)"

output_path = replace(data_path, "raw_data/" => "results/")
isdir(output_path) || mkpath(output_path)

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


function plt_specloc_vs_mu(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap, # choose :gap or :invariant
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

    # Choose the column to plot
    y_col = type == :gap ? :specloc_gap : :specloc_invariant
    
    plt = plot(
        filtered_rows.mu,
        filtered_rows[!, y_col],
        xlabel="mu",
        ylabel=string(y_col),
        title="Spectral Localiser $(y_col) vs mu"
    )

    # display(plt)
    savefig(plt, filename)

end

# function plt_specloc_mu_QCSC_ratio(
#     df::DataFrame,
#     filename::String;
#     atol::Float64=1e-8,
#     rtol::Float64=1e-5,
#     verbose::Bool=true,
#     fixed_values...
# )

#     # approx-aware filtering on row-level context columns
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
#         error("No gaps match the specified fixed values.")
#     end

#     plt = plot(
#         filtered_rows.mu,
#         filtered_rows.specloc_gap ./ filtered_rows.qcsc_gap,
#         xlabel="mu",
#         ylabel="specloc_gap / qcsc_gap",
#         title="Spectral Localiser Gap / QCSC Gap vs mu"
#     )

#     # display(plt)
#     savefig(plt, filename)

# end

function plt_specloc_butterfly(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap, # choose :gap or :invariant
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
    slopes = unique(getindex.(filtered_rows.sequence_id, 1))
    sort!(slopes)
    mu_values = unique(filtered_rows.mu)
    sort!(mu_values)
    
    # Choose the column to plot
    z_col = type == :gap ? :specloc_gap : :specloc_invariant
    
    # Create a dictionary for fast lookup: (slope, mu) -> z_value
    # Extract slope from sequence_id and pair with mu for each row
    group_dict = Dict()
    for row in eachrow(filtered_rows)
        slope = getindex(row.sequence_id, 1)
        mu = row.mu
        key = (slope, mu)
        if !haskey(group_dict, key)
            group_dict[key] = row[z_col]
        end
    end
    
    # Create matrix for heatmap indexed by (slope, mu)
    z_matrix = fill(NaN, length(slopes), length(mu_values))
    
    for (i, slope) in enumerate(slopes)
        for (j, mu) in enumerate(mu_values)
            key = (slope, mu)
            if haskey(group_dict, key)
                val = group_dict[key]
                z_matrix[i, j] = type == :gap ? log(val) : val
            end
        end
    end
    
    # Determine title and colorbar label based on type
    title_str = type == :gap ? "Spectral Localiser Gap Butterfly" : "Spectral Localiser Invariant Butterfly"
    colorbar_title_str = type == :gap ? "log(Gap)" : "Invariant"
    
    plt = heatmap(
        mu_values,
        slopes,
        z_matrix,
        xlabel="mu",
        ylabel="Slope",
        title=title_str,
        colorbar_title=colorbar_title_str
    )

    # display(plt)
    savefig(plt, filename)

end

function plt_specloc_mu_vs_N(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap, # choose :gap or :invariant
    atol::Float64=1e-8,
    rtol::Float64=1e-5,
    verbose::Bool=true,
    log_N::Bool=true,
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

    # Extract unique Ns and mu values, sort them
    Ns = unique(filtered_rows.N)
    sort!(Ns)
    mu_values = unique(filtered_rows.mu)
    sort!(mu_values)
    
    # Choose the column to plot
    z_col = type == :gap ? :specloc_gap : :specloc_invariant
    
    # Create a dictionary for fast lookup: (N, mu) -> z_value
    grouped = combine(groupby(filtered_rows, [:N, :mu]), z_col => first => :z_value)
    group_dict = Dict((row.N, row.mu) => row.z_value for row in eachrow(grouped))

    # Create matrix for heatmap indexed by (N, mu)
    z_matrix = fill(NaN, length(Ns), length(mu_values))

    for (i, N) in enumerate(Ns)
        for (j, mu) in enumerate(mu_values)
            key = (N, mu)
            if haskey(group_dict, key)
                val = group_dict[key]
                z_matrix[i, j] = type == :gap ? log(val) : val
            end
        end
    end

    if log_N
        Ns = log.(Ns)
    end
    
    # Determine title and colorbar label based on type
    title_str = type == :gap ? "Spectral Localiser Gap Butterfly" : "Spectral Localiser Invariant Butterfly"
    colorbar_title_str = type == :gap ? "log(Gap)" : "Invariant"
    
    plt = heatmap(
        mu_values,
        Ns,
        z_matrix,
        xlabel="mu",
        ylabel="N",
        title=title_str,
        colorbar_title=colorbar_title_str
    )

    # display(plt)
    savefig(plt, filename)

end

function plt_specloc_mu_vs_kappa(
    df::DataFrame,
    filename::String;
    type::Symbol=:gap, # choose :gap or :invariant
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

    # Extract unique kappas and mu values, sort them
    kappas = unique(filtered_rows.kappa)
    sort!(kappas)
    mu_values = unique(filtered_rows.mu)
    sort!(mu_values)

    # Choose the column to plot
    z_col = type == :gap ? :specloc_gap : :specloc_invariant
    
    # Create a dictionary for fast lookup: (kappa, mu) -> z_value
    grouped = combine(groupby(filtered_rows, [:kappa, :mu]), z_col => first => :z_value)
    group_dict = Dict((row.kappa, row.mu) => row.z_value for row in eachrow(grouped))

    # Create matrix for heatmap indexed by (kappa, mu)
    z_matrix = fill(NaN, length(kappas), length(mu_values))

    for (i, kappa) in enumerate(kappas)
        for (j, mu) in enumerate(mu_values)
            key = (kappa, mu)
            if haskey(group_dict, key)
                val = group_dict[key]
                z_matrix[i, j] = type == :gap ? log(val) : val
            end
        end
    end

    # Determine title and colorbar label based on type
    title_str = type == :gap ? "Spectral Localiser Gap vs Mu" : "Spectral Localiser Invariant vs Mu"
    colorbar_title_str = type == :gap ? "log(Gap)" : "Invariant"

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




kappas_in_df = unique(df.kappa)
# println("Unique kappa values in DataFrame: ", kappas_in_df)

plt_types = [:gap, :invariant]
kappas = kappas_in_df # [1e-2, 1e-1, 1e0]
Ns = unique(df.N)
x0s = unique(df.x0)
t_ns = unique(df.t_n)
t2s = [x[2] for x in t_ns]
rhos = t2s

function add_rho_to_df!(df::DataFrame)
    df.rho = map(x -> x[2] / x[1], df.t_n)
end

add_rho_to_df!(df)
println("type of df.rho: ", typeof(df.rho))
println("type of df: ", typeof(df))

# # show dataframe 
# show(df)
# # print keyanmes
# println("column names in df: ", names(df))


# # plt specloc_gap and specloc_invariant vs mu
# for plt_type in plt_types
#     plt_specloc_vs_mu(
#         df,
#         joinpath(output_path, "specloc_$(plt_type)_vs_mu.png"),
#         type=plt_type
#     )
# end


# for plt_type in plt_types
#     for kappa in kappas
#         for N in Ns
#             plt_specloc_butterfly(
#                 df, 
#                 joinpath(output_path, "specloc_$(plt_type)_butterfly_kappa$(kappa)_N$(N).png"), 
#                 type=plt_type,
#                 kappa=kappa,
#                 N=N
#             )
#         end
#     end
# end

# for plt_type in plt_types
#     for kappa in kappas
#         plt_specloc_mu_vs_N(
#             df, 
#             joinpath(output_path, "specloc_$(plt_type)_vs_N_$(minimum(Ns))-$(maximum(Ns))_kappa$(kappa).png"), 
#             type=plt_type,
#             kappa=kappa,
#         )
#     end
# end


# for plt_type in plt_types
#     for x0 in x0s
#         plt_specloc_mu_vs_kappa(
#             df, 
#             joinpath(output_path, "specloc_$(plt_type)_mu_vs_kappa_N$(Ns[1])_x0$(x0).png"), 
#             type=plt_type,
#             x0=x0
#         )
#     end
# end

function avg_over_x0s_for_kappa_plot(df::DataFrame)
    group_cols = [:mu, :kappa, :rho]
    grouped = groupby(df, group_cols)
    avg_df = combine(grouped, :specloc_gap => mean => :specloc_gap, :specloc_invariant => mean => :specloc_invariant)
    return avg_df
end

averaged_df = avg_over_x0s_for_kappa_plot(df)
for plt_type in plt_types
    for rho in rhos
        plt_specloc_mu_vs_kappa(
            averaged_df, 
            joinpath(output_path, "specloc_$(plt_type)_mu_vs_kappa_N$(Ns[1])_rho$(rho)_avg_over_x0s.png"), 
            log_kappa=true,
            type=plt_type,
            rho=rho
        )
    end
end

plt_type = :gap
kappa=1e-2

plt_specloc_butterfly(
    df, 
    joinpath(output_path, "specloc_$(plt_type)_butterfly_kappa$(kappa).png"), 
    type=plt_type,
    kappa=kappa
)


plt_type = :invariant
kappa=1e-2

plt_specloc_butterfly(
    df, 
    joinpath(output_path, "specloc_$(plt_type)_butterfly_kappa$(kappa).png"), 
    type=plt_type,
    kappa=kappa
)