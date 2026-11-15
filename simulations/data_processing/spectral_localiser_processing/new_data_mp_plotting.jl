using DataFrames
using JLD2
using Glob
using Plots
using ProgressMeter
using Statistics
using BSON


data_path_new = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N(5-10000-8)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-601)_Delta(0p05-0p05-1)"
data_path_old = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/N_range_data/hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1-0.0_nbins1000_npb1_N(100-1000-10)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"

output_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N(5-10000-8)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-601)_Delta(0p05-0p05-1)"
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


df_new = process_jld2_files_fluid(data_path_new)
df_old = process_bson_files(data_path_old)


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
# restricted_Ns = filter(N -> N <= 1000, Ns)  # Restrict to N <= 1000 for plotting
# # for N in restricted_Ns
# #     plt_mp_vs_mu(df, joinpath(output_path, "mp_vs_mu_$N.png"), N=N)
# # end


function plt_mp_vs_mu_multiline_L(
    df::DataFrame,
    filename::String;
    L_col::Symbol=:N,
    restrict_Ns::Bool=false,
    ignore_Ns::Vector{Int}=Int[],
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
    if restrict_Ns
        L_values = filter(L -> !(L in ignore_Ns), L_values)  # Exclude ignored Ns for plotting
    end
    sort!(L_values)

    if isempty(L_values)
        error("No values found for $(L_col) after filtering.")
    end

    nL = length(L_values)
    tvals = nL == 1 ? [0.5] : collect(range(0.0, 1.0; length=nL))
    line_colors = [get(cgrad(:Blues), t) for t in tvals]

    p = plot(
        xlabel="", #"μ",
        ylabel="", #"MP",
        title="",#"Majorana Polarization vs Chemical Potential (by $(String(L_col)))",
        legend=false, #:outerright,
        linewidth=2,
        grid=false,
        ylims=(-1.01, 0.01),
        xlims=(0.0, 3.0), #(2.35,2.41)
        axis=false,
        size=(1000, 200)
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
            color=line_colors[i],
            linewidth=2
        )
    end

    savefig(p, filename)

    if verbose
        println("Plot saved to: $filename")
    end

    return p
end

function mp_size_trend(df::DataFrame; min_N::Int=100, max_N::Int=1000, min_points::Int=3)
    size_rows = filter(row -> min_N <= row.N <= max_N, df)
    if isempty(size_rows)
        error("No rows found with N between $min_N and $max_N.")
    end

    size_rows = combine(groupby(size_rows, [:mu, :N]), :mp => mean => :mp)
    sort!(size_rows, [:mu, :N])
    trends = DataFrame(
        mu=Float64[],
        robust_slope=Float64[],
        kendall_tau=Float64[],
        direction_agreement=Float64[],
        n_sizes=Int[],
        min_N=Int[],
        max_N=Int[],
        phase=Union{Missing, Int}[]
    )

    for mu_rows in groupby(size_rows, :mu)
        Ns = Float64.(mu_rows.N)
        MPs = Float64.(mu_rows.mp)
        pair_slopes = Float64[]
        concordant = 0
        discordant = 0
        for i in 1:(length(Ns) - 1), j in (i + 1):length(Ns)
            push!(pair_slopes, (MPs[j] - MPs[i]) / (Ns[j] - Ns[i]))
            if MPs[j] > MPs[i]
                concordant += 1
            elseif MPs[j] < MPs[i]
                discordant += 1
            end
        end

        slope = length(Ns) >= min_points ? median(pair_slopes) : NaN
        comparable_pairs = concordant + discordant
        tau = comparable_pairs == 0 ? NaN : (concordant - discordant) / comparable_pairs
        agreement = comparable_pairs == 0 || !isfinite(slope) || iszero(slope) ? NaN :
                (slope > 0 ? concordant : discordant) / comparable_pairs
        phase = !isfinite(slope) || iszero(slope) ? missing : (slope > 0 ? 0 : 1)
        push!(trends, (Float64(mu_rows.mu[1]), slope, tau, agreement,
                       length(Ns), Int(minimum(Ns)), Int(maximum(Ns)), phase))
    end

    return size_rows, trends
end

function plot_mp_with_size_invariant(
    df::DataFrame,
    filename::String;
    min_N::Int=100,
    max_N::Int=1000,
    diagnostics_filename::String=replace(filename, r"\.[^.]+$" => "_diagnostics.csv")
)
    size_rows, trends = mp_size_trend(df; min_N=min_N, max_N=max_N)
    mu_values = trends.mu
    y_limits = (-1.01, 0.01)
    p = plot(
        xlabel="Chemical potential, μ",
        ylabel="Majorana polarization, MP",
        title="MP vs μ: green = phase 0 (dMP/dN > 0), red = phase 1 (dMP/dN < 0); N=$(min_N)–$(max_N)",
        legend=false,
        grid=false,
        ylims=y_limits,
        xlims=extrema(mu_values),
        size=(1100, 360)
    )

    if length(mu_values) > 1
        boundaries = [(mu_values[i] + mu_values[i + 1]) / 2 for i in 1:(length(mu_values) - 1)]
        left_edges = vcat(mu_values[1] - (boundaries[1] - mu_values[1]), boundaries)
        right_edges = vcat(boundaries, mu_values[end] + (mu_values[end] - boundaries[end]))
        for i in eachindex(mu_values)
            phase = trends.phase[i]
            ismissing(phase) && continue
            shade = phase == 0 ? :seagreen : :indianred
            shape = Shape([left_edges[i], right_edges[i], right_edges[i], left_edges[i]],
                          [y_limits[1], y_limits[1], y_limits[2], y_limits[2]])
            plot!(p, shape; fillcolor=shade, fillalpha=0.18, linealpha=0, label="")
        end
    end

    for N in sort(unique(size_rows.N))
        N_rows = filter(:N => ==(N), size_rows)
        sort!(N_rows, :mu)
        plot!(p, N_rows.mu, N_rows.mp; linewidth=1.5, label="", color=:steelblue)
    end

    open(diagnostics_filename, "w") do io
        println(io, "mu,robust_slope_dMP_dN,kendall_tau,direction_agreement,n_sizes,min_N,max_N,phase")
        for row in eachrow(trends)
            phase = ismissing(row.phase) ? "" : string(row.phase)
            println(io, "$(row.mu),$(row.robust_slope),$(row.kendall_tau),$(row.direction_agreement),$(row.n_sizes),$(row.min_N),$(row.max_N),$phase")
        end
    end

    savefig(p, filename)
    println("Invariant plot saved to: $filename")
    println("Per-mu trend diagnostics saved to: $diagnostics_filename")
    println("Phase assignment: positive robust slope = 0; negative robust slope = 1; zero/insufficient data = unclassified.")
    println("Direction agreement is the fraction of non-tied size pairs matching the Theil-Sen slope sign.")
    return p, trends
end

# plt_mp_vs_mu_multiline_L(df, joinpath(output_path, "mp_vs_mu_all_L_restricted.png"); restrict_Ns=true, max_N=1000)

## combine dfs df_new and df_old filtering for unique N
Ns_new = unique(df_new.N)
df_old_filtered = filter(row -> !(row.N in Ns_new), df_old)
# combine gracefully if one df is missing columns the other one has
df_combined = vcat(df_new, df_old_filtered; cols=:union)

Ns_combined = unique(df_combined.N)
ignore_Ns = [5, 10, 50, 800, 5000, 10000]

plt_mp_vs_mu_multiline_L(df_combined, joinpath(output_path, "mp_vs_mu_all_L_all_restricted.svg")
; restrict_Ns=true, ignore_Ns=ignore_Ns)

invariant_plot_path = joinpath(output_path, "mp_vs_mu_N100_1000_theilsen_invariant.svg")
invariant_diagnostics_path = joinpath(output_path, "mp_vs_mu_N100_1000_theilsen_diagnostics.csv")
invariant_plot, invariant_trends = plot_mp_with_size_invariant(
    df_combined,
    invariant_plot_path;
    min_N=100,
    max_N=1000,
    diagnostics_filename=invariant_diagnostics_path
)
