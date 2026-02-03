module MPGapSplitUnpacking

using BSON
using DataFrames
using PrettyTables
using Printf
using ProgressMeter

include("../../bson_unpacker.jl")

function load_raw_dataset(folder_path_hof)
    println("Loading raw Hofstadter data from folder: ", folder_path_hof)
    df_base = process_bson_files(folder_path_hof)
    isempty(df_base) && error("No data unpacked — aborting.")
    println("Calculating norms and extracting parameters...")
    df_base = calc_norms_df!(df_base)
    df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
    df_base[!, :phason] = getindex.(df_base.sequence_id, 2)
    println("Raw dataset loaded. Rows: $(nrow(df_base)).")
    return df_base
end

# Helper to count gaps (runs of zeros)
function count_zero_runs(vec::AbstractVector{<:Integer})
    cnt = 0
    runlen = 0
    min_gap_length = 1
    for v in vec
        if v == 0
            runlen += 1
        else
            if runlen >= min_gap_length
                cnt += 1
            end
            runlen = 0
        end
    end
    if runlen >= min_gap_length
        cnt += 1
    end
    return cnt
end

# Eigenvalue analysis functions moved to processing.jl

function process_tolerance_sweeps(
    df_base::DataFrame,
    output_dir::String,
    mp_tol_range::Vector{Float64};
    mp_targ::Float64=-1.0,
    plt_N_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_phason_range_arg=nothing
)
    isdir(output_dir) || mkpath(output_dir)
    println("Processing tolerance sweeps. Output root: ", output_dir)

    # 1. Group by the "Folder" parameters
    grouped_folders = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phason])
    
    combo_to_file = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64, Float64}, String}()

    p = Progress(length(grouped_folders), desc="Processing folders: ")

    for g_folder in grouped_folders
        N_val = g_folder.N[1]
        Delta_val = g_folder.Delta[1]
        t_n_vec = g_folder.t_n[1]
        phason_val = g_folder.phason[1]
        t_n_tuple = Tuple(t_n_vec)

        # --- Filtering (Folder Level) ---
        if !isnothing(plt_N_range_arg) && !(N_val in plt_N_range_arg); next!(p); continue; end
        if !isnothing(plt_Delta_range_arg) && !(Delta_val in plt_Delta_range_arg); next!(p); continue; end
        if !isnothing(plt_phason_range_arg) && !(phason_val in plt_phason_range_arg); next!(p); continue; end
        if !isnothing(plt_tn_range_arg)
            tn_range_tuples = [Tuple(x) for x in plt_tn_range_arg]
            if !(t_n_tuple in tn_range_tuples); next!(p); continue; end
        end

        # --- Determine Mu Range String ---
        all_mus = g_folder.mu
        mu_min = minimum(all_mus)
        mu_max = maximum(all_mus)
        mu_len = length(unique(all_mus))
        mu_str = @sprintf("mu(%.1f-%.1f-%d)", mu_min, mu_max, mu_len)

        # --- Create Folder ---
        folder_name = @sprintf("N%d_Delta%.3f_t1%.3f_t2%.3f_phason%.3f_%s", 
                               N_val, Delta_val, t_n_tuple[1], t_n_tuple[2], phason_val, mu_str)
        full_folder_path = joinpath(output_dir, folder_name)
        isdir(full_folder_path) || mkpath(full_folder_path)

        # --- Process Phis inside this folder ---
        grouped_phi = DataFrames.groupby(g_folder, :phi)

        for g_phi in grouped_phi
            phi_val = g_phi.phi[1]

            if !isnothing(plt_phi_range_arg) && !(phi_val in plt_phi_range_arg); continue; end

            # Sort by mu to ensure the vector is ordered
            g_phi_sorted = sort(g_phi, :mu)
            
            # Extract the full vectors from the columns
            mu_vals = g_phi_sorted.mu
            mp_vals = g_phi_sorted.mp
            eig_vals = g_phi_sorted.eigenvalues
            # println("type of eig_vals: ", typeof(eig_vals))
            max_eig = nothing
            if :eigenvalues in propertynames(g_phi_sorted)
                mu0_row = findfirst(x -> isapprox(x, 0.0; atol=1e-8), mu_vals)
                if !isnothing(mu0_row)
                    mu0_eig = g_phi_sorted.eigenvalues[mu0_row]
                    if !(ismissing(mu0_eig) || isempty(mu0_eig))
                        max_eig = maximum(abs.(mu0_eig))
                    else
                        @warn "No valid eigenvalues at mu=0.0 for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val. max_eig will be nothing."
                    end
                else
                    @warn "No mu=0.0 found for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val. max_eig will be nothing."
                end
            else
                @warn "No eigenvalues column found for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val. max_eig will be nothing."
            end

            # println("Max eig: ", max_eig)
            
            # Build Matrix
            n_mu = length(mu_vals)
            n_tols = length(mp_tol_range)
            disc_mp_matrix = zeros(Int8, n_mu, n_tols)
            mp_real = real.(mp_vals) 

            for (i, tol) in enumerate(mp_tol_range)
                disc_mp_matrix[:, i] .= Int8.(abs.(mp_real .- mp_targ) .<= tol)
            end

            # Save BSON
            fname = @sprintf("sweep_phi%.15f.bson", phi_val)
            fpath = joinpath(full_folder_path, fname)

            data_dict = Dict(
                :mu => mu_vals,
                :mp_raw => mp_vals,
                :mp_tol_range => mp_tol_range,
                :disc_mp_matrix => disc_mp_matrix,
                :params => (N=N_val, Delta=Delta_val, t_n=t_n_tuple, phi=phi_val, phason=phason_val, mp_targ=mp_targ),
                :max_eig => max_eig
            )
            BSON.bson(fpath, data_dict)
            
            combo_to_file[(N_val, Delta_val, t_n_tuple, phi_val, phason_val)] = fpath
        end
        next!(p)
    end
    
    BSON.bson(joinpath(output_dir, "sweep_index.bson"), Dict(:combo_to_file => combo_to_file))
    
    return combo_to_file
end

end  # module MPGapSplitUnpacking