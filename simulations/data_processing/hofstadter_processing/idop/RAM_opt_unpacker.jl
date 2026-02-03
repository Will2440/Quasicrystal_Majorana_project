module IDOPUnpacking

using BSON
using DataFrames
using PrettyTables
using Printf
using ProgressMeter
using Glob

include("../../bson_unpacker.jl")

# --- NEW HELPER FUNCTION ---
function _extract_min_eigs(evals, mode::Symbol)
    if ismissing(evals) || isempty(evals)
        return missing, missing
    end

    if mode == :maj
        # :maj mode -> assume saved eigs are central ones (usually 2)
        if length(evals) == 2
            # Assuming sorted: [neg, pos]
            return evals[1], evals[2]
        else
            # Fallback for :maj if length != 2
            p = missing; n = missing
            for x in evals
                if x > 0
                    if ismissing(p) || x < p; p = x; end
                elseif x < 0
                    if ismissing(n) || x > n; n = x; end
                end
            end
            return n, p
        end
    elseif mode == :all
        # :all mode -> use central indices (assuming sorted) for reliability
        len = length(evals)
        if len >= 2
            mid = len ÷ 2
            return evals[mid], evals[mid+1]
        else
            return evals[1], evals[end]
        end
    else
        # Default fallback (value based)
        p = missing; n = missing
        for x in evals
            if x > 0
                if ismissing(p) || x < p; p = x; end
            elseif x < 0
                if ismissing(n) || x > n; n = x; end
            end
        end
        return n, p
    end
end
# ---------------------------

function _mirror_data_if_needed(df_sorted::DataFrame, mirror::Bool)
    if mirror && !isempty(df_sorted)
        # Check start value
        mu_start = df_sorted.mu[1]
        
        # Only mirror if we start near 0 (positive side)
        # Tolerance 0.1 allows for mu starting slightly above 0
        if mu_start > -1e-6 && mu_start < 0.1
            # If start is exactly 0 (within tight tolerance), we skip it to avoid double counting 0
            # If start is > 0 (e.g. 0.005), we include it in the mirror
            start_idx = abs(mu_start) < 1e-6 ? 2 : 1
            
            if start_idx <= nrow(df_sorted)
                mirrored_part = df_sorted[start_idx:end, :]
                reverse!(mirrored_part)
                mirrored_part.mu = -mirrored_part.mu # Negate mu
                
                # @info "Mirrored data" original_start=mu_start new_start=mirrored_part.mu[1] new_len=(nrow(mirrored_part)+nrow(df_sorted))
                return vcat(mirrored_part, df_sorted)
            end
        end
    end
    return df_sorted
end

function unpack_and_slice_by_slope(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing,
    mirror_data_about_mu0::Bool=false
)
    println("Unpacking Hofstadter data from folder: ", folder_path_hof)

    # Unpack data
    df_base = process_bson_files(folder_path_hof)
    isempty(df_base) && error("No data unpacked — aborting.")

    # Remove eigenvalues to save RAM, as they are not needed for IDOP processing
    select!(df_base, Not(:eigenvalues))

    df_base = calc_norms_df!(df_base)
    df_base = discretise_mp!(df_base, mp_targ, mp_tol)
    df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
    df_base[!, :phason] = getindex.(df_base.sequence_id, 2)
    println("DataFrame keynames: $(names(df_base))")

    # Get unique ranges
    N_range = sort(unique(df_base.N))
    t_n_range = sort(unique(df_base.t_n))
    Delta_range = sort(unique(df_base.Delta))
    mu_range = sort(unique(df_base.mu))
    phi_range = sort(unique(df_base.phi))
    phason_range = sort(unique(df_base.phason))
    println("N_range: ", N_range)
    println("t_n_range: ", t_n_range)
    println("Delta_range: ", Delta_range)
    println("mu_range (length, max, min): ", (length(mu_range), maximum(mu_range), minimum(mu_range)))
    println("phi_range (length, max, min): ", (length(phi_range), maximum(phi_range), minimum(phi_range)))
    println("phason_range (length, max, min): ", (length(phason_range), maximum(phason_range), minimum(phason_range)))

    # Set plotting ranges (or use full ranges if not provided)
    plt_N_range = isnothing(plt_N_range_arg) ? N_range : plt_N_range_arg
    plt_phi_range = isnothing(plt_phi_range_arg) ? phi_range : plt_phi_range_arg
    plt_Delta_range = isnothing(plt_Delta_range_arg) ? Delta_range : plt_Delta_range_arg
    plt_tn_range = isnothing(plt_tn_range_arg) ? t_n_range : plt_tn_range_arg
    plt_phason_range = isnothing(plt_phason_range_arg) ? phason_range : plt_phason_range_arg

    println("plt_N_range: ", plt_N_range)
    println("plt_phi_range: ", plt_phi_range)
    println("plt_Delta_range: ", plt_Delta_range)
    println("plt_tn_range: ", plt_tn_range)
    println("plt_phason_range: ", plt_phason_range)

    # Ensure base_slice_dir exists
    isdir(base_slice_dir) || mkpath(base_slice_dir)
    println("Saving base slices to folder: ", base_slice_dir)

    # Dictionary to map combos to file paths (key is (N, Delta, t_n, phason))
    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64}, String}()

    println("Caching base slices grouped by (N, Delta, t_n, phason)...")
    grouped = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phason])
    total_groups = length(grouped)

    @showprogress for (i, g) in enumerate(grouped)
        N_val = g.N[1]
        N_val in plt_N_range || continue

        Delta_val = g.Delta[1]
        Delta_val in plt_Delta_range || continue

        t_n_vec = g.t_n[1]
        t_n_tuple = Tuple(t_n_vec)
        phason_val = g.phason[1]
        phason_val in plt_phason_range || continue

        # Filter the group for plt_phi_range
        g_filtered = filter(row -> row.phi in plt_phi_range, g)
        isempty(g_filtered) && continue

        # Create slice DataFrame: for each phi, collect mu_mp as vector
        slice_rows = []
        for phi_sub in DataFrames.groupby(g_filtered, :phi)
            # --- FIX: Sort by mu ---
            phi_sub_sorted = sort(phi_sub, :mu)
            
            phi_sub_sorted = _mirror_data_if_needed(phi_sub_sorted, mirror_data_about_mu0)

            phi_val = phi_sub_sorted.phi[1]
            mu_range = [row.mu for row in eachrow(phi_sub_sorted)]
            mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phi_sub_sorted)]
            
            # --- UPDATED EXTRACTION LOGIC ---
            min_pos_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phi_sub_sorted))
            min_neg_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phi_sub_sorted))
            
            for (k, row) in enumerate(eachrow(phi_sub_sorted))
                n, p = _extract_min_eigs(row.eigenvalues, eig_save_mode)
                min_neg_eigs[k] = n
                min_pos_eigs[k] = p
            end
            # --------------------------------

            push!(slice_rows, (phi=phi_val, mu_range=mu_range, mu_mp=mu_mp, min_pos_eigs=min_pos_eigs, min_neg_eigs=min_neg_eigs))
        end
        slice_df = DataFrame(slice_rows)

        fname = @sprintf("slice_N%d_Delta%.5f_t1%.5f_t2%.5f_phason%.5f.bson",
                         N_val, Delta_val, t_n_tuple[1], t_n_tuple[2], phason_val)
        fpath = joinpath(base_slice_dir, fname)
        BSON.bson(fpath, Dict(:df_slice => slice_df))
        combo_to_files[(N_val, Delta_val, t_n_tuple, phason_val)] = fpath
    end

    df_base = nothing
    GC.gc()
    println("Stored $(length(combo_to_files)) sliced combos.")

    # Save combo_to_files for reuse in processing script
    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")

    return combo_to_files
end

function unpack_and_slice_by_phason(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing,
    eig_save_mode::Symbol=:all,  # <--- NEW ARGUMENT
    mirror_data_about_mu0::Bool=false
)
    println("Unpacking Hofstadter data from folder: ", folder_path_hof)

    # Unpack data
    df_base = process_bson_files(folder_path_hof)
    isempty(df_base) && error("No data unpacked — aborting.")

    df_base = calc_norms_df!(df_base)
    df_base = discretise_mp!(df_base, mp_targ, mp_tol)
    df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
    df_base[!, :phason] = getindex.(df_base.sequence_id, 2)
    println("DataFrame keynames: $(names(df_base))")

    # Get unique ranges
    N_range = sort(unique(df_base.N))
    t_n_range = sort(unique(df_base.t_n))
    Delta_range = sort(unique(df_base.Delta))
    mu_range = sort(unique(df_base.mu))
    phi_range = sort(unique(df_base.phi))
    phason_range = sort(unique(df_base.phason))
    println("N_range: ", N_range)
    println("t_n_range: ", t_n_range)
    println("Delta_range: ", Delta_range)
    println("mu_range (length, max, min): ", (length(mu_range), maximum(mu_range), minimum(mu_range)))
    println("phi_range (length, max, min): ", (length(phi_range), maximum(phi_range), minimum(phi_range)))
    println("phason_range (length, max, min): ", (length(phason_range), maximum(phason_range), minimum(phason_range)))

    # Set plotting ranges (or use full ranges if not provided)
    plt_N_range = isnothing(plt_N_range_arg) ? N_range : plt_N_range_arg
    plt_phi_range = isnothing(plt_phi_range_arg) ? phi_range : plt_phi_range_arg
    plt_Delta_range = isnothing(plt_Delta_range_arg) ? Delta_range : plt_Delta_range_arg
    plt_tn_range = isnothing(plt_tn_range_arg) ? t_n_range : plt_tn_range_arg
    plt_phason_range = isnothing(plt_phason_range_arg) ? phason_range : plt_phason_range_arg

    println("plt_N_range: ", plt_N_range)
    println("plt_phi_range: ", plt_phi_range)
    println("plt_Delta_range: ", plt_Delta_range)
    println("plt_tn_range: ", plt_tn_range)
    println("plt_phason_range: ", plt_phason_range)

    # Ensure base_slice_dir exists
    isdir(base_slice_dir) || mkpath(base_slice_dir)
    println("Saving base slices to folder: ", base_slice_dir)

    # Dictionary to map combos to file paths (key is (N, Delta, t_n, phi))
    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64}, String}()

    println("Caching base slices grouped by (N, Delta, t_n, phi)...")
    grouped = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phi])
    total_groups = length(grouped)

    @showprogress for (i, g) in enumerate(grouped)
        N_val = g.N[1]
        N_val in plt_N_range || continue

        Delta_val = g.Delta[1]
        Delta_val in plt_Delta_range || continue

        t_n_vec = g.t_n[1]
        t_n_tuple = Tuple(t_n_vec)
        phi_val = g.phi[1]
        phi_val in plt_phi_range || continue

        # Filter the group for plt_phason_range
        g_filtered = filter(row -> row.phason in plt_phason_range, g)
        isempty(g_filtered) && continue

        # Create slice DataFrame: for each phason, collect mu_mp as vector
        slice_rows = []
        for phason_sub in DataFrames.groupby(g_filtered, :phason)
            # --- FIX: Sort by mu to ensure vector order is monotonic ---
            phason_sub_sorted = sort(phason_sub, :mu)
            
            phason_sub_sorted = _mirror_data_if_needed(phason_sub_sorted, mirror_data_about_mu0)

            phason_val = phason_sub_sorted.phason[1]
            mu_range = [row.mu for row in eachrow(phason_sub_sorted)]
            mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phason_sub_sorted)]
            raw_mp = [row.mp for row in eachrow(phason_sub_sorted)]

            # --- UPDATED EXTRACTION LOGIC ---
            min_pos_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phason_sub_sorted))
            min_neg_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phason_sub_sorted))
            
            for (k, row) in enumerate(eachrow(phason_sub_sorted))
                n, p = _extract_min_eigs(row.eigenvalues, eig_save_mode)
                min_neg_eigs[k] = n
                min_pos_eigs[k] = p
            end
            # --------------------------------

            push!(slice_rows, (phason=phason_val, mu_range=mu_range, mu_mp=mu_mp, raw_mp=raw_mp, min_pos_eigs=min_pos_eigs, min_neg_eigs=min_neg_eigs))
        end
        slice_df = DataFrame(slice_rows)

        fname = "slice_N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phi$(phi_val).bson"
        fpath = joinpath(base_slice_dir, fname)
        BSON.bson(fpath, Dict(:df_slice => slice_df))
        combo_to_files[(N_val, Delta_val, t_n_tuple, phi_val)] = fpath
    end

    df_base = nothing
    GC.gc()
    println("Stored $(length(combo_to_files)) sliced combos.")

    # Save combo_to_files for reuse in processing script
    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")

    return combo_to_files
end

#################################################
########## No MP Processing Unpacker ############
#################################################
function unpack_and_slice_by_slope_raw_mp(
    folder_path_hof,
    base_slice_dir;
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing,
    mirror_data_about_mu0::Bool=false
)
    println("Unpacking Hofstadter data from folder: ", folder_path_hof)

    # Unpack data
    df_base = process_bson_files(folder_path_hof)
    isempty(df_base) && error("No data unpacked — aborting.")

    # Remove eigenvalues to save RAM, as they are not needed for IDOP processing
    select!(df_base, Not(:eigenvalues))

    df_base = calc_norms_df!(df_base)
    # DO NOT discretise mp here!
    df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
    df_base[!, :phason] = getindex.(df_base.sequence_id, 2)
    println("DataFrame keynames: $(names(df_base))")

    # Get unique ranges
    N_range = sort(unique(df_base.N))
    t_n_range = sort(unique(df_base.t_n))
    Delta_range = sort(unique(df_base.Delta))
    mu_range = sort(unique(df_base.mu))
    phi_range = sort(unique(df_base.phi))
    phason_range = sort(unique(df_base.phason))
    println("N_range: ", N_range)
    println("t_n_range: ", t_n_range)
    println("Delta_range: ", Delta_range)
    println("mu_range (length, max, min): ", (length(mu_range), maximum(mu_range), minimum(mu_range)))
    println("phi_range (length, max, min): ", (length(phi_range), maximum(phi_range), minimum(phi_range)))
    println("phason_range (length, max, min): ", (length(phason_range), maximum(phason_range), minimum(phason_range)))

    # Set plotting ranges (or use full ranges if not provided)
    plt_N_range = isnothing(plt_N_range_arg) ? N_range : plt_N_range_arg
    plt_phi_range = isnothing(plt_phi_range_arg) ? phi_range : plt_phi_range_arg
    plt_Delta_range = isnothing(plt_Delta_range_arg) ? Delta_range : plt_Delta_range_arg
    plt_tn_range = isnothing(plt_tn_range_arg) ? t_n_range : plt_tn_range_arg
    plt_phason_range = isnothing(plt_phason_range_arg) ? phason_range : plt_phason_range_arg

    println("plt_N_range: ", plt_N_range)
    println("plt_phi_range: ", plt_phi_range)
    println("plt_Delta_range: ", plt_Delta_range)
    println("plt_tn_range: ", plt_tn_range)
    println("plt_phason_range: ", plt_phason_range)

    # Ensure base_slice_dir exists
    isdir(base_slice_dir) || mkpath(base_slice_dir)
    println("Saving base slices to folder: ", base_slice_dir)

    # Dictionary to map combos to file paths (key is (N, Delta, t_n, phason))
    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64}, String}()

    println("Caching base slices grouped by (N, Delta, t_n, phason)...")
    grouped = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phason])
    total_groups = length(grouped)

    @showprogress for (i, g) in enumerate(grouped)
        N_val = g.N[1]
        N_val in plt_N_range || continue

        Delta_val = g.Delta[1]
        Delta_val in plt_Delta_range || continue

        t_n_vec = g.t_n[1]
        t_n_tuple = Tuple(t_n_vec)
        phason_val = g.phason[1]
        phason_val in plt_phason_range || continue

        # Filter the group for plt_phi_range
        g_filtered = filter(row -> row.phi in plt_phi_range, g)
        isempty(g_filtered) && continue

        # Create slice DataFrame: for each phi, collect mu_mp as vector
        slice_rows = []
        for phi_sub in DataFrames.groupby(g_filtered, :phi)
            # --- FIX: Sort by mu ---
            phi_sub_sorted = sort(phi_sub, :mu)
            
            phi_sub_sorted = _mirror_data_if_needed(phi_sub_sorted, mirror_data_about_mu0)

            phi_val = phi_sub_sorted.phi[1]
            mu_range = [row.mu for row in eachrow(phi_sub_sorted)]
            mp_raw = [row.mp for row in eachrow(phi_sub_sorted)]
            push!(slice_rows, (phi=phi_val, mu_range=mu_range, mp_raw=mp_raw))
        end
        slice_df = DataFrame(slice_rows)

        fname = @sprintf("slice_N%d_Delta%.5f_t1%.5f_t2%.5f_phason%.5f.bson",
                         N_val, Delta_val, t_n_tuple[1], t_n_tuple[2], phason_val)
        fpath = joinpath(base_slice_dir, fname)
        BSON.bson(fpath, Dict(:df_slice => slice_df))
        combo_to_files[(N_val, Delta_val, t_n_tuple, phason_val)] = fpath
    end

    df_base = nothing
    GC.gc()
    println("Stored $(length(combo_to_files)) sliced combos.")

    # Save combo_to_files for reuse in processing script
    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")

    return combo_to_files
end



#################################################
########### Partitioned unpackers ##############
#################################################

_in_range(val, rng) = isnothing(rng) ? true : (val in rng)

function _partition_file_paths(file_paths::Vector{String}, n_partitions::Int)
    n_partitions = max(1, n_partitions)
    total = length(file_paths)
    chunk = total == 0 ? 0 : ceil(Int, total / n_partitions)
    partitions = Vector{Vector{String}}()
    chunk == 0 && return partitions
    for i in 1:n_partitions
        start_idx = (i-1) * chunk + 1
        start_idx > total && break
        end_idx = min(i * chunk, total)
        push!(partitions, file_paths[start_idx:end_idx])
    end
    return partitions
end

# function _merge_slice!(fpath::String, slice_df::DataFrame, key_col::Symbol)
#     if isfile(fpath)
#         existing = BSON.load(fpath)[:df_slice]
#         combined = vcat(existing, slice_df; cols=:union)
#         combined = combined[end:-1:1, :]
#         unique!(combined, key_col)
#         slice_df = combined[end:-1:1, :]
#     end
#     BSON.bson(fpath, Dict(:df_slice => slice_df))
# end

function _merge_slice!(fpath::String, new_df::DataFrame, key_col::Symbol)
    # If file doesn't exist, just save the new dataframe
    if !isfile(fpath)
        BSON.bson(fpath, Dict(:df_slice => new_df))
        return
    end

    # Load existing data
    existing_df = BSON.load(fpath)[:df_slice]
    
    # We need to merge 'new_df' into 'existing_df'.
    # 1. Identify rows with keys that are already in 'existing_df'
    # 2. For those rows, concatenate the vector columns and re-sort by mu
    # 3. For completely new rows, validly append them
    
    # Map key -> index in existing_df
    existing_map = Dict(r[key_col] => i for (i, r) in enumerate(eachrow(existing_df)))
    
    rows_to_add = []
    
    for new_row in eachrow(new_df)
        k = new_row[key_col]
        if haskey(existing_map, k)
            # --- MERGE LOGIC ---
            idx = existing_map[k]
            old_row = existing_df[idx, :]
            
            # We assume 'mu_range' exists and is the sorting key
            # (Note: In 'unpack_by_slope', 'mu_range' might not exist in older versions, 
            #  but 'unpack_by_phason' creates it. If missing, we can't safely sort.)
            if :mu_range in propertynames(existing_df)
                mu_old = old_row.mu_range
                mu_new = new_row.mu_range
                
                # Combine mu to find the sort permutation
                mu_comb = vcat(mu_old, mu_new)
                perm = sortperm(mu_comb)
                
                # Update all vector columns
                for col_name in names(existing_df)
                    val_old = old_row[col_name]
                    val_new = new_row[col_name]
                    
                    # Merge if both are vectors and it's not the key column
                    if isa(val_old, Vector) && isa(val_new, Vector) && col_name != string(key_col)
                        comb = vcat(val_old, val_new)
                        # Apply sort permutation if lengths match
                        if length(comb) == length(perm)
                            existing_df[idx, col_name] = comb[perm]
                        end
                    end
                end
            else
                # Fallback purely for safety if mu_range is missing (e.g. legacy slope mode)
                # Just append vectors without sorting (risky but better than overwrite)
                for col_name in names(existing_df)
                    val_old = old_row[col_name]
                    val_new = new_row[col_name]
                    if isa(val_old, Vector) && isa(val_new, Vector) && col_name != string(key_col)
                        existing_df[idx, col_name] = vcat(val_old, val_new)
                    end
                end
            end
        else
            # New unique row
            push!(rows_to_add, new_row)
        end
    end
    
    if !isempty(rows_to_add)
        to_append = DataFrame(rows_to_add)
        # Ensure column order matches before appending
        select!(to_append, names(existing_df))
        append!(existing_df, to_append)
    end

    BSON.bson(fpath, Dict(:df_slice => existing_df))
end


function unpack_and_slice_by_slope_partitioned(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing,
    n_partitions::Int=1,
    mirror_data_about_mu0::Bool=false
)
    println("Unpacking Hofstadter data (partitioned slope mode) from folder: ", folder_path_hof)

    file_paths = sort(glob("*.bson", folder_path_hof))
    partitions = _partition_file_paths(file_paths, n_partitions)
    isempty(partitions) && error("No BSON files found to unpack.")
    isdir(base_slice_dir) || mkpath(base_slice_dir)

    allowed_N = plt_N_range_arg
    allowed_phi = plt_phi_range_arg
    allowed_Delta = plt_Delta_range_arg
    allowed_tn = plt_tn_range_arg
    allowed_phason = plt_phason_range_arg

    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64}, String}()
    failed_overall = String[]
    seen_N = Set{Int}()
    seen_phi = Set{Float64}()
    seen_phason = Set{Float64}()

    for (idx, subset) in enumerate(partitions)
        isempty(subset) && continue
        desc = "partition $(idx)/$(n_partitions)"
        println("Processing $desc with $(length(subset)) files…")
        res = process_bson_file_list(subset; desc=desc)
        append!(failed_overall, res.failed_files)
        df_chunk = res.df
        isempty(df_chunk) && continue

        if :eigenvalues in names(df_chunk)
            select!(df_chunk, Not(:eigenvalues))
        end

        df_chunk = calc_norms_df!(df_chunk)
        df_chunk = discretise_mp!(df_chunk, mp_targ, mp_tol)
        df_chunk[!, :phi] = getindex.(df_chunk.sequence_id, 1)
        df_chunk[!, :phason] = getindex.(df_chunk.sequence_id, 2)

        union!(seen_N, unique(df_chunk.N))
        union!(seen_phi, unique(df_chunk.phi))
        union!(seen_phason, unique(df_chunk.phason))

        grouped = DataFrames.groupby(df_chunk, [:N, :Delta, :t_n, :phason])
        for g in grouped
            N_val = g.N[1]
            Delta_val = g.Delta[1]
            t_n_vec = g.t_n[1]; t_n_tuple = Tuple(t_n_vec)
            phason_val = g.phason[1]
            _in_range(N_val, allowed_N) || continue
            _in_range(Delta_val, allowed_Delta) || continue
            _in_range(t_n_tuple, allowed_tn) || continue
            _in_range(phason_val, allowed_phason) || continue

            g_filtered = isnothing(allowed_phi) ? g : filter(row -> row.phi in allowed_phi, g)
            isempty(g_filtered) && continue

            slice_rows = []
            for phi_sub in DataFrames.groupby(g_filtered, :phi)
                # --- FIX: Sort by mu ---
                phi_sub_sorted = sort(phi_sub, :mu)
                
                phi_sub_sorted = _mirror_data_if_needed(phi_sub_sorted, mirror_data_about_mu0)

                phi_val = phi_sub_sorted.phi[1]
                mu_range = [row.mu for row in eachrow(phi_sub_sorted)]
                mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phi_sub_sorted)]
                push!(slice_rows, (phi=phi_val, mu_range=mu_range, mu_mp=mu_mp))
            end
            isempty(slice_rows) && continue
            slice_df = DataFrame(slice_rows)

            fname = @sprintf("slice_N%d_Delta%.5f_t1%.5f_t2%.5f_phason%.5f.bson",
                             N_val, Delta_val, t_n_tuple[1], t_n_tuple[2], phason_val)
            fpath = joinpath(base_slice_dir, fname)
            _merge_slice!(fpath, slice_df, :phi)
            combo_to_files[(N_val, Delta_val, t_n_tuple, phason_val)] = fpath
        end

        df_chunk = nothing
        GC.gc()
    end

    isempty(combo_to_files) && error("No data unpacked — aborting.")
    if !isempty(failed_overall)
        @warn "Unreadable BSON files detected during partitioned slope unpack" failed_count=length(failed_overall)
    end

    println("Seen N values: ", sort!(collect(seen_N)))
    println("Seen phi values: ", sort!(collect(seen_phi)))
    println("Seen phason values: ", sort!(collect(seen_phason)))

    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")
    return combo_to_files
end

function unpack_and_slice_by_phason_partitioned(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing,
    n_partitions::Int=1,
    eig_save_mode::Symbol=:all,  # <--- NEW ARGUMENT
    mirror_data_about_mu0::Bool=false
)
    println("Unpacking Hofstadter data (partitioned phason mode) from folder: ", folder_path_hof)

    file_paths = sort(glob("*.bson", folder_path_hof))
    partitions = _partition_file_paths(file_paths, n_partitions)
    isempty(partitions) && error("No BSON files found to unpack.")
    isdir(base_slice_dir) || mkpath(base_slice_dir)

    allowed_N = plt_N_range_arg
    allowed_phi = plt_phi_range_arg
    allowed_Delta = plt_Delta_range_arg
    allowed_tn = plt_tn_range_arg
    allowed_phason = plt_phason_range_arg

    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64}, String}()
    failed_overall = String[]
    seen_N = Set{Int}()
    seen_phi = Set{Float64}()
    seen_phason = Set{Float64}()

    for (idx, subset) in enumerate(partitions)
        isempty(subset) && continue
        desc = "partition $(idx)/$(n_partitions)"
        println("Processing $desc with $(length(subset)) files…")
        res = process_bson_file_list(subset; desc=desc)
        append!(failed_overall, res.failed_files)
        df_chunk = res.df
        isempty(df_chunk) && continue

        df_chunk = calc_norms_df!(df_chunk)
        df_chunk = discretise_mp!(df_chunk, mp_targ, mp_tol)
        df_chunk[!, :phi] = getindex.(df_chunk.sequence_id, 1)
        df_chunk[!, :phason] = getindex.(df_chunk.sequence_id, 2)

        union!(seen_N, unique(df_chunk.N))
        union!(seen_phi, unique(df_chunk.phi))
        union!(seen_phason, unique(df_chunk.phason))

        grouped = DataFrames.groupby(df_chunk, [:N, :Delta, :t_n, :phi])
        for g in grouped
            N_val = g.N[1]
            Delta_val = g.Delta[1]
            t_n_vec = g.t_n[1]; t_n_tuple = Tuple(t_n_vec)
            phi_val = g.phi[1]
            _in_range(N_val, allowed_N) || continue
            _in_range(Delta_val, allowed_Delta) || continue
            _in_range(t_n_tuple, allowed_tn) || continue
            _in_range(phi_val, allowed_phi) || continue

            g_filtered = isnothing(allowed_phason) ? g : filter(row -> row.phason in allowed_phason, g)
            isempty(g_filtered) && continue

            slice_rows = []
            for phason_sub in DataFrames.groupby(g_filtered, :phason)
                # --- FIX: Sort by mu ---
                phason_sub_sorted = sort(phason_sub, :mu)
                
                phason_sub_sorted = _mirror_data_if_needed(phason_sub_sorted, mirror_data_about_mu0)

                phason_val = phason_sub_sorted.phason[1]
                mu_range = [row.mu for row in eachrow(phason_sub_sorted)]
                mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phason_sub_sorted)]
                
                # --- UPDATED EXTRACTION LOGIC ---
                min_pos_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phason_sub_sorted))
                min_neg_eigs = Vector{Union{Float64, Missing}}(undef, nrow(phason_sub_sorted))
                
                for (k, row) in enumerate(eachrow(phason_sub_sorted))
                    n, p = _extract_min_eigs(row.eigenvalues, eig_save_mode)
                    min_neg_eigs[k] = n
                    min_pos_eigs[k] = p
                end
                # --------------------------------

                push!(slice_rows, (phason=phason_val, mu_range=mu_range, mu_mp=mu_mp,
                                   min_pos_eigs=min_pos_eigs, min_neg_eigs=min_neg_eigs))
            end
            isempty(slice_rows) && continue
            slice_df = DataFrame(slice_rows)

            fname = "slice_N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phi$(phi_val).bson"
            fpath = joinpath(base_slice_dir, fname)
            _merge_slice!(fpath, slice_df, :phason)
            combo_to_files[(N_val, Delta_val, t_n_tuple, phi_val)] = fpath
        end

        df_chunk = nothing
        GC.gc()
    end

    isempty(combo_to_files) && error("No data unpacked — aborting.")
    if !isempty(failed_overall)
        @warn "Unreadable BSON files detected during partitioned phason unpack" failed_count=length(failed_overall)
    end

    println("Seen N values: ", sort!(collect(seen_N)))
    println("Seen phi values: ", sort!(collect(seen_phi)))
    println("Seen phason values: ", sort!(collect(seen_phason)))

    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")
    return combo_to_files
end

end  # IDOPUnpacking