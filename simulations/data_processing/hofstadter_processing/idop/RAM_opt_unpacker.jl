module IDOPUnpacking

using BSON
using DataFrames
using PrettyTables
using Printf
using ProgressMeter
using Glob

include("../../bson_unpacker.jl")

function unpack_and_slice_by_slope(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_phi_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phason_range_arg=nothing
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
            phi_val = phi_sub.phi[1]
            mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phi_sub)]
            push!(slice_rows, (phi=phi_val, mu_mp=mu_mp))
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
    plt_phason_range_arg=nothing
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
            phason_val = phason_sub.phason[1]
            mu_range = [row.mu for row in eachrow(phason_sub)]
            mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phason_sub)]
            min_pos_eigs = [minimum(eig for eig in row.eigenvalues if eig > 0) for row in eachrow(phason_sub)]
            min_neg_eigs = [maximum(eig for eig in row.eigenvalues if eig < 0) for row in eachrow(phason_sub)]
            push!(slice_rows, (phason=phason_val, mu_range=mu_range, mu_mp=mu_mp, min_pos_eigs=min_pos_eigs, min_neg_eigs=min_neg_eigs))
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
########### Partintioned unpackers ##############
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

function _merge_slice!(fpath::String, slice_df::DataFrame, key_col::Symbol)
    if isfile(fpath)
        existing = BSON.load(fpath)[:df_slice]
        combined = vcat(existing, slice_df; cols=:union)
        combined = combined[end:-1:1, :]
        unique!(combined, key_col)
        slice_df = combined[end:-1:1, :]
    end
    BSON.bson(fpath, Dict(:df_slice => slice_df))
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
    n_partitions::Int=1
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
                phi_val = phi_sub.phi[1]
                mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phi_sub)]
                push!(slice_rows, (phi=phi_val, mu_mp=mu_mp))
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
    n_partitions::Int=1
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
                phason_val = phason_sub.phason[1]
                mu_range = [row.mu for row in eachrow(phason_sub)]
                mu_mp = [row.disc_mp == 1 ? row.mu : missing for row in eachrow(phason_sub)]
                min_pos_eigs = [minimum(eig for eig in row.eigenvalues if eig > 0; init=missing) for row in eachrow(phason_sub)]
                min_neg_eigs = [maximum(eig for eig in row.eigenvalues if eig < 0; init=missing) for row in eachrow(phason_sub)]
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