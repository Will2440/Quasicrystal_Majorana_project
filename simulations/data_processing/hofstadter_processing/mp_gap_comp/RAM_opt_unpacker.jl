module MPGapUnpacking

using BSON
using DataFrames
using PrettyTables
using Printf
using ProgressMeter

include("../../bson_unpacker.jl")

function unpack_and_slice(
    folder_path_hof,
    base_slice_dir;
    mp_targ=-1.0,
    mp_tol=5e-2,
    plt_N_range_arg=nothing,
    plt_Delta_range_arg=nothing,
    plt_tn_range_arg=nothing,
    plt_phi_range_arg=nothing,
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
    plt_Delta_range = isnothing(plt_Delta_range_arg) ? Delta_range : plt_Delta_range_arg
    plt_tn_range = isnothing(plt_tn_range_arg) ? t_n_range : plt_tn_range_arg
    plt_phi_range = isnothing(plt_phi_range_arg) ? phi_range : plt_phi_range_arg
    plt_phason_range = isnothing(plt_phason_range_arg) ? phason_range : plt_phason_range_arg

    println("plt_N_range: ", plt_N_range)
    println("plt_Delta_range: ", plt_Delta_range)
    println("plt_tn_range: ", plt_tn_range)
    println("plt_phi_range: ", plt_phi_range)
    println("plt_phason_range: ", plt_phason_range)

    # Ensure base_slice_dir exists
    isdir(base_slice_dir) || mkpath(base_slice_dir)
    println("Saving base slices to folder: ", base_slice_dir)

    # Dictionary to map combos to file paths
    combo_to_files = Dict{Tuple{Int, Float64, Tuple{Vararg{Float64}}, Float64, Float64}, String}()

    println("Caching base slices grouped by (N, Delta, t_n, phi, phason)...")
    grouped = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phi, :phason])
    total_groups = length(grouped)

    @showprogress for (i, g) in enumerate(grouped)
        N_val = g.N[1]
        N_val in plt_N_range || continue

        Delta_val = g.Delta[1]
        Delta_val in plt_Delta_range || continue

        t_n_vec = g.t_n[1]
        t_n_vec in plt_tn_range || continue

        phi_val = g.phi[1]
        phi_val in plt_phi_range || continue

        phason_val = g.phason[1]
        phason_val in plt_phason_range || continue

        # Create slice DataFrame: rows for each mu in the group
        slice_df = DataFrame(phi=phi_val, phason=phason_val, mu=g.mu, eigenvalues=g.eigenvalues, mp=g.mp, disc_mp=g.disc_mp)

        bytes = Base.summarysize(slice_df)
        row_cnt = nrow(slice_df)
        println("Group $i/$total_groups -> rows=$row_cnt, size≈$(round(bytes/1e6, digits=2)) MB (N=$N_val, Delta=$Delta_val, t_n=$t_n_vec, phi=$phi_val, phason=$phason_val)")

        t_n_tuple = Tuple(t_n_vec)
        fname = @sprintf("slice_N%d_Delta%.5f_t1%.5f_t2%.5f_phi%.5f_phason%.5f.bson",
                         N_val, Delta_val, t_n_tuple[1], t_n_tuple[2], phi_val, phason_val)
        fpath = joinpath(base_slice_dir, fname)
        BSON.bson(fpath, Dict(:df_slice => slice_df))
        combo_to_files[(N_val, Delta_val, t_n_tuple, phi_val, phason_val)] = fpath
    end

    df_base = nothing
    GC.gc()
    println("Stored $(length(combo_to_files)) sliced combos.")

    # Save combo_to_files for reuse in processing script
    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")

    return combo_to_files
end


end  # module MPGapUnpacking