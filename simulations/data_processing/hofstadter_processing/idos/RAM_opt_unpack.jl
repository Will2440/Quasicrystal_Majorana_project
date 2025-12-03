module IDOSUnpacking

using BSON
using DataFrames
using PrettyTables
using Printf
using ProgressMeter

include("../../bson_unpacker.jl")

function unpack_and_slice_by_slope(
    folder_path_hof, 
    base_slice_dir; 
    plt_mu_range_arg=nothing, 
    plt_tn_range_arg=nothing,
    plt_Delta_range_arg=nothing, 
    plt_phason_range_arg=nothing
)
   println("Unpacking Hofstadter data from folder: ", folder_path_hof)

    # df_base = process_bson_files_threaded(folder_path_hof)
    df_base = process_bson_files(folder_path_hof)
    isempty(df_base) && error("No data unpacked — aborting.")

    df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
    df_base[!, :phason] = getindex.(df_base.sequence_id, 2)
    println("DataFrame keynames: $(names(df_base))")

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

    plt_tn_range = isnothing(plt_tn_range_arg) ? t_n_range : plt_tn_range_arg
    plt_mu_range = isnothing(plt_mu_range_arg) ? mu_range : plt_mu_range_arg
    plt_Delta_range = isnothing(plt_Delta_range_arg) ? Delta_range : plt_Delta_range_arg
    plt_phason_range = isnothing(plt_phason_range_arg) ? phason_range : plt_phason_range_arg

    println("plt_tn_range: ", plt_tn_range)
    println("plt_mu_range: ", plt_mu_range)
    println("plt_Delta_range: ", plt_Delta_range)
    println("plt_phason_range: ", plt_phason_range)

    # base_slice_dir = normpath(joinpath(base_results_dir, "idos_base_slices"))
    isdir(base_slice_dir) || mkpath(base_slice_dir)
    println("Saving base slices to folder: ", base_slice_dir)

    combo_to_files = Dict{Tuple{Float64, Tuple{Vararg{Float64}}, Float64, Float64}, String}()
    println("Caching base slices grouped by (Delta, t_n, phason, mu)...")
    grouped = DataFrames.groupby(df_base, [:Delta, :t_n, :phason, :mu])
    total_groups = length(grouped)

    @showprogress for (i, g) in enumerate(grouped)
        mu_val = g.mu[1]
        mu_val in plt_mu_range || continue

        slice_df = DataFrame(g)
        Delta = slice_df.Delta[1]
        t_n_vec = slice_df.t_n[1]
        t_n_tuple = Tuple(t_n_vec)
        phason = slice_df.phason[1]

        bytes = Base.summarysize(slice_df)
        row_cnt = nrow(slice_df)
        eig_len = length(slice_df.eigenvalues[1])
        # println("Group $i/$total_groups -> rows=$row_cnt, eig_len=$eig_len, size≈$(round(bytes/1e6, digits=2)) MB (Delta=$Delta, t_n=$t_n_vec, phason=$phason, mu=$mu_val)")

        fname = @sprintf("slice_Delta%.5f_t1%.5f_t2%.5f_phason%.5f_mu%.5f.bson",
                         Delta, t_n_tuple[1], t_n_tuple[2], phason, mu_val)
        fpath = joinpath(base_slice_dir, fname)
        BSON.bson(fpath, Dict(:df_slice => slice_df))
        combo_to_files[(Delta, t_n_tuple, phason, mu_val)] = fpath
    end

    df_base = nothing
    GC.gc()
    println("Stored $(length(combo_to_files)) sliced combos.")

    # Save combo_to_files for reuse in processing script
    BSON.bson(joinpath(base_slice_dir, "combo_to_files.bson"), Dict(:combo_to_files => combo_to_files))
    println("combo_to_files saved to $(joinpath(base_slice_dir, "combo_to_files.bson"))")

    return combo_to_files
end

# Run if executed as a script
if abspath(PROGRAM_FILE) == @__FILE__
    unpack_and_slice(folder_path_hof, base_slice_dir)
end

end  # IDOSUnpacking