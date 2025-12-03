using PrettyTables

include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")

using .IDOSPlotting
using .IDOSProcessing
using .BSONUnpacker


data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_algorithmic_generalised_N(200-200-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-0.0-1)_Delta(0.0-0.0-1)"#"sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturm_range_full" #"sturmian_sweep_t1-t2_const_mapping"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "np", run_set_name, data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)


# using Glob
# println("Checking folder: ", folder_path_hof)
# println("Does folder exist? ", isdir(folder_path_hof))
# bson_files = glob("*.bson", folder_path_hof)
# println("BSON files found: ", bson_files)
# if !isempty(bson_files)
#     # Check first file's keys
#     using BSON
#     raw = BSON.load(bson_files[1])
#     println("Keys in first BSON: ", collect(keys(raw)))
# end


df = BSONUnpacker.unpack_bson_hofstadter(folder_path_hof)
## Rename columns
# rename!(df, :sequence_id => :phi)
sequence_ids = collect(df.sequence_id)
df[!, :phi] = getindex.(df.sequence_id, 1)
df[!, :phason] = getindex.(df.sequence_id, 2)
df = IDOSProcessing.compute_eigs_norm!(df; eigs_col=:eigenvalues, out_col=:eigs_norm)
df = IDOSProcessing.compute_eigs_inner_norm!(df; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
println("DataFrame keynames after renaming: $(names(df))")

N_range = unique(df.N)
t_n_range = unique(df.t_n)
Delta_range = unique(df.Delta)
mu_range = unique(df.mu)
phi_range = unique(df.phi)
phason_range = unique(df.phason)
println("N_range: ", N_range)
println("t_n_range: ", t_n_range)
println("Delta_range: ", Delta_range)
println("mu_range (length, max, min): ", (length(mu_range), maximum(mu_range), minimum(mu_range)))
println("phi_range (length, max, min): ", (length(phi_range), maximum(phi_range), minimum(phi_range)))
println("phason_range (length, max, min): ", (length(phason_range), maximum(phason_range), minimum(phason_range)))

N = N_range[1]
Delta = Delta_range[1]
t_n = t_n_range[1]

outdir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "np", run_set_name, data_folder))#, "idos"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)

plt_tn_range = t_n_range
plt_mu_range = mu_range #[0.0, 1.0, 2.0, 3.0]
plt_Delta_range = Delta_range
plt_phason_range = phason_range #[0.0]

# slope_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/sturm_grad_sets/sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"
slope_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/sturm_seq_sets/sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_algorithmic_generalised.bson"

println("plt_tn_range: ", plt_tn_range)
println("plt_mu_range: ", plt_mu_range)
println("plt_Delta_range: ", plt_Delta_range)
println("plt_phason_range: ", plt_phason_range)


for t_n in plt_tn_range
    for phason in plt_phason_range
        for Delta in plt_Delta_range
            for mu in plt_mu_range
                IDOSPlotting.plt_eval_projections_fully_labelled(
                    df,
                    :eigenvalues,
                    :phi,
                    joinpath(outdir, "eval_projections_comps_mu$(mu)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
                    # colour_rats = false, #comment these in for using plt_eval_projections
                    # colour_comps = true,
                    colour_all = true,
                    transform_y = false,
                    grad_filepath = slope_file,
                    mu = mu,
                    # N = N,
                    Delta = Delta,
                    t_n = t_n,
                    phason = phason
                )

                # IDOSPlotting.plt_eval_projections(
                #     df,
                #     :eigs_norm,
                #     :phi,
                #     joinpath(outdir, "norm_eval_projections_comps_mu$(mu)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
                #     colour_rats = false,
                #     colour_comps = true,
                #     # colour_all = true,
                #     grad_filepath = slope_file,
                #     mu = mu,
                #     # N = N,
                #     Delta = Delta,
                #     t_n = t_n,
                #     phason = phason
                # )

                # IDOSPlotting.plt_eval_projections(
                #     df,
                #     :eigs_inner_norm,
                #     :phi,
                #     joinpath(outdir, "inner_norm_eval_projections_comps_mu$(mu)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
                #     colour_rats = false,
                #     colour_comps = true,
                #     # colour_all = true,
                #     grad_filepath = slope_file,
                #     mu = mu,
                #     # N = N,
                #     Delta = Delta,
                #     # t_n = t_n,
                #     phason = phason
                # )
            end
        end
    end
end

# df = compute_idos_df!(df)

# ################################################
# ## Gap-labelled plateau finding and plotting ###
# ################################################
# threshold = 0.05
# df = compute_plateaus_on_idos_df!(df; threshold=threshold)
# n = 5
# p_range = collect(-n:n)
# q_max = n
# df = compute_gap_labels_qlim!(df; p_range=p_range, q_max=q_max)

# # ## Optional plot single IDOS to check (not mu_range safe!)
# # row_idx = 4500
# # row = df[row_idx, :]
# # energies = sort(real(row.eigenvalues))
# # idos = row.idos
# # gap_labels = row.gap_labels

# # plot_idos_with_gaps(energies, idos, gap_labels)
# for phason in plt_phason_range
#     for Delta in plt_Delta_range
#         for mu in plt_mu_range
#             IDOSPlotting.plt_plateaus_vs_phi_coloured_legend(
#                 df, 
#                 joinpath(outdir, "idos_plateaus_vs_phi_coloured_legend_mu$(mu)_Delta$(Delta)_phason$(phason)_thresh$(threshold)_prange$(p_range)_qmax$(q_max).png"); 
#                 mu=mu,
#                 Delta=Delta,
#                 phason=phason
#             )
#         end
#     end
# end

# ##########################################
# ### Plotting gap-areas in energy space ###
# ##########################################

# # @show names(df)
# # @show typeof(df.gap_labels)
# # @show eltype(df.gap_labels)
# # if !isempty(df.gap_labels)
# #     @show typeof(first(df.gap_labels))
# #     @show first(df.gap_labels)
# # else
# #     println(":gap_labels column is empty")
# # end
# # pretty_table(first(df_gaps, min(10, nrow(df_gaps))))
# # println("Total number of gaps found: ", nrow(df_gaps))

# extreme_cmap = ["#000033", "#FFFF66", "#CC0000"]

# for phason in plt_phason_range
#     for Delta in plt_Delta_range
#         for mu in plt_mu_range

#             # IDOSPlotting.plt_coloured_gaps_energy_vs_phi(
#             #     df,
#             #     joinpath(outdir, "gaps_energy_vs_phi_coloured_legend_mu$(mu)_Delta$(Delta)_phason$(phason)_thresh$(threshold)_prange$(p_range)_qmax$(q_max).png");
#             #     cmap=:viridis, #extreme_cmap, # :viridis
#             #     mu=mu,
#             #     Delta=Delta,
#             #     phason=phason,
#             #     atol=1e-8,
#             #     rtol=1e-6,
#             #     verbose=false
#             # )

#             IDOSPlotting.plt_qled_coloured_gaps_energy_vs_phi(
#                 df,
#                 joinpath(outdir, "abs_qled_gaps_energy_vs_phi_coloured_legend_mu$(mu)_Delta$(Delta)_phason$(phason)_thresh$(threshold)_prange$(p_range)_qmax$(q_max).png");
#                 cmap=:RdBu, #extreme_cmap, # :viridis
#                 mu=mu,
#                 Delta=Delta,
#                 phason=phason,
#                 atol=1e-8,
#                 rtol=1e-6,
#                 verbose=false
#             )
#         end
#     end
# end