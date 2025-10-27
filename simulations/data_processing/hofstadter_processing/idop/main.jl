include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")

data_folder = "Sturmian_K12_L6_balanced_bins500_mpb1_N(100-100-1)_t1(1.0-1.0-1_t2(1.5-1.5-1)_mu(0.0-3.0-301)_Delta(0.05-0.05-1)"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "np", "sturmian_sweep", data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)

df = unpack_bson_hofstadter(folder_path_hof; mp_targ=-1.0, mp_tol=1e-1, eval_targ=0.0, eval_tol=1e-2)
rename!(df, :sequence_id => :phi)
println("DataFrame keynames after renaming: $(names(df))")

N_range = sort(unique(df.N))
t_n_range = sort(unique(df.t_n))
Delta_range = sort(unique(df.Delta))
mu_range = sort(unique(df.mu))
phi_range = sort(unique(df.phi))
println("N_range: ", N_range)
println("t_n_range: ", t_n_range)
println("Delta_range: ", Delta_range)
println("mu_range (length, max, min): ", (length(mu_range), maximum(mu_range), minimum(mu_range)))
println("phi_range (length, max, min): ", (length(phi_range), maximum(phi_range), minimum(phi_range)))

N = N_range[1]
Delta = Delta_range[1]
t_n = t_n_range[1]

outdir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "sturmian_sweeps", data_folder))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)


plt_discrete_phase_projections(
    df,
    :mu,
    :phi,
    :disc_mp,
    joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi.png");
    # N = N,
    # Delta = Delta,
    # t_n = [t_n]
)

df_mp_idop = prep_df_for_IDOP(df; N=N, Delta=Delta, t_n=t_n)
df_mp_idop = compute_idop_df!(df_mp_idop; disc_variable=:disc_mp)

# ## Compute and check just a few phi values
# test_idx = 1:20:length(df_mp_idop.phi)
# for i in test_idx
#     phi = df_mp_idop.phi[i]

#     plateaus_mp = find_idop_plateaus_from_idop(
#         df_mp_idop;
#         delta_idop_thresh = 1e-4,
#         min_mu_span = 0.0,
#         min_samples = 1,
#         phi = phi
#     )
#     plt_idop_plateaus_check(
#         df_mp_idop, plateaus_mp, joinpath(outdir, "idop_plateaus_test/idop_plateaus_mp_phi$(phi).png");
#         N = N, Delta = 0.05, t_n = [1.0, 1.5], phi = phi
#     )
# end

df_mp_idop = compute_idop_plateaus_all!(df_mp_idop; delta_idop_thresh=1e-3, min_mu_span=0.0, min_samples=5)

plt_plateaus_vs_phi_idop(
    df_mp_idop,
    joinpath(outdir, "sturmian_sweep_idop_plateaus_vs_phi_N$(N)_Delta$(Delta)_t_n$(t_n).png")
)