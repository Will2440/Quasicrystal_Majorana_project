include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")


data_folder = "Sturmian_K12_L6_balanced_bins500_mpb1_N(200-200-1)_t1(1.0-1.0-1_t2(1.5-1.5-1)_mu(0.0-3.0-151)_Delta(0.05-0.05-1)"
run_set_name = "sturmian_sweep"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "np", run_set_name, data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)


mp_tol = 1e-2
df = unpack_bson_hofstadter(folder_path_hof; mp_targ=-1.0, mp_tol=mp_tol, eval_targ=0.0, eval_tol=1e-2)
rename!(df, :sequence_id => :phi)
println("DataFrame keynames after renaming: $(names(df))")
# sequence_ids = collect(df.sequence_id)
# df[!, :phi] = getindex.(df.sequence_id, 1)
# df[!, :phason] = getindex.(df.sequence_id, 2)
# println("DataFrame keynames after renaming: $(names(df))")

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

plt_tn_range = t_n_range
plt_mu_range = mu_range #[1.0, 2.0, 3.0]
plt_Delta_range = Delta_range
# plt_phason_range = phason_range #[0.0]

# # for t_n in plt_tn_range
# # for phason in plt_phason_range
# for Delta in plt_Delta_range
#     plt_discrete_phase_projections(
#         df,
#         :mu,
#         :phi,
#         :disc_mp,
#         joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_.png");#phason$(phason).png");
#         # N = N,
#         Delta = Delta,
#         # phason = phason,
#         # t_n = [t_n]
#     )
# end
# # end
# # end

df_mp_preidop = prep_df_for_IDOP(df; N=N, Delta=Delta, t_n=t_n)

## Brute force normalisation approaches which do not account for missing edge values:
# df_mp_preidop = compute_phase_norm!(df_mp_preidop; mu_col=:mu_values, disc_col=:disc_mp, mu_norm_col=:mu_values_norm, disc_norm_col=:mp_disc_norm)
# df_mp_preidop = compute_phase_norm!(df_mp_preidop; mu_col=:mu_values, disc_col=:disc_mp, mu_mp_col=:mu_mp, mu_mp_norm_col=:mu_mp_norm)


## Fit-based normalisation approached and plotting
# Fit the regression line or polynomial
line_data = fit_final_mp_line(df_mp_preidop, ignore_n_largest_phi=5)
poly_data = fit_final_mp_poly(df_mp_preidop, order=2, ignore_n_largest_phi=5)
println("Polyfit coefficients: ", poly_data.poly)

# Normalise using the polynomial
df_mp_preidop_poly = compute_phase_norm_with_poly!(
    deepcopy(df_mp_preidop);
    mu_col=:mu_values,
    disc_col=:disc_mp,
    phi_col=:phi,
    mu_mp_col=:mu_mp,
    mu_mp_norm_col=:mu_mp_norm,
    disc_norm_col=:mp_disc_norm,
    poly=poly_data.poly
)

# Normalise using the regression line
df_mp_preidop_line = compute_phase_norm_with_line!(
    deepcopy(df_mp_preidop);
    mu_col=:mu_values,
    disc_col=:disc_mp,
    phi_col=:phi,
    mu_mp_col=:mu_mp,
    mu_mp_norm_col=:mu_mp_norm,
    disc_norm_col=:mp_disc_norm,
    line_coeffs=line_data
)

# # Plot with the regression line overlay
# plt_discrete_phase_projections(
#     df,
#     :mu,
#     :phi,
#     :disc_mp,
#     joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_with_line.png");#phason$(phason).png");
#     final_line_data=line_data
# )

# # Plot with the polynomial overlay
# plt_discrete_phase_projections(
#     df,
#     :mu,
#     :phi,
#     :disc_mp,
#     joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_with_poly.png");#phason$(phason).png");
#     final_poly_data=poly_data
# )

# # Plot the linear normalised data
# plt_discrete_phase_projections(
#     df_mp_preidop_line,
#     :mu_mp_norm,
#     :phi,
#     :mp_disc_norm,
#     joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_line_normalised.png");#phason$(phason).png");
# )

# # Plot the polynomial normalised data
# plt_discrete_phase_projections(
#     df_mp_preidop_poly,
#     :mu_mp_norm,
#     :phi,
#     :mp_disc_norm,
#     joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_poly_normalised.png");#phason$(phason).png");
# )


df_idop_norm = compute_idop_df!(df_mp_preidop_poly; disc_variable=:mp_disc_norm)
df_idop = compute_idop_df!(df_mp_preidop_poly; disc_variable=:disc_mp)

# ## Compute and check just a few phi values
# sub_outdir = joinpath(outdir, "idop_plateaus_test")
# isdir(sub_outdir) || mkpath(sub_outdir)
# test_idx = 1:20:length(df_idop_norm.phi)
# for i in test_idx
#     phi = df_idop_norm.phi[i]

#     plateaus_mp = find_idop_plateaus_from_idop(
#         df_idop_norm;
#         delta_idop_thresh = 1e-4,
#         min_mu_span = 0.0,
#         min_samples = 1,
#         phi = phi
#     )
#     plt_idop_plateaus_check(
#         df_idop_norm, plateaus_mp, joinpath(sub_outdir, "idop_plateaus_mp_phi$(phi)_mptol$(mp_tol).png");
#         # N = N, 
#         # Delta = 0.05, 
#         # t_n = [1.0, 1.5], 
#         phi = phi
#     )
# end

df_idop_norm_plateaus = compute_idop_plateaus_all!(df_idop_norm; delta_idop_thresh=1e-4, min_mu_span=0.0, min_samples=1)

plt_plateaus_vs_phi_idop(
    df_idop_norm_plateaus,
    joinpath(outdir, "sturmian_sweep_idop_plateaus_vs_phi_N$(N)_Delta$(Delta)_t_n$(t_n)_mptol$(mp_tol).png")
)