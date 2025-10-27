include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")

data_folder = "Sturmian_K12_L6_balanced_bins500_mpb1_N(100-100-1)_t1(1.0-1.0-1_t2(1.5-1.5-1)_mu(0.0-3.0-301)_Delta(0.05-0.05-1)"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "np", "sturmian_sweep", data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)

df = unpack_bson_hofstadter(folder_path_hof)
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

plt_mu_range = [0.0, 1.0, 2.0, 3.0]
for mu in plt_mu_range
    plt_eval_projections(
        df,
        :phi,
        joinpath(outdir, "sturmian_sweep_eval_projections_mu$(mu).png");
    mu = mu,
    # N = N,
    # Delta = Delta,
    # t_n = [t_n]
    )
end

df_idos = compute_idos_df!(df)
df_idos_plateaus = compute_plateaus_df!(df_idos; threshold=0.1)

for mu in plt_mu_range
    plt_plateaus_vs_phi(
        df_idos_plateaus,
        joinpath(outdir, "sturmian_sweep_idos_plateaus_vs_phi_mu$(mu).png");
        mu = mu
    )
end

