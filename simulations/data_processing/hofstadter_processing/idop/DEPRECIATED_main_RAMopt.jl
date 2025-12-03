using BSON
using DataFrames
include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")

data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturmian_sweep_t1-t2_const_mapping"
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

mp_targ = -1.0
mp_tol_range = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1]
eval_targ = 0.0
eval_tol = 1e-2

label_max = 5  # for gap labelling p,q ranges

println("Unpacking Hofstadter data from: ", folder_path_hof)
df_base = process_bson_files(folder_path_hof)
if isempty(df_base)
    error("No data unpacked — aborting.")
end

df_base = calc_norms_df!(df_base)
df_base = discretise_evals!(df_base, eval_targ, eval_tol)
df_base[!, :phi] = getindex.(df_base.sequence_id, 1)
df_base[!, :phason] = getindex.(df_base.sequence_id, 2)

base_slice_dir = normpath(joinpath(base_results_dir, "base_slices"))
isdir(base_slice_dir) || mkpath(base_slice_dir)

combo_to_files = Dict{Tuple{Float64, Float64, Float64}, Vector{String}}()
println("Caching base slices...")
grouped = DataFrames.groupby(df_base, [:N, :Delta, :t_n, :phason])
for g in grouped
    N = Float64(first(g.N))
    Delta = Float64(first(g.Delta))
    t_n = Float64(first(g.t_n))
    phason = Float64(first(g.phason))
    fname = @sprintf("slice_N%.0f_Delta%.4f_tn%.4f_phason%.4f.bson", N, Delta, t_n, phason)
    fpath = joinpath(base_slice_dir, fname)
    BSON.bson(fpath, Dict(:df_slice => g))
    push!(get!(combo_to_files, (Delta, t_n, phason), String[]), fpath)
end
df_base = nothing
GC.gc()
println("Stored $(length(combo_to_files)) base combos.")

for mp_tol in mp_tol_range
    println("\n=== Processing mp_tol = $(mp_tol) ===")
    outdir = normpath(joinpath(base_results_dir, "mp_tol$(mp_tol)"))
    isdir(outdir) || mkpath(outdir)

    # ---- Phase A: basic discretised phase projections per (Delta, t_n, phason) ----
    for ((Delta, t_n, phason), files) in combo_to_files
        df_combo = DataFrame()
        for f in files
            data = BSON.load(f)
            df_slice = data[:df_slice]
            df_slice = discretise_mp!(df_slice, mp_targ, mp_tol)
            append!(df_combo, df_slice)
            df_slice = nothing
            GC.gc()
        end
        if isempty(df_combo)
            @warn "No data for Delta=$Delta, t_n=$t_n, phason=$phason at mp_tol=$mp_tol (basic plot skipped)."
            continue
        end
        try
            plt_discrete_phase_projections(
                df_combo,
                :mu,
                :phi,
                :disc_mp,
                joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
                Delta=Delta,
                t_n=t_n,
                phason=phason
            )
        catch e
            @warn "Basic projection failed for Delta=$Delta, t_n=$t_n, phason=$phason at mp_tol=$mp_tol" error=e
        end
        df_combo = nothing
        GC.gc()
    end

    # ---- Phase B: IDOP, poly fits, gap labelling per (N, Delta, t_n, phason) slice ----
    for files in values(combo_to_files)
        for f in files
            data = BSON.load(f)
            df_slice = data[:df_slice]
            Delta = Float64(first(df_slice.Delta))
            t_n = Float64(first(df_slice.t_n))
            phason = Float64(first(df_slice.phason))
            N = Int(first(df_slice.N))
            try
                df_slice = discretise_mp!(df_slice, mp_targ, mp_tol)
                df_mp_preidop = prep_df_for_IDOP(df_slice; N=N, Delta=Delta, t_n=t_n, phason=phason)
                if isempty(df_mp_preidop)
                    @warn "Empty pre-IDOP slice for N=$N, Delta=$Delta, t_n=$t_n, phason=$phason at mp_tol=$mp_tol"
                    continue
                end

                poly_data = fit_final_mp_poly(df_mp_preidop, order=2, ignore_n_largest_phi=5)
                compute_phase_norm_with_poly_inplace!(
                    df_mp_preidop;
                    mu_col=:mu_values,
                    disc_col=:disc_mp,
                    phi_col=:phi,
                    mu_mp_col=:mu_mp,
                    mu_mp_norm_col=:mu_mp_norm,
                    disc_norm_col=:mp_disc_norm,
                    poly=poly_data.poly
                )

                plt_discrete_phase_projections(
                    df_slice,
                    :mu,
                    :phi,
                    :disc_mp,
                    joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason)_with_poly_N$(N).png");
                    final_poly_data=poly_data,
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    N=N
                )

                plt_discrete_phase_projections(
                    df_mp_preidop,
                    :mu_mp_norm,
                    :phi,
                    :mp_disc_norm,
                    joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason)_poly_normalised_N$(N).png");
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    N=N
                )

                df_idop_norm = compute_idop_df!(df_mp_preidop; disc_variable=:mp_disc_norm)
                df_idop_norm_plateaus = compute_idop_plateaus_all!(df_idop_norm; delta_idop_thresh=1e-4, min_mu_span=0.0, min_samples=1)
                filter_idop_plateaus_at_1!(df_idop_norm_plateaus; atol=1e-6)

                plt_plateaus_vs_phi_idop(
                    df_idop_norm_plateaus,
                    joinpath(outdir, "sturmian_sweep_idop_plateaus_vs_phi_N$(N)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol).png");
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    N=N
                )

                p_range = collect(-label_max:label_max)
                q_max = label_max
                df_idop_norm_plateaus = compute_gap_labels_qlim!(df_idop_norm_plateaus; p_range=p_range, q_max=q_max)

                plt_plateaus_vs_phi_coloured_legend(
                    df_idop_norm_plateaus,
                    joinpath(outdir, "idop_plateaus_vs_phi_gap_labeled_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol)_qmax$(q_max)_prange$(maximum(abs.(p_range)))_N$(N).png");
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    N=N
                )

                plt_qled_coloured_gaps_mu_vs_phi(
                    df_idop_norm_plateaus,
                    joinpath(outdir, "qled_gaps_mu_vs_phi_axismu_values_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol)_N$(N).png");
                    mu_axis=:mu_values,
                    cmap=:RdBu,
                    atol=1e-8,
                    rtol=1e-6,
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    N=N
                )
            catch e
                @warn "IDOP pipeline failed for N=$N, Delta=$Delta, t_n=$t_n, phason=$phason, mp_tol=$mp_tol" error=e
            end

            df_slice = nothing
            data = nothing
            GC.gc()
        end
    end
end

println("\nAll mp_tol runs complete. Base slices left on disk at: ", base_slice_dir)