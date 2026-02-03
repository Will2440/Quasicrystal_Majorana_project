include("plotting.jl")
include("processing.jl")
include("../../bson_unpacker.jl")

using .IDOPPlotting
using .IDOPProcessing


data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturmian_sweep_t1-t2_const_mapping"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)


# mp_tol = 5e-2
mp_tol_range = [5e-2] # [1e-3, 5e-3, 1e-2, 5e-2, 1e-1]

for mp_tol in mp_tol_range

    #############################################
    ########### Unpacking Data ##################
    #############################################

    println("Using mp_tol = ", mp_tol)

    df = unpack_bson_hofstadter(folder_path_hof; mp_targ=-1.0, mp_tol=mp_tol, eval_targ=0.0, eval_tol=1e-2)
    # rename!(df, :sequence_id => :phi)
    # println("DataFrame keynames after renaming: $(names(df))")
    sequence_ids = collect(df.sequence_id)
    df[!, :phi] = getindex.(df.sequence_id, 1)
    df[!, :phason] = getindex.(df.sequence_id, 2)
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

    outdir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder, "mp_tol$(mp_tol)"))
    isdir(outdir) || mkpath(outdir)
    println("Saving results to folder: ", outdir)

    #############################################
    ############# Basic Plotting ################
    #############################################
  
    plt_N_range = N_range
    plt_tn_range = t_n_range
    # plt_mu_range = mu_range #[1.0, 2.0, 3.0]
    plt_Delta_range = Delta_range
    plt_phason_range = phason_range #[0.0]


    # for t_n in plt_tn_range
    #     for phason in plt_phason_range
    #         for Delta in plt_Delta_range
    #             try
    #                 IDOPPlotting.smallData_plt_discrete_phase_projections(
    #                     df,
    #                     :mu,
    #                     :phi,
    #                     :disc_mp,
    #                     joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
    #                     Delta=Delta,
    #                     phason=phason,
    #                     t_n=t_n
    #                 )
    #             catch e
    #                 @warn "Error in basic plotting for Delta=$Delta, t_n=$t_n, phason=$phason: $e. Skipping."
    #             end
    #         end
    #     end
    # end


    ###############################################################
    ####### IDOP and gap-label processing and plotting ############
    ###############################################################
    for N in plt_N_range
        for Delta in plt_Delta_range
            for t_n in plt_tn_range
                for phason in plt_phason_range
                    try
                        df_mp_preidop = IDOPProcessing.prep_df_for_IDOP(df; N=N, Delta=Delta, t_n=t_n, phason=phason)
                        
                        if isempty(df_mp_preidop)
                            @warn "df_mp_preidop is empty for N=$N, Delta=$Delta, t_n=$t_n, phason=$phason. Skipping."
                            continue
                        end
                        
                        println("DataFrame keynames after prepping for IDOP: $(names(df_mp_preidop))")
                        
                        ## Fit-based normalisation approached and plotting
                        poly_data = IDOPProcessing.fit_final_mp_poly(df_mp_preidop, order=2, ignore_n_largest_phi=5)
                        println("Polyfit coefficients: ", poly_data.poly)
                        
                        # Normalise using the polynomial
                        df_mp_preidop_poly = IDOPProcessing.compute_phase_norm_with_poly!(
                            deepcopy(df_mp_preidop);
                            mu_col=:mu_values,
                            disc_col=:disc_mp,
                            phi_col=:phi,
                            mu_mp_col=:mu_mp,
                            mu_mp_norm_col=:mu_mp_norm,
                            disc_norm_col=:mp_disc_norm,
                            poly=poly_data.poly
                        )
                        
                        ######################################################
                        ###### Plotting phase projection normalisations ######
                        ######################################################
                        # Plot with the polynomial overlay
                        IDOPPlotting.smallData_plt_discrete_phase_projections(
                            df,
                            :mu,
                            :phi,
                            :disc_mp,
                            joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason)_with_poly.png");
                            final_poly_data=poly_data,
                            Delta=Delta,
                            phason=phason,
                            t_n=t_n
                        )
                        
                        # Plot the polynomial normalised data
                        IDOPPlotting.smallData_plt_discrete_phase_projections(
                            df_mp_preidop_poly,
                            :mu_mp_norm,
                            :phi,
                            :mp_disc_norm,
                            joinpath(outdir, "sturmian_sweep_disc_mp_mu_phi_mptol$(mp_tol)_Delta$(Delta)_tn$(t_n)_phason$(phason)_poly_normalised.png");
                            Delta=Delta,
                            phason=phason,
                            t_n=t_n
                        )
                        
                        # ######################################################
                        # ############# Calculate and check IDOP ###############
                        # ######################################################
                        # df_idop_norm = IDOPProcessing.compute_idop_df!(df_mp_preidop_poly; disc_variable=:mp_disc_norm)
                        # df_idop = IDOPProcessing.compute_idop_df!(df_mp_preidop; disc_variable=:disc_mp)
                        
                        # ######################################################
                        # ################# Find IDOP Plateaus #################
                        # ######################################################
                        # df_idop_norm_plateaus = IDOPProcessing.compute_idop_plateaus_all!(df_idop_norm; delta_idop_thresh=1e-4, min_mu_span=0.0, min_samples=1)
                        # IDOPProcessing.filter_idop_plateaus_at_1!(df_idop_norm_plateaus; atol=1e-6)
                        
                        # IDOPPlotting.plt_plateaus_vs_phi_idop(
                        #     df_idop_norm_plateaus,
                        #     joinpath(outdir, "sturmian_sweep_idop_plateaus_vs_phi_N$(N)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol).png");
                        #     Delta=Delta,
                        #     phason=phason,
                        #     t_n=t_n
                        # )
                        
                        # # Compute gap labels
                        # n = 5
                        # p_range = collect(-n:n)
                        # q_max = n
                        # df_idop_norm_plateaus = IDOPProcessing.compute_gap_labels_qlim!(df_idop_norm_plateaus; p_range=p_range, q_max=q_max)
                        
                        # # plot gap-labelled idop plateaus vs phi
                        # IDOPPlotting.plt_plateaus_vs_phi_coloured_legend(
                        #     df_idop_norm_plateaus,
                        #     joinpath(outdir, "idop_plateaus_vs_phi_gap_labeled_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol)_qmax$(q_max)_prange$(maximum(abs.(p_range))).png");
                        #     Delta=Delta,
                        #     phason=phason,
                        #     t_n=t_n
                        # )
                        
                        # ## plot gap-labelled mp gaps from idop plateaus (mu vs phi -- choose mu norm or not)
                        # mu_axis = :mu_values
                        
                        # IDOPPlotting.plt_qled_coloured_gaps_mu_vs_phi(
                        #     df_idop_norm_plateaus,
                        #     joinpath(outdir, "qled_gaps_mu_vs_phi_axis$(mu_axis)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_mptol$(mp_tol).png");
                        #     mu_axis=mu_axis,
                        #     cmap=:RdBu,
                        #     atol=1e-8,
                        #     rtol=1e-6,
                        #     verbose=false,
                        #     Delta=Delta,
                        #     phason=phason,
                        #     t_n=t_n
                        # )
                        
                    catch e
                        @warn "Error in IDOP processing for N=$N, Delta=$Delta, t_n=$t_n, phason=$phason: $e. Skipping this combination."
                    end
                end # end of phason loop
            end # end of Delta loop
        end # end of t_n loop
    end # end of N loop


end # end of mp_tol loop