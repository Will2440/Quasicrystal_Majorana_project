using CSV
using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter

# Set logger level from JULIA_LOG_LEVEL environment variable, default to Info if not set
level_str = uppercase(get(ENV, "JULIA_LOG_LEVEL", "info"))
level = if level_str == "DEBUG"
    Logging.Debug
elseif level_str == "INFO"
    Logging.Info
elseif level_str == "WARN"
    Logging.Warn
elseif level_str == "ERROR"
    Logging.Error
else
    Logging.Info  # fallback
end
global_logger(ConsoleLogger(stderr, level))

include("plotting.jl")
include("processing.jl")
include("RAM_opt_unpack.jl")

using .IDOSUnpacking
using .IDOSPlotting
using .IDOSProcessing

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "hof_style_slopes_N1000_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)"
run_set_name = "MBData_final"
# data_folder = "hof_style_slopes_N500_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-5p0-7)_mu(0p0-0p0-1)_Delta(0p05-0p5-6)"
# run_set_name = "hof_style_slopes" 
# data_folder = "hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"
# run_set_name = "QC-SC_comp_small_L_check"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))


######################################
############### Unpack? ##############
######################################

# Define ranges for repacking or plotting -- otherwise the full range will be.
repack_mu_range = nothing #[0.0, 0.5, 1.0, 2.0]

## Defines base_slice_dir for accessing slices using this script
base_slice_dir = normpath(joinpath(base_results_dir, "idos_base_slices", "mu_range$(repack_mu_range)"))

unpack_mode = :by_phason #:by_slope or :by_phason

unpack_data = true  # Set to true to unpack data, false to process existing slices
fixed_param = nothing
if unpack_data
    # Run unpacking and slicing
    if unpack_mode == :by_phason
        combo_to_files = IDOSUnpacking.unpack_and_slice_by_phason(
            folder_path_hof, 
            base_slice_dir; 
            plt_mu_range_arg=repack_mu_range
        )
        fixed_param = :phason
    elseif unpack_mode == :by_slope
        combo_to_files = IDOSUnpacking.unpack_and_slice_by_slope(
            folder_path_hof, 
            base_slice_dir; 
            plt_mu_range_arg=repack_mu_range
        )
        fixed_param = :slope
    else
        error("Invalid unpack_mode: $unpack_mode. Use :by_slope or :by_phason.")
    end
else
    # Processing existing slices
    combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
    combo_to_files = combo_data[:combo_to_files]
end

## Dictionary to hold phason-independent spectra
phason_indep_spectra = Dict{NamedTuple, Vector{Float64}}()


######################################
########### Extra variables ##########
######################################

slope_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/sturm_grad_sets/sturmian_slopes_K100_L4_balanced_bins5000_mpb2_r50_comp-false_tailoption-0.bson"
idos_plateau_width_threshold = 0.001
N = 500
idos_plateau_height_threshold = 2/N
max_gap_label_value = 20
mp_tol = 0.1

normalise_type = :outer # :inner, :outer, :both
plot_norm_or_not = :eigenvalues # :eigenvalues or :eigs_interp_norm

## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "idos", "mu_range$(repack_mu_range)", "plat_thresh$(idos_plateau_width_threshold),$(idos_plateau_height_threshold)_maxGapLabel$(max_gap_label_value)_mptol$(mp_tol)", "norm_type$(normalise_type)_plotNorm$(plot_norm_or_not)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)


######################################
########## Plotting Execution ########
######################################

function q_anchor_colormap(qs::Vector{Int})
    lerp(c1::RGB, c2::RGB, t) = RGB((1-t)*c1.r + t*c2.r,
                                    (1-t)*c1.g + t*c2.g,
                                    (1-t)*c1.b + t*c2.b)

    qs_sorted = sort(qs)
    
    # --- SENTINEL HANDLING ---
    SC_GAP_Q = 9999
    has_sentinel = SC_GAP_Q in qs_sorted
    
    # Filter out sentinel for gradient calculation
    qs_phys = filter(!=(SC_GAP_Q), qs_sorted)
    # -------------------------

    cmap = fill(RGB(0.7,0.7,0.7), length(qs_phys))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs_phys)
    idx_pos1     = findfirst(==(  1), qs_phys)
    idx_pos_high = length(qs_phys)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    cmap[idx_neg_high] = col_neg_high
    
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end

    if !isnothing(idx_neg1) && idx_neg1 > idx_neg_high
        span = idx_neg1 - idx_neg_high
        for (j,i) in enumerate(idx_neg_high+1:idx_neg1-1)
            t = j / span
            if t < 1/3
                cmap[i] = lerp(col_neg_high, col_neg_purp, t/(1/3))
            elseif t < 2/3
                cmap[i] = lerp(col_neg_purp, col_neg_dark, (t-1/3)/(1/3))
            else
                cmap[i] = lerp(col_neg_dark, col_neg_light, (t-2/3)/(1/3))
            end
        end
    end

    if !isnothing(idx_pos1) && idx_pos_high > idx_pos1
        span = idx_pos_high - idx_pos1
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high)
            t = j / span
            if t < 0.33
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.33)
            elseif t < 0.66
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.33)/0.33)
            else
                cmap[i] = lerp(col_pos_yell, col_pos_high, (t-0.66)/0.34)
            end
        end
    end

    if !isnothing(idx_neg1) && !isnothing(idx_pos1) && idx_pos1 - idx_neg1 > 1
        span = idx_pos1 - idx_neg1
        for (j,i) in enumerate(idx_neg1+1:idx_pos1-1)
            t = j / span
            cmap[i] = lerp(col_neg_light, RGB(1.0,0.0,0.0), t)
        end
    end

    # --- RECONSTRUCT MAP WITH SENTINEL ---
    # Map q -> color
    color_dict = Dict(q => c for (q,c) in zip(qs_phys, cmap))
    
    if has_sentinel
        color_dict[SC_GAP_Q] = RGB(0.85, 0.85, 0.85) # Light Grey for SC gaps
    end
    
    # Return colors in the order of the original sorted input
    return [color_dict[q] for q in qs_sorted]
end

plt_t_n = [(1.0, 2.0)]
plt_Delta = 0.05
plt_slope = nothing #0.6180048661800487
# plt_phason = [0.0]

@showprogress for ((Delta, t_n_tuple, param_val, mu_val), f) in combo_to_files
    if fixed_param == :slope
        phi_val = param_val
        phason_val = "ALL" # or specific value if extracted from rows later, but "ALL" implies integration
        loop_fixed_desc = "phi=$(phi_val)"
    elseif fixed_param == :phason
        phason_val = param_val
        phi_val = "ALL"
        loop_fixed_desc = "phason=$(phason_val)"
    else
        # Fallback for safety
        phi_val = param_val 
        phason_val = 0.0 #"UNKNOWN"
        loop_fixed_desc = "param=$(param_val)"
    end

    if !isnothing(plt_slope) && !(abs(phi_val - plt_slope) < 1e-4)
        continue
    end

    # if !isnothing(plt_phason) && !(phason_val in plt_phason)
    #     continue
    # end

    if !isnothing(repack_mu_range) && !(mu_val in repack_mu_range)
        continue
    end

    if !isnothing(plt_t_n) && !(t_n_tuple in plt_t_n)
        continue
    end

    if !isnothing(plt_Delta) && !(abs(Delta - plt_Delta) < 1e-10)
        continue
    end

    p_range = collect(-max_gap_label_value:max_gap_label_value)
    q_max = max_gap_label_value

    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)


    ##########################################
    ############ Processing Section ##########
    ##########################################

    # ## 1:eigenvalues,
    #     poly_deg=4, # Adjust degree if needed (2 for parabola, 4 for more complex dome)
    #     coeffs_out_file=coeffs_out_file
    # )


    # ## 2) Compute phason-independent spectrum
    # tol=0.01
    # grid_spacing=0.003 #(approx bandwidth/2N)
    # phason_indep_spec = IDOSProcessing.compute_phason_independent_spectrum_gridTol(
    #     df_slice; 
    #     eigs_col=:eigenvalues,
    #     tolerance=tol,
    #     grid_spacing=grid_spacing,
    #     verbose=false
    # )

    # # Extract metadata for key to save spectrum
    # try
    #     phi_val = first(df_slice.phi)
        
    #     row1 = first(df_slice)
    #     N_val = hasproperty(row1, :N) ? row1.N : length(row1.eigenvalues)
        
    #     spec_key = (slope=phi_val, N=N_val, Delta=Delta, t_n=t_n_tuple, mu=mu_val)
    #     phason_indep_spectra[spec_key] = phason_indep_spec
    # catch e
    #     @warn "Could not extract metadata for phason_indep_spectra key" error=e
    # end


    # ## plot old vs new spectra for comparison
    # try
    #     spec_outdir = joinpath(outdir, "phason_indep_spectra_comparison", "tolerance_$(tol)_gridSpacing_$(grid_spacing)")
    #     isdir(spec_outdir) || mkpath(spec_outdir)

    #     raw_specs_col = :eigenvalues
    #     IDOSPlotting.plt_phason_indep_spectrum_comparison(
    #         df_slice,
    #         raw_specs_col,
    #         phason_indep_spec,
    #         joinpath(spec_outdir, "phason_indep_spectrum_comparison_$(loop_fixed_desc)_Delta$(Delta)_t1$(t_n[1])_t2$(t_n[2])_mu$(mu_val).png");
    #         slope=phi_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         mu=mu_val
    #     )
    # catch e
    #     @warn "plt_phason_indep_spectrum_comparison failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val error=e
    # end


    ## 3) Compute IDOS and plateaus
    df_slice = IDOSProcessing.compute_idos_df!(df_slice)
    df_slice = IDOSProcessing.compute_plateaus_on_idos_df_extra_thresh!(
        df_slice; 
        threshold=idos_plateau_width_threshold, 
        eigs_col=plot_norm_or_not
    )

    ## 4) Compute gap labels
    try
        # df_slice = IDOSProcessing.compute_gap_labels_qlim!(df_slice; p_range=p_range, q_max=q_max)
        df_slice = IDOSProcessing.compute_gap_labels_parsimonious!(
            df_slice; 
            p_range=p_range, 
            q_max=q_max, 
            tol=1.5/N,
            merge_threshold=idos_plateau_height_threshold,
            check_mp=true,
            mp_col=:mp,
            mp_threshold=mp_tol,
            always_mask_central=true
        )
    catch e
        @warn "Gap label computation failed for mu=$mu_val, Delta=$Delta, t_n=$t_n, phason=$phason_val (likely due to Inf values): $e"
        # Skip to next iteration or continue without gap labels
        df_slice = nothing
        data = nothing
        GC.gc()
        continue
    end

    ## filter out missing gap label vals
    gap_filtered_df = filter(row -> all(label -> !ismissing(label.p) && !ismissing(label.q), row.gap_labels), df_slice)

    ## filter out rationals for some plotting
    n=11 # max denominator
    rationals = Set{Float64}()
    for r in 1:n
        push!(rationals, 1.0/r)
        push!(rationals, 1.0-1.0/r)
    end
    rational_filtered_df = filter(row -> !any(isapprox(row.phi, rat; atol=1e-9) for rat in rationals), gap_filtered_df)


    # # --- DEBUG: summarise slice contents before plotting ---
    # @info "Slice summary" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val n_rows=nrow(gap_filtered_df)
    # if nrow(gap_filtered_df) > 0
    #     # unique values actually present in this slice
    #     mus_present     = sort(unique(gap_filtered_df.mu))
    #     Deltas_present  = sort(unique(gap_filtered_df.Delta))
    #     phis_present    = sort(unique(gap_filtered_df.phi))
    #     phasons_present = :phason in names(gap_filtered_df) ? sort(unique(gap_filtered_df.phason)) : Float64[]
    #     ngaps_per_row   = [length(r.gap_labels) for r in eachrow(gap_filtered_df)]

    #     @info "Slice value ranges" mus_present=mus_present Deltas_present=Deltas_present phis_present=phis_present phasons_present=phasons_present
    #     @info "Gap stats" min_gaps_per_row = minimum(ngaps_per_row) max_gaps_per_row = maximum(ngaps_per_row)
    # else
    #     @warn "ga) Eigenvalue normalisation

        # a) Direct normalisation by bandwidth
    # df_slice = IDOSProcessing.compute_eigs_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_norm)
    # df_slice = IDOSProcessing.compute_eigs_inner_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
    df_slice = IDOSProcessing.compute_eigs_bandwidth_norm!(
        df_slice; 
        eigs_col=:eigenvalues, 
        out_col=:eigs_outer_norm,
        min_out_col=:global_min,
        max_out_col=:global_max
    )
    
    # # Create norm_save_df: Keep only parameter columns + global_min/global_max, remove all others
    # parameter_cols = [:phi, :phason, :mu, :Delta, :t_n, :N]  # Core parameter columns (adjust if your df_slice has different names)
    # norm_save_df = select(df_slice, parameter_cols..., :global_min, :global_max)
    # norm_save_path = joinpath(outdir, "eigs_bandwidth_norms_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_$(loop_fixed_desc).csv")
    # CSV.write(norm_save_path, norm_save_df)

    
        ## b) By linear interpolation
    # -- OLD Linear Interpolation --
    # df_slice = IDOSProcessing.compute_interpolated_normalisation!(
    #     df_slice; 
    #     out_col=:eigs_interp_norm, 
    #     mode=normalise_type, 
    #     ref_indices=(5, 1)
    # )

    #     ## c) By polynomial (saving coefficients)
    # coeffs_out_file = joinpath(outdir, "smooth_norm_coeffs_Delta$(Delta)_t1$(t_n[1])_t2$(t_n[2])_$(loop_fixed_desc)_thresh$(idos_plateau_width_threshold).txt")

    # # -- NEW Smooth Polynomial Envelope --
    # df_slice = IDOSProcessing.compute_smooth_normalisation!(
    #     df_slice;
    #     out_col=:eigs_interp_norm, # reusing this name for compatibility with downstream plotting
    #     eigs_col=p_filtered_df is empty after filtering missing p/q" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val
    # end
    # # --- END DEBUG ---


    ##########################################
    ############# Plotting Section ###########
    ##########################################

    # ## 1) Plot raw eigenvalues at fixed mu vs phi
    # try
    #     evl_outdir = joinpath(outdir, "eigenvalue_projections")
    #     isdir(evl_outdir) || mkpath(evl_outdir)
        
    #     IDOSPlotting.plt_eval_projections(
    #         df_slice,
    #         :eigenvalues,
    #         :phi,
    #         joinpath(evl_outdir, "eval_projections_comps_mu$(mu_val)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
    #         colour_rats=false,
    #         colour_comps=false,
    #         grad_filepath=nothing, #slope_file,
    #         mu=mu_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         phason=phason_val
    #     )
    # catch e
    #     @warn "plt_eval_projections (eigenvalues) failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val error=e
    # end

    # ## 2) Plot normalised eigenvalues at fixed mu vs phi
    # try
    #     norm_evl_outdir = joinpath(outdir, "norm_eigenvalue_projections")
    #     isdir(norm_evl_outdir) || mkpath(norm_evl_outdir)
        
    #     IDOSPlotting.plt_eval_projections(
    #         df_slice,
    #         :eigs_interp_norm, # Use :eigs_norm with old snap normalisation, use :eigs_inner_norm with new normalisation method
    #         :phi,
    #         joinpath(norm_evl_outdir, "norm_eval_projections_comps_mu$(mu_val)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
    #         colour_rats = false,
    #         colour_comps = false,
    #         # colour_all = true,
    #         grad_filepath = nothing, #slope_file,
    #         mu = mu_val,
    #         Delta = Delta,
    #         t_n = t_n,
    #         phason = phason_val
    #     )
    # catch e
    #     @warn "plt_eval_projections (eigs_norm) failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val error=e
    # end

    # ## 3) Plot IDOS plateaus vs phi coloured by gap labels
    # try
    #     IDOSPlotting.plt_plateaus_vs_phi_coloured_legend(
    #         df_slice,
    #         joinpath(outdir, "idos_plateaus_vs_phi_coloured_legend_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason_val)_thresh$(idos_plateau_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png");
    #         mu=mu_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         phason=phason_val
    #     )
    # catch e
    #     @warn "plt_plateaus_vs_phi_coloured_legend failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason error=e
    # end


        # IDOSPlotting.plt_qled_coloured_gaps_energy_vs_phi(
        #     gap_filtered_df,
        #     joinpath(outdir, "abs_qled_gaps_energy_vs_phi_coloured_legend_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_thresh$(idos_plateau_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png");
        #     cmap=:RdBu,
        #     mu=mu_val,
        #     Delta=Delta,
        #     t_n=t_n,
        #     phason=phason,
        #     atol=1e-8,
        #     rtol=1e-6,
        #     verbose=false
        # )


    ## 4) Plot coloured gaps energy vs phi with different q optionality colour modes
        # colour_modes = [:qled_p, :qled_notp, :absq_p, :absq_notp]
        # colour_schemes = [:RdBu, :Reds, :RdBu, :Reds]

        unique_qs = sort(unique(g.q for g in Iterators.flatten(rational_filtered_df.gap_labels)))
        q_signed_custom = q_anchor_colormap(unique_qs)

        colour_modes = [:qled_notp]
        colour_schemes = [q_signed_custom]
        for (i, colour_mode) in enumerate(colour_modes)

            colour_plot_outdir = joinpath(outdir, "q_colouring", "colourmode_$(colour_mode)")
            isdir(colour_plot_outdir) || mkpath(colour_plot_outdir)

            # try 
            #     png_path = joinpath(
            #         colour_plot_outdir,
            #         "coloured_gaps_norats_energy_vs_phi_colormode_$(colour_mode)_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_thresh$(idos_plateau_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png"
            #     )

            #     IDOSPlotting.plt_coloured_gaps_q_optionality(
            #         rational_filtered_df, #gap_filtered_df,
            #         png_path;
            #         cmap=colour_schemes[i],
            #         mu=mu_val,
            #         Delta=Delta,
            #         t_n=t_n,
            #         phason=phason,
            #         atol=1e-8,
            #         rtol=1e-6,
            #         verbose=true,
            #         normalise_energy=false,
            #         colour_mode=colour_mode # :qled_p, :qled_notp, :absq_p, :absq_notp
            #     )
            # catch e
            #     @warn "plt_coloured_gaps_q_optionality failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason colour_mode=colour_mode png_path=png_path error=e
            # end

            png_path = joinpath(
                colour_plot_outdir,
                "coloured_gaps_normalised_energy_vs_phi_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason_val)_thresh$(idos_plateau_width_threshold),$(idos_plateau_height_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png"
            )

            try 
                IDOSPlotting.plt_coloured_gaps_q_optionality(
                    rational_filtered_df, #gap_filtered_df,
                    png_path;
                    cmap=colour_schemes[i],
                    mu=mu_val,
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason_val,
                    atol=1e-8,
                    rtol=1e-6,
                    verbose=true,
                    normalise_energy=false,
                    colour_mode=colour_mode, # :qled_p, :qled_notp, :absq_p, :absq_notp
                    line_width=1,
                    plt_xlims=nothing, #(-1.0, 1.0),
                    plt_ylims=(0.0, 1.0)
                )
            catch e
                @warn "plt_coloured_gaps_q_optionality failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason_val colour_mode=colour_mode png_path=png_path error=e
            end
        end

    # ## 5) Plot eigenvalues verus phason at fixed phi (and mu, Delta, t_n)
    # try
    #     evl_phason_outdir = joinpath(outdir, "eigenvalue_projections_vs_phason")
    #     isdir(evl_phason_outdir) || mkpath(evl_phason_outdir)
    #     IDOSPlotting.plt_eval_projections_vs_phason(
    #         df_slice,
    #         :eigenvalues,
    #         :phason,
    #         joinpath(evl_phason_outdir, "LINE_eval_projections_vs_phason_comps_mu$(mu_val)_Delta$(Delta)_tn$(t_n)_phi$(phi_val).png");
    #         mu=mu_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         phi=phi_val,
    #         lims=:half,
    #         reverse_axes=true,
    #         plt_type=:line,
    #         size=(1000,400)
    #     )
    # catch e
    #     @warn "plt_eval_projections_vs_phason (eigenvalues) failed" mu_val=mu_val Delta=Delta t_n=t_n phi=phi_val error=e
    # end
    

    df_slice = nothing
    data = nothing
    # rm(f)
    GC.gc()
end


# Save the phason independent spectra
spec_savepath = joinpath(outdir, "phason_independent_spectra.bson")
BSON.bson(spec_savepath, Dict(:phason_indep_spectra => phason_indep_spectra))
println("Saved phason independent spectra to: ", spec_savepath)



# Save as CSV for easier access (convert Dict to DataFrame)
df_spectra = DataFrame(
    key = collect(keys(phason_indep_spectra)),
    spectrum = collect(values(phason_indep_spectra))
)
csv_savepath = joinpath(outdir, "phason_independent_spectra.csv")
CSV.write(csv_savepath, df_spectra)
println("Also saved as CSV to: ", csv_savepath)


println("\nAll IDOS processing complete. Results in: ", outdir)