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
data_folder = "sturmian_slopes_K100_L4_balanced_bins5000_mpb2_r50_comp-false_tailoption-0_N500_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-0p0-1)_Delta(0p0-0p1-2)"#"sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturmian_sweep_t1-t2_const_mapping"

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

unpack_data = false  # Set to true to unpack data, false to process existing slices
if unpack_data
    # Run unpacking and slicing
    combo_to_files = IDOSUnpacking.unpack_and_slice(folder_path_hof, base_slice_dir; 
        plt_mu_range_arg=repack_mu_range
    )
else
    # Processing existing slices
    combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
    combo_to_files = combo_data[:combo_to_files]
end



######################################
########### Extra variables ##########
######################################

slope_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/sturm_grad_sets/sturmian_slopes_K100_L4_balanced_bins5000_mpb2_r50_comp-false_tailoption-0.bson"
idos_plateau_threshold = 0.005
max_gap_label_value = 10

## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "idos", "mu_range$(repack_mu_range)", "plat_thresh$(idos_plateau_threshold)_maxGapLabel$(max_gap_label_value)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)


######################################
########## Plotting Execution ########
######################################

## Colour Colorscheme
function q_anchor_colormap(qs::Vector{Int})
    lerp(c1::RGB, c2::RGB, t) = RGB((1-t)*c1.r + t*c2.r,
                                    (1-t)*c1.g + t*c2.g,
                                    (1-t)*c1.b + t*c2.b)

    qs_sorted = sort(qs)
    cmap = fill(RGB(0.7,0.7,0.7), length(qs_sorted))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs_sorted)
    idx_pos1     = findfirst(==(  1), qs_sorted)
    idx_pos_high = length(qs_sorted)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    cmap[idx_neg_high] = col_neg_high
    cmap[idx_pos_high] = col_pos_high
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
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high-1)
            t = j / span
            if t < 0.5
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.5)
            else
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.5)/0.5)
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

    cmap
end

@showprogress for ((Delta, t_n_tuple, phason, mu_val), f) in combo_to_files
    if !isnothing(repack_mu_range) && !(mu_val in repack_mu_range)
        continue
    end
    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)

    df_slice = IDOSProcessing.compute_eigs_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_norm)
    df_slice = IDOSProcessing.compute_eigs_inner_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
    df_slice = IDOSProcessing.compute_idos_df!(df_slice)
    df_slice = IDOSProcessing.compute_plateaus_on_idos_df!(df_slice; threshold=idos_plateau_threshold)
    p_range = collect(-max_gap_label_value:max_gap_label_value)
    q_max = max_gap_label_value
    # df_slice = compute_gap_labels_qlim!(df_slice; p_range=p_range, q_max=q_max)
    try
        df_slice = IDOSProcessing.compute_gap_labels_qlim!(df_slice; p_range=p_range, q_max=q_max)
    catch e
        @warn "Gap label computation failed for mu=$mu_val, Delta=$Delta, t_n=$t_n, phason=$phason (likely due to Inf values): $e"
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


    # --- DEBUG: summarise slice contents before plotting ---
    @info "Slice summary" mu_val=mu_val Delta=Delta t_n=t_n phason=phason n_rows=nrow(gap_filtered_df)
    if nrow(gap_filtered_df) > 0
        # unique values actually present in this slice
        mus_present     = sort(unique(gap_filtered_df.mu))
        Deltas_present  = sort(unique(gap_filtered_df.Delta))
        phis_present    = sort(unique(gap_filtered_df.phi))
        phasons_present = :phason in names(gap_filtered_df) ? sort(unique(gap_filtered_df.phason)) : Float64[]
        ngaps_per_row   = [length(r.gap_labels) for r in eachrow(gap_filtered_df)]

        @info "Slice value ranges" mus_present=mus_present Deltas_present=Deltas_present phis_present=phis_present phasons_present=phasons_present
        @info "Gap stats" min_gaps_per_row = minimum(ngaps_per_row) max_gaps_per_row = maximum(ngaps_per_row)
    else
        @warn "gap_filtered_df is empty after filtering missing p/q" mu_val=mu_val Delta=Delta t_n=t_n phason=phason
    end
    # --- END DEBUG ---

    # ## 1) raw eigenvalues at fixed mu vs phi
    # try
    #     evl_outdir = joinpath(outdir, "eigenvalue_projections")
    #     isdir(evl_outdir) || mkpath(evl_outdir)
        
    #     IDOSPlotting.plt_eval_projections(
    #         df_slice,
    #         :eigenvalues,
    #         :phi,
    #         joinpath(evl_outdir, "eval_projections_comps_mu$(mu_val)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
    #         colour_rats=false,
    #         colour_comps=true,
    #         grad_filepath=slope_file,
    #         mu=mu_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         phason=phason
    #     )
    # catch e
    #     @warn "plt_eval_projections (eigenvalues) failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason error=e
    # end

    # ## 2) normalised eigenvalues at fixed mu vs phi
    # try
    #     norm_evl_outdir = joinpath(outdir, "norm_eigenvalue_projections")
    #     isdir(norm_evl_outdir) || mkpath(norm_evl_outdir)
        
    #     IDOSPlotting.plt_eval_projections(
    #         df_slice,
    #         :eigs_norm,
    #         :phi,
    #         joinpath(norm_evl_outdir, "norm_eval_projections_comps_mu$(mu_val)_Delta$(Delta)_tn$(t_n)_phason$(phason).png");
    #         colour_rats = false,
    #         colour_comps = true,
    #         # colour_all = true,
    #         grad_filepath = slope_file,
    #         mu = mu_val,
    #         Delta = Delta,
    #         t_n = t_n,
    #         phason = phason
    #     )
    # catch e
    #     @warn "plt_eval_projections (eigs_norm) failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason error=e
    # end

    # ## 3) IDOS plateaus vs phi coloured by gap labels
    # try
    #     IDOSPlotting.plt_plateaus_vs_phi_coloured_legend(
    #         df_slice,
    #         joinpath(outdir, "idos_plateaus_vs_phi_coloured_legend_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_thresh$(idos_plateau_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png");
    #         mu=mu_val,
    #         Delta=Delta,
    #         t_n=t_n,
    #         phason=phason
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


    ## 4) coloured gaps energy vs phi with different q optionality colour modes
        # colour_modes = [:qled_p, :qled_notp, :absq_p, :absq_notp]
        # colour_schemes = [:RdBu, :Reds, :RdBu, :Reds]

        unique_qs = sort(unique(g.q for g in Iterators.flatten(gap_filtered_df.gap_labels)))
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

            try 
                png_path = joinpath(
                    colour_plot_outdir,
                    "coloured_gaps_normalised_norats_energy_vs_phi_colormode_$(colour_mode)_mu$(mu_val)_Delta$(Delta)_t_n$(t_n)_phason$(phason)_thresh$(idos_plateau_threshold)_prange$(maximum(abs.(p_range)))_qmax$(q_max).png"
                )

                IDOSPlotting.plt_coloured_gaps_q_optionality(
                    rational_filtered_df, #gap_filtered_df,
                    png_path;
                    cmap=colour_schemes[i],
                    mu=mu_val,
                    Delta=Delta,
                    t_n=t_n,
                    phason=phason,
                    atol=1e-8,
                    rtol=1e-6,
                    verbose=true,
                    normalise_energy=true,
                    colour_mode=colour_mode, # :qled_p, :qled_notp, :absq_p, :absq_notp
                    line_width=1
                )
            catch e
                @warn "plt_coloured_gaps_q_optionality failed" mu_val=mu_val Delta=Delta t_n=t_n phason=phason colour_mode=colour_mode png_path=png_path error=e
            end
        end

    df_slice = nothing
    data = nothing
    # rm(f)
    GC.gc()
end

# rm(base_slice_dir, recursive=true)

println("\nAll IDOS processing complete. Results in: ", outdir)