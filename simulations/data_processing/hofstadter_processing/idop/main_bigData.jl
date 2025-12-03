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
include("RAM_opt_unpacker.jl")

using .IDOPUnpacking
using .IDOPPlotting
using .IDOPProcessing

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "sturmian_slopes_K100_L4_balanced_bins5000_mpb2_r50_comp-false_tailoption-0_N500_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(-3p0-3p0-601)_Delta(0p05-0p05-1)"#"sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "full_mu_range_negative_included"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define the unpacking mode: :unpack_by_phason, :unpack_by_slope, :load_phason, or :load_slope
# - :unpack_by_phason: Unpack and slice data with slices per phason (varying phi within slices)
# - :unpack_by_slope: Unpack and slice data with slices per phi (slope), containing range of phason
# - :load_phason / :load_slope: Load existing combo_to_files.bson without unpacking
# -------- IMPORTANT --------
unpack_mode = :load_slope  # Default to loading existing phason-based slices
partitioned_unpack = true  # Use partitioned unpacking to reduce RAM usage
n_partitions = 4
# ---------------------------

# Define ranges for repacking or plotting -- otherwise the full range will be used.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phason_range = nothing
repack_phi_range = nothing

# MP discretisation parameters
mp_targ = -1.0
mp_tol = 5e-2

# IDOP plateau parameters        
delta_idop_thresh = 1e-4
min_mu_span = 0.0
min_samples = 1

# Gap labelling parameters
n = 5
gap_label_p_range = collect(-n:n)
gap_label_q_max = n

# Determine base_slice_dir based on unpack_mode for consistent file locations
slice_suffix = if unpack_mode in [:unpack_by_phason, :load_phason]
    "_phason"
elseif unpack_mode in [:unpack_by_slope, :load_slope]
    "_slope"
else
    error("Invalid unpack_mode: $unpack_mode. Must be :unpack_by_phason, :unpack_by_slope, :load_phason, or :load_slope.")
end
base_slice_dir = normpath(joinpath(base_results_dir, "idop_base_slices", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)$(slice_suffix)"))

if unpack_mode == :unpack_by_phason
    if partitioned_unpack
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_phason_partitioned(folder_path_hof, base_slice_dir;
            n_partitions=n_partitions,
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range
        )
    else
        # Run unpacking and slicing per phason (fixed phi, varying phason)
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_phason(folder_path_hof, base_slice_dir;
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range
        )
    end
elseif unpack_mode == :unpack_by_slope
    # Run unpacking and slicing per phi (fixed phason, varying phi)
    if partitioned_unpack
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope_partitioned(folder_path_hof, base_slice_dir;
            n_partitions=n_partitions,
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phason_range_arg=repack_phason_range
    )
    else
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope(folder_path_hof, base_slice_dir;
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range
        )
    end
elseif unpack_mode in [:load_phason, :load_slope]
    # Load existing slices
    combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
    combo_to_files = combo_data[:combo_to_files]
    println("keys of comb_to_files.bson: ", keys(combo_data))
else
    error("Unexpected unpack_mode: $unpack_mode")
end


######################################
########### Extra variables ##########
######################################

## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "idop", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)"))
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

# Determine fixed param based on mode for consistent handling
if unpack_mode in [:unpack_by_phason, :load_phason]
    fixed_param = :phi
    fixed_param_name = "phi"
    varying_param = :phason
    varying_param_name = "phason"
    x_axis_param = varying_param  # For plotting: x = phason, y = mu_mp
else  # :unpack_by_slope or :load_slope
    fixed_param = :phason
    fixed_param_name = "phason"
    varying_param = :phi
    varying_param_name = "phi"
    x_axis_param = fixed_param  # For plotting: x = phason (fixed per slice), but accumulate across slices
end


@showprogress for ((N_val, Delta_val, t_n_tuple, fixed_val), f) in combo_to_files

    slice_suffix_plotname = "$(fixed_param_name)$(fixed_val)"

    # Apply filtering based on fixed_param
    if !isnothing(repack_N_range) && !(N_val in repack_N_range)
        continue
    end
    if !isnothing(repack_Delta_range) && !(Delta_val in repack_Delta_range)
        continue
    end
    if !isnothing(repack_tn_range) && !(t_n_tuple in repack_tn_range)
        continue
    end
    # Filter based on fixed_param
    range_to_check = if fixed_param == :phason
        repack_phason_range
    elseif fixed_param == :phi
        repack_phi_range
    else
        nothing
    end
    if !isnothing(range_to_check) && !(fixed_val in range_to_check)
        continue
    end

    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)

    # --- DEBUG: summarise slice contents ---
    @info "Slice summary" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val n_rows=nrow(df_slice)
    if nrow(df_slice) > 0
        # unique values actually present in this slice
        varying_present = sort(unique(df_slice[!, varying_param]))
        mu_mp_lengths = length.(df_slice.mu_mp)
        min_mu_mp_len = minimum(mu_mp_lengths)
        max_mu_mp_len = maximum(mu_mp_lengths)
        total_mu_points = sum(mu_mp_lengths)

        @info "Slice value ranges" $varying_param_name=varying_present min_mu_mp_len=min_mu_mp_len max_mu_mp_len=max_mu_mp_len total_mu_points=total_mu_points
    else
        @warn "df_slice is empty" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val
    end
    # --- END DEBUG ---

    # IDOP processing: Prepare the DataFrame for IDOP (expand mu_mp into vectors if needed, or use as-is)
    # Assuming df_slice has :phi and :mu_mp (vector of mu where disc_mp == 1)
    # For full IDOP, you might call prep_df_for_IDOP or other functions here
    # For now, keep it minimal

    # ## 1) Plot discrete phase projections
    # try
    #     phase_proj_outdir = joinpath(outdir, "phase_projections")
    #     isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)
        
    #     IDOPPlotting.plt_discrete_phase_projections(
    #         df_slice,
    #         :mu_mp,
    #         :phi,
    #         joinpath(phase_proj_outdir, "phase_proj_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$(slice_suffix).png");
    #         # N=N_val,
    #         # Delta=Delta_val,
    #         # t_n=t_n,
    #         # phason=phason_val
    #     )

    # catch e
    #     @warn "plt_discrete_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    # ## 2) Plot mu_mp vs phason for fixed N, Delta, t_n and slope (phi)
    # try
    #     mu_mp_phason_outdir = joinpath(outdir, "mu_mp_vs_phason")
    #     isdir(mu_mp_phason_outdir) || mkpath(mu_mp_phason_outdir)
        
    #     IDOPPlotting.plt_mu_mp_vs_phason(
    #         df_slice,
    #         joinpath(mu_mp_phason_outdir, "mu_mp_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # $fixed_param => fixed_val  # Filter by fixed param (phi for by_phason mode)
    #     )
    # catch e
    #     @warn "plt_mu_mp_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    # ## 3) Plot lowest eigs vs phason
    # try
    #     eigs_phason_outdir = joinpath(outdir, "lowest_eigs_vs_phason")
    #     isdir(eigs_phason_outdir) || mkpath(eigs_phason_outdir)
        
    #     IDOPPlotting.plt_lowest_eigs_vs_phason(
    #         df_slice,
    #         joinpath(eigs_phason_outdir, "lowest_eigs_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # $fixed_param => fixed_val
    #     )
    # catch e
    #     @warn "plt_lowest_eigs_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end


    # ## 4) Plot gap eigs vs phason
    # # Process gap eigenvalues
    # df_slice = IDOPProcessing.process_gap_eigenvalues(df_slice)

    # try
    #     gap_eigs_phason_outdir = joinpath(outdir, "gap_eigs_vs_phason")
    #     isdir(gap_eigs_phason_outdir) || mkpath(gap_eigs_phason_outdir)
        
    #     IDOPPlotting.plt_gap_eigs_vs_phason(
    #         df_slice,
    #         joinpath(gap_eigs_phason_outdir, "gap_eigs_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png")
    #     )
    # catch e
    #     @warn "plt_gap_eigs_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    ## 5) IDOP and gap labelling processing
    df_slice = IDOPProcessing.compute_idop_df!(df_slice; disc_variable=:mu_mp, rescale_idop=false)
    
    ## currently manually add mu_values column
    mu_values = collect(range(-3.0, 3.0, 601))
    df_slice[!, :mu_values] = [mu_values for _ in 1:nrow(df_slice)]

    df_slice = IDOPProcessing.compute_idop_plateaus_all!(
        df_slice; 
        delta_idop_thresh=delta_idop_thresh, 
        min_mu_span=min_mu_span, 
        min_samples=min_samples
    )

    IDOPProcessing.filter_idop_plateaus_at_1!(df_slice; atol=1e-6)
    IDOPProcessing.filter_idop_plateaus_at_0!(df_slice; atol=1e-6)

    df_slice = IDOPProcessing.compute_gap_labels_qlim_new!(
        df_slice;
        p_range=gap_label_p_range,
        q_max=gap_label_q_max
    )


    ## 6) Gap labelled plots
    
    # Define colourscheme
    unique_qs = sort(unique(g.q for g in Iterators.flatten(df_slice.gap_labels)))
    unique_qs = Int.(unique_qs)
    # println("typeof(unique_qs) = ", typeof(unique_qs), ", unique_qs = ", unique_qs)
    q_signed_custom = q_anchor_colormap(unique_qs)
    
    try
        coloured_gaps_outdir = joinpath(outdir, "coloured_gaps")
        isdir(coloured_gaps_outdir) || mkpath(coloured_gaps_outdir)

        colour_mode = :qled_notp  # Options: :qled_p, :qled_notp, :absq_p, :absq_notp

        IDOPPlotting.plt_coloured_gaps_q_optionality(
            df_slice,
            joinpath(coloured_gaps_outdir, "qled_gaps_mu_vs_phi_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
            # p_range=nothing,
            # q_max=nothing,
            colour_mode=colour_mode,
            cmap=q_signed_custom, #:RdBu,
            atol=1e-8,
            rtol=1e-6,
            verbose=false,
            # Delta=Delta,
            # phason=phason,
            # t_n=t_n
        )
    catch e
        @warn "plt_coloured_gaps_q_optionality failed"  N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    end

    df_slice = nothing
    data = nothing
    GC.gc()
end

println("\nAll IDOP processing complete. Results in: ", outdir)