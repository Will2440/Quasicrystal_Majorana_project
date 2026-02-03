using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter

## data variable for found mp_tol ranges
struct MPTolResult
    N::Int64
    Delta::Float64
    t_n::Tuple{Float64, Float64}
    phi::Float64
    phason::Float64
    qc_region_count::Int64
    mp_tol_ranges::Union{Nothing, Tuple{Float64, Float64}}
end

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
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)"
run_set_name = "MBData_final"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define the unpacking mode: :unpack_by_slope or :load_slope
# - :unpack_by_slope: Unpack and slice data with slices per phi (slope), containing range of phason
# - :load_slope: Load existing combo_to_files.bson without unpacking
# -------- IMPORTANT --------
unpack_mode = :load_slope
# ---------------------------

eig_save_mode = :all        # Options: :all (full spectrum), :maj (central 2 only)
mirror_data_about_mu0 = true  # Whether to mirror data about mu=0 for unpacking

# Define ranges for repacking or plotting -- otherwise the full range will be used.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phason_range = nothing
repack_phi_range = nothing #[0.6180048661800487]

# MP discretisation parameters
mp_tol_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_tol_sweeps/mptarg-1.0_mptolrange1.0e-10-1.0-1001_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_rangenothing/found_mp_tol_ranges.bson"
mp_tol_df = IDOPProcessing.extract_mp_tol(mp_tol_filepath)
mp_targ = -1.0

slice_suffix = "_slope"
base_slice_dir = normpath(joinpath(base_results_dir, "idop_base_slices", "mptarg$(mp_targ)_mptolVariable_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)$(slice_suffix)_muduplicated-$(mirror_data_about_mu0)"))


if unpack_mode == :unpack_by_slope
    combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope_raw_mp(
        folder_path_hof, 
        base_slice_dir;
        plt_N_range_arg=repack_N_range,
        plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phason_range_arg=repack_phason_range,
        mirror_data_about_mu0=mirror_data_about_mu0
    )
elseif unpack_mode == :load_slope
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
outdir = normpath(joinpath(base_results_dir, "idop", "mptarg$(mp_targ)_mptolVariable_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)_mpduplicated-$(mirror_data_about_mu0)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)


######################################
########## Plotting Execution ########
######################################


# Determine fixed param based on mode for consistent handling
fixed_param = :phason
fixed_param_name = "phason"
varying_param = :phi
varying_param_name = "phi"
x_axis_param = fixed_param  # For plotting: x = phason (fixed per slice), but accumulate across slices



@showprogress for ((N_val, Delta_val, t_n_tuple, fixed_val), f) in combo_to_files

    slice_suffix_plotname = "$(fixed_param_name)$(fixed_val)"

    # # One-off filter for specific phi
    # if fixed_param == :phi && abs(fixed_val - 0.6180048661800487) > 1e-1
    #     continue
    # end

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
    mu_values = Float64[]
    phason_vals = Float64[]
    phason_range = missing
    target_phason_actual = missing
    idop_labels_formatted = Vector{Tuple{Float64, Float64, Int}}()
    phi_vals = Float64[]
    phi_value_for_record = fixed_param == :phi ? fixed_val : missing
    phi_range = missing

    if !(fixed_param in propertynames(df_slice))
        df_slice[!, fixed_param] .= fixed_val
    end

    @info "mp_tol range from mp_tol_df (min, max): ", maximum(mp_tol_df.mp_tol), minimum(mp_tol_df.mp_tol)

    ## 1) Add mp_disc column to df_slice
    IDOPProcessing.add_mp_disc_column!(
        df_slice,
        mp_tol_df;
        mp_targ=mp_targ,
        raw_col=:mp_raw,
    )

    ## 2) Add mu_mp column to df_slice
    IDOPProcessing.add_mu_mp_column!(
        df_slice,
    )

    

    # --- FIX: Ensure fixed parameter column exists in DataFrame ---
    # When loading slices by phason, 'phi' is fixed and might not be in the DF columns.
    # When loading slices by slope, 'phason' is fixed and might not be in the DF columns.
    if !(fixed_param in propertynames(df_slice))
        df_slice[!, fixed_param] .= fixed_val
    end

    # --- DEBUG: summarise slice contents ---
    @info "Slice summary" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val n_rows=nrow(df_slice)
    if nrow(df_slice) > 0
        # unique phi values present in this slice
        varying_present = sort(unique(df_slice[!, varying_param]))
        @info "Unique $varying_param_name values in slice:" varying_present

        # mu_range summary
        if :mu_range in propertynames(df_slice)
            mu_lengths = length.(df_slice.mu_range)
            min_mu_len = minimum(mu_lengths)
            max_mu_len = maximum(mu_lengths)
            @info "mu_range vector lengths: min=$min_mu_len, max=$max_mu_len"
            @info "First mu_range: min=$(minimum(df_slice.mu_range[1])), max=$(maximum(df_slice.mu_range[1]))"
        end

        # mp_raw summary
        if :mp_raw in propertynames(df_slice)
            mp_lengths = length.(df_slice.mp_raw)
            min_mp_len = minimum(mp_lengths)
            max_mp_len = maximum(mp_lengths)
            @info "mp_raw vector lengths: min=$min_mp_len, max=$max_mp_len"
            # Print a sample of mp_raw values for the first phi
            @info "Sample mp_raw (first phi):" df_slice.mp_raw[1][1:min(5, length(df_slice.mp_raw[1]))]
        end

        # mu_mp summary (after processing)
        if :mu_mp in propertynames(df_slice)
            mu_mp_lengths = length.(df_slice.mu_mp)
            min_mu_mp_len = minimum(mu_mp_lengths)
            max_mu_mp_len = maximum(mu_mp_lengths)
            total_mu_mp_points = sum(mu_mp_lengths)
            @info "mu_mp vector lengths: min=$min_mu_mp_len, max=$max_mu_mp_len, total=$total_mu_mp_points"
            # Print a sample of mu_mp values for the first phi
            @info "Sample mu_mp (first phi):" df_slice.mu_mp[1][1:min(5, length(df_slice.mu_mp[1]))]
        end

        # Print phi and phason ranges
        if :phi in propertynames(df_slice)
            phi_vals = [df_slice.phi[i] for i in 1:nrow(df_slice) if !ismissing(df_slice.phi[i])]
            if !isempty(phi_vals)
                @info "phi range: min=$(minimum(phi_vals)), max=$(maximum(phi_vals))"
            end
        end
        if :phason in propertynames(df_slice)
            phason_vals = [df_slice.phason[i] for i in 1:nrow(df_slice) if !ismissing(df_slice.phason[i])]
            if !isempty(phason_vals)
                @info "phason range: min=$(minimum(phason_vals)), max=$(maximum(phason_vals))"
            end
        end
    else
        @warn "df_slice is empty" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val
    end
    # --- END DEBUG ---

    # coeff_filename = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/hof_style_slopes/hof_style_slopes_N500_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-5p0-7)_mu(0p0-0p0-1)_Delta(0p05-0p5-6)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_interp_norm/smooth_norm_coeffs_Delta0.05_t11.0_t22.0_phason0.0_thresh0.001,0.004_prange20_qmax20.txt"
    # df_slice = IDOPProcessing.apply_smooth_normalisation!(
    #     df_slice,
    #     coeff_filename;
    #     target_col=:mu_mp,
    #     out_col=:mu_mp_norm
    # )

    norm_filename = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_outer_norm/eigs_bandwidth_norms_mu0.0_Delta0.05_t_n[1.0, 2.0]_phason=0.0.csv"
    df_slice = IDOPProcessing.apply_outer_normalisation!(
        df_slice,
        norm_filename;
        target_col=:mu_mp,
        out_col=:mu_mp_norm
    )

    # filter out phi=0.0 and 1.0 if present (these are trivial cases)
    df_slice = filter(row -> !(row.phi == 0.0 || row.phi == 1.0), df_slice)

    ## 2a) Plot discrete phase projections
    try
        phase_proj_outdir = joinpath(outdir, "phase_projections")
        isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)
        
        IDOPPlotting.plt_discrete_phase_projections(
            df_slice,
            :mu_mp,
            :phi,
            joinpath(phase_proj_outdir, "TRANSPCOL_phase_proj_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
            # plt_xlims=(-1.0,1.0),
            # plt_ylims=(0.0,1.0)
            # N=N_val,
            # Delta=Delta_val,
            # t_n=t_n,
            # phason=phason_val
        )

    catch e
        @warn "plt_discrete_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    end

    ## 2b) Plot normalised discrete phase projections
    try
        phase_proj_outdir = joinpath(outdir, "phase_projections")
        isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)
        
        IDOPPlotting.plt_discrete_phase_projections(
            df_slice,
            :mu_mp_norm,
            :phi,
            joinpath(phase_proj_outdir, "TRANSPCOL_norm_phase_proj_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
            plt_xlims=(-1.0,1.0),
            plt_ylims=(0.0,1.0)
            # N=N_val,
            # Delta=Delta_val,
            # t_n=t_n,
            # phason=phason_val
        )

    catch e
        @warn "plt_discrete_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    end



    df_slice = nothing
    data = nothing
    GC.gc()
end


println("\nAll IDOP processing complete. Results in: ", outdir)