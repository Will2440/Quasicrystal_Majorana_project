"""
This script is for making the duel plots of labelled energy spectrum + raw mp vs mu
"""

using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter
using LaTeXStrings
using CSV

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
include("../idop/processing.jl")
include("../idos/processing.jl")
include("RAM_opt_unpacker.jl")
include("processing.jl")

using .MPGapUnpacking
using .MPGapPlotting
using .MPGapProcessing
using .IDOPProcessing
using .IDOSProcessing


######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
# data_folder = "hof_style_slopes_N1000_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)"
# run_set_name = "MBData_final"
# data_folder = "hof_style_slopes_N500_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(-3p0-3p0-1201)_Delta(0p05-0p05-1)"
# run_set_name = "hof_style_slopes"
data_folder = "hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"
run_set_name = "QC-SC_comp_small_L_check"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define ranges for repacking or plotting -- otherwise the full range will be used.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phi_range = nothing #[0.6180048661800487]
repack_phason_range = [0.53] #nothing
mu_plot_min, mu_plot_max = 0.0, 3.0 #nothing, nothing

mp_targ = -1.0
mp_tol = 1.2e-1

qc_region_counts_log = []

## Defines base_slice_dir for accessing slices using this script
base_slice_dir = normpath(joinpath(base_results_dir, "mp_gap_comp_base_slices", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)"))

## Defines directory where sweep data (from main_splitUnpack.jl) is stored
# Match parameters from main_splitUnpack.jl
sweep_mp_tol_min = 1e-10
sweep_mp_tol_max = 1.0
sweep_mp_tol_len = 1001
sweep_results_dir = normpath(joinpath(base_results_dir, "mp_tol_sweeps", "mptarg$(mp_targ)_mptolrange$(sweep_mp_tol_min)-$(sweep_mp_tol_max)-$(sweep_mp_tol_len)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)"))

# --- LOAD SWEEP LOOKUP INDEX ---
sweep_lookup = Dict()
sweep_index_path = joinpath(sweep_results_dir, "sweep_index.bson")
sweep_index_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/QC-SC_comp_small_L_check/hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)/mp_tol_sweeps/mptarg-1.0_mptolrange1.0e-10-1.0-1001_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.53]/sweep_index.bson"
if isfile(sweep_index_path)
    println("Loading sweep index from: $sweep_index_path")
    sweep_lookup = BSON.load(sweep_index_path)[:combo_to_file]
else
    @warn "No sweep index found at $sweep_index_path. Sweep plots will be empty."
end
# -------------------------------


unpack_data = false  # Set to true to unpack data, false to process existing slices
if unpack_data
    # Run unpacking and slicing
    combo_to_files = MPGapUnpacking.unpack_and_slice(folder_path_hof, base_slice_dir;
        mp_targ=mp_targ,
        mp_tol=mp_tol,
        plt_N_range_arg=repack_N_range,
        plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phi_range_arg=repack_phi_range,
        plt_phason_range_arg=repack_phason_range
    )
else
    # Processing existing slices
    combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
    combo_to_files = combo_data[:combo_to_files]
end

######################################
########### Extra variables ##########
######################################

## for IDOS
Delta=0.05
idos_plateau_width_threshold = 0.025 #Delta #0.001
N = 200
idos_plateau_height_threshold = 2.5/N
idos_label_err_tol = 1.5/N
idos_mp_tol = 0.1
max_gap_label_value = 20
max_idos_gap_label_value = max_gap_label_value

normalise_type = :outer # :inner, :outer, :both
plot_norm_or_not = :eigenvalues # :eigenvalues or :eigs_interp_norm


## for IDOP
delta_idop_thresh = 1e-4
min_mu_span = 0.0
min_samples = 2
max_idop_gap_label_value = max_gap_label_value

## for plot
max_plot_gap_label_value = max_gap_label_value



## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "mp_gap_comp", "mptarg$(mp_targ)_mptol$(mp_tol)_idos$(idos_plateau_width_threshold),$(idos_plateau_height_threshold),$(max_idos_gap_label_value)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phi_range$(repack_phi_range)_phason_range$(repack_phason_range)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)

######################################
########## Processing Execution ########
######################################

test_index = nothing
if !isnothing(test_index)
    keys_list = collect(keys(combo_to_files))
    if test_index < 1 || test_index > length(keys_list)
        error("test_index $test_index is out of range (1 to $(length(keys_list)))")
    end
    selected_key = keys_list[test_index]
    combos_to_process = [(selected_key, combo_to_files[selected_key])]
else
    combos_to_process = collect(combo_to_files)
end


@showprogress for ((N_val, Delta_val, t_n_tuple, phi_val, phason_val), f) in combos_to_process
    # Apply any additional filtering if needed
    if !isnothing(repack_N_range) && !(N_val in repack_N_range)
        continue
    end
    if !isnothing(repack_Delta_range) && !(Delta_val in repack_Delta_range)
        continue
    end
    if !isnothing(repack_tn_range) && !(t_n_tuple in repack_tn_range)
        continue
    end
    if !isnothing(repack_phi_range) && !(phi_val in repack_phi_range)
        continue
    end
    if !isnothing(repack_phason_range) && !(phason_val in repack_phason_range)
        continue
    end

    println("-----------------------------")
    println("phi_cal: ", phi_val)
    println("-----------------------------")

    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)

    # discretise_evals!(df_slice, 0.0, 1e-2)

    # ## 1) Processing on energy structure: normalise, IDOS and gap labelling
    # ## Using functions from idos/processing.jl
    # df_eigs = df_slice
    # try
    #     # df_eigs = IDOSProcessing.compute_eigs_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_norm)
    #     # df_eigs = IDOSProcessing.compute_eigs_inner_norm!(df_eigs; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
    #     df_eigs = IDOSProcessing.compute_interpolated_normalisation!(
    #         df_slice; 
    #         out_col=:eigs_interp_norm, 
    #         mode=normalise_type, 
    #         ref_indices=(5, 1)
    #     )
    #     df_eigs = IDOSProcessing.compute_idos_df!(df_eigs)
    #     # df_eigs = IDOSProcessing.compute_plateaus_on_idos_df!(df_eigs; threshold=idos_plateau_threshold)
    #     df_eigs = IDOSProcessing.compute_plateaus_on_idos_df_extra_thresh!(
    #         df_slice; 
    #         threshold=idos_plateau_width_threshold, 
    #         # idos_threshold=nothing,
    #         eigs_col=plot_norm_or_not
    #     )
    #     p_range = collect(-max_idos_gap_label_value:max_idos_gap_label_value)
    #     q_max = max_idos_gap_label_value
    #     # df_eigs = IDOSProcessing.compute_gap_labels_qlim!(df_eigs; p_range=p_range, q_max=q_max)
    #     df_slice = IDOSProcessing.compute_gap_labels_parsimonious!(
    #         df_slice; 
    #         p_range=p_range, 
    #         q_max=q_max, 
    #         tol=idos_label_err_tol, 
    #         merge_threshold=idos_plateau_height_threshold,
    #         check_mp=false,
    #         mp_col=:mp,
    #         mp_threshold=idos_mp_tol,
    #         always_mask_central=true
    #     )        
    # catch e
    #     @warn "Gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    # end

    # ## 2) Processing on phase transition
    # ## Using functions from idop/processing.jl
    # df_phase = DataFrame()
    # try
    #     df_mp_preidop = IDOPProcessing.prep_df_for_IDOP_single(df_slice, phi_val, phason_val, N_val, Delta_val, t_n)
    #     df_idop = IDOPProcessing.compute_idop_df!(df_mp_preidop; disc_variable=:disc_mp)
    #     df_idop_plateaus = IDOPProcessing.compute_idop_plateaus_all!(df_idop; delta_idop_thresh=delta_idop_thresh, min_mu_span=min_mu_span, min_samples=min_samples)
    #     p_range = collect(-max_idop_gap_label_value:max_idop_gap_label_value)
    #     q_max = max_idop_gap_label_value
    #     df_phase = IDOPProcessing.compute_gap_labels_qlim_new!(df_idop_plateaus; p_range=p_range, q_max=q_max)
    # catch e
    #     @warn "IDOP gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    # end

    # ## 3) Load Sweep Data (MP vs Tolerance)
    # mp_sweep_matrix = nothing
    # mp_sweep_tols = nothing
    # try
    #     # Use the lookup dict instead of constructing the path manually
    #     # Key: (N, Delta, t_n_tuple, phi, phason)
    #     lookup_key = (N_val, Delta_val, t_n_tuple, phi_val, phason_val)

    #     override_path = nothing #"/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/QC-SC_comp_small_L_check/hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)/mp_tol_sweeps/mptarg-1.0_mptolrange1.0e-15-1.0-1501_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_range[0.53]/sweep_index.bson"

    #     if !isnothing(override_path)
    #         sweep_data = BSON.load(override_path)
    #         mp_sweep_matrix = sweep_data[:disc_mp_matrix]
    #         mp_sweep_tols = sweep_data[:mp_tol_range]
    #     else
    #         if haskey(sweep_lookup, lookup_key)
    #             sweep_path = sweep_lookup[lookup_key]
    #             if isfile(sweep_path)
    #                 sweep_data = BSON.load(sweep_path)
    #                 mp_sweep_matrix = sweep_data[:disc_mp_matrix]
    #                 mp_sweep_tols = sweep_data[:mp_tol_range]
    #                 # println("Loaded sweep data from: $sweep_path")
    #                 # println("mp_sweep_matrix size: ", size(mp_sweep_matrix))
    #                 # println("mp+sweep_matrix contents (first 5 rows, first 5 cols):")
    #                 # println(mp_sweep_matrix[1:min(5,end), 1:min(5,end)])
    #             else
    #                 @warn "Sweep file indexed but not found on disk: $sweep_path"
    #             end
    #         else
    #             # Silent fail or debug log if needed, as some combos might not have sweeps
    #             # @debug "No sweep entry for $lookup_key"
    #         end
    #     end
    # catch e
    #     @warn "Failed to load sweep data for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple: $e"
    # end

    # ## 4) Generate Plots
    # try
    #     mp_gap_comp_outdir = joinpath(outdir, "N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))")
    #     isdir(mp_gap_comp_outdir) || mkpath(mp_gap_comp_outdir)

        
    #     if !isnothing(mu_plot_max) && !isnothing(mu_plot_min)
    #         # --- Manual Filtering for Plotting ---
            
    #         # 1. Filter Slice (Energy & MP)
    #         # Ensure sorted by mu for consistency
    #         sort!(df_slice, :mu)
    #         mask_mu = (df_slice.mu .>= mu_plot_min) .& (df_slice.mu .<= mu_plot_max)
    #         df_slice_plot = df_slice[mask_mu, :]
            
    #         # 2. Filter Sweep Matrix
    #         mp_sweep_matrix_plot = nothing
    #         if !isnothing(mp_sweep_matrix)
    #             if size(mp_sweep_matrix, 1) == nrow(df_slice)
    #                 mp_sweep_matrix_plot = mp_sweep_matrix[mask_mu, :]
    #             else
    #                 @warn "Sweep matrix rows $(size(mp_sweep_matrix, 1)) != df_slice rows $(nrow(df_slice)). Skipping sweep plot."
    #             end
    #         end
            
    #         # 3. Filter IDOP Gaps (Optional, but cleans up "messy" labels)
    #         # Assuming df_phase is the "packed" format (1 row, col :gap_labels -> Vector{NamedTuple})
    #         df_phase_plot = df_phase
    #         if nrow(df_phase) == 1 && "gap_labels" in names(df_phase)
    #             gaps = df_phase[1, :gap_labels]
    #             if isa(gaps, Vector)
    #                 # Filter gaps that overlap with the plot range
    #                 filtered_gaps = filter(g -> (g.E_high >= mu_plot_min) && (g.E_low <= mu_plot_max), gaps)
    #                 df_phase_plot = DataFrame(gap_labels = [filtered_gaps])
    #             end
    #         end
    #         @info "Filtered for mu range [$mu_plot_min, $mu_plot_max]."
    #     end

    #     MPGapPlotting.plt_mp_gap_comparison_full(
    #         df_slice_plot,  # Use filtered slice
    #         joinpath(mp_gap_comp_outdir, "mp_gap_comp_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phi$(phi_val)_phason$(phason_val).svg");
    #         gap_segments=df_slice_plot,  # Use filtered slice (contains gap labels from step 1)
    #         mp_gap_segments=df_phase_plot,  # Use filtered IDOP gaps
    #         full_range=false,
    #         mu_range_lim=(mu_plot_min, mu_plot_max),
    #         q_range=collect(-max_plot_gap_label_value:max_plot_gap_label_value),
    #         plot_disc_mp=false,
    #         mp_sweep_matrix=mp_sweep_matrix_plot, # Use filtered matrix
    #         mp_sweep_tols=mp_sweep_tols,
    #         raw_mp_tol=mp_tol,
    #         raw_mp_target=mp_targ
    #     )

    #     println("Final plot saved to: ", joinpath(mp_gap_comp_outdir, "mp_gap_comp_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phi$(phi_val)_phason$(phason_val).svg"))
    # catch e
    #     @warn "Plotting failed for  phi=$phi_val, phason=$phason_val, N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple: $e"
    # end


    # #########################################################################################################
    # #########################################################################################################
    # ########################################## TESTING QC/SC CRITERION ######################################
    # #########################################################################################################
    # ## 5) IDOS computation
    # eigs_norm_col = :eigs_norm
    # coeffs_out_file = joinpath(outdir, "idos_norm_coeffs_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phi$(phi_val)_phason$(phason_val).csv")
    # df_slice = IDOSProcessing.compute_smooth_normalisation!(
    #     df_slice;
    #     out_col=:eigs_norm,
    #     eigs_col=:eigenvalues,
    #     poly_deg=4, 
    #     coeffs_out_file=coeffs_out_file
    # )
    # df_slice = IDOSProcessing.compute_idos_df!(df_slice)
    # df_slice = IDOSProcessing.compute_plateaus_on_idos_df_extra_thresh!(
    #     df_slice; 
    #     threshold=idos_plateau_width_threshold,
    #     eigs_col=:eigenvalues #eigs_norm_col
    # )
    # df_slice = IDOSProcessing.compute_gap_labels_parsimonious!(
    #     df_slice; 
    #     p_range=collect(-max_idos_gap_label_value:max_idos_gap_label_value), 
    #     q_max=max_idos_gap_label_value, 
    #     tol=idos_label_err_tol, 
    #     merge_threshold=idos_plateau_height_threshold,
    #     check_mp=true,
    #     mp_col=:mp,
    #     mp_threshold=mp_tol,
    #     always_mask_central=true
    # )
    # ## df_slice now has column :eigs_norm, :idos, :plateaus_idxd, :plateaus, :gap_labels


    # ## 5) Bulk vs QC energy gap processing
    #     ## a) Compute bulk gaps
    #         #  Ignoring the MBS zero energy state when degenerate within tol
    #         #  calculate the energy difference between two lowest energy states,
    # # zero_energy_tol=1e-8
    
    # tol_range = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    # for zero_energy_tol in tol_range
    #     df_slice, bulk_gaps = MPGapProcessing.compute_bulk_gaps!(
    #         df_slice;
    #         eig_col=:eigenvalues,
    #         tol=zero_energy_tol,
    #     )

    #     mid_gaps, outer_gaps, robust_bulk_gaps = MPGapProcessing.compute_bulk_gaps_robust(df_slice)

    #         ##  b) Compute QC gaps at mu=0.0
    #             # Over energy_range==mu_range calculate the gap size at each point
    #             # 1) Compute directly on the eigenvalues
    #     E_range, qc_gaps_direct = MPGapProcessing.compute_qc_gaps_direct(
    #         df_slice;
    #         mu_val=0.0
    #     )
    #             # 2) Compute on the IDOS plateaus
    #     E_range, qc_gaps_idos = MPGapProcessing.compute_qc_gaps_idos(
    #         df_slice;
    #         mu_val=0.0,
    #         plateaus_col=:gap_labels
    #     )
    #             # 3) Compute robust QC gaps directly (skipping one eval each time)
    #     E_range, qc_gaps_direct_robust, merged_count = MPGapProcessing.compute_qc_gaps_direct_robust(
    #         df_slice;
    #         mu_val=0.0
    #     )

    #     println("merged_count of small gaps in robust QC gap computation: ", merged_count)


    #     ## 6) Bulk vs QC energy gap plotting
    #     bulk_vs_qc_gap_plt_outdir = joinpath(outdir, "bulk_vs_qc_gaps")
    #     isdir(bulk_vs_qc_gap_plt_outdir) || mkpath(bulk_vs_qc_gap_plt_outdir)

    #     xlim = (0.0, maximum(df_slice.mu))
    #     logscale = false

    #     ## a) Plot bulk vs QC gaps from eigs direct
    #     try 
    #         filename = "bulk_vs_qc_gaps_from_eigs_direct_tol$(zero_energy_tol)_logscale$(logscale)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_bulk_vs_qc_gaps(
    #             df_slice.mu,
    #             bulk_gaps,
    #             E_range,
    #             qc_gaps_direct,
    #             savepath;
    #             plot_title="Bulk vs QC Gaps from Eigs Direct",
    #             xlims=xlim,
    #             logyscale=logscale
    #         )
    #     catch e
    #         @warn "Plotting bulk vs QC gaps from eigs direct failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end

    #     ## b) Plot mid_gaps, outer_gaps and bulk_gaps together
    #     try
    #         filename = "bulk_mid_outer_gaps_tol$(zero_energy_tol)_logscale$(logscale)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_all_gap_components(
    #             df_slice.mu,
    #             bulk_gaps,
    #             E_range,
    #             qc_gaps_direct,
    #             savepath;
    #             plot_title="Bulk, Mid, and Outer Gaps with QC gaps from Eigs Direct",
    #             xlims=xlim,
    #             logyscale=logscale,
    #             mid_gaps=mid_gaps,
    #             outer_gaps=outer_gaps
    #         )
    #     catch e
    #         @warn "Plotting bulk, mid, outer gaps failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end

    #     ## c) Plot bulk vs QC gaps from IDOS plateaus
    #     try 
    #         filename = "bulk_vs_qc_gaps_from_idos_tol$(zero_energy_tol)_logscale$(logscale)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_bulk_vs_qc_gaps(
    #             df_slice.mu,
    #             bulk_gaps,
    #             E_range,
    #             qc_gaps_idos,
    #             savepath;
    #             plot_title="Bulk vs QC Gaps from IDOS Plateaus",
    #             xlims=xlim,
    #             logyscale=logscale
    #         )
    #     catch e
    #         @warn "Plotting bulk vs QC gaps from IDOS plateaus failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end

    #     ## d) Plot bulk vs QC gap difference from direct eigs
    #     try
    #         filename = "bulk_vs_qc_gaps_direct_diff_tol$(zero_energy_tol)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_bulk_vs_qc_diff(
    #             df_slice.mu,
    #             bulk_gaps,
    #             E_range,
    #             qc_gaps_direct,
    #             savepath;
    #             plot_title="Bulk - QC Gap Difference from Eigs Direct",
    #             xlims=xlim
    #         )
    #     catch e
    #         @warn "Plotting bulk vs QC gap difference from IDOS plateaus failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end

    #     ## e) Plot robust bulk gaps vs QC gaps robust from eigs robust direct
    #     try
    #         filename = "bulk_vs_qc_gaps_robust_from_eigs_robust_direct_tol$(zero_energy_tol)_logscale$(logscale)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_bulk_vs_qc_gaps(
    #             df_slice.mu,
    #             robust_bulk_gaps,   
    #             E_range,
    #             qc_gaps_direct_robust,
    #             savepath;
    #             plot_title="Bulk vs QC Gaps Robust from Eigs Direct",
    #             xlims=xlim,
    #             logyscale=logscale
    #         )
    #     catch e
    #         @warn "Plotting bulk vs QC gaps robust from eigs direct failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end

    #     ## f) Plot bulk vs QC gap difference from robust bulk and QC gaps
    #     try
    #         filename = "bulk_vs_qc_gaps_robust_direct_diff_tol$(zero_energy_tol)"
    #         savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
    #         MPGapPlotting.plt_bulk_vs_qc_diff(
    #             df_slice.mu,
    #             robust_bulk_gaps,
    #             E_range,
    #             qc_gaps_direct_robust,
    #             savepath;
    #             plot_title="Bulk - QC Gap Difference from Eigs Direct Robust",
    #             xlims=xlim
    #         )
    #     catch e
    #         @warn "Plotting bulk vs QC gap difference from robust direct failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    #     end
    # end
    # #########################################################################################################
    # #########################################################################################################





    #########################################################################################################
    ### ACTUAL QC/SC PROCESSING NEEDED
    #########################################################################################################
    # println("Testing what evals data is saved in df_slice:")
    # for (i, evals) in enumerate(df_slice.eigenvalues[1:min(10, end)])
    #     println("Row $i: typeof=", typeof(evals), ", length=", ismissing(evals) ? "missing" : length(evals), ", values=", evals)
    # end
    # println("the middle two evals are: ", df_slice.eigenvalues[1, div(end,2)], ", ", df_slice.eigenvalues[1, div(end,2)+1])

    ## a) Compute robust bulk gaps
    mid_gaps, outer_gaps, robust_bulk_gaps = MPGapProcessing.compute_bulk_gaps_robust(
        df_slice
    )

    # println("got here 1")

    ## b) Load phason-independent spectra from IDOS results
    # This assumes the structure generated by idos/main_bigData.jl
    # phason_spec_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_interp_norm/phason_independent_spectra.bson"
    phason_spec_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/QC-SC_comp_small_L_check/hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_outer_norm/phason_independent_spectra.bson"

    if !(@isdefined PHASON_INDEP_DICT)
        if isfile(phason_spec_path)
            @info "Loading phason independent spectra from: $phason_spec_path"
            PHASON_INDEP_DICT = BSON.load(phason_spec_path)[:phason_indep_spectra]
        else
            @warn "Phason independent spectra file not found at: $phason_spec_path"
            PHASON_INDEP_DICT = Dict()
        end
    end

    # Find the spectrum for this combo at mu=0.0 (approximate match for slope and mu)
    closest_key = nothing
    min_dist = Inf
    for k in keys(PHASON_INDEP_DICT)
        # Require exact match for integer/fixed fields
        if k.N != N_val || k.Delta != Delta_val || k.t_n != t_n_tuple
            continue
        end
        # Compute distance in floating-point fields (slope and mu)
        dist = sqrt((k.slope - phi_val)^2 + (k.mu - 0.0)^2)
        if dist < min_dist
            min_dist = dist
            closest_key = k
        end
    end

    spectrum = if !isnothing(closest_key)
        if min_dist > 1e-4
            @warn "Closest phason independent spectrum found but distance $min_dist > 1e-4 for target: slope=$phi_val, N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, mu=0.0. Using it anyway."
        end
        PHASON_INDEP_DICT[closest_key]
    else
        @warn "No phason independent spectrum found at all for target: slope=$phi_val, N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, mu=0.0"
        Float64[]
    end

    # Attach to df_slice: spectrum only for mu=0.0 row, empty elsewhere
    phason_indep_col = Vector{Vector{Float64}}(undef, nrow(df_slice))
    for (i, row) in enumerate(eachrow(df_slice))
        phason_indep_col[i] = (row.mu == 0.0) ? spectrum : Float64[]
    end
    df_slice[!, :phason_indep_eigs] = phason_indep_col
    # --------------------------------------------------


    ## manually remove some phason states from eigenvalues at mu=0.0

    # energies at which to remove states
    remove_evals_nearest = [0.21, 0.25, 1.35, 1.39, 1.97, 1.99]
    df_slice[!, :manual_eigenvalues] = copy(df_slice.eigenvalues)

    # remove nearest eval to each target energy
    mu0_row_idx = findfirst(df_slice.mu .== 0.0)
    if !isnothing(mu0_row_idx)
        evals_at_mu0 = df_slice.eigenvalues[mu0_row_idx]
        for target_eval in remove_evals_nearest
            if !isempty(evals_at_mu0)
                diffs = abs.(evals_at_mu0 .- target_eval)
                nearest_idx = argmin(diffs)
                eval_to_remove = evals_at_mu0[nearest_idx]
                # create new df column called manual_eigenvalues with all evals except removed ones
                new_evals = filter(x -> x != eval_to_remove, evals_at_mu0)
                df_slice.manual_eigenvalues[mu0_row_idx] = new_evals
                @info "Removed eval $eval_to_remove near target $target_eval at mu=0.0 for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val"
                # update evals_at_mu0 for next iteration
                evals_at_mu0 = new_evals
            end
        end
    end

    # manually remove all evals within small windows
    remove_evals_windows = [(1.1, 1.2)]
    for (E_low, E_high) in remove_evals_windows
        mu0_row_idx = findfirst(df_slice.mu .== 0.0)
        if !isnothing(mu0_row_idx)
            evals_at_mu0 = df_slice.manual_eigenvalues[mu0_row_idx]
            new_evals = filter(x -> x < E_low || x > E_high, evals_at_mu0)
            df_slice.manual_eigenvalues[mu0_row_idx] = new_evals
            @info "Removed evals in window [$E_low, $E_high] at mu=0.0 for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val"
        end
    end

    # ## c) Compute QC gaps at mu=0.0 (CHOOSE ONE METHOD)
    # #     ## i) skip-merging method
    # E_range, qc_gaps_direct_robust, merged_count = MPGapProcessing.compute_qc_gaps_direct_robust(
    #     df_slice;
    #     mu_val=0.0,
    #     eig_col=:phason_indep_eigs,
    #     hit_eig_tol=1e-9
    # )
    merged_count =  NaN

    # println("got here 2")
    E_range, qc_gaps_direct_robust = MPGapProcessing.compute_qc_gaps_direct(
        df_slice;
        mu_val=0.0,
        eig_col=:manual_eigenvalues,
        # eig_tol=1e-9,
        N_val=N_val
    )


    #     ## ii) standard direct method on phason indpendent spectra
    # E_range, qc_gaps_direct_robust = MPGapProcessing.compute_qc_gaps_direct(
    #     df_slice;
    #     mu_val=0.0,
    #     eig_col=:phason_indep_eigs
    # )

    # d) Count QC regions using robust bulk and QC gaps
    xlim = (0.0, 3.0) #maximum(df_slice.mu))

    qc_region_count = MPGapProcessing.count_qc_beats_sc_regions(
        df_slice.mu,
        robust_bulk_gaps,
        E_range,
        qc_gaps_direct_robust,
        xlims=xlim
    )

    ## e) Log results
    push!(qc_region_counts_log, (
        N=N_val,
        Delta=Delta_val,
        t_n=t_n_tuple,
        phi=phi_val,
        phason=phason_val,
        merged_count=merged_count,
        qc_region_count=qc_region_count
    ))

    ##########################
    #### Optional Plotting ####
    ###########################
    bulk_vs_qc_gap_plt_outdir = joinpath(outdir, "bulk_vs_qc_gaps")
    isdir(bulk_vs_qc_gap_plt_outdir) || mkpath(bulk_vs_qc_gap_plt_outdir)
    ## a) Plot robust bulk gaps vs QC gaps robust from eigs robust direct
    try
        println("got here 3")
        logscale=false
        filename = "MANUAL_bulk_vs_qc_gaps_robust_from_eigs_robust_direct_logscale$(logscale)_phi$(phi_val)_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phason$(phason_val)"
        savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
        println("Plotting to: ", savepath)
        MPGapPlotting.plt_bulk_vs_qc_gaps(
            df_slice.mu,
            robust_bulk_gaps,   
            E_range,
            qc_gaps_direct_robust,
            savepath;
            plot_title="Bulk vs QC Gaps Robust from Eigs Direct",
            xlims=xlim,
            logyscale=logscale
        )
    catch e
        @warn "Plotting bulk vs QC gaps robust from eigs direct failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    end

    ## b) Plot bulk vs QC gap difference from robust bulk and QC gaps
    try
        println("got here 4")
        filename = "MANUAL_bulk_vs_qc_gaps_robust_direct_diff_phi$(phi_val)_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phason$(phason_val)"
        savepath = joinpath(bulk_vs_qc_gap_plt_outdir, "$(filename).svg")
        MPGapPlotting.plt_bulk_vs_qc_diff(
            df_slice.mu,
            robust_bulk_gaps,
            E_range,
            qc_gaps_direct_robust,
            savepath;
            plot_title="Bulk - QC Gap Difference from Eigs Direct Robust",
            xlims=xlim
        )
    catch e
        @warn "Plotting bulk vs QC gap difference from robust direct failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    end


    df_slice = nothing
    data = nothing
    GC.gc()
end

# Save qc_region_counts_log to a CSV for later analysis
qc_log_df = DataFrame(qc_region_counts_log)
qc_log_savepath = joinpath(outdir, "qc_region_counts_log.csv")
CSV.write(qc_log_savepath, qc_log_df)

println("\nAll MP gap comp processing complete. Results in: ", outdir)