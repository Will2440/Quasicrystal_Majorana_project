# using BSON
# using DataFrames
# using PrettyTables
# using Printf
# using Logging
# using Colors, ColorSchemes, Plots
# using ProgressMeter

# # Set logger level from JULIA_LOG_LEVEL environment variable, default to Info if not set
# level_str = uppercase(get(ENV, "JULIA_LOG_LEVEL", "info"))
# level = if level_str == "DEBUG"
#     Logging.Debug
# elseif level_str == "INFO"
#     Logging.Info
# elseif level_str == "WARN"
#     Logging.Warn
# elseif level_str == "ERROR"
#     Logging.Error
# else
#     Logging.Info  # fallback
# end
# global_logger(ConsoleLogger(stderr, level))

# include("plotting.jl")
# include("../idop/processing.jl")
# include("../idos/processing.jl")
# include("RAM_opt_unpacker.jl")

# using .MPGapUnpacking
# using .MPGapPlotting
# using .IDOPProcessing
# using .IDOSProcessing

# ######################################
# ########## Path Determination ########
# ######################################

# ## Determine the data folder for saving, should be consistent with raw data location
# data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
# run_set_name = "sturmian_sweep_t1-t2_const_mapping"

# ## Ensure data is being lifted from the correct raw data folder
# folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

# ## Defines results directory for saving slices
# base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

# ######################################
# ############### Unpack? ##############
# ######################################

# # Define ranges for repacking or plotting -- otherwise the full range will be.
# repack_N_range = nothing
# repack_Delta_range = [0.05] #nothing
# repack_tn_range = [(1.0,1.5)] #nothing
# repack_phi_range = nothing
# repack_phason_range = nothing

# mp_targ = -1.0
# mp_tol = 5e-2

# ## Defines base_slice_dir for accessing slices using this script
# base_slice_dir = normpath(joinpath(base_results_dir, "mp_gap_comp_base_slices", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)"))


# unpack_data = true  # Set to true to unpack data, false to process existing slices
# if unpack_data
#     # Run unpacking and slicing
#     combo_to_files = MPGapUnpacking.unpack_and_slice(folder_path_hof, base_slice_dir;
#         mp_targ=mp_targ,
#         mp_tol=mp_tol,
#         plt_N_range_arg=repack_N_range,
#         plt_Delta_range_arg=repack_Delta_range,
#         plt_tn_range_arg=repack_tn_range,
#         plt_phi_range_arg=repack_phi_range,
#         plt_phason_range_arg=repack_phason_range
#     )
# else
#     # Processing existing slices
#     combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
#     combo_to_files = combo_data[:combo_to_files]
# end

# ######################################
# ########### Extra variables ##########
# ######################################

# ## for IDOS
# idos_plateau_threshold = 0.05
# max_idos_gap_label_value = 5

# ## for IDOP
# delta_idop_thresh = 1e-4
# min_mu_span = 0.0
# min_samples = 1
# max_idop_gap_label_value = 5



# ## Defines outdir for saving plots in this script
# outdir = normpath(joinpath(base_results_dir, "mp_gap_comp", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phi_range$(repack_phi_range)_phason_range$(repack_phason_range)"))
# isdir(outdir) || mkpath(outdir)
# println("Saving results to folder: ", outdir)

# ######################################
# ########## Processing Execution ########
# ######################################

# test_index = 102
# if !isnothing(test_index)
#     keys_list = collect(keys(combo_to_files))
#     if test_index < 1 || test_index > length(keys_list)
#         error("test_index $test_index is out of range (1 to $(length(keys_list)))")
#     end
#     selected_key = keys_list[test_index]
#     combos_to_process = [(selected_key, combo_to_files[selected_key])]
# else
#     combos_to_process = collect(combo_to_files)
# end

# # @showprogress for ((N_val, Delta_val, t_n_tuple, phi_val, phason_val), f) in combo_to_files
# @showprogress for ((N_val, Delta_val, t_n_tuple, phi_val, phason_val), f) in combos_to_process
#     # Apply any additional filtering if needed
#     if !isnothing(repack_N_range) && !(N_val in repack_N_range)
#         continue
#     end
#     if !isnothing(repack_Delta_range) && !(Delta_val in repack_Delta_range)
#         continue
#     end
#     if !isnothing(repack_tn_range) && !(t_n_tuple in repack_tn_range)
#         continue
#     end
#     if !isnothing(repack_phi_range) && !(phi_val in repack_phi_range)
#         continue
#     end
#     if !isnothing(repack_phason_range) && !(phason_val in repack_phason_range)
#         continue
#     end

#     data = BSON.load(f)
#     df_slice = data[:df_slice]
#     t_n = collect(t_n_tuple)

#     # discretise_evals!(df_slice, 0.0, 1e-2)

#     ## 1) Processing on energy structure: normalise, IDOS and gap labelling
#     ## Using functions from idos/processing.jl
#     df_eigs = df_slice
#     try
#         df_eigs = IDOSProcessing.compute_eigs_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_norm)
#         df_eigs = IDOSProcessing.compute_eigs_inner_norm!(df_eigs; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
#         df_eigs = IDOSProcessing.compute_idos_df!(df_eigs)
#         df_eigs = IDOSProcessing.compute_plateaus_on_idos_df!(df_eigs; threshold=idos_plateau_threshold)
#         p_range = collect(-max_idos_gap_label_value:max_idos_gap_label_value)
#         q_max = max_idos_gap_label_value
#         df_eigs = IDOSProcessing.compute_gap_labels_qlim!(df_eigs; p_range=p_range, q_max=q_max)
#     catch e
#         @warn "Gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
#     end

#     ## 2) Processing on phase transition
#     ## Using functions from idop/processing.jl
#     df_phase = DataFrame()
#     try
#         df_mp_preidop = IDOPProcessing.prep_df_for_IDOP_single(df_slice, phi_val, phason_val, N_val, Delta_val, t_n)
#         df_idop = IDOPProcessing.compute_idop_df!(df_mp_preidop; disc_variable=:disc_mp)
#         df_idop_plateaus = IDOPProcessing.compute_idop_plateaus_all!(df_idop; delta_idop_thresh=delta_idop_thresh, min_mu_span=min_mu_span, min_samples=min_samples)
#         p_range = collect(-max_idop_gap_label_value:max_idop_gap_label_value)
#         q_max = max_idop_gap_label_value
#         df_phase = IDOPProcessing.compute_gap_labels_qlim_new!(df_idop_plateaus; p_range=p_range, q_max=q_max)
#     catch e
#         @warn "IDOP gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
#     end

#     # println("df_eigs: ", df_eigs)
#     # println("df_phase: ", df_phase)
#     # println("gap labels column in df_eigs: ", df_eigs.gap_labels)
#     # println("gap labels column in df_pahase: ", df_phase.gap_labels)

#     try
#         mp_gap_comp_outdir = joinpath(outdir, "test")#"N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))")
#         isdir(mp_gap_comp_outdir) || mkpath(mp_gap_comp_outdir)

#         MPGapPlotting.plt_mp_gap_comparison(
#             df_slice,  # Raw slice DataFrame with :mu, :eigenvalues, :mp, :disc_mp
#             joinpath(mp_gap_comp_outdir, "mp_gap_comp_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phi$(phi_val)_phason$(phason_val).png");
#             gap_segments=df_eigs,  # Processed DataFrame with gap columns (:gap_low, :gap_high, :q)
#             mp_gap_segments=df_phase,  # Processed DataFrame with IDOP gap columns (:gap_low, :gap_high, :q)
#             full_range=false,
#             q_range=collect(-5:5),
#             plot_disc_mp=false
#         )
#     catch e
#         @warn "Plotting failed for  phi=$phi_val, phason=$phason_val, N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple: $e"
#     end

#     df_slice = nothing
#     data = nothing
#     GC.gc()
# end

# println("\nAll MP gap comp processing complete. Results in: ", outdir)




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
include("../idop/processing.jl")
include("../idos/processing.jl")
include("RAM_opt_unpacker.jl")

using .MPGapUnpacking
using .MPGapPlotting
using .IDOPProcessing
using .IDOSProcessing

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturmian_sweep_t1-t2_const_mapping"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define ranges for repacking or plotting -- otherwise the full range will be.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phi_range = nothing
repack_phason_range = nothing

mp_targ = -1.0
mp_tol = 5e-2

## Defines base_slice_dir for accessing slices using this script
base_slice_dir = normpath(joinpath(base_results_dir, "mp_gap_comp_base_slices", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)"))

# ## Defines directory where sweep data (from main_splitUnpack.jl) is stored
# sweep_results_dir = normpath(joinpath(base_results_dir, "mp_tol_sweeps", "mptarg$(mp_targ)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)"))

# # --- LOAD SWEEP LOOKUP INDEX ---
# sweep_lookup = Dict()
# sweep_index_path = joinpath(sweep_results_dir, "sweep_index.bson")
# if isfile(sweep_index_path)
#     println("Loading sweep index from: $sweep_index_path")
#     sweep_lookup = BSON.load(sweep_index_path)[:combo_to_file]
# else
#     @warn "No sweep index found at $sweep_index_path. Sweep plots will be empty."
# end
# # -------------------------------


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
idos_plateau_threshold = 0.05
max_idos_gap_label_value = 5

## for IDOP
delta_idop_thresh = 1e-4
min_mu_span = 0.0
min_samples = 1
max_idop_gap_label_value = 5



## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "mp_gap_comp", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phi_range$(repack_phi_range)_phason_range$(repack_phason_range)"))
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

# @showprogress for ((N_val, Delta_val, t_n_tuple, phi_val, phason_val), f) in combo_to_files
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

    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)

    # discretise_evals!(df_slice, 0.0, 1e-2)

    ## 1) Processing on energy structure: normalise, IDOS and gap labelling
    ## Using functions from idos/processing.jl
    df_eigs = df_slice
    try
        df_eigs = IDOSProcessing.compute_eigs_norm!(df_slice; eigs_col=:eigenvalues, out_col=:eigs_norm)
        df_eigs = IDOSProcessing.compute_eigs_inner_norm!(df_eigs; eigs_col=:eigenvalues, out_col=:eigs_inner_norm)
        df_eigs = IDOSProcessing.compute_idos_df!(df_eigs)
        df_eigs = IDOSProcessing.compute_plateaus_on_idos_df!(df_eigs; threshold=idos_plateau_threshold)
        p_range = collect(-max_idos_gap_label_value:max_idos_gap_label_value)
        q_max = max_idos_gap_label_value
        df_eigs = IDOSProcessing.compute_gap_labels_qlim!(df_eigs; p_range=p_range, q_max=q_max)
    catch e
        @warn "Gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    end

    ## 2) Processing on phase transition
    ## Using functions from idop/processing.jl
    df_phase = DataFrame()
    try
        df_mp_preidop = IDOPProcessing.prep_df_for_IDOP_single(df_slice, phi_val, phason_val, N_val, Delta_val, t_n)
        df_idop = IDOPProcessing.compute_idop_df!(df_mp_preidop; disc_variable=:disc_mp)
        df_idop_plateaus = IDOPProcessing.compute_idop_plateaus_all!(df_idop; delta_idop_thresh=delta_idop_thresh, min_mu_span=min_mu_span, min_samples=min_samples)
        p_range = collect(-max_idop_gap_label_value:max_idop_gap_label_value)
        q_max = max_idop_gap_label_value
        df_phase = IDOPProcessing.compute_gap_labels_qlim_new!(df_idop_plateaus; p_range=p_range, q_max=q_max)
    catch e
        @warn "IDOP gap label computation failed for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple, phi=$phi_val, phason=$phason_val: $e"
    end

    # ## 3) Load Sweep Data (MP vs Tolerance)
    # mp_sweep_matrix = nothing
    # mp_sweep_tols = nothing
    # try
    #     # Use the lookup dict instead of constructing the path manually
    #     # Key: (N, Delta, t_n_tuple, phi, phason)
    #     lookup_key = (N_val, Delta_val, t_n_tuple, phi_val, phason_val)
        
    #     if haskey(sweep_lookup, lookup_key)
    #         sweep_path = sweep_lookup[lookup_key]
    #         if isfile(sweep_path)
    #             sweep_data = BSON.load(sweep_path)
    #             mp_sweep_matrix = sweep_data[:disc_mp_matrix]
    #             mp_sweep_tols = sweep_data[:mp_tol_range]
    #         else
    #             @warn "Sweep file indexed but not found on disk: $sweep_path"
    #         end
    #     else
    #         # Silent fail or debug log if needed, as some combos might not have sweeps
    #         # @debug "No sweep entry for $lookup_key"
    #     end
    # catch e
    #     @warn "Failed to load sweep data for N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple: $e"
    # end

    try
        mp_gap_comp_outdir = joinpath(outdir, "N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))")
        isdir(mp_gap_comp_outdir) || mkpath(mp_gap_comp_outdir)

        MPGapPlotting.plt_mp_gap_comparison(
            df_slice,  # Raw slice DataFrame with :mu, :eigenvalues, :mp, :disc_mp
            joinpath(mp_gap_comp_outdir, "mp_gap_comp_N$(N_val)_Delta$(Delta_val)_tn$(join(t_n, "_"))_phi$(phi_val)_phason$(phason_val).png");
            gap_segments=df_eigs,  # Processed DataFrame with gap columns (:gap_low, :gap_high, :q)
            mp_gap_segments=df_phase,  # Processed DataFrame with IDOP gap columns (:gap_low, :gap_high, :q)
            full_range=false,
            q_range=collect(-5:5),
            plot_disc_mp=false,
            # Pass the new sweep data to the plotting function
            # mp_sweep_matrix=mp_sweep_matrix,
            # mp_sweep_tols=mp_sweep_tols
        )
    catch e
        @warn "Plotting failed for  phi=$phi_val, phason=$phason_val, N=$N_val, Delta=$Delta_val, t_n=$t_n_tuple: $e"
    end

    df_slice = nothing
    data = nothing
    GC.gc()
end

println("\nAll MP gap comp processing complete. Results in: ", outdir)