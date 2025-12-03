using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using ProgressMeter
using Plots # Required for plotting

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

include("split_unpacker.jl")
include("plotting.jl") # Include plotting module

using .MPGapSplitUnpacking
using .MPGapPlotting # Use plotting module

######################################
########## Path Determination ########
######################################

data_folder = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_const_mapping_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-2p5-3)_mu(0p0-3p0-601)_Delta(0p0-0p2-5)"
run_set_name = "sturmian_sweep_t1-t2_const_mapping"

folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack ###############
######################################

# Define ranges for filtering
repack_N_range = nothing
repack_Delta_range = [0.01, 0.05, 0.1] #nothing
repack_tn_range = [(1.0, 1.5)] #nothing
repack_phi_range = nothing
repack_phason_range = nothing

mp_targ = -1.0

# Define the full range of tolerances to sweep
mp_tol_range = 10 .^ range(log10(1e-10), log10(1e-0), length=1001)
println("mp_tol_range (min, max, length): ", (minimum(mp_tol_range), maximum(mp_tol_range), length(mp_tol_range)))

# Output directory for the sweep data
sweep_results_dir = normpath(joinpath(base_results_dir, "mp_tol_sweeps", "mptarg$(mp_targ)_mptolrange$(minimum(mp_tol_range))-$(maximum(mp_tol_range))-$(length(mp_tol_range))_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)"))

unpack_data = false

if unpack_data
    # 1. Load raw data ONCE
    df_base_raw = MPGapSplitUnpacking.load_raw_dataset(folder_path_hof)

    # 2. Process all tolerances and save aggregated BSONs
    combo_to_files = MPGapSplitUnpacking.process_tolerance_sweeps(
        df_base_raw, 
        sweep_results_dir, 
        mp_tol_range;
        mp_targ=mp_targ,
        plt_N_range_arg=repack_N_range,
        plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phi_range_arg=repack_phi_range,
        plt_phason_range_arg=repack_phason_range
    )
    
    # Clean up
    df_base_raw = nothing
    GC.gc()
else
    # Load existing index if not unpacking
    index_path = joinpath(sweep_results_dir, "sweep_index.bson")
    if isfile(index_path)
        combo_data = BSON.load(index_path)
        combo_to_files = combo_data[:combo_to_file]
    else
        error("No sweep index found at $index_path")
    end
end

println("\nSweep processing complete. Data in: ", sweep_results_dir)

######################################
########## Plotting Checks ###########
######################################

println("\nGenerating MP Tolerance Check Plots...")

# Create a directory for the check plots
check_plots_dir = joinpath(sweep_results_dir, "check_plots")
isdir(check_plots_dir) || mkpath(check_plots_dir)

# Iterate through the generated files and plot
p = Progress(length(combo_to_files), desc="Plotting checks: ")

for ((N_val, Delta_val, t_n_tuple, phi_val, phason_val), fpath) in combo_to_files
    
    # Construct a filename for the plot
    plt_folder = joinpath(check_plots_dir, "N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phason$(phason_val)")
    isdir(plt_folder) || mkpath(plt_folder)
    plot_fname = "check_plt_phi$(phi_val).png"
    save_path = joinpath(plt_folder, plot_fname)

    try
        MPGapPlotting.plt_mu_vs_mptol_check(
            fpath; # Pass the BSON file path directly
            N=N_val,
            mode=:full, # Highlight ranges matching expected gaps
            savepath=save_path
        )
    catch e
        @warn "Plotting failed for $fpath: $e"
    end
    
    next!(p)
end

println("Check plots saved to: ", check_plots_dir)