"""
This script is for making the mp_tol vs mu plots
"""

using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using ProgressMeter
using Plots # Required for plotting
using DelimitedFiles
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

include("split_unpacker.jl")
include("processing.jl") # Include processing module
include("plotting.jl") # Include plotting module

using .MPGapSplitUnpacking
using .MPGapProcessing # Use processing module
using .MPGapPlotting # Use plotting module

######################################
########## Path Determination ########
######################################

# data_folder = "hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1001-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p0-801)_Delta(0p05-0p05-1)"
# run_set_name = "fig_2_final"
# data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)"
# run_set_name = "MBData_final"
# data_folder = "hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)"
# run_set_name = "QC-SC_comp_small_L_check"
data_folder = "hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(1p0-2p0-7)_mu(0p0-3p0-1201)_Delta(0p01-0p5-7)"
run_set_name = "QC-SC_comp_N500"


folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))
isdir(base_results_dir) || mkpath(base_results_dir)

######################################
############### Unpack ###############
######################################

# Define ranges for filtering
repack_N_range = nothing
repack_Delta_range = [0.01] #nothing
repack_tn_range = nothing
repack_phi_range = nothing #[0.6180048661800487] 
repack_phason_range = nothing

mp_targ = -1.0

data_mu_range_full = false
plot_mu_range_full = false


# Define the full range of tolerances to sweep
mp_tol_range = 10 .^ range(log10(1e-15), log10(1e-0), length=1501)
println("mp_tol_range (min, max, length): ", (minimum(mp_tol_range), maximum(mp_tol_range), length(mp_tol_range)))

# Output directory for the sweep data
sweep_results_dir = normpath(joinpath(base_results_dir, "mp_tol_sweeps", "mptarg$(mp_targ)_mptolrange$(minimum(mp_tol_range))-$(maximum(mp_tol_range))-$(length(mp_tol_range))_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phi_range$(repack_phi_range)_phason_range$(repack_phason_range)"))
isdir(sweep_results_dir) || mkpath(sweep_results_dir)
# println("Sweep results directory: ", sweep_results_dir)

unpack_data = false 

# Load raw Hofstadter dataset for subsequent plotting steps
df_base_raw = MPGapSplitUnpacking.load_raw_dataset(folder_path_hof)

if unpack_data
    println("got here 1")
    # Process tolerances and persist sweep BSONs
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
else
    # Reuse precomputed sweep index
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
######## Fetch Expected Gaps #########
######################################

# expected_gaps_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/QC-SC_comp_small_L_check/hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-601)_Delta(0p05-0p05-1)/mp_gap_comp/mptarg-1.0_mptol0.01_idos0.05,0.004,99_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_rangenothing/qc_region_counts_log.csv"
# data_matrix, header_row = readdlm(expected_gaps_path, ',', header=true)
# expected_gaps_df = DataFrame(data_matrix, vec(header_row))
# println("Expected gaps data: ", names(expected_gaps_df))

# ## data variable for found mp_tol ranges
# struct MPTolResult
#     N::Int64
#     Delta::Float64
#     t_n::Tuple{Float64, Float64}
#     phi::Float64
#     phason::Float64
#     qc_region_count::Int64
#     mp_tol_ranges::Union{Nothing, Tuple{Float64, Float64}}
# end

# found_mp_tol_ranges = MPTolResult[]


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
    
    ## 0) Load the sweep data
    sweep_data = nothing
    try
        sweep_data = BSON.load(fpath)
    catch e
        @warn "Failed to load sweep data" fpath exception=e
        next!(p)
        continue
    end
    # println("type of sweep_data: ", typeof(sweep_data))
    # println("keys of sweep_data: ", keys(sweep_data))

    ## 1) Find mp_tol ranges that yield expected gap counts
        
    # # a) Load the row of expected_gaps_df matching the current parameters
    # expected_gaps_row = expected_gaps_df[(expected_gaps_df.N .== N_val) .&
    #                                      (abs.(expected_gaps_df.Delta .- Delta_val) .< 1e-10) .&
    #                                      (abs.(expected_gaps_df.phi .- phi_val) .< 1e-10) .&
    #                                      (abs.(expected_gaps_df.phason .- phason_val) .< 1e-10), :]
    # @info "Expected gaps row for N$(N_val), Delta$(Delta_val), phi$(phi_val), phason$(phason_val):" expected_gaps_row

    # qc_region_count = expected_gaps_row[1, :qc_region_count]
    # qc_region_count = Int(qc_region_count)

    qc_region_count = 5

    # ## b) Define mu_range limits by eigenvalues
    # max_eigenvalue = get(sweep_data, :max_eig, nothing)
    # mu_values = get(sweep_data, :mu, nothing)
    # min_mu = isnothing(mu_values) ? nothing : minimum(mu_values)
    # max_mu = isnothing(mu_values) ? nothing : maximum(mu_values)
    # if min_mu === nothing || max_eigenvalue === nothing
    #     @warn "min_mu or max_eigenvalue is nothing, skipping mu range restriction" fpath
    # elseif isapprox(min_mu, 0.0)
    #     mu_count_range = (0.0, max_eigenvalue)
    # else
    #     # mu_count_range = (-max_eigenvalue, max_eigenvalue)
    #     mu_count_range = (0.0, max_eigenvalue)
    # end
    # @info "Using mu range for gap count restriction: " mu_count_range
    
    # ## c) Analyze the sweep data to find mp_tol ranges yielding the expected gap count
    # highest_tol_range = nothing
    # try
    #     highest_tol_range = MPGapProcessing.find_mptol_ranges_for_expected_gaps(
    #         sweep_data,
    #         mu_count_range,
    #         qc_region_count
    #     )
    #     @info "mp_tol ranges yielding expected gap count $(qc_region_count):" highest_tol_range
    # catch e
    #     @warn "Failed to find mp_tol ranges for expected gaps" fpath exception=e
    # end

    # ## d) Save the mp_tol range with all parameters to 
    # push!(found_mp_tol_ranges, MPTolResult(
    #     N_val, 
    #     Delta_val, 
    #     t_n_tuple, 
    #     phi_val, 
    #     phason_val, 
    #     qc_region_count, 
    #     highest_tol_range
    # ))


    
    ## 2) Make mp_tol chec plots with computed ranges
    ## a) Construct a filename for the plot
    # plt_folder = joinpath(check_plots_dir) #, "N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phason$(phason_val)")
    plt_folder = check_plots_dir
    isdir(plt_folder) || mkpath(plt_folder)

    ## b) Generate mp_tol check heatmap OLD -- don't parse in precomputed mp_tol ranges
    # check_path = joinpath(plt_folder, check_fname)
    try
        name = "check_plt_phi$(phi_val)_N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phason$(phason_val).svg"
        folder = joinpath(plt_folder, "raw")
        isdir(folder) || mkpath(folder)
        MPGapPlotting.plt_mu_vs_mptol_check(
            sweep_data;
            N=N_val,
            mode=:full,
            savepath=joinpath(folder, name)
        )
    catch e
        @warn "Plotting failed" exception=e
    end

    try
        name = "boundary_curve_phi$(phi_val)_N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phason$(phason_val).svg"
        folder = joinpath(plt_folder, "boundary_curve")
        isdir(folder) || mkpath(folder)
        MPGapPlotting.plt_boundary_curve_with_peak_heights(
            sweep_data;
            N=N_val,
            savepath=joinpath(folder, name)
        )
    catch e
        @warn "Boundary curve plotting failed" exception=e
    end

    try
        name = "normalized_check_plt_phi$(phi_val)_N$(N_val)_Delta$(Delta_val)_t1$(t_n_tuple[1])_t2$(t_n_tuple[2])_phason$(phason_val).svg"
        folder = joinpath(plt_folder, "norm")
        isdir(folder) || mkpath(folder)
        MPGapPlotting.plt_mu_vs_mptol_check_normalized(
            sweep_data;
            N=N_val,
            mode=:full,
            savepath=joinpath(folder, name)
        )
    catch e
        @warn "Normalized check plotting failed" exception=e
    end

    # ## c) Generate mp_tol check heatmap NEW -- parse in precomputed mp_tol ranges
    # check_fname = "check_plt_phi$(phi_val).svg"
    # check_path = joinpath(plt_folder, check_fname)
    # try
    #     MPGapPlotting.plt_mu_vs_mptol_check_with_precomputed_ranges(
    #         sweep_data,
    #         highest_tol_range,
    #         mu_count_range,
    #         check_path
    #     )
    # catch e
    #     @warn "Plotting failed" check_path exception=e
    # end





    # # Export gap counts per tolerance for inspection
    # try
    #     mp_tols = get(sweep_data, :mp_tol_range, nothing)
    #     analysis = get(sweep_data, :analysis, nothing)
    #     gap_counts = (!isnothing(analysis) && haskey(analysis, :gap_counts)) ? analysis[:gap_counts] : nothing
    #     if !isnothing(mp_tols) && !isnothing(gap_counts)
    #         csv_fname = "gaps_vs_tol_phi$(phi_val).csv"
    #         csv_path = joinpath(plt_folder, csv_fname)
    #         open(csv_path, "w") do io
    #             println(io, "mp_tol,num_gaps")
    #             writedlm(io, hcat(mp_tols, gap_counts), ',')
    #         end
    #     else
    #         @warn "CSV skipped due to missing tolerance data" fpath
    #     end
    # catch e
    #     @warn "CSV generation failed" fpath exception=e
    # end

    # # Build dataframe slice for the requested parameter tuple
    # slice_mask = (df_base_raw[!, :N] .== N_val) .&
    #              (abs.(df_base_raw[!, :Delta] .- Delta_val) .< 1e-10) .&
    #              (abs.(df_base_raw[!, :phi] .- phi_val) .< 1e-10) .&
    #              (abs.(df_base_raw[!, :phason] .- phason_val) .< 1e-10)
    # slice_mask = BitVector(slice_mask)
    # tn_mask = map(x -> Tuple(x) == t_n_tuple, df_base_raw[!, :t_n])
    # slice_mask .&= BitVector(tn_mask)
    # df_slice = df_base_raw[slice_mask, :]

    # if nrow(df_slice) == 0
    #     @warn "No dataframe rows matched combination" N_val Delta_val phi_val phason_val t_n_tuple
    #     next!(p)
    #     continue
    # end

    # # Process eigenvalue analysis
    # try
    #     MPGapProcessing.process_eigenvalue_analysis!(df_slice, sweep_data; plot_mu_range_full=plot_mu_range_full)
    # catch e
    #     @warn "Eigenvalue analysis failed" fpath exception=e
    #     next!(p)
    #     continue
    # end

    # df_slice, sweep_data = MPGapProcessing.adjust_mu_range(df_slice, sweep_data, data_mu_range_full, plot_mu_range_full)
    # analysis = get(sweep_data, :analysis, nothing)
    # energy_gap_ranges_for_plot = (analysis !== nothing && haskey(analysis, :energy_gap_ranges)) ? analysis[:energy_gap_ranges] : nothing

    # projection_fname = "gap_projection_phi$(phi_val).svg"
    # projection_path = joinpath(plt_folder, "Delta_gap_criterion", projection_fname)
    # isdir(dirname(projection_path)) || mkpath(dirname(projection_path))
    # @debug "Generating energy-gap projection plot, saving to: " projection_path
    # try
    #     MPGapPlotting.plt_eigs_and_mptol_with_sc_gap_condition(
    #         df_slice,
    #         sweep_data,
    #         projection_path;
    #         plot_full_y_range=plot_mu_range_full,
    #         energy_gap_ranges_override=energy_gap_ranges_for_plot
    #     )
    # catch e
    #     @warn "Energy-gap projection plotting failed" projection_path exception=e
    # end

    # if isfile(projection_path)
    #     @debug "Saved gap projection plot" projection_path
    # else
    #     @warn "Gap projection plot not found after save attempt" projection_path
    # end

    next!(p)
end

df_base_raw = nothing
GC.gc()

# # Save found_mp_tol_ranges to a CSV for later analysis
# found_mp_tol_ranges_df = DataFrame(found_mp_tol_ranges)
# found_mp_tol_ranges_savepath = joinpath(sweep_results_dir, "found_mp_tol_ranges.csv")
# println("Saving found mp_tol ranges to CSV: ", found_mp_tol_ranges_savepath)
# CSV.write(found_mp_tol_ranges_savepath, found_mp_tol_ranges_df)

# # Also save to BSON for structured data preservation
# bson_savepath = joinpath(sweep_results_dir, "found_mp_tol_ranges.bson")
# println("Saving found mp_tol ranges to BSON: ", bson_savepath)
# BSON.bson(bson_savepath, Dict(:found_mp_tol_ranges => found_mp_tol_ranges))

println("Check plots saved to: ", check_plots_dir)