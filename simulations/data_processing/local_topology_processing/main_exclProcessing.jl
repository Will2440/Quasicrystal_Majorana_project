using DataFrames
# using BSON
using JLD2
using Base.Threads
using ProgressMeter
using Logging

# include("../bson_unpacker.jl")
include("../jld2_unpacker.jl")
include("local_mp_scales.jl")
include("spectral_function.jl")
include("plotting.jl")

using .LocMPScales
using .SpecFunc
using .Plotting


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


println("--------------------------------------------------------------")
println("Starting data processing for local topology analysis...")

## (1) Unpack data
data_folder = "slopes_N1000_literal_0.61803_phason_0.61803-1-0.61803_N(610-610-1)_t1(1p0-1p0-1)_t2(0p0-1p0-101)_mu(0p01-0p01-1)_Delta(0p01-0p01-1)_prec512"
precision_str = "hp"
run_set_name = "local_mp_data"
# folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "simulations", "raw_data", run_set_name, precision_str, data_folder))
# println("Unpacking Hofstadter data from folder: ", folder_path_hof)

folder_path_hof = "/Volumes/External Drive/Quasicrystal_Majorana_Project/PRB_test_results/simulations/raw_data/local_mp_data/hp/slopes_N1000_literal_0.61803_phason_0.61803-1-0.61803_N(610-610-1)_t1(1p0-1p0-1)_t2(0p0-1p0-101)_mu(0p01-0p01-1)_Delta(0p01-0p01-1)_prec512"
println("Unpacking Hofstadter data from external drive path: ", folder_path_hof)


raw_df = process_jld2_files_fluid(folder_path_hof)

omega_range = collect(-3.0:0.1:3.0)
println("WARNING: check <map_to> in SpecFunc.regen_sequence matches what you used in cut_and_project.jl to generate for solving!")
println("WARNING: check the omega_range ($(omega_range[1]):$(round(omega_range[2]-omega_range[1], digits=2)):$(omega_range[end])) for spectral function calculation in main.jl match what you use here for processing/plotting!")

## (4) Generate plots
println("Generating plots...")
plot_type = "spectral_and_mp_all_scales_precalc_test"
plt_filepath = normpath(joinpath(@__DIR__, "..", "..", "..", "simulations", "results", run_set_name, precision_str, data_folder, plot_type))
isdir(plt_filepath) || mkpath(plt_filepath)

spec_and_mp_df = select(raw_df, Not([:eigenvalues, :eigenvectors]))
## TEMP: rename df columns to match what plotting function expects
# colmap = Dict(
#     "spec_func" => :spectral_function,
#     :spec_func => :spectral_function,
#     "loc_scales_mp" => :mp_local_scales,
#     :loc_scales_mp => :mp_local_scales
# )

# for old in keys(colmap)
#     if old in names(spec_and_mp_df) || old in propertynames(spec_and_mp_df)
#         rename!(spec_and_mp_df, old => colmap[old])
#     end
# end

# println("NEW DATAFRAME NAMES: ", names(spec_and_mp_df))

Delta_range = unique(spec_and_mp_df[!, :Delta])
t_n_range = unique(spec_and_mp_df[!, :t_n])
mu_range = unique(spec_and_mp_df[!, :mu])
N_range = unique(spec_and_mp_df[!, :N])
slope_range = unique(getindex.(spec_and_mp_df[!, :sequence_id], 1))
phason_range = unique(getindex.(spec_and_mp_df[!, :sequence_id], 2))

println("Parameter ranges for plotting:")
println("Delta_range: ", Delta_range)
println("t_n_range: ", t_n_range)
println("mu_range: ", mu_range)
println("N_range: ", N_range)
println("slope_range: ", slope_range)
println("phason_range: ", phason_range)

total_plots = length(Delta_range) * length(t_n_range) * length(mu_range) * length(N_range) * length(slope_range) * length(phason_range)
println("Total plots to generate: ", total_plots)

@showprogress for Delta in Delta_range, t_n in t_n_range, mu in mu_range, N in N_range, slope in slope_range, phason in phason_range
    plot_filename = "spectral_and_mp_scales_Delta$(Delta)_tn$(t_n)_mu$(mu)_N$(N)_slope$(slope)_phason$(phason).pdf"
    Plotting.plt_spectral_and_mp_scales(
        spec_and_mp_df,
        omega_range,
        joinpath(plt_filepath, plot_filename);
        Delta=Delta, 
        t_n=t_n, 
        mu=mu, 
        N=N, 
        sequence_id=(slope, phason)
    )
end

println("Plotting complete. Plots saved to: ", plt_filepath)


println("Re-arranging plots into subfolders...")
for Delta in Delta_range
    delta_folder = joinpath(plt_filepath, "Delta_$(Delta)")
    isdir(delta_folder) || mkpath(delta_folder)
    for file in filter(f -> occursin("Delta$(Delta)", f), readdir(plt_filepath))
        mv(joinpath(plt_filepath, file), joinpath(delta_folder, file))
    end
end


## (5) Plot IPR 