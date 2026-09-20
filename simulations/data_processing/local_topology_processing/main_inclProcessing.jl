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
data_folder = "slopes_N1000_literal_0.61803_phason_0.61803-1-0.61803_mapto[2,1]_N(89-89-1)_t1(0.0-1.0-101)_t2(1.0-1.0-1)_mu(0.01-0.01-1)_Delta(0.001-0.05-4)"
precision_str = "hp"
run_set_name = "local_mp_data"
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "simulations", "raw_data", run_set_name, precision_str, data_folder))
println("Unpacking Hofstadter data from folder: ", folder_path_hof)

raw_df = process_jld2_files_fluid(folder_path_hof)

# println(raw_df)


## (2) Calculate local MP scales
println("Calculating local MP scales...")
local_mp_results = Vector{Vector{NamedTuple}}(undef, nrow(raw_df))
@showprogress for i in 1:nrow(raw_df)
    eigenvec = raw_df[i, :eigenvectors]
    # local_mp_results[i] = LocMPScales.calc_local_scales_mp_optimised(eigenvec)
    local_mp_results[i] = LocMPScales.calc_local_scales_mp_allscales_optimised(eigenvec)
end
# mp_df = deepcopy(raw_df)
# mp_df[!, :mp_local_scales] = local_mp_results
println("Local MP scales calculation complete.")


## (3) Calculate Spectral function
println("Calculating spectral functions...")
println("WARNING: check <map_to> in SpecFunc.regen_sequence matches what you used in cut_and_project.jl to generate for solving!")
eta = 0.1
omega_range = collect(-3.0:0.1:3.0)
spec_func_df = SpecFunc.calc_spectral_func_for_all_optimised!(raw_df, eta, omega_range)
spec_func_results = spec_func_df[!, :spectral_function]
println("Spectral function calculation complete.")


## (4) Generate plots
println("Generating plots...")
plot_type = "spectral_and_mp_all_scales_precalc_test"
plt_filepath = normpath(joinpath(@__DIR__, "..", "..", "..", "simulations", "results", run_set_name, precision_str, data_folder, plot_type))
isdir(plt_filepath) || mkpath(plt_filepath)

spec_and_mp_df = select(raw_df, Not([:eigenvalues, :eigenvectors]))
spec_and_mp_df[!, :mp_local_scales] = local_mp_results
spec_and_mp_df[!, :spectral_function] = spec_func_results

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