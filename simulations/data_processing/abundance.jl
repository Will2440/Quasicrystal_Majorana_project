project_root = @__DIR__
include(joinpath(project_root, "bson_unpacker.jl"))

using DataFrames
using CSV
using LaTeXStrings
using Measures
using Plots
using Unzip
using Base.Threads
using ProgressMeter

###########################################################
######### Sec 1: Auxilliary Integral Functions ############
###########################################################

function trapz_2d(x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64})
    """
    Performs 2D trapezoidal integration over x and y with respect to z,
    normalized by the total range of x and y.

    Parameters:
    x -- Vector of x values (independent variable for the first dimension).
    y -- Vector of y values (independent variable for the second dimension).
    z -- Matrix of z values where integration is performed over both x and y.

    Returns:
    Float64 -- The normalized numerical integral of z with respect to x and y.
    """
    dx = x[end] - x[1]  # Total range of x
    dy = y[end] - y[1]  # Total range of y

    integral = 0.0  # Initialize the integral value

    # Perform 2D integration using the trapezoidal rule
    for i in 1:(length(x) - 1)
        for j in 1:(length(y) - 1)
            # Average values in the grid cell (trapezoidal rule)
            cell_value = (z[j, i] + z[j+1, i] + z[j, i+1] + z[j+1, i+1]) / 4
            # Multiply by cell area (dx * dy)
            cell_area = (x[i+1] - x[i]) * (y[j+1] - y[j])
            integral += cell_value * cell_area
        end
    end

    # Normalize the integral by the full parameter space
    return integral #/ (dx * dy)
end

function calc_integral_under_curve(df::DataFrame, value_col::Symbol, x_variable::Symbol, y_variable::Symbol; fixed_values...)
    """
    Calculates the integral under the curve defined by two parameters (e.g., mu_t and rho) with values of another column (e.g., mp_disc).

    Parameters:
    df           -- DataFrame containing the data.
    value_col    -- Symbol representing the column whose values will be integrated.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.

    Returns:
    Float64      -- The computed integral value.
    """

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)

        if !isnothing(x_index) && !isnothing(y_index)
            value_data = row[value_col]
            value_matrix[y_index, x_index] = real(value_data)
        end
    end

    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    # Compute the integral using the trapezoidal rule
    integral_value = trapz_2d(x_values, y_values, value_matrix)

    return integral_value
end

function calc_percentage_under_curve(df::DataFrame, condition_col::Symbol, x_variable::Symbol, y_variable::Symbol; fixed_values...)
    """
    Calculates the percentage of points satisfying a condition (e.g., mp_disc == -1.0) over the total points in the specified parameter space.

    Parameters:
    df           -- DataFrame containing the data.
    condition_col -- Symbol representing the column to check the condition against.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.

    Returns:
    Float64      -- The computed percentage value.
    """

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    total_points = nrow(filtered_df)
    if total_points == 0
        return 0.0
    end

    condition_points = nrow(filter(row -> row[condition_col] == -1.0, filtered_df))

    percentage_value = (condition_points / total_points) * 100.0

    return percentage_value
end

function save_integrals_to_csv(filepath::String, mp_integral::Float64, ipr_integral::Float64, mbs_integral::Float64)
    df = DataFrame(
        Integral_Type = ["mp_integral", "ipr_integral", "mbs_integral"],
        Value = [mp_integral, ipr_integral, mbs_integral]
    )
    
    CSV.write(filepath, df)
end

function save_percentages_to_csv(filepath::String, mp_percentage::Float64, ipr_percentage::Float64, mbs_percentage::Float64)
    df = DataFrame(
        Percentage_Type = ["mp_percentage", "ipr_percentage", "mbs_percentage"],
        Value = [mp_percentage, ipr_percentage, mbs_percentage]
    )
    
    CSV.write(filepath, df)
end



###########################################################
########### Sec 2: Main Abundance Processing ##############
###########################################################

function standard_abundance_processing(
    df::DataFrame,
    N::Int,
    Delta::Float64,
    seq_name::String,
    save_filepath::String;
)
    Delta_safe = replace(string(Delta), "." => "p")

    mp_integral = calc_integral_under_curve(df, :mp_disc, :mu_t, :rho, N=N, Delta_t=Delta)
    ipr_integral = calc_integral_under_curve(df, :ipr_masked, :mu_t, :rho, N=N, Delta_t=Delta)
    mbs_integral = calc_integral_under_curve(df, :maj_gap_masked, :mu_t, :rho, N=N, Delta_t=Delta)

    filename = "$(save_filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_integrals.csv"
    save_integrals_to_csv(filename, mp_integral, ipr_integral, mbs_integral)

end

function standard_abundance_processing_PQC(
    df::DataFrame,
    N::Int,
    Delta::Float64,
    seq_name::String,
    save_filepath::String;
    sigma::Float64 = 2.0
)
    Delta_safe = replace(string(Delta), "." => "p")
    sigma_safe = replace(string(sigma), "." => "p")

    mp_integral = calc_integral_under_curve(df, :mp_disc, :mu_t, :rho, N=N, Delta_t=Delta, sigma=sigma)
    ipr_integral = calc_integral_under_curve(df, :ipr_masked, :mu_t, :rho, N=N, Delta_t=Delta, sigma=sigma)
    mbs_integral = calc_integral_under_curve(df, :maj_gap_masked, :mu_t, :rho, N=N, Delta_t=Delta, sigma=sigma)

    filename = "$(save_filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_sigma$(sigma_safe)_integrals.csv"
    save_integrals_to_csv(filename, mp_integral, ipr_integral, mbs_integral)

end

function percentage_abundance_processing(
    df::DataFrame,
    N::Int,
    Delta::Float64,
    seq_name::String,
    save_filepath::String;
)
    Delta_safe = replace(string(Delta), "." => "p")

    mp_percentage = calc_percentage_under_curve(df, :mp_disc, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    ipr_percentage = calc_percentage_under_curve(df, :ipr_masked, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    mbs_percentage = calc_percentage_under_curve(df, :maj_gap_masked, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    filename = "$(save_filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_percentage.csv"
    save_percentages_to_csv(filename, mp_percentage, ipr_percentage, mbs_percentage)

end



###########################################################
#################### Sec 3: Plotting ######################
###########################################################

function combine_integrals(folder_path::String)
    file_paths = glob("*integrals.csv", folder_path)
    combined_rows = DataFrame(Delta=Float64[], mp_integral=Float64[], ipr_integral=Float64[], mbs_integral=Float64[])

    for file_path in file_paths
        filename = split(file_path, "/") |> last
        delta_match = match(r"(\d+p\d+)", filename)
        if delta_match === nothing
            println("Could not parse Delta from filename: $filename")
            continue
        end
        delta_str = replace(delta_match.match, "p" => ".")
        delta = parse(Float64, delta_str)

        df = CSV.read(file_path, DataFrame; delim=",", header=true)
        row = Dict(row.Integral_Type => row.Value for row in eachrow(df))

        push!(combined_rows, (
            Delta = delta,
            mp_integral = abs(row["mp_integral"]),
            ipr_integral = row["ipr_integral"],
            mbs_integral = row["mbs_integral"]
        ))
    end

    sort!(combined_rows, :Delta)

    return combined_rows
end

function combine_percentages(folder_path::String)
    file_paths = glob("*percentage.csv", folder_path)
    combined_rows = DataFrame(Delta=Float64[], mp_percentage=Float64[], ipr_percentage=Float64[], mbs_percentage=Float64[])

    for file_path in file_paths
        filename = split(file_path, "/") |> last
        delta_match = match(r"(\d+p\d+)", filename)
        if delta_match === nothing
            println("Could not parse Delta from filename: $filename")
            continue
        end
        delta_str = replace(delta_match.match, "p" => ".")
        delta = parse(Float64, delta_str)

        df = CSV.read(file_path, DataFrame; delim=",", header=true)
        row = Dict(row.Percentage_Type => row.Value for row in eachrow(df))

        push!(combined_rows, (
            Delta = delta,
            mp_percentage = abs(row["mp_percentage"]),
            ipr_percentage = row["ipr_percentage"],
            mbs_percentage = row["mbs_percentage"]
        ))
    end

    sort!(combined_rows, :Delta)

    return combined_rows
end

function combine_integrals_PQC(folder_path::String)
    file_paths = glob("*integrals.csv", folder_path)
    combined_rows = DataFrame(Delta=Float64[], sigma=Float64[], mp_integral=Float64[], ipr_integral=Float64[], mbs_integral=Float64[])

    for file_path in file_paths
        filename = split(file_path, "/") |> last
        delta_match = match(r"Delta(\d+p\d+)", filename)
        sigma_match = match(r"sigma(\d+p\d+)", filename)
        if delta_match === nothing || sigma_match === nothing
            println("Could not parse Delta or sigma from filename: $filename")
            continue
        end
        delta_str = replace(delta_match.captures[1], "p" => ".")
        sigma_str = replace(sigma_match.captures[1], "p" => ".")
        delta = parse(Float64, delta_str)
        sigma = parse(Float64, sigma_str)

        df = CSV.read(file_path, DataFrame; delim=",", header=true)
        row = Dict(row.Integral_Type => row.Value for row in eachrow(df))

        push!(combined_rows, (
            Delta = delta,
            sigma = sigma,
            mp_integral = abs(row["mp_integral"]),
            ipr_integral = row["ipr_integral"],
            mbs_integral = row["mbs_integral"]
        ))
    end

    sort!(combined_rows, [:Delta, :sigma])

    return combined_rows
end

function plot_multiple_mp_integrals(dfs::Vector{DataFrame}, labels::Vector{LaTeXString}, save_path::String, data_type::Symbol;
                                    markers::Vector{Symbol} = repeat([:o], length(dfs)),
                                    colors::Vector{Symbol} = nothing,
                                    N::Int = 0,
                                    window::String = ""
)

    plt = Plots.plot(
        xlabel = L"\Delta / t_A",
        ylabel = L"\textnormal{Abundance}",
        title = L"\textnormal{Abundance for } N =" * string(N) * ", " * L"\textnormal{ window range: }" * window,
        legend = :bottomright, #:outerbottomleft,
        grid = false,
        size = (800, 600),
        labelfontsize = 20,
        tickfontsize = 12,
        legendfontsize = 20,
        left_margin = 5mm,
        bottom_margin = 5mm,
        # xlims = (0.0, 2.0),
        # ylims = (-0.01, 35.0)
    )

    if data_type == :integral
        for (i, df) in enumerate(dfs)
            filtered = [(Δ, mp) for (Δ, mp) in zip(df.Delta, df.mp_integral) if Δ <= 2.0]
            delta_filtered, mp_filtered = unzip(filtered)

            Plots.plot!(
                plt,
                delta_filtered, #df.Delta,
                mp_filtered, #df.mp_integral,
                # df.mbs_integral,
                label = labels[i],
                lw = 2,
                marker = markers[i],
                # color = colors === nothing ? nothing : colors[i]
            )
        end
    elseif data_type == :percentage
        for (i, df) in enumerate(dfs)
            filtered = [(Δ, mp) for (Δ, mp) in zip(df.Delta, df.mp_percentage) if Δ <= 2.0]
            delta_filtered, mp_filtered = unzip(filtered)

            Plots.plot!(
                plt,
                delta_filtered, #df.Delta,
                mp_filtered, #df.mp_integral,
                # df.mbs_integral,
                label = labels[i],
                lw = 2,
                marker = markers[i],
                # color = colors === nothing ? nothing : colors[i]
            )
        end
    end

    Plots.savefig(save_path)
end


function plot_multiple_sigmas(df::DataFrame, save_path::String;
                              metric::Symbol = :mp,
                              x_max::Float64 = 2.0,
                              xlabel = L"\Delta / t_A",
                              ylabel = L"\textnormal{Abundance}",
                              title = L"\textnormal{Abundance vs } \Delta",
                              markers = nothing,
                              colors = nothing,
                              lw::Real = 2,
                              legendpos = :bottomright
)

    metric_map = Dict(
        :mp  => :mp_integral,
        :ipr => :ipr_integral,
        :mbs => :mbs_integral
    )
    haskey(metric_map, metric) || error("metric must be one of $(collect(keys(metric_map)))")

    col = metric_map[metric]

    unique_sigmas = sort(unique(df.sigma))
    n = length(unique_sigmas)

    if markers === nothing
        default_markers = [:o, :s, :d, :x, :utriangle, :star5, :hex, :+, :circle, :square]
        markers = [default_markers[mod1(i, length(default_markers))] for i in 1:n]
    end
    if colors === nothing
        palette = Plots.palette(:viridis, n)
        colors = palette
    end

    plt = Plots.plot(
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        legend = legendpos,
        grid = false,
        size = (1200, 600),
        labelfontsize = 20,
        tickfontsize = 12,
        legendfontsize = 16,
        left_margin = 5mm,
        bottom_margin = 5mm
    )

    for (i, σ) in enumerate(unique_sigmas)
        sub = df[df.sigma .== σ, :]
        sub = sub[sub.Delta .<= x_max, :]
        sort!(sub, :Delta)
        Plots.plot!(
            plt,
            sub.Delta,
            sub[!, col],
            label = L"\sigma = " * string(σ),
            lw = lw,
            marker = markers[i],
            color = colors[i]
        )
    end

    Plots.savefig(save_path)
end



##########################################################
############## Sec 4: Run Calculations ###################
##########################################################

 # Run basic processing once on single folder of data
data_folder = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance_data/PQC_N(50-50-1)_t1(1.0-1.0-1_t2(0.0-10.0-101)_t3(0.001-100.0-17)_mu(0.0-10.0-101)_Delta(0.1-2.0-20)" 
# data_folder = joinpath(project_root, "../../simulations/raw_data/np/test_run/GQC_N(50-50-1)_t1(1.0-1.0-1_t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(2.0-2.0-1)/")

# mp_tol = 0.01
# df = unpack_bason_standard(data_folder; mp_tol=mp_tol)

# # # fudge mp_disc back to -1.0 after topological region was cut in restricted solver
# # df[findall(row -> row[:mu] <= 1.0, eachrow(df)), :mp_disc] .= -1.0

# N = 50
# # Delta = 2.0
# seq_name = "PQC"

# save_integral_filepath = joinpath(project_root, "../../simulations/results/mu_rho_integrals/GQC_N50_mu0-40_rho0-40_mptol0.1/")
save_integral_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/extreme_sigma_test/PQC_N(50-50-1)_t1(1.0-1.0-1_t2(0.0-10.0-101)_t3(0.001-100.0-17)_mu(0.0-10.0-101)_Delta(0.1-2.0-20)_mptol0.01/"
isdir(save_integral_filepath) || mkpath(save_integral_filepath)

# # standard_abundance_processing(df, N, Delta, seq_name, save_integral_filepath)
# sigma_range = [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 100.0]
# Delta_range = collect(range(0.1, 2.0, 20))
# @showprogress Threads.@threads for Delta in Delta_range
#     for sigma in sigma_range  
#         standard_abundance_processing_PQC(df, N, Delta, seq_name, save_integral_filepath; sigma=sigma)
#     end
# end

all_sig_df = combine_integrals_PQC(save_integral_filepath)

plot_multiple_sigmas(all_sig_df, joinpath(save_integral_filepath, "PQC_N$(N)_abundance_vs_Delta_all_sigmas.png");
                     metric=:mp,
                     x_max=2.0,
                     xlabel=L"\Delta / t_A",
                     ylabel=L"\textnormal{Abundance}",
                     title=L"\textnormal{Abundance vs } \Delta \textnormal{ for PQC, } N=50",
                     legendpos=:outerbottomleft
)






# # Run basic processing on multiple folders of data
# mp_tol = 0.1
# N=100

# root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/"

# NC_data_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/NC_N(100-100-1)_t1(1.0-1.0-1)_t2(0.0-5.0-51)_mu(0.0-5.0-51)_Delta(0.1-2.0-20)"
# df_NC = unpack_bason_standard(NC_data_filepath; mp_tol=mp_tol)
# GQC_data_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/GQC_N(100-100-1)_t1(1.0-1.0-1)_t2(0.0-5.0-51)_mu(0.0-5.0-51)_Delta(0.1-2.0-20)"
# df_GQC = unpack_bason_standard(GQC_data_filepath; mp_tol=mp_tol)
# SQC_data_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/SQC_N(100-100-1)_t1(1.0-1.0-1)_t2(0.0-5.0-51)_mu(0.0-5.0-51)_Delta(0.1-2.0-20)"
# df_SQC = unpack_bason_standard(SQC_data_filepath; mp_tol=mp_tol)
# TMQC_data_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/TMQC_N(100-100-1)_t1(1.0-1.0-1)_t2(0.0-5.0-51)_mu(0.0-5.0-51)_Delta(0.1-2.0-20)"
# df_TMQC = unpack_bason_standard(TMQC_data_filepath; mp_tol=mp_tol)
# PQC_sig3_data_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/PQC_N(100-100-1)_t1(1.0-1.0-1)_t2(0.0-5.0-51)_mu(0.0-5.0-51)_Delta(0.1-2.0-20)"
# df_PQC_sig3 = unpack_bason_standard(PQC_sig3_data_filepath; mp_tol=mp_tol)

# raw_dfs = [df_NC, df_GQC, df_SQC, df_TMQC, df_PQC_sig3]
# seq_names = ["NC", "GQC", "SQC", "TMQC", "PQC_sig3"]
# window = "5"
# save_integral_filepaths = [
#     "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N$(N)/windowrange_$(window)/NC_N100_mu0-10_rho0-10/",
#     "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N$(N)/windowrange_$(window)/GQC_N100_mu0-10_rho0-10/",
#     "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N$(N)/windowrange_$(window)/SQC_N100_mu0-10_rho0-10/",
#     "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N$(N)/windowrange_$(window)/TMQC_N100_mu0-10_rho0-10/",
#     "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N$(N)/windowrange_$(window)/PQC_N100_sig3_mu0-10_rho0-10/"
# ]

# for (df, seq_name, save_integral_filepath) in zip(raw_dfs, seq_names, save_integral_filepaths)
#     isdir(save_integral_filepath) || mkpath(save_integral_filepath)

#     Delta_range = collect(range(0.1, 2.0, 20))
#     @showprogress Threads.@threads for Delta in Delta_range    
#         standard_abundance_processing(df, N, Delta, seq_name, save_integral_filepath)
#         percentage_abundance_processing(df, N, Delta, seq_name, save_integral_filepath)
#     end
# end


# ###########################################################
# ################## Sec 4: Run Plotting ####################
# ###########################################################
# # NC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/NC_N50_mu0-10_rho0-10/"
# df_NC_p = combine_percentages(save_integral_filepaths[1])
# # GQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/GQC_N50_mu0-10_rho0-10/"
# df_GQC_p = combine_percentages(save_integral_filepaths[2])
# # SQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/SQC_N50_mu0-10_rho0-10/"
# df_SQC_p = combine_percentages(save_integral_filepaths[3])
# # TMQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/TMQC_N50_mu0-10_rho0-10/"
# df_TMQC_p = combine_percentages(save_integral_filepaths[4])
# # PQC_sig3_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/PQC_N50_sig3_mu0-10_rho0-10/"
# df_PQC_sig3_p = combine_percentages(save_integral_filepaths[5])

# # NC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/NC_N50_mu0-10_rho0-10/"
# df_NC_i = combine_integrals(save_integral_filepaths[1])
# # GQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/GQC_N50_mu0-10_rho0-10/"
# df_GQC_i = combine_integrals(save_integral_filepaths[2])
# # SQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/SQC_N50_mu0-10_rho0-10/"
# df_SQC_i = combine_integrals(save_integral_filepaths[3])
# # TMQC_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/TMQC_N50_mu0-10_rho0-10/"
# df_TMQC_i = combine_integrals(save_integral_filepaths[4])
# # PQC_sig3_integrals_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/PQC_N50_sig3_mu0-10_rho0-10/"
# df_PQC_sig3_i = combine_integrals(save_integral_filepaths[5])


# # window = "5"
# final_plots_savepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/mu_rho_integrals/abundance_plots/N100/windowrange_$(window)/"

# plot_multiple_mp_integrals([df], #[df_NC_p, df_GQC_p, df_SQC_p, df_TMQC_p, df_PQC_sig3_p], 
#     [L"\textnormal{NC}", L"\textnormal{GQC}", L"\textnormal{SQC}", L"\textnormal{TMQC}", L"\textnormal{PQC} \ \sigma=3.0"],# L"\textnormal{PQC} \ \sigma=4.0"], 
#     "$(final_plots_savepath)all_QCs_percentage_abundance_N$(N)_range$(window)_mptol$(mp_tol).png",
#     :percentage;
#     markers = [:o, :s, :d, :o, :x], 
#     colors = [:blue, :red, :green, :orange, :purple],
#     N = N,
#     window = window
# )



# plot_multiple_mp_integrals([df_NC_i, df_GQC_i, df_SQC_i, df_TMQC_i, df_PQC_sig3_i], 
#     [L"\textnormal{NC}", L"\textnormal{GQC}", L"\textnormal{SQC}", L"\textnormal{TMQC}", L"\textnormal{PQC} \ \sigma=3.0"],# L"\textnormal{PQC} \ \sigma=4.0"], 
#     "$(final_plots_savepath)all_QCs_integral_abundance_N$(N)_range$(window)_mptol$(mp_tol).png",
#     :integral;
#     markers = [:o, :s, :d, :o, :x], 
#     colors = [:blue, :red, :green, :orange, :purple],
#     N = N,
#     window = window
# )

