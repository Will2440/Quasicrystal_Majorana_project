"""
    file name:   bson_unpacker.jl
    created:     24/09/2025
    last edited: 25/09/2025

    overview:


    structure:
        - Section 1: .bson upacker into dataframe -- ready to use on the data fodlers provided
        - Section 2: Generalised plotting functions, can be used individualy, or within the following function in Sec 3.
        - Section 3: Additional code to calculate the abundance integral and ssave this data
        - Section 4: Contains a single function which can be run on the dataframe created in Sec 1 to generate all the standard plots seen in the data folders
"""


using BSON
using DataFrames
using Glob
using Plots
using CSV
using LsqFit
# using PyPlot
using LaTeXStrings
using LinearAlgebra
using PyPlot
using ProgressMeter
using Measures
using Statistics
using Base.Threads


##############################################################################
# # Section 1 # #
##############################################################################
# # helpful functions to unpack the .bson files contained within a data folder, formatting into a dataframe

function calc_rho(t_n::Vector{Float64})
    rho = t_n[2] / t_n[1]
    return rho
end

function normalise_to_t(t_n::Vector{Float64})
    norm = t_n[1]
    return norm
end

function bson_to_dataframe(file_path::String)

    bson_data = BSON.load(file_path)
    dataframe = DataFrame(bson_data)
    
    return dataframe
end

function process_bson_files(folder_path::String)
    # Initialize an empty DataFrame with the specified structure
    combined_dataframe = DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        # sequence_name = String[],
        mp = Float64[],
        ipr = Float64[],
        maj_gap = Float64[],
        eigenvalues = Vector{Float64}[]
        # mbs_gap = Float64[],
        # loc_len = Float64[]
    )
    
    file_paths = glob("*.bson", folder_path)
    files_read = 0
    for file_path in file_paths
        bson_data = BSON.load(file_path)
        if haskey(bson_data, :results_df)
            # Extract the DataFrame stored in the `results_df` key
            results_df = DataFrame(bson_data[:results_df])
            
            append!(combined_dataframe, results_df)
            files_read += 1
        end
    end
    println("number of files read: $files_read")
    
    return combined_dataframe
end

function calc_norms_df!(df::DataFrame)
    df.rho = [row.t_n[2] / row.t_n[1] for row in eachrow(df)]
    df.mu_t = [row.mu / row.t_n[1] for row in eachrow(df)]
    df.Delta_t = [row.Delta / row.t_n[1] for row in eachrow(df)]

    return df
end

function calc_mp_disc!(df::DataFrame, tol::Float64)
    df.mp_disc = [isapprox(row.mp, -1.0, atol=tol) ? -1.0 : 0.0 for row in eachrow(df)]
    return df
end

function mask_df!(df::DataFrame)
    df.mbs_gap_disc = [row.mp_disc == -1.0 ? row.mbs_gap : 0.0 for row in eachrow(df)]
    # df.ipr_disc = [row.mp_disc == -1.0 ? row.ipr : 0.0 for row in eachrow(df)]
    # df.loc_len_disc = [row.loc_len_disc == -1.0 ? row.loc_len : 0.0 for row in eachrow(df)]
    
    return df
end



# # # Example usage
folder_path = "/Users/Will/Documents/FINAL_PROJECT/simulations/phase_transition_fractality/data/hp/GQC_N55-610-6_g8-13_mu0.0-2.5-501_rho1.5-1.5-1_Delta0.05"

# for file in glob("*.bson", folder_path)
#     df = DataFrame(BSON.load(file))
#     println(names(df))
# end

combined_df = process_bson_files(folder_path)
norm_df = calc_norms_df!(combined_df)
tol = 1e-1
disc_mp_df = calc_mp_disc!(norm_df, tol)
# df = mask_df!(disc_mp_df)
df = disc_mp_df

# # df is now the final dataframe containing all raw data

# test_file = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/NC_N50_mu0.0-5.0-501_rho0.0-5.0-501_Delta0.1,0.2,0.3,0.5/NC_N50_tn1.0_0.1_Delta0.1_mu0.0-5.0_501.bson"
# # test_df = bson_to_dataframe(test_file)
# # println(stat(test_file).size)

# bson_data = BSON.load(test_file)
# println("Keys in file ", test_file, ": ", keys(bson_data))




##############################################################################
# # Section 2 # #
##############################################################################
# # Following functions make plots of the raw data, they are designed to be used for any variable name which are the key names of the dataframe

function plt_eigenvalues_generalised(df::DataFrame, x_variable::Symbol; fixed_values...)
    """
    Plots eigenvalues as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    plt = Plots.scatter(
        legend=false, 
        xlabel=string(x_variable), 
        ylabel="E", 
        grid=true,
        size=(800,600)
    )

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = filtered_df[:, x_variable]
    all_eigenvalues = filtered_df[:, :eigenvalues]

    num_eigenvalues = length(all_eigenvalues[1])

    for i in 1:num_eigenvalues
        eigenvalue_line = [eigenvalues[i] for eigenvalues in all_eigenvalues]
        plot!(x_values, eigenvalue_line)
    end

    Plots.display(plt)
end

function plt_mp_generalised(df::DataFrame, x_variable::Symbol, filename::String; fixed_values...)
    """
    Plots Majorana Polarization (MP) as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    plt = Plots.plot(
        legend=false, 
        xlabel=string(x_variable), 
        ylabel="MP", 
        grid=true,
        size=(800,600)
    )

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = filtered_df[:, x_variable]
    mp_values = [real(mp[floor(Int, length(mp) / 2) + 1]) for mp in filtered_df[:, :mp]]

    sort_indices = sortperm(x_values)
    x_values = x_values[sort_indices]
    mp_values = mp_values[sort_indices]

    Plots.plot!(x_values, mp_values, line=:solid) #marker=:circle,

    Plots.savefig(filename)
end

function plt_mp_heatmap(df::DataFrame, colour_var::Symbol, x_variable::Symbol, y_variable::Symbol, filename::String; fixed_values...)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
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
            value = row[colour_var]  # Already Vector{ComplexF64}
    
            # # Ensure it's non-empty and valid
            # if isa(mp_data, Vector{ComplexF64}) && !isempty(mp_data)
            #     mid_index = div(length(mp_data), 2) + 1
            #     mp_matrix[y_index, x_index] = real(mp_data[mid_index])
            # else
            #     mp_matrix[y_index, x_index] = NaN
            # end
            value_matrix[y_index, x_index] = real(value)
        end
    end

    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    plt = Plots.heatmap(
        x_values, 
        y_values, 
        value_matrix, 
        xlabel=string(x_variable), 
        ylabel=string(y_variable), 
        color=:viridis,
        clims=(minimum(value_matrix), maximum(value_matrix)), # (0.0, 1.0),
        colorbar=true,
        grid=false,
        size=(800,600)
    )
    
    # display(plt)
    Plots.savefig(filename)
end

function plt_mp_heatmap_with_contour(
    df::DataFrame, 
    colour_var::Symbol, 
    x_variable::Symbol, 
    y_variable::Symbol, 
    filename::String,
    tol::Float64;
    fixed_values...
)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables,
    and adds a single contour line where the value crosses from < -1 + tol to > -1 + tol.

    Parameters:
    df           -- DataFrame containing the data.
    colour_var   -- Symbol for the variable used for heatmap coloring.
    x_variable   -- Symbol for the variable to plot on the x-axis.
    y_variable   -- Symbol for the variable to plot on the y-axis.
    filename     -- String: Path to save the plot.
    tol          -- Float64: Tolerance used for detecting crossings (default 1e-2).
    fixed_values -- Keyword arguments specifying fixed values for filtering the DataFrame.
    """
    # # Filter DataFrame for fixed values
    # filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    
    filtered_df = filter(
    row -> all(getproperty(row, key) == value for (key, value) in fixed_values) && getproperty(row, x_variable) != 0.0, df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # Extract and sort unique x and y values
    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    # Initialize value matrix
    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    # Populate the value matrix
    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)
        if !isnothing(x_index) && !isnothing(y_index)
            value_matrix[y_index, x_index] = real(row[colour_var])
        end
    end

    # Ensure data is sorted properly
    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    # Define the level for the contour line
    contour_level = -1 + tol

    # Create the heatmap
    heatmap_plot = Plots.heatmap(
        x_values, 
        y_values, 
        value_matrix, 
        xlabel=L"\mu/t_A", #string(x_variable), 
        ylabel=L"\rho", #string(y_variable), 
        color=:viridis,
        labelfontsize = 35,
        tickfontsize = 24,
        left_margin=10mm,
        bottom_margin=5mm,
        right_margin=10mm,
        top_margin=5mm,
        clims=(minimum(value_matrix), maximum(value_matrix)),
        colorbar=true,
        grid=false,
        size=(1600, 1200)
    )

    # Add the contour line for the level
    Plots.contour!(
        x_values, 
        y_values, 
        value_matrix, 
        levels=[contour_level], 
        color=:red, 
        lw=2.0, 
        legend=false
    )

    # Save the plot to the specified file
    Plots.savefig(heatmap_plot, filename)
end

function plt_mp_heatmap_mbs_disc_final(
    df::DataFrame, 
    colour_var::Symbol, 
    x_variable::Symbol, 
    y_variable::Symbol, 
    filename::String; 
    fixed_values...
)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables.
    Values equal to 0.0 are shown as grey, while non-zero values are plotted in color.

    Parameters:
    df           -- DataFrame containing the data.
    colour_var   -- Symbol representing the variable used for coloring the heatmap.
    x_variable   -- Symbol for the variable to plot on the x-axis.
    y_variable   -- Symbol for the variable to plot on the y-axis.
    filename     -- String: Path to save the heatmap plot.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    # Filter the DataFrame for the fixed values
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # Extract unique x and y values, sorted
    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    # Initialize the value matrix
    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)
        if !isnothing(x_index) && !isnothing(y_index)
            value = real(row[colour_var])
            # Set 0.0 values to NaN to exclude them from the color map
            value_matrix[y_index, x_index] = value == 0.0 ? NaN : value
        end
    end

    # Ensure y_values and value_matrix are sorted consistently
    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    # Define the colormap and the background color for NaN (grey)
    heatmap_plot = Plots.heatmap(
        x_values,
        y_values,
        value_matrix,
        xlabel=L"\mu/t_A",
        ylabel=L"\rho",
        color=:plasma,
        nan_color=:grey,           # Set NaN (0.0 values) to be displayed as grey
        bg_inside=:grey,
        labelfontsize=35,
        tickfontsize=24,
        left_margin=10mm,
        bottom_margin=5mm,
        right_margin=10mm,
        top_margin=5mm,
        clims=(0.0, 2.0), #(minimum(skipmissing(value_matrix)), maximum(skipmissing(value_matrix))),
        colorbar=true,
        grid=false,
        size=(1600, 1200)
    )

    # Save the plot
    Plots.savefig(heatmap_plot, filename)
end





##############################################################################
# # Section 3 # #
##############################################################################
# # This contains functions for further processing of the data: abundance processing and quantity saving.

function trapz_2d(x::Vector{Float64}, y::Matrix{Float64})
    """
    Performs trapezoidal integration on a matrix y with respect to x.

    Parameters:
    x -- Vector of x values (independent variable)
    y -- Matrix of y values where integration is performed along the first dimension.

    Returns:
    integral -- The numerical integral of y with respect to x.
    """
    integral = zeros(size(y, 2))
    for j in 1:size(y, 2)
        for i in 1:(length(x) - 1)
            integral[j] += (x[i+1] - x[i]) * (y[i+1, j] + y[i, j]) / 2
        end
    end
    return sum(integral)
end

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

function exponential_decay(x::Vector{Float64}, p::Vector{Float64})
    return p[1] .* exp.(-p[2] .* x) .+ p[3]
end

function one_over_N(x::Vector{Float64}, p::Vector{Float64})
    return p[1] ./ (x .+ p[2]) .+ p[3]
end

function fit_curve(x_values::Vector{Int}, y_values::Vector{Float64}, fit_type::String)
    x_values_float = Float64.(x_values)
    if fit_type == "exp"
        model = exponential_decay
        p0 = [maximum(y_values), 0.1, minimum(y_values)]
    elseif fit_type == "1/N"
        model = one_over_N
        p0 = [maximum(y_values), 1.0, minimum(y_values)]
    else
        return nothing
    end

    fit = LsqFit.curve_fit(model, x_values_float, y_values, p0)
    return fit, model
end

function hp_ipr_plot_multi_mu_fit(
    df::DataFrame,
    fit_types::Dict{Float64, String},  # fit_types: μ => "none" | "exp" | "1/N"
    x_variable::Symbol,  # Variable to plot on x-axis (e.g., :N)
    dpi::Int,
    filename::String;
    fixed_values...
)
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    mu_groups = Dict{Float64, Vector{Tuple{Int64, Float64}}}()
    for row in eachrow(filtered_df)
        N = row[:N]
        mu = row[:mu]
        ipr_value = row[x_variable]  # Directly assign the Float64 value

        if !haskey(mu_groups, mu)
            mu_groups[mu] = Vector{Tuple{Int64, Float64}}()
        end
        push!(mu_groups[mu], (N, ipr_value))
    end

    fig, ax = subplots(dpi=dpi)
    # ax.set_title("(QGC) High Precision IPR Values vs N for Different μ")
    ax.set_xlabel(L"L")
    ax.set_ylabel("IPR (a.u)")
    ax.grid(false)

    marker_styles = ["o", "s", "D", "^", "v", "p", "*", "h", "x", "+"]  # Define a set of marker styles
    color_map = Dict()  # Store colors for each `μ`
    
    # Loop through the μ groups
    for (idx, (mu, entries)) in enumerate(sort(collect(mu_groups), by = x -> x[1]))

        sorted_entries = sort(entries, by = x -> x[1])
        sorted_N_values = [entry[1] for entry in sorted_entries]
        sorted_ipr_values = [entry[2] for entry in sorted_entries]
    
        # Select the marker style based on the index
        marker_style = marker_styles[mod1(idx, length(marker_styles))]
    
        # Scatter plot with the selected marker style and store the color
        scatter_handle = ax.scatter(sorted_N_values, sorted_ipr_values, label="μ = $(round(mu, digits = 3))", marker=marker_style)
        color_map[mu] = scatter_handle.get_edgecolor()  # Extract and store the color of the scatter points
    
        # Fit curve if applicable
        fit_type = get(fit_types, mu, "none")
        if fit_type != "none"
            fit_result, model = fit_curve(sorted_N_values, sorted_ipr_values, fit_type)
            if fit_result !== nothing
                fit_params = fit_result.param
                min_N, max_N = minimum(sorted_N_values), maximum(sorted_N_values)
                fit_x_values = collect(range(min_N, max_N, length=5 * length(sorted_N_values)))
                fitted_curve = model(fit_x_values, fit_params)
                # formatted_params = map(x -> @sprintf("%.5f", x), fit_params)
    
                # Use the stored color for the fit line to match the scatter points
                ax.plot(fit_x_values, fitted_curve, label="Fit ($(fit_type))", linewidth=1.5, color=color_map[mu])
            end
        end
    end

    ax.legend()
    fig.tight_layout()
    fig.savefig(filename, dpi=dpi)
end

function save_integrals_to_csv(filepath::String, mp_integral::Float64, ipr_integral::Float64, mbs_integral::Float64)
    df = DataFrame(
        Integral_Type = ["mp_integral", "ipr_integral", "mbs_integral"],
        Value = [mp_integral, ipr_integral, mbs_integral]
    )
    
    CSV.write(filepath, df)
end




##############################################################################
# # Section 4 # #
##############################################################################
# # This makes use of the above plotting and integral calc functions to produce the standard set of plots to be extracted from the raw data


function standard_process_and_plot(
    df::DataFrame,
    N::Int,
    rho::Float64,
    Delta::Float64,
    seq_name::String,
    filepath::String
)
    # rho_safe = replace(string(rho), "." => "p")
    Delta_safe = replace(string(Delta), "." => "p")

    # mp_integral = calc_integral_under_curve(df, :mp_disc, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    # ipr_integral = calc_integral_under_curve(df, :ipr_disc, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    # mbs_integral = calc_integral_under_curve(df, :mbs_gap_disc, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_integrals.csv"
    # save_integrals_to_csv(filename, mp_integral, ipr_integral, mbs_integral)


    # filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
    # plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)#, sequence_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp"
    plt_mp_heatmap(df, :mp, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp_disc"
    plt_mp_heatmap(df, :mp_disc, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr"
    # plt_mp_heatmap(df, :ipr, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr_disc"
    # plt_mp_heatmap(df, :ipr_disc, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap"
    # plt_mp_heatmap(df, :mbs_gap, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap_disc"
    # plt_mp_heatmap_mbs_disc_final(df, :mbs_gap_disc, :mu_t, :rho, filename, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    return
end



# # Example usage of individual plotting functions from above

filepath = "/Users/Will/Documents/FINAL_PROJECT/simulations/phase_transition_fractality/data/hp/GQC_N55-610-6_g8-13_mu0.0-2.5-501_rho1.5-1.5-1_Delta0.05/"
N = 89
rho = 0.0
seq_name = "GQC"
# # standard_process_and_plot(df, N, rho, Delta, seq_name, filepath)
# Delta_range = [0.1, 0.2, 0.3, 0.5] #[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0] #[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
# # for Delta in Delta_range
# #     standard_process_and_plot(df, N, rho, Delta, seq_name, filepath)
# # end

########################
Delta=0.05
standard_process_and_plot(df, N, rho, Delta, seq_name, filepath)
#######################
N_range = [55, 89, 144, 233, 377, 610] #[5, 8, 13, 21, 34, 55, 89, 144]
# Delta = 1.0
# rho = 1.5

for N in N_range
    standard_process_and_plot(df, N, rho, Delta, seq_name, filepath)
end

# rho_safe = replace(string(rho), "." => "p")
# Delta_safe = replace(string(Delta), "." => "p")
# filepath = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_fractality_results/"
# filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"

# for N in N_range
#     filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
#     plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)
# end



# ########## Use just this for the mp fractality of phase transition data #########
# # plt_eigenvalues_generalised(df, :mu_t)#, :rho, :Delta) # # N.B. currently not saved to .bsons
# filename = "/Users/Will/Documents/FINAL_PROJECT/simulations/phase_transition_fractality/data/mp_plt"
# plt_mp_generalised(df, :mu_t, filename)

