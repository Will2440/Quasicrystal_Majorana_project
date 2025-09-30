include("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_processing/bson_unpacker.jl")

using Plots
using LaTeXStrings
using Measures


###########################################################
########### Sec 1: Generic Plotting Functions #############
###########################################################
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

function plt_mp_heatmap(df::DataFrame, colour_var::Symbol, x_variable::Symbol, y_variable::Symbol, filename::String, title::String; fixed_values...)
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
        # xlabel=L"\mu/t_A",
        # ylabel=L"\rho",
        title=title,
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


###########################################################
########### Sec 2: Standard Plotting Processes ############
###########################################################


function standard_plotting(
    df::DataFrame,
    seq_name::String,
    mp_tol::Float64,
    filepath::String;
    N::Int64=50,
    rho::Float64=1.5,
    Delta::Float64=0.1,
)
    rho_safe = replace(string(rho), "." => "p")
    Delta_safe = replace(string(Delta), "." => "p")

    # filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
    # plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp"
    plt_mp_heatmap(df, :mp, :mu_t, :rho, filename, "MP for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp_disc"
    plt_mp_heatmap(df, :mp_disc, :mu_t, :rho, filename, "MP discretised (tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr"
    plt_mp_heatmap(df, :ipr, :mu_t, :rho, filename, "IPR for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr_disc"
    plt_mp_heatmap(df, :ipr_masked, :mu_t, :rho, filename, "IPR phase masked (MP tol=$mp_tol) for $(seq_name) N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap"
    plt_mp_heatmap(df, :maj_gap, :mu_t, :rho, filename, "MBS gap for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap_disc"
    # plt_mp_heatmap_mbs_disc_final(df, :maj_gap_masked, :mu_t, :rho, filename, "MBS gap masked (MP tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    return
end

function standard_plotting_PQC(
    df::DataFrame,
    seq_name::String,
    mp_tol::Float64,
    filepath::String;
    N::Int64=50,
    rho::Float64=1.5,
    sigma::Float64=2.0,
    Delta::Float64=0.1,
)
    rho_safe = replace(string(rho), "." => "p")
    sig_safe = replace(string(sigma), "." => "p")
    Delta_safe = replace(string(Delta), "." => "p")

    # filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
    # plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mp"
    plt_mp_heatmap(df, :mp, :mu_t, :rho, filename, "MP for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mp_disc"
    plt_mp_heatmap(df, :mp_disc, :mu_t, :rho, filename, "MP discretised (tol=$mp_tol) for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_ipr"
    plt_mp_heatmap(df, :ipr, :mu_t, :rho, filename, "IPR for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_ipr_disc"
    plt_mp_heatmap(df, :ipr_masked, :mu_t, :rho, filename, "IPR phase masked (MP tol=$mp_tol) for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mbs_gap"
    plt_mp_heatmap(df, :maj_gap, :mu_t, :rho, filename, "MBS gap for $(seq_name), N=$N, sigma=$sigma, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap_disc"
    # plt_mp_heatmap_mbs_disc_final(df, :maj_gap_masked, :mu_t, :rho, filename, "MBS gap masked (MP tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    return
end


###########################################################
###################### Sec 3: Run #########################
###########################################################

filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/all_crystal_grad_testruns/restricted_mu_vs_rho_mp_heatmaps/GQC_N(50-50-1)_t1(1.0-1.0-101__t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(0.0-2.0-21)/"
mp_tol = 0.1
df = unpack_bason_standard(filepath; mp_tol=mp_tol)


N = 50
seq_name = "GQC"
Delta=0.5

# # One instance
# standard_plotting(df, seq_name, mp_tol, filepath; Delta=Delta)

# Instances of all Delta
Delta_range = collect(range(0.0, 2.0, 21))
for Delta in Delta_range
    standard_plotting_PQC(df, seq_name, mp_tol, filepath; Delta=Delta)
end


# # Instances of Delta and sigma for PQC
# Delta_range = collect(range(0.0, 2.0, 21))
# sigma_range = [2.0, 3.0, 4.0]
# seq_name = "PQC"
# for sigma in sigma_range
#     for Delta in Delta_range
#         standard_plotting_PQC(df, seq_name, mp_tol, filepath; Delta=Delta, sigma=sigma)
#     end
# end