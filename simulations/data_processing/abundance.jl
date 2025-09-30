
include("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_processing/bson_unpacker.jl")

using DataFrames
using CSV

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

function save_integrals_to_csv(filepath::String, mp_integral::Float64, ipr_integral::Float64, mbs_integral::Float64)
    df = DataFrame(
        Integral_Type = ["mp_integral", "ipr_integral", "mbs_integral"],
        Value = [mp_integral, ipr_integral, mbs_integral]
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
    filepath::String;
)
    Delta_safe = replace(string(Delta), "." => "p")

    mp_integral = calc_integral_under_curve(df, :mp_disc, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    ipr_integral = calc_integral_under_curve(df, :ipr_masked, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)
    mbs_integral = calc_integral_under_curve(df, :maj_gap_masked, :mu_t, :rho, N=N, Delta_t=Delta)#, sequence_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_integrals.csv"
    save_integrals_to_csv(filename, mp_integral, ipr_integral, mbs_integral)

end


###########################################################
###################### Sec 3: Run #########################
###########################################################

filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/mu_vs_rho_mp_heatmaps/GQC_N(200-200-1)_t1(1.0-1.0-51__t2(0.0-10.0-51)_mu(0.0-10.0-51)_Delta(0.1-0.1-1)/"
mp_tol = 0.1
df = unpack_bason_standard(filepath; mp_tol=mp_tol)

N = 200
Delta = 0.1
seq_name = "GQC"

standard_abundance_processing(df, N, Delta, seq_name, filepath)