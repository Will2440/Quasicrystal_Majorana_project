using DataFrames
using Plots

function plt_eval_projections(
    df::DataFrame,
    y_variable::Symbol,
    savepath::String;
    fixed_values...
)
    """
    Plots projected eigenvalues along x-axis as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Dictionary specifying fixed values for all other DataFrame columns.
    """

    plt = Plots.scatter(
        legend=false, 
        xlabel="Eigenvalues",
        ylabel=string(y_variable),
        grid=true,
        size=(800,600)
    )

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    y_values = filtered_df[:, y_variable]
    all_eigenvalues = filtered_df[:, :eigenvalues]
    
    # For each row, plot all eigenvalues at the corresponding y_value
    for (i, y) in enumerate(y_values)
        eigenvalues = all_eigenvalues[i]
        # Plot each eigenvalue as a point at y
        Plots.scatter!(plt, eigenvalues, fill(y, length(eigenvalues)), markersize=1, markerstrokewidth=0)
    end

    Plots.savefig(savepath)
    # Plots.display(plt)
end


function plt_idos_vs_energy(df::DataFrame, savepath::String; fixed_values...)
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plt = Plots.plot(
        xlabel="Energy",
        ylabel="IDOS",
        legend=:topright,
        grid=true,
        size=(800,600),
        title="IDOS vs Energy"
    )

    for (i, row) in enumerate(eachrow(filtered_df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
        Plots.plot!(plt, energies, idos, label=label_str)
    end

    Plots.savefig(savepath)
    # Plots.display(plt)
end

function plt_idos_plateaus_check(df::DataFrame, plateaus::Vector{Vector{Float64}}, savepath::String; fixed_values...)
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plt = Plots.plot(
        xlabel="Energy",
        ylabel="IDOS",
        legend=:topright,
        grid=true,
        size=(800,600),
        title="IDOS vs Energy with Plateaus to check threshold performance"
    )

    for (i, row) in enumerate(eachrow(filtered_df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
        Plots.plot!(plt, energies, idos, label=label_str)
        # Overlay horizontal dotted lines at each provided plateau value
        for plateau in plateaus[i]
            Plots.hline!(plt, [plateau], linestyle=:dot, color=:black, label=false)
        end
    end

    Plots.savefig(savepath)
    # Plots.display(plt)
end

function plt_plateaus_vs_phi(df::DataFrame, savepath::String; fixed_values...)
    # Filter by fixed values if needed
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # Sort by phi
    sorted_df = sort(filtered_df, :phi)

    x_points = Float64[]
    y_points = Float64[]
    for row in eachrow(sorted_df)
        phi = row.phi
        for plateau in row.plateaus
            push!(x_points, plateau)
            push!(y_points, phi)
        end
    end

    plt = Plots.scatter(
        x_points, y_points,
        markersize=2.5,
        markerstrokewidth=0,
        legend=false,
        grid=true,
        xlabel="IDOS Plateau Value",
        ylabel="phi",
        size=(900,600),
        title="IDOS Plateaus vs phi"
    )
    Plots.savefig(savepath)
    # Plots.display(plt)
end