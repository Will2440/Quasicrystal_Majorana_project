using Plots
using DataFrames

# function plt_discrete_phase_projections(
#     df::DataFrame,
#     x_variable::Symbol,
#     y_variable::Symbol,
#     disc_variable::Symbol,
#     savepath::String;
#     fixed_values...
# )
#     """
#     Plot rows where a discretised indicator column marks a match.

#     - x_variable, y_variable : Symbols for the x and y axes (e.g. :mu, :phi, :mu_values_norm).
#         * Can be scalar (Float64) or vector (Vector{Float64}) per row.
#     - disc_variable : Symbol of the discretised column to use (:disc_eigenvalues, :disc_mp, or :mp_disc_norm).
#         * Can be scalar (Int/Float) or vector (Vector{Int}) per row.
#         * For vectors, plots points where the element == 1, using corresponding indices from x_variable/y_variable if they are vectors.
#     - fixed_values... filters the DataFrame the same way other plotting helpers do.
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     x_points = Float64[]
#     y_points = Float64[]

#     for row in eachrow(filtered_df)
#         # read axis values (use proper try/catch blocks to avoid parse errors)
#         x_val = nothing
#         y_val = nothing
#         try
#             x_val = getproperty(row, x_variable)
#         catch
#             continue
#         end
#         try
#             y_val = getproperty(row, y_variable)
#         catch
#             continue
#         end

#         disc = get(row, disc_variable, nothing)
#         if disc === nothing
#             # debug help: uncomment to see what's wrong
#             @show disc_variable, hasproperty(row, disc_variable), names(df)
#             continue
#         end

#         # Handle disc as vector or scalar
#         if disc isa AbstractVector
#             # Vector case: iterate and plot where disc[i] == 1
#             for i in 1:length(disc)
#                 if disc[i] == 1
#                     x_p = x_val isa AbstractVector ? x_val[i] : x_val
#                     y_p = y_val isa AbstractVector ? y_val[i] : y_val
#                     push!(x_points, x_p)
#                     push!(y_points, y_p)
#                 end
#             end
#         else
#             # Scalar case
#             if disc === missing
#                 continue
#             end
#             d = Int(round(disc))  # coerce numeric to Int (1/0)
#             if d == 1
#                 push!(x_points, x_val)
#                 push!(y_points, y_val)
#             end
#         end
#     end

#     println(length(x_points), " points to plot for discretised variable ", string(disc_variable))
#     println(length(y_points), " points to plot for discretised variable ", string(disc_variable))

#     plt = Plots.scatter(
#         x_points, y_points;
#         markersize=1.5,
#         markerstrokewidth=0,
#         legend=false,
#         xlabel=string(x_variable),
#         ylabel=string(y_variable),
#         grid=true,
#         size=(800,600)
#     )

#     Plots.savefig(savepath)
#     # Plots.display(plt)
# end


function plt_discrete_phase_projections(
    df::DataFrame,
    x_variable::Symbol,
    y_variable::Symbol,
    disc_variable::Symbol,
    savepath::String;
    final_line_data=nothing,
    final_poly_data=nothing,
    fixed_values...
)
    """
    Plot rows where a discretised indicator column marks a match.

    - x_variable, y_variable : Symbols for the x and y axes (e.g. :mu, :phi, :mu_mp_norm).
        * Can be scalar (Float64) or vector (Vector{Union{Missing, Float64}}) per row.
    - disc_variable : Symbol of the discretised column to use (:disc_eigenvalues, :disc_mp, or :mp_disc_norm).
        * Can be scalar (Int/Float) or vector (Vector{Int}) per row.
        * For vectors, plots points where the element == 1, using corresponding indices from x_variable/y_variable if they are vectors.
        * Skips points where x_variable or y_variable is missing.
    - fixed_values... filters the DataFrame the same way other plotting helpers do.
    """
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_points = Float64[]
    y_points = Float64[]

    for row in eachrow(filtered_df)
        # read axis values (use proper try/catch blocks to avoid parse errors)
        x_val = nothing
        y_val = nothing
        try
            x_val = getproperty(row, x_variable)
        catch
            continue
        end
        try
            y_val = getproperty(row, y_variable)
        catch
            continue
        end

        disc = get(row, disc_variable, nothing)
        if disc === nothing
            # debug help: uncomment to see what's wrong
            @show disc_variable, hasproperty(row, disc_variable), names(df)
            continue
        end

        # Handle disc as vector or scalar
        if disc isa AbstractVector
            # Vector case: iterate and plot where disc[i] == 1
            for i in 1:length(disc)
                if disc[i] == 1
                    x_p = x_val isa AbstractVector ? x_val[i] : x_val
                    y_p = y_val isa AbstractVector ? y_val[i] : y_val
                    # Skip if missing
                    if !ismissing(x_p) && !ismissing(y_p)
                        push!(x_points, x_p)
                        push!(y_points, y_p)
                    end
                end
            end
        else
            # Scalar case
            if disc === missing
                continue
            end
            d = Int(round(disc))  # coerce numeric to Int (1/0)
            if d == 1
                # Skip if missing
                if !ismissing(x_val) && !ismissing(y_val)
                    push!(x_points, x_val)
                    push!(y_points, y_val)
                end
            end
        end
    end

    println(length(x_points), " points to plot for discretised variable ", string(disc_variable))
    println(length(y_points), " points to plot for discretised variable ", string(disc_variable))

    plt = Plots.scatter(
        x_points, y_points;
        markersize=1.5,
        markerstrokewidth=0,
        legend=false,
        xlabel=string(x_variable),
        ylabel=string(y_variable),
        grid=true,
        size=(800,600)
    )

    # Overlay final regression line if provided
    if final_line_data !== nothing
        a, b = final_line_data.a, final_line_data.b
        # Use the actual phi range from the data used to fit the line
        phi_line = sort(final_line_data.phi_vals)
        mu_line = a .* phi_line .+ b
        # Plot the line in μ-φ space (or normalise if needed for your x-axis)
        Plots.plot!(mu_line, phi_line; color=:red, linewidth=2, label="Final MP regression")
    end

    if final_poly_data !== nothing
        p = final_poly_data.poly
        phi_line = range(minimum(final_poly_data.phi_vals), maximum(final_poly_data.phi_vals), length=200)
        mu_line = p.(phi_line)
        Plots.plot!(mu_line, phi_line; color=:red, linewidth=2, label="Final MP poly fit")
    end

    Plots.savefig(savepath)
    # Plots.display(plt)
end


function plt_idop_vs_mu(df::DataFrame, savepath::String; fixed_values...)
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plt = Plots.plot(
        xlabel = "mu",
        ylabel = "IDOP",
        legend = :topright,
        grid = true,
        size = (800, 600),
        title = "IDOP vs mu"
    )

    for (i, row) in enumerate(eachrow(filtered_df))
        mu_vals = row.mu_values       # Vector{Float64}
        idop_vals = row.idop          # Vector{Float64}
        label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
        Plots.plot!(plt, mu_vals, idop_vals, label = label_str)
    end

    Plots.savefig(savepath)
    Plots.display(plt)
end

function plt_idop_plateaus_check(df::DataFrame, plateaus::Vector{Vector{Float64}}, savepath::String; fixed_values...)
    """Plot IDOP vs mu for filtered rows and overlay horizontal lines at provided plateau values (per-row)."""
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end
    if length(plateaus) != nrow(filtered_df)
        error("Length of plateaus vector must match number of filtered rows.")
    end

    plt = Plots.plot(
        xlabel = "mu",
        ylabel = "IDOP",
        legend = :topright,
        grid = true,
        size = (900, 600),
        title = "IDOP vs mu with detected plateaus"
    )

    for (i, row) in enumerate(eachrow(filtered_df))
        mu_vals = row.mu_values
        idop_vals = row.idop
        label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
        Plots.plot!(plt, mu_vals, idop_vals, label = label_str)
        for p in plateaus[i]
            Plots.hline!(plt, [p], linestyle = :dot, color = :black, label = false)
        end
    end

    Plots.savefig(savepath)
    # Plots.display(plt)
end

function plt_plateaus_vs_phi_idop(df::DataFrame, savepath::String; fixed_values...)
    """
    Scatter plot of IDOP plateau values vs phi.
    - df must contain columns :phi and :plateaus (Vector{Float64} per row).
    - fixed_values... filters rows before plotting (same pattern as other helpers).
    """
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    sorted_df = sort(filtered_df, :phi)

    x_points = Float64[]
    y_points = Float64[]

    for row in eachrow(sorted_df)
        phi = row.phi
        for p in row.plateaus
            push!(x_points, p)
            push!(y_points, phi)
        end
    end

    plt = Plots.scatter(
        x_points, y_points;
        markersize = 2.5,
        markerstrokewidth = 0,
        legend = false,
        grid = true,
        xlabel = "IDOP Plateau Value",
        ylabel = "phi",
        size = (900, 600),
        title = "IDOP Plateaus vs phi"
    )

    Plots.savefig(savepath)
    # Plots.display(plt)
end