using Plots
using DataFrames

function plt_discrete_phase_projections(
    df::DataFrame,
    x_variable::Symbol,
    y_variable::Symbol,
    disc_variable::Symbol,
    savepath::String;
    fixed_values...
)
    """
    Plot rows where a discretised indicator column marks a match.

    - x_variable, y_variable : Symbols for the x and y axes (e.g. :mu, :phi).
    - disc_variable : Symbol of the discretised column to use (:disc_eigenvalues or :disc_mp).
        * :disc_eigenvalues is expected to be a Vector{Int} of length >= 1 (take first element).
        * :disc_mp is expected to be a scalar Int/Float indicating match (no indexing).
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

        # accept both Vector (disc_eigenvalues) and scalar (disc_mp)
        if disc_variable == :disc_eigenvalues || string(disc_variable) == "disc_eigenvalues"
            if !(disc isa AbstractVector) || length(disc) == 0
                continue
            end
            d = disc[1]
        elseif disc_variable == :disc_mp || string(disc_variable) == "disc_mp"
            # allow missing/Nullable values to be treated as non-match
            if disc === missing
                continue
            end
            d = Int(round(disc))  # coerce numeric to Int (1/0)
        else
            @warn "Unknown discretised variable" disc_variable
            continue
        end

        if d != 1
            continue
        end

        push!(x_points, x_val)
        push!(y_points, y_val)
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
    Plots.display(plt)
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