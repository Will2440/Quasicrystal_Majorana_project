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