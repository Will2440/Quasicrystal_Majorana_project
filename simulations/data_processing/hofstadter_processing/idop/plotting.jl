module IDOPPlotting

using Plots
using DataFrames
using Measures

## Updated to not error on empty sets
function smallData_plt_discrete_phase_projections(
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
    # Filter by fixed values
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_rows = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_rows)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return  # Skip instead of erroring
    end

    x_points = Float64[]
    y_points = Float64[]

    for row in eachrow(filtered_rows)
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

## new function to work directly with the bigData handling method
function plt_discrete_phase_projections(
    df::DataFrame,
    x_variable::Symbol,
    y_variable::Symbol,
    savepath::String;
    final_line_data=nothing,
    final_poly_data=nothing,
    fixed_values...
)
    """
    Plot mu_mp vs phi, where mu_mp is a vector of Union{Missing, Float64} per row,
    plotting only non-missing points.

    - fixed_values... filters the DataFrame the same way other plotting helpers do.
    """
    # Filter by fixed values
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_rows = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_rows)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return  # Skip instead of erroring
    end

    x_points = Float64[]
    y_points = Float64[]

    for row in eachrow(filtered_rows)
        # read axis values
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

        # Assume x_val is Vector{Union{Missing, Float64}}, y_val is Float64
        if x_val isa AbstractVector && !ismissing(y_val)
            for i in 1:length(x_val)
                x_p = x_val[i]
                if !ismissing(x_p)
                    push!(x_points, x_p)
                    push!(y_points, y_val)
                end
            end
        end
    end

    println(length(x_points), " points to plot")

    plt = Plots.scatter(
        x_points, y_points;
        markersize=0.1,
        markerstrokewidth=0,
        legend=false,
        xlabel="mu_mp",
        ylabel="phi",
        grid=true,
        size=(800,600)
    )

    # Overlay final regression line if provided
    if final_line_data !== nothing
        a, b = final_line_data.a, final_line_data.b
        # Use the actual phi range from the data used to fit the line
        phi_line = sort(final_line_data.phi_vals)
        mu_line = a .* phi_line .+ b
        # Plot the line in μ-φ space
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
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
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
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
    end
    if length(plateaus) != nrow(filtered_df)
        @warn "Length of plateaus vector must match number of filtered rows. Skipping plot."
        return
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
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
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

###########################################
########## Gap-labelled plotting ##########
###########################################

"""
Plot IDOP plateaus vs. phi, colored by gap labels (p/q).
"""
function plt_plateaus_vs_phi_coloured_legend(
    df::DataFrame,
    savepath::AbstractString;
    fixed_values...
)
    # Filter by fixed values (flexible, like other plotting functions)
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_df)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
    end
    
    x_points = Float64[]
    y_points = Float64[]
    colors = String[]
    
    for row in eachrow(filtered_df)
        phi = row.phi
        for label in row.gap_labels
            push!(x_points, label.plateau)
            push!(y_points, phi)
            push!(colors, "$(label.p)/$(label.q)")  # Color by p/q label
        end
    end
    
    plt = Plots.scatter(
        x_points, y_points,
        group=colors,
        markersize=2.5,
        markerstrokewidth=0,
        legend=:outerright,
        grid=true,
        xlabel="IDOP Plateau Value",
        ylabel="phi",
        size=(900,600),
        title="IDOP Plateaus vs phi (Gap-Labeled)"
    )
    Plots.savefig(savepath)
end

function plt_qled_coloured_gaps_mu_vs_phi(
    df::DataFrame,
    savepath::AbstractString;
    mu_axis::Symbol = :mu_mp_norm,  # Choose :mu_values or :mu_mp_norm for x-axis
    cmap=:RdBu,
    atol=1e-8,
    rtol=1e-6,
    verbose=false,
    fixed_values...
)
    # Filter by fixed values (flexible, like other plotting functions)
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_df)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
            for k in keys(Dict(fixed_values))
                if k in names(df)
                    vals = unique(df[:, k])
                    @info "Unique $k" n=length(vals) sample=first(vals, min(10, length(vals)))
                end
            end
        end
        @warn "No gaps match the specified fixed values. Skipping plot."
        return
    end

    # Collect all gaps
    mu_starts = Float64[]
    mu_ends = Float64[]
    phi_vals = Float64[]
    ps = Int[]
    qs = Int[]
    pq_set = Set{Tuple{Int,Int}}()
    q_to_ps = Dict{Int, Vector{Int}}()

    for row in eachrow(filtered_df)
        phi = row.phi
        gap_labels = row.gap_labels  # Vector{NamedTuple} with :plateau, :p, :q, :err
        idop_vals = row.idop  # Vector of IDOP values
        mu_vals = row[mu_axis]  # Vector of μ values (raw or normalized), may contain Missing
        
        for label in gap_labels
            plateau_val = label.plateau
            p = label.p
            q = label.q
            
            # Find μ ranges where IDOP ≈ plateau_val (flat regions)
            ranges = find_idop_ranges(idop_vals, plateau_val, atol)
            for (start_idx, end_idx) in ranges
                mu_start = mu_vals[start_idx]
                mu_end = mu_vals[end_idx]
                
                # Skip if any μ value is missing
                if ismissing(mu_start) || ismissing(mu_end)
                    continue
                end
                
                push!(mu_starts, mu_start)
                push!(mu_ends, mu_end)
                push!(phi_vals, phi)
                push!(ps, p)
                push!(qs, q)
                push!(pq_set, (p, q))
                q_to_ps[q] = get(q_to_ps, q, Int[])
                if !(p in q_to_ps[q])
                    push!(q_to_ps[q], p)
                end
            end
        end
    end

    if isempty(mu_starts)
        @warn "No valid gap labels to plot (possibly due to missing μ values)."
        return
    end

    # Determine unique q ordering (stable/sorted)
    q_list = sort(collect(keys(q_to_ps)))
    n_q = length(q_list)

    # Build base colours for q from cmap
    base_colors = nothing
    if isa(cmap, Symbol) || isa(cmap, String)
        base_colors = cgrad(cmap, n_q; categorical=true).colors
    elseif isa(cmap, AbstractVector)
        if length(cmap) >= n_q
            base_colors = cmap[1:n_q]
        else
            base_colors = cgrad(cmap, n_q; categorical=true).colors
        end
    else
        base_colors = cgrad(:viridis, n_q; categorical=true).colors
    end

    # Build final color for each (p,q): same hue per q, vary value by index of p within that q
    pq_to_color = Dict{Tuple{Int,Int}, Colorant}()
    for (qi, q) in enumerate(q_list)
        ps_for_q = sort(q_to_ps[q])   # deterministic ordering of p within q
        n_pq = length(ps_for_q)
        base_hsv = HSV(base_colors[qi])
        for (pi, pval) in enumerate(ps_for_q)
            t = n_pq == 1 ? 0.5 : (pi - 1) / (n_pq - 1)   # 0..1
            v_new = clamp(0.45 + 0.50 * t, 0.0, 1.0)
            s_new = clamp(0.6 + 0.35 * (1 - t), 0.0, 1.0)
            color_pq = RGB(HSV(base_hsv.h, s_new, v_new))
            pq_to_color[(pval, q)] = color_pq
        end
    end

    # Axes ranges
    xmin = minimum(mu_starts); xmax = maximum(mu_ends)
    ymin = minimum(phi_vals); ymax = maximum(phi_vals)
    x_span = xmax - xmin; x_span = x_span == 0 ? 1.0 : x_span
    y_span = ymax - ymin; y_span = y_span == 0 ? 1.0 : y_span

    # Band half-height in phi units
    h = max(1e-6, 0.5 * min(0.02 * y_span, (0.9 * y_span) / max(1, length(unique(phi_vals)))))

    plt = plot(xlabel=mu_axis == :mu_mp_norm ? "Normalized mu" : "mu", ylabel="phi", legend=false, grid=true, size=(900,600), title="Q-Led Colored Gaps in $(mu_axis == :mu_mp_norm ? "Normalized " : "")mu vs phi (IDOP)")

    # Draw each gap rectangle using pq_to_color mapping
    for i in eachindex(mu_starts)
        col = pq_to_color[(ps[i], qs[i])]
        phi = phi_vals[i]
        xs = [mu_starts[i], mu_ends[i], mu_ends[i], mu_starts[i], mu_starts[i]]
        ys = [phi - h, phi - h, phi + h, phi + h, phi - h]
        plot!(plt, xs, ys; seriestype=:shape, fillcolor=col, fillalpha=0.85, linecolor=:transparent, label=false)
    end

    # Legend-like annotations (same layout as IDOS)
    pq_list = sort(collect(pq_set), by = x -> (x[2], x[1]))  # sort by q then p
    right_margin = 0.40 * x_span
    plot!(plt, xlim=(xmin, xmax + right_margin))

    x_col1 = xmax + 0.03 * x_span
    x_col2 = xmax + 0.20 * x_span

    # Column 1: pq_list order
    N1 = length(pq_list)
    y_start1 = ymax - 0.02 * y_span
    y_step1 = min(0.03 * y_span, (0.9 * y_span) / max(1, N1))
    ys1 = [y_start1 - (i-1) * y_step1 for i in 1:N1]
    for (i, pq) in enumerate(pq_list)
        col = pq_to_color[pq]
        scatter!(plt, [x_col1 - 0.01 * x_span], [ys1[i]];
                 markersize=6, markerstrokewidth=0, color=col, label=false)
        annotate!(plt, (x_col1, ys1[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    # Column 2: ordered by q
    pq_by_q = sort(pq_list, by = x -> x[2])
    N2 = length(pq_by_q)
    y_start2 = ymax - 0.02 * y_span
    y_step2 = min(0.03 * y_span, (0.9 * y_span) / max(1, N2))
    ys2 = [y_start2 - (i-1) * y_step2 for i in 1:N2]
    for (i, pq) in enumerate(pq_by_q)
        col = pq_to_color[pq]
        scatter!(plt, [x_col2 - 0.01 * x_span], [ys2[i]];
                 markersize=6, markerstrokewidth=0, color=col, label=false)
        annotate!(plt, (x_col2, ys2[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    Plots.savefig(savepath)
end

# Helper function to find ranges where IDOP ≈ plateau_val
function find_idop_ranges(idop_vals::Vector{Float64}, plateau_val::Float64, atol::Float64)
    ranges = Tuple{Int, Int}[]
    in_range = false
    start_idx = 0
    for (i, val) in enumerate(idop_vals)
        if isapprox(val, plateau_val; atol=atol)
            if !in_range
                start_idx = i
                in_range = true
            end
        else
            if in_range
                push!(ranges, (start_idx, i-1))
                in_range = false
            end
        end
    end
    if in_range
        push!(ranges, (start_idx, length(idop_vals)))
    end
    return ranges
end


function plt_coloured_gaps_q_optionality(
    df::DataFrame,
    savepath::String;
    p_range::Union{AbstractVector{<:Integer},Nothing}=nothing,
    q_max::Union{Int,Nothing}=nothing,
    cmap::Any = :viridis,
    xlabel::String="µ",
    ylabel::String="φ",
    title::Union{String,Nothing}=nothing,
    verbose::Bool=false,
    atol::Real=1e-9,
    rtol::Real=0.0,
    colour_mode::Symbol = :qled_p,
    line_width::Int = 4,
    fixed_values...
)
    # approx-aware filtering on row-level context columns
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)

    filtered_rows = isempty(fixed_values) ? df :
                    filter(row -> _row_matches(row, fixed_values), df)

    if isempty(filtered_rows)
        if verbose
            @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
        end
        error("No gaps match the specified fixed values.")
    end

    # Collect plateaus (using plateau as the "energy" value for compatibility)
    lows  = Float64[]; highs = Float64[]; phis = Float64[]
    ps    = Int[];     qs    = Int[]
    pq_set = Set{Tuple{Int,Int}}()
    q_to_ps = Dict{Int, Vector{Int}}()

    row_idx = 0
    for row in eachrow(filtered_rows)
        row_idx += 1
        gls = getproperty(row, :gap_labels)
        (gls === nothing || isempty(gls)) && continue
        phi = Float64(row.phi)
        for g in gls
            p = Int(g.p); q = Int(g.q)
            if p_range !== nothing && !(p in p_range); continue; end
            if q_max   !== nothing && abs(q) > q_max; continue; end
            plateau = Float64(g.N_gap)
            mu_low = Float64(g.E_low)
            mu_high = Float64(g.E_high)
            push!(lows,  mu_low)
            push!(highs, mu_high)
            # push!(lows,  plateau)
            # push!(highs, plateau)  # Set highs to plateau for line plotting
            push!(phis,  phi)
            push!(ps,    p)
            push!(qs,    q)
            push!(pq_set, (p,q))
            q_to_ps[q] = get(q_to_ps, q, Int[])
            if !(p in q_to_ps[q])
                push!(q_to_ps[q], p)
            end
        end
    end

    isempty(lows) && error("No gaps remain after filtering.")

    # Keys for coloring: either q or |q|
    key_list = Int[]
    if colour_mode in [:qled_p, :qled_notp]
        key_list = sort(collect(keys(q_to_ps)))
    elseif colour_mode in [:absq_p, :absq_notp]
        key_set = Set{Int}()
        for q in keys(q_to_ps)
            push!(key_set, abs(q))
        end
        key_list = sort(collect(key_set))
    else
        error("Invalid colour_mode=$colour_mode (use :qled_p, :qled_notp, :absq_p, :absq_notp)")
    end
    n_key = length(key_list)

    # Base colours
    base_colors = if isa(cmap, Symbol) || isa(cmap, String)
        cgrad(cmap, n_key; categorical=true).colors
    elseif isa(cmap, AbstractVector)
        length(cmap) >= n_key ? cmap[1:n_key] : cgrad(cmap, n_key; categorical=true).colors
    else
        cgrad(:viridis, n_key; categorical=true).colors
    end

    # Map each (p,q) -> colour and legend label
    pq_to_color = Dict{Tuple{Int,Int}, Colorant}()
    pq_to_label = Dict{Tuple{Int,Int}, String}()

    # Helper: legend text per mode
    legend_label(p,q) = begin
        if colour_mode == :qled_p
            "(p=$p, q=$q)"
        elseif colour_mode == :qled_notp
            "q = $q"
        elseif colour_mode == :absq_p
            "(p=$p, q=$q, |q|=$(abs(q)))"
        elseif colour_mode == :absq_notp
            "|q| = $(abs(q))"
        else
            ""
        end
    end

    for (ki, key) in enumerate(key_list)
        if colour_mode in [:qled_p, :absq_p]
            # Vary colour by p within key
            ps_for_key = Int[]
            if colour_mode == :qled_p
                ps_for_key = sort(q_to_ps[key])
            elseif colour_mode == :absq_p
                for q in keys(q_to_ps)
                    if abs(q) == key
                        append!(ps_for_key, q_to_ps[q])
                    end
                end
                ps_for_key = sort(unique(ps_for_key))
            end
            n_p = length(ps_for_key)
            base_hsv = HSV(base_colors[ki])

            for (pi, pval) in enumerate(ps_for_key)
                t = n_p == 1 ? 0.5 : (pi - 1) / (n_p - 1)
                v_new = clamp(0.45 + 0.50 * t, 0.0, 1.0)
                s_new = clamp(0.6 + 0.35 * (1 - t), 0.0, 1.0)
                color_p = RGB(HSV(base_hsv.h, s_new, v_new))

                if colour_mode == :qled_p
                    pq = (pval, key)
                    pq_to_color[pq] = color_p
                    pq_to_label[pq] = legend_label(pq[1], pq[2])
                else # :absq_p
                    for q in keys(q_to_ps)
                        if abs(q) == key && pval in q_to_ps[q]
                            pq = (pval, q)
                            pq_to_color[pq] = color_p
                            pq_to_label[pq] = legend_label(pq[1], pq[2])
                        end
                    end
                end
            end

        else
            # :qled_notp or :absq_notp – same colour for all p in key
            color_key = base_colors[ki]
            if colour_mode == :qled_notp
                q = key
                for p in q_to_ps[q]
                    pq = (p, q)
                    pq_to_color[pq] = color_key
                    pq_to_label[pq] = legend_label(pq[1], pq[2])  # label is "q = q"
                end
            else # :absq_notp
                for q in keys(q_to_ps)
                    if abs(q) == key
                        for p in q_to_ps[q]
                            pq = (p, q)
                            pq_to_color[pq] = color_key
                            pq_to_label[pq] = legend_label(pq[1], pq[2])  # label is "|q| = |q|"
                        end
                    end
                end
            end
        end
    end

    # Axes ranges
    xmin = minimum(lows); xmax = maximum(highs)
    ymin = minimum(phis); ymax = maximum(phis)
    x_span = xmax - xmin; x_span = x_span == 0 ? 1.0 : x_span
    y_span = ymax - ymin; y_span = y_span == 0 ? 1.0 : y_span

    # Band half-height in φ units (adjusted for line plotting)
    h = max(1e-6, 0.5 * min(0.01 * y_span, (0.9 * y_span) / max(1, length(unique(phis)))))

    # Default titles by mode if not provided
    default_title = begin
        if colour_mode == :qled_p
            "IDOP Plateaus vs φ (colour by p and q)"
        elseif colour_mode == :qled_notp
            "IDOP Plateaus vs φ (colour by q)"
        elseif colour_mode == :absq_p
            "IDOP Plateaus vs φ (colour by |q| and p)"
        elseif colour_mode == :absq_notp
            "IDOP Plateaus vs φ (colour by |q|)"
        else
            "IDOP Plateaus vs φ"
        end
    end

    plt = plot(
        xlabel = xlabel,
        ylabel = ylabel,
        legend = :outerright,
        grid   = true,
        size   = (800,600),
        title  = isnothing(title) ? default_title : title,
    )

    # Group indices by (p,q) so we can plot each group once for legend
    # Manual grouping to avoid dependency if not wanted:
    pq_to_indices = Dict{Tuple{Int,Int}, Vector{Int}}()
    for i in eachindex(lows)
        pq = (ps[i], qs[i])
        push!(get!(pq_to_indices, pq, Int[]), i)
    end

    # # Plot each group with one legend entry
    if colour_mode in (:qled_notp, :absq_notp)
        key_to_indices = Dict{Int, Vector{Int}}()
        key_to_color   = Dict{Int, Colorant}()

        for (pq, idxs) in pq_to_indices
            key = colour_mode == :qled_notp ? pq[2] : abs(pq[2])
            push!(get!(key_to_indices, key, Int[]), idxs...)
            key_to_color[key] = pq_to_color[pq]
        end

        for key in sort(collect(keys(key_to_indices)))
            idxs = key_to_indices[key]
            col  = key_to_color[key]
            lbl  = colour_mode == :qled_notp ? "q = $key" : "|q| = $key"

            xs_all = Float64[]; ys_all = Float64[]
            for i in idxs
                phi = phis[i]
                # Plot as horizontal lines at plateau
                append!(xs_all, [lows[i], highs[i], NaN])
                append!(ys_all, [phi, phi, NaN])
            end

            plot!(plt, xs_all, ys_all; seriestype=:path, linecolor=col, linewidth=line_width, label=lbl)
        end
    else
        # Default: one legend entry per (p,q)
        for (pq, idxs) in sort(collect(pq_to_indices); by=x->(x[1][2], x[1][1]))
            col = pq_to_color[pq]
            lbl = pq_to_label[pq]
            xs_all = Float64[]; ys_all = Float64[]
            for i in idxs
                phi = phis[i]
                # Plot as horizontal lines at plateau
                append!(xs_all, [lows[i], highs[i], NaN])
                append!(ys_all, [phi, phi, NaN])
            end
            plot!(plt, xs_all, ys_all;
                seriestype=:path,
                linecolor=col,
                linewidth=line_width,
                label=lbl)
        end
    end

    # Add some right margin for legend
    plot!(plt, xlim=(xmin, xmax + 0.10 * x_span))

    savefig(plt, savepath)
    return plt
end

###########################################
#### Phason winding in MP gaps plotting ###
############################################

function plt_mu_mp_vs_phason(
    df::DataFrame,
    savepath::String;
    fixed_values...
)
    """
    Plot phason vs mu_mp for fixed values (e.g., phi).

    - mu_mp is a vector of Union{Missing, Float64} per row, plotting only non-missing points (no normalization).
    - Highlights missing mu ranges (0.0-0.6, 1.4-1.8, 2.2-3.0) with shaded gray areas.
    - X-axis (mu_mp) limited to 0.0-3.0.
    """
    # Filter by fixed values
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_rows = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_rows)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
    end

    x_points = Float64[]  # Now mu_mp
    y_points = Float64[]  # Now phason

    for row in eachrow(filtered_rows)
        # read phason values
        phason_val = nothing
        mu_mp_val = nothing
        try
            phason_val = getproperty(row, :phason)
        catch
            continue
        end
        try
            mu_mp_val = getproperty(row, :mu_mp)
        catch
            continue
        end

        # Assume mu_mp_val is Vector{Union{Missing, Float64}}, phason_val is Float64
        if mu_mp_val isa AbstractVector && !ismissing(phason_val)
            for mu_p in mu_mp_val
                if !ismissing(mu_p)
                    push!(x_points, mu_p)  # mu_mp as x
                    push!(y_points, phason_val)  # phason as y
                end
            end
        end
    end

    println(length(x_points), " points to plot")

    plt = Plots.scatter(
        x_points, y_points;
        markersize=1.5,
        markerstrokewidth=0,
        legend=false,
        xlabel="mu_mp",
        ylabel="phason",
        grid=true,
        size=(800,600),
        xlim=(0.0, 3.0)  # Set x-axis limits (mu_mp)
    )

    # Add shaded regions for missing data (vertical spans on x-axis)
    Plots.vspan!([0.0, 0.6], color=:gray, alpha=0.3, label="Missing data")
    Plots.vspan!([1.4, 1.8], color=:gray, alpha=0.3, label="")
    Plots.vspan!([2.2, 3.0], color=:gray, alpha=0.3, label="")

    Plots.savefig(savepath)
end

function plt_lowest_eigs_vs_phason(
    df::DataFrame,
    savepath::String;
    fixed_values...
)
    """
    Plot lowest eigenvalues vs phason (y) and mu (x), with two side-by-side subplots.
    Points are colored by the eigenvalue value.
    - Highlights missing mu ranges (0.0-0.6, 1.4-1.8, 2.2-3.0) with shaded gray areas.
    - X-axis limited to 0.0-3.0.
    - Colorbar shows the range of +E and -E values with explicit ticks.
    """
    # Filter by fixed values
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_rows = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_rows)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
    end

    # Prepare data for plotting (filter out missing eigenvalues)
    mu_pos = Float64[]
    phason_pos = Float64[]
    pos_eigs = Float64[]
    mu_neg = Float64[]
    phason_neg = Float64[]
    neg_eigs = Float64[]

    for row in eachrow(filtered_rows)
        phason_val = getproperty(row, :phason)
        mu_range = getproperty(row, :mu_range)
        min_pos_eigs = getproperty(row, :min_pos_eigs)
        min_neg_eigs = getproperty(row, :min_neg_eigs)

        for (mu, pos_eig, neg_eig) in zip(mu_range, min_pos_eigs, min_neg_eigs)
            if !ismissing(pos_eig)
                push!(mu_pos, mu)
                push!(phason_pos, phason_val)
                push!(pos_eigs, pos_eig)
            end
            if !ismissing(neg_eig)
                push!(mu_neg, mu)
                push!(phason_neg, phason_val)
                push!(neg_eigs, neg_eig)
            end
        end
    end

    # println(length(mu_pos), " positive points, ", length(mu_neg), " negative points to plot")

    # Create side-by-side plots
    plt = Plots.plot(
        layout=(1,2),
        size=(1200,600),
        grid=true,
        left_margin = 8mm,
        right_margin = 10mm,
        top_margin = 5mm,
        bottom_margin = 10mm
    )

    # Left: Positive eigenvalues
    if !isempty(pos_eigs)
        min_pos = minimum(pos_eigs)
        max_pos = maximum(pos_eigs)
        # println("Positive eigenvalue range: [$min_pos, $max_pos]")
        ticks_pos = range(min_pos, max_pos, length=5)  # 5 ticks for range
        Plots.scatter!(
            plt[1],
            mu_pos, phason_pos;
            marker_z=pos_eigs,
            label="min(+E)",
            markersize=2,
            markerstrokewidth=0,
            color=:viridis,
            xlabel="mu",
            ylabel="phason",
            title="Lowest Positive Eigenvalues",
            colorbar=true,
            colorbar_title="+E",
            # colorbar_ticks=ticks_pos,
            clims=(min_pos, max_pos),
            xlim=(0.0, 3.0),
            ylims=(0.0, 1.0)
        )
    end
    # Shaded missing regions
    Plots.vspan!(plt[1], [0.0, 0.6], color=:gray, alpha=0.3, label="Missing data")
    Plots.vspan!(plt[1], [1.4, 1.8], color=:gray, alpha=0.3, label="")
    Plots.vspan!(plt[1], [2.2, 3.0], color=:gray, alpha=0.3, label="")

    # Right: Negative eigenvalues
    if !isempty(neg_eigs)
        min_neg = minimum(neg_eigs)
        max_neg = maximum(neg_eigs)
        # println("Negative eigenvalue range: [$min_neg, $max_neg]")
        ticks_neg = range(min_neg, max_neg, length=5)  # 5 ticks for range
        Plots.scatter!(
            plt[2],
            mu_neg, phason_neg;
            marker_z=neg_eigs,
            label="min(-E)",
            markersize=2,
            markerstrokewidth=0,
            color=:viridis,
            xlabel="mu",
            ylabel="phason",
            title="Lowest Negative Eigenvalues",
            colorbar_title="-E",
            colorbar=:right,
            # colorbar_ticks=ticks_neg,
            clims=(min_neg, max_neg),
            xlim=(0.0, 3.0),
            ylims=(0.0,1.0)
        )
    end
    # Shaded missing regions
    Plots.vspan!(plt[2], [0.0, 0.6], color=:gray, alpha=0.3, label="Missing data")
    Plots.vspan!(plt[2], [1.4, 1.8], color=:gray, alpha=0.3, label="")
    Plots.vspan!(plt[2], [2.2, 3.0], color=:gray, alpha=0.3, label="")

    Plots.savefig(plt, savepath)
end

function plt_gap_eigs_vs_phason(
    df::DataFrame,
    savepath::String;
    fixed_values...
)
    """
    Plot gap eigenvalues vs phason (y) and mu (x), with two side-by-side subplots:
    - Left: Positive gap eigenvalues (closest to zero in gaps).
    - Right: Negative gap eigenvalues.
    - Highlights missing mu ranges (0.0-0.6, 1.4-1.8, 2.2-3.0) with shaded gray areas.
    - X-axis limited to 0.0-3.0.
    """
    # Filter by fixed values
    _kv_match(row, key, val) = begin
        rv = getproperty(row, key)
        (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
    end
    _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
    filtered_rows = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
    if isempty(filtered_rows)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return
    end

    # Prepare data for plotting
    mu_all = Float64[]
    phason_all = Float64[]
    pos_gap_eigs_all = Float64[]
    neg_gap_eigs_all = Float64[]

    for row in eachrow(filtered_rows)
        phason_val = getproperty(row, :phason)
        mu_range = getproperty(row, :mu_range)
        gap_pos_eigs = getproperty(row, :gap_pos_eigs)
        gap_neg_eigs = getproperty(row, :gap_neg_eigs)

        for (mu, pos_eig, neg_eig) in zip(mu_range, gap_pos_eigs, gap_neg_eigs)
            push!(mu_all, mu)
            push!(phason_all, phason_val)
            push!(pos_gap_eigs_all, pos_eig)
            push!(neg_gap_eigs_all, neg_eig)
        end
    end

    println(length(mu_all), " points to plot")

    # Create side-by-side plots
    plt = Plots.plot(
        layout=(1,2),
        size=(1200,600),
        grid=true
    )

    # Left: Positive gap eigenvalues
    Plots.scatter!(
        plt[1],
        mu_all, phason_all;
        z=pos_gap_eigs_all,
        markersize=2,
        markerstrokewidth=0,
        color=:viridis,
        xlabel="mu",
        ylabel="phason",
        title="Positive Gap Eigenvalues",
        colorbar=true,
        xlim=(0.0, 3.0)
    )
    # Shaded missing regions
    Plots.vspan!(plt[1], [0.0, 0.6], color=:gray, alpha=0.3, label="Missing data")
    Plots.vspan!(plt[1], [1.4, 1.8], color=:gray, alpha=0.3, label="")
    Plots.vspan!(plt[1], [2.2, 3.0], color=:gray, alpha=0.3, label="")

    # Right: Negative gap eigenvalues
    Plots.scatter!(
        plt[2],
        mu_all, phason_all;
        z=neg_gap_eigs_all,
        markersize=2,
        markerstrokewidth=0,
        color=:viridis,
        xlabel="mu",
        ylabel="phason",
        title="Negative Gap Eigenvalues",
        colorbar=true,
        xlim=(0.0, 3.0)
    )
    # Shaded missing regions
    Plots.vspan!(plt[2], [0.0, 0.6], color=:gray, alpha=0.3, label="Missing data")
    Plots.vspan!(plt[2], [1.4, 1.8], color=:gray, alpha=0.3, label="")
    Plots.vspan!(plt[2], [2.2, 3.0], color=:gray, alpha=0.3, label="")

    Plots.savefig(plt, savepath)
end



end # module IDOPPlotting