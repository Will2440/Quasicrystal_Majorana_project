using DataFrames
using Plots
using BSON
using Statistics
using Colors

size = (1200, 900)

function load_sturmian_grad_bson(path::AbstractString)
    @assert isfile(path) "BSON file not found: $path"
    raw = BSON.load(path)

    # normalize keys to Symbol -> value mapping
    norm = Dict{Symbol,Any}()
    for (k,v) in pairs(raw)
        # keys from BSON.load are usually strings; convert to Symbol safely
        kn = try
            Symbol(k)
        catch
            Symbol(string(k))
        end
        norm[kn] = v
    end

    # 1) If any value is already a DataFrame, return it
    for v in values(norm)
        if v isa DataFrame
            return deepcopy(v)   # return a copy to avoid accidental mutation of loaded object
        end
    end

    # 2) If top-level has prefix & slope arrays, build a DataFrame
    if (haskey(norm, :prefix) || haskey(norm, Symbol("prefix"))) &&
       (haskey(norm, :slope)  || haskey(norm, Symbol("slope")))
        pref = haskey(norm, :prefix) ? norm[:prefix] : norm[Symbol("prefix")]
        slp  = haskey(norm, :slope)  ? norm[:slope]  : norm[Symbol("slope")]
        return DataFrame(prefix = pref, slope = slp)
    end

    # 3) If the file contains a single value which is a Vector{NamedTuple}, convert it
    if length(norm) == 1
        v = first(values(norm))
        if isa(v, Vector) && !isempty(v) && v[1] isa NamedTuple
            return DataFrame(v)
        end
    end

    error("No DataFrame or recognizable (prefix,slope) data found in BSON file: $path. Keys: $(collect(keys(norm)))")
end

function plt_eval_projections(
    df::DataFrame,
    evals_col::Symbol,
    y_variable::Symbol,
    savepath::String;
    colour_rats::Bool=false,
    colour_comps::Bool=false,
    grad_filepath::Union{String,Nothing}=nothing,
    fixed_values...
)
    """
    Plots projected eigenvalues along x-axis as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Dictionary specifying fixed values for all other DataFrame columns.
    """

    x_label = evals_col == :eigenvalues ? "Eigenvalues" :
              evals_col == :eigs_norm  ? "Normalised Eigenvalues" :
              string(evals_col)

    # show legend only when colouring modes request it
    legend_setting = (colour_comps || colour_rats) ? :topright : false

    plt = Plots.scatter(
        legend=legend_setting,
        xlabel=x_label,
        ylabel=string(y_variable),
        grid=true,
        size=size
    )

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    y_values = filtered_df[:, y_variable]
    all_eigenvalues = filtered_df[:, evals_col]

    # helper for float-safe comparison
    matches_any = function(val, arr; tol=1e-8)
        v = Float64(val)
        for a in arr
            if abs(v - Float64(a)) < tol
                return true
            end
        end
        return false
    end

    if colour_comps
        # grad_filepath must be provided
        if grad_filepath === nothing
            error("colour_comps=true requires grad_filepath to be provided.")
        end

        grad_data = load_sturmian_grad_bson(grad_filepath)

        # extract a vector of phis from whatever load returned
        raw_phis = Float64[]
        if grad_data isa DataFrame
            if :slope in names(grad_data)
                raw_phis = Float64.(grad_data[:, :slope])
            elseif :phi in names(grad_data)
                raw_phis = Float64.(grad_data[:, :phi])
            else
                # pick first numeric column
                found = false
                for nm in names(grad_data)
                    col = grad_data[:, nm]
                    if eltype(col) <: Number || all(x->x isa Number, col)
                        raw_phis = Float64.(col)
                        found = true
                        break
                    end
                end
                if !found
                    error("Could not find numeric phi/slope column in grad file DataFrame.")
                end
            end
        elseif isa(grad_data, AbstractVector)
            raw_phis = Float64.(grad_data)
        else
            error("Unexpected grad file format returned by load_sturmian_grad_bson.")
        end

        # classification buckets
        raw_x = Float64[]; raw_y = Float64[]
        comp_x = Float64[]; comp_y = Float64[]
        other_x = Float64[]; other_y = Float64[]

        for (i, y) in enumerate(y_values)
            eigenvalues = real.(all_eigenvalues[i])
            yval = Float64(y)
            if matches_any(yval, raw_phis)
                append!(raw_x, eigenvalues); append!(raw_y, fill(yval, length(eigenvalues)))
            elseif matches_any(yval, 1 .- raw_phis)
                append!(comp_x, eigenvalues); append!(comp_y, fill(yval, length(eigenvalues)))
            else
                append!(other_x, eigenvalues); append!(other_y, fill(yval, length(eigenvalues)))
            end
        end

        # plot order: other (background), comps (semi-transparent), raws (top)
        if !isempty(other_x)
            Plots.scatter!(plt, other_x, other_y; markersize=1, markerstrokewidth=0, color=:gray, alpha=0.95, label="rationals")
        end
        if !isempty(comp_x)
            Plots.scatter!(plt, comp_x, comp_y; markersize=1, markerstrokewidth=0, color=:red, alpha=0.95, label="(1-φ)")
        end
        if !isempty(raw_x)
            Plots.scatter!(plt, raw_x, raw_y; markersize=1, markerstrokewidth=0, color=:blue, alpha=0.35, label="φ")
        end

    elseif colour_rats
        # helper: detect values equal to 1/r or 1 - 1/r within tol (r up to maxr)
        is_1_over_r = function(phi; tol=1e-8, maxr=15)
            # ensure numeric
            ph = Float64(phi)
            for r in 1:maxr
                if abs(ph - 1/r) < tol || abs(ph - (1 - 1/r)) < tol
                    return true
                end
            end
            return false
        end

        rat_x = Float64[]
        rat_y = Float64[]
        irr_x = Float64[]
        irr_y = Float64[]

        for (i, y) in enumerate(y_values)
            eigenvalues = real.(all_eigenvalues[i])
            yval = Float64(y)
            if is_1_over_r(yval)
                append!(rat_x, eigenvalues)
                append!(rat_y, fill(yval, length(eigenvalues)))
            else
                append!(irr_x, eigenvalues)
                append!(irr_y, fill(yval, length(eigenvalues)))
            end
        end

        # plot irrationals first (background colour), rationals on top
        if !isempty(irr_x)
            Plots.scatter!(plt, irr_x, irr_y; markersize=1, markerstrokewidth=0, color=:gray, alpha=0.35, label="irrational φ")
        end
        if !isempty(rat_x)
            Plots.scatter!(plt, rat_x, rat_y; markersize=1, markerstrokewidth=0, color=:blue, alpha=0.95, label="1/r or 1-1/r φ")
        end
    else
        # default behaviour: plot all eigenvalues the same way (no legend)
        for (i, y) in enumerate(y_values)
            eigenvalues = real.(all_eigenvalues[i])
            Plots.scatter!(plt, eigenvalues, fill(Float64(y), length(eigenvalues)); markersize=1, markerstrokewidth=0, label=false)
        end
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
        size=size,
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
        size=size,
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


#######################################
############ Gap Labelling ############
#######################################

function plot_idos_with_gaps(energies::Vector{Float64},
                             idos::Vector{Float64},
                             gap_labels::Vector{NamedTuple};
                             linewidth=1.5,
                             markersize=3,
                             label_fontsize=8)

    x_min, x_max = extrema(energies)

    plt = plot(energies, idos;
        seriestype = :line,
        linewidth  = linewidth,
        markersize = markersize,
        label      = "IDOS",
        xlabel     = "Energy",
        ylabel     = "Integrated DOS",
        size = (1200, 800)
    )

    drawn = Set{Float64}()

    for gap in gap_labels
        N_gap = gap.N_gap
        p     = gap.p
        q     = gap.q

        if N_gap in drawn
            continue
        end
        push!(drawn, N_gap)

        plot!(plt, [x_min, x_max], [N_gap, N_gap];
              linestyle = :dash, label = "")

        x_label = x_max + 0.02 * (x_max - x_min)
        annotate!(plt, x_label, N_gap, text("($p,$q)", label_fontsize))
    end

    return plt
end

function plt_plateaus_vs_phi_coloured_legend(df::DataFrame, savepath::String;
                                      fixed_values...)
    # 1. Filter rows
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # 2. Sort by phi
    sorted_df = sort(filtered_df, :phi)

    # 3. Collect unique (p,q)
    pq_set = Set{Tuple{Int,Int}}()
    for row in eachrow(sorted_df)
        for gap in row.gap_labels
            push!(pq_set, (gap.p, gap.q))
        end
    end
    pq_list = collect(pq_set)
    pq_to_idx = Dict{Tuple{Int,Int},Int}((pq, i) for (i, pq) in enumerate(pq_list))

    # 4. Build plot data
    x_points = Float64[]
    y_points = Float64[]
    c_points_idx = Int[]   # integer index into palette

    for row in eachrow(sorted_df)
        phi = row.phi
        for gap in row.gap_labels
            push!(x_points, gap.N_gap)
            push!(y_points, phi)
            push!(c_points_idx, pq_to_idx[(gap.p, gap.q)])
        end
    end

    if isempty(x_points)
        error("No plateau points to plot.")
    end

    # determine z-range for consistent mapping
    minz = minimum(c_points_idx)
    maxz = maximum(c_points_idx)
    denom = maxz - minz
    denom = denom == 0 ? 1 : denom

    # 5. Main scatter: use marker_z and fix zlims so further scatter! calls use identical mapping
    plt = scatter(
        x_points, y_points;
        marker_z = c_points_idx,
        zlims = (minz, maxz),
        markersize = 3,
        markerstrokewidth = 0,
        legend = false,
        colorbar = false,
        grid = true,
        xlabel = "IDOS Plateau Value (N_gap)",
        ylabel = "phi",
        size = (1000, 700),
        title = "Gap-labelled IDOS Plateaus vs phi"
    )

    # compute ranges and reserve room on the right for annotation columns
    xmin, xmax = extrema(x_points)
    ymin, ymax = extrema(y_points)
    x_span = xmax - xmin
    y_span = ymax - ymin
    x_span = x_span == 0.0 ? max(1.0, abs(xmin)) : x_span
    y_span = y_span == 0.0 ? max(1.0, abs(ymin)+1.0) : y_span

    right_margin = 0.40 * x_span     # reserve 40% of span on right for annotations
    plot!(plt, xlim = (xmin, xmax + right_margin))

    # annotation column x positions
    x_col1 = xmax + 0.03 * x_span
    x_col2 = xmax + 0.20 * x_span

    # --- Column 1: preserve original pq_list order (top-down) ---
    N1 = length(pq_list)
    y_start1 = ymax - 0.02 * y_span
    # compute spacing to fit entries within ~90% of vertical span
    y_step1 = min(0.03 * y_span, (0.9 * y_span) / max(1, N1))
    ys1 = [y_start1 - (i-1) * y_step1 for i in 1:N1]

    for (i, pq) in enumerate(pq_list)
        idx = pq_to_idx[pq]
        # use marker_z for the small marker so Plots uses the same colormap mapping
        scatter!(plt, [x_col1 - 0.01 * x_span], [ys1[i]];
                 marker_z = [idx],
                 zlims = (minz, maxz),
                 markersize = 6,
                 markerstrokewidth = 0,
                 label = false)
        annotate!(plt, (x_col1, ys1[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    # --- Column 2: ordered by q ---
    pq_by_q = sort(pq_list, by = x -> x[2])
    N2 = length(pq_by_q)
    y_start2 = ymax - 0.02 * y_span
    y_step2 = min(0.03 * y_span, (0.9 * y_span) / max(1, N2))
    ys2 = [y_start2 - (i-1) * y_step2 for i in 1:N2]

    for (i, pq) in enumerate(pq_by_q)
        idx = pq_to_idx[pq]
        scatter!(plt, [x_col2 - 0.01 * x_span], [ys2[i]];
                 marker_z = [idx],
                 zlims = (minz, maxz),
                 markersize = 6,
                 markerstrokewidth = 0,
                 label = false)
        annotate!(plt, (x_col2, ys2[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    savefig(plt, savepath)
    # display(plt)
end


# ## for plotting the gap-labelled areas in energy space
function plt_coloured_gaps_energy_vs_phi(
    df::DataFrame,
    savepath::String;
    p_range::Union{AbstractVector{<:Integer},Nothing}=nothing,
    q_max::Union{Int,Nothing}=nothing,
    cmap::Any = :viridis,
    xlabel::String="Energy",
    ylabel::String="phi",
    title::String="Gaps in Energy vs phi (coloured by (p,q))",
    verbose::Bool=false,
    atol::Real=1e-9,
    rtol::Real=0.0,
    fixed_values...
)
    # @assert :gap_labels in names(df) "df must have a :gap_labels column"
    # @assert :phi in names(df) "df must have a :phi column"

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
            for k in keys(Dict(fixed_values))
                if k in names(df)
                    vals = unique(df[:, k])
                    @info "Unique $k" n=length(vals) sample=first(vals, min(10, length(vals)))
                end
            end
        end
        error("No gaps match the specified fixed values.")
    end

    println("Number of rows after fixed_value filtering: ", nrow(filtered_rows))

    # Collect all gaps after per-gap filters, preserving row order
    lows  = Float64[]; highs = Float64[]; phis = Float64[]
    ps    = Int[];     qs    = Int[]
    pq_set = Set{Tuple{Int,Int}}()

    println("Collected gaps after pre-filtering.")

    for row in eachrow(filtered_rows)
        gls = getproperty(row, :gap_labels)
        (gls === nothing || isempty(gls)) && continue
        phi = Float64(row.phi)
        for g in gls
            p = Int(g.p); q = Int(g.q)
            if p_range !== nothing && !(p in p_range); continue; end
            if q_max !== nothing && abs(q) > q_max; continue; end
            push!(lows,  Float64(g.E_low))
            push!(highs, Float64(g.E_high))
            push!(phis,  phi)
            push!(ps,    p)
            push!(qs,    q)
            push!(pq_set, (p,q))
        end
    end

    if isempty(lows)
        error("No gaps remain after applying p_range/q_max/fixed_value filtering.")
    end

    println("Total gaps to plot: ", length(lows))

    # Color mapping by (p,q)
    pq_list = sort(collect(pq_set))
    pq_to_idx = Dict{Tuple{Int,Int},Int}((pq, i) for (i, pq) in enumerate(pq_list))
    # palette = cgrad(:viridis, length(pq_list); categorical=true).colors

    # Build palette from cmap argument
    if cmap === nothing
        palette = cgrad(:viridis, length(pq_list); categorical=true).colors
    elseif isa(cmap, Symbol) || isa(cmap, String)
        palette = cgrad(cmap, length(pq_list); categorical=true).colors
    elseif isa(cmap, AbstractVector)
        # if vector long enough, use entries directly; otherwise build gradient from stops
        if length(cmap) >= length(pq_list)
            palette = cmap[1:length(pq_list)]
        else
            palette = cgrad(cmap, length(pq_list); categorical=true).colors
        end
    else
        # fallback: try to build a gradient from cmap
        palette = cgrad(cmap, length(pq_list); categorical=true).colors
    end

    # Axes ranges
    xmin = minimum(lows); xmax = maximum(highs)
    ymin = minimum(phis); ymax = maximum(phis)
    x_span = xmax - xmin; x_span = x_span == 0 ? 1.0 : x_span
    y_span = ymax - ymin; y_span = y_span == 0 ? 1.0 : y_span

    # Band half-height in phi units
    h = max(1e-6, 0.5 * min(0.02 * y_span, (0.9 * y_span) / max(1, length(unique(phis)))))

    plt = plot(xlabel=xlabel, ylabel=ylabel, legend=false, grid=true, size=size, title=title)

    println("Plotting gaps...")

    # Draw each gap rectangle
    for i in eachindex(lows)
        idx = pq_to_idx[(ps[i], qs[i])]
        col = palette[idx]
        phi = phis[i]
        xs = [lows[i], highs[i], highs[i], lows[i], lows[i]]
        ys = [phi - h, phi - h, phi + h, phi + h, phi - h]
        plot!(plt, xs, ys; seriestype=:shape, fillcolor=col, fillalpha=0.85, linecolor=:transparent, label=false)
    end

    println("Adding legend annotations...")

    # Legend-like annotations (same layout as plateaus legend)
    right_margin = 0.40 * x_span
    plot!(plt, xlim=(xmin, xmax + right_margin))

    x_col1 = xmax + 0.03 * x_span
    x_col2 = xmax + 0.20 * x_span

    N1 = length(pq_list)
    y_start1 = ymax - 0.02 * y_span
    y_step1 = min(0.03 * y_span, (0.9 * y_span) / max(1, N1))
    ys1 = [y_start1 - (i-1) * y_step1 for i in 1:N1]
    for (i, pq) in enumerate(pq_list)
        idx = pq_to_idx[pq]
        scatter!(plt, [x_col1 - 0.01 * x_span], [ys1[i]];
                 markersize=6, markerstrokewidth=0, color=palette[idx], label=false)
        annotate!(plt, (x_col1, ys1[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    pq_by_q = sort(pq_list, by = x -> x[2])
    N2 = length(pq_by_q)
    y_start2 = ymax - 0.02 * y_span
    y_step2 = min(0.03 * y_span, (0.9 * y_span) / max(1, N2))
    ys2 = [y_start2 - (i-1) * y_step2 for i in 1:N2]
    for (i, pq) in enumerate(pq_by_q)
        idx = pq_to_idx[pq]
        scatter!(plt, [x_col2 - 0.01 * x_span], [ys2[i]];
                 markersize=6, markerstrokewidth=0, color=palette[idx], label=false)
        annotate!(plt, (x_col2, ys2[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
    end

    println("Saving figure to ", savepath)

    Plots.savefig(plt, savepath)
end

function plt_qled_coloured_gaps_energy_vs_phi(
    df::DataFrame,
    savepath::String;
    p_range::Union{AbstractVector{<:Integer},Nothing}=nothing,
    q_max::Union{Int,Nothing}=nothing,
    cmap::Any = :viridis,
    xlabel::String="Energy",
    ylabel::String="phi",
    title::String="Gaps in Energy vs phi (coloured by (p,q))",
    verbose::Bool=false,
    atol::Real=1e-9,
    rtol::Real=0.0,
    fixed_values...
)
    # @assert :gap_labels in names(df) "df must have a :gap_labels column"
    # @assert :phi in names(df) "df must have a :phi column"

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
            for k in keys(Dict(fixed_values))
                if k in names(df)
                    vals = unique(df[:, k])
                    @info "Unique $k" n=length(vals) sample=first(vals, min(10, length(vals)))
                end
            end
        end
        error("No gaps match the specified fixed values.")
    end

    # Collect all gaps after per-gap filters, preserving row order
    lows  = Float64[]; highs = Float64[]; phis = Float64[]
    ps    = Int[];     qs    = Int[]
    pq_set = Set{Tuple{Int,Int}}()
    q_to_ps = Dict{Int, Vector{Int}}()  # keep track of p values per q (order)

    for row in eachrow(filtered_rows)
        gls = getproperty(row, :gap_labels)
        (gls === nothing || isempty(gls)) && continue
        phi = Float64(row.phi)
        for g in gls
            p = Int(g.p); q = Int(g.q)
            if p_range !== nothing && !(p in p_range); continue; end
            if q_max !== nothing && abs(q) > q_max; continue; end
            push!(lows,  Float64(g.E_low))
            push!(highs, Float64(g.E_high))
            push!(phis,  phi)
            push!(ps,    p)
            push!(qs,    q)
            push!(pq_set, (p,q))
            q_to_ps[q] = get(q_to_ps, q, Int[])    # ensure vector exists
            if !(p in q_to_ps[q])
                push!(q_to_ps[q], p)
            end
        end
    end

    if isempty(lows)
        error("No gaps remain after applying p_range/q_max/fixed_value filtering.")
    end

    # Determine unique q ordering (stable/sorted)
    q_list = sort(collect(keys(q_to_ps)))
    n_q = length(q_list)

    # Build base colours for q from cmap
    base_colors = nothing
    if isa(cmap, Symbol) || isa(cmap, String)
        base_colors = cgrad(cmap, n_q; categorical=true).colors
    elseif isa(cmap, AbstractVector)
        # if user provided a vector of stops or colors
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
        # convert base to HSV
        base_hsv = HSV(base_colors[qi])
        # generate brightness/saturation variations for p. distribute values in [0.45,0.95]
        for (pi, pval) in enumerate(ps_for_q)
            t = n_pq == 1 ? 0.5 : (pi - 1) / (n_pq - 1)   # 0..1
            v_new = clamp(0.45 + 0.50 * t, 0.0, 1.0)
            s_new = clamp(0.6 + 0.35 * (1 - t), 0.0, 1.0)
            color_pq = RGB(HSV(base_hsv.h, s_new, v_new))
            pq_to_color[(pval, q)] = color_pq
        end
    end

    # Axes ranges
    xmin = minimum(lows); xmax = maximum(highs)
    ymin = minimum(phis); ymax = maximum(phis)
    x_span = xmax - xmin; x_span = x_span == 0 ? 1.0 : x_span
    y_span = ymax - ymin; y_span = y_span == 0 ? 1.0 : y_span

    # Band half-height in phi units
    h = max(1e-6, 0.5 * min(0.02 * y_span, (0.9 * y_span) / max(1, length(unique(phis)))))

    plt = plot(xlabel=xlabel, ylabel=ylabel, legend=false, grid=true, size=size, title=title)

    # Draw each gap rectangle using pq_to_color mapping
    for i in eachindex(lows)
        col = pq_to_color[(ps[i], qs[i])]
        phi = phis[i]
        xs = [lows[i], highs[i], highs[i], lows[i], lows[i]]
        ys = [phi - h, phi - h, phi + h, phi + h, phi - h]
        plot!(plt, xs, ys; seriestype=:shape, fillcolor=col, fillalpha=0.85, linecolor=:transparent, label=false)
    end

    # Legend-like annotations (same layout as plateaus legend) - use q-led ordering but show (p,q)
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

    savefig(plt, savepath)
    return plt
end