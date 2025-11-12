# using Plots
# using DataFrames


# function plt_discrete_phase_projections(
#     df::DataFrame,
#     x_variable::Symbol,
#     y_variable::Symbol,
#     disc_variable::Symbol,
#     savepath::String;
#     final_line_data=nothing,
#     final_poly_data=nothing,
#     fixed_values...
# )
#     """
#     Plot rows where a discretised indicator column marks a match.

#     - x_variable, y_variable : Symbols for the x and y axes (e.g. :mu, :phi, :mu_mp_norm).
#         * Can be scalar (Float64) or vector (Vector{Union{Missing, Float64}}) per row.
#     - disc_variable : Symbol of the discretised column to use (:disc_eigenvalues, :disc_mp, or :mp_disc_norm).
#         * Can be scalar (Int/Float) or vector (Vector{Int}) per row.
#         * For vectors, plots points where the element == 1, using corresponding indices from x_variable/y_variable if they are vectors.
#         * Skips points where x_variable or y_variable is missing.
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
#                     # Skip if missing
#                     if !ismissing(x_p) && !ismissing(y_p)
#                         push!(x_points, x_p)
#                         push!(y_points, y_p)
#                     end
#                 end
#             end
#         else
#             # Scalar case
#             if disc === missing
#                 continue
#             end
#             d = Int(round(disc))  # coerce numeric to Int (1/0)
#             if d == 1
#                 # Skip if missing
#                 if !ismissing(x_val) && !ismissing(y_val)
#                     push!(x_points, x_val)
#                     push!(y_points, y_val)
#                 end
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

#     # Overlay final regression line if provided
#     if final_line_data !== nothing
#         a, b = final_line_data.a, final_line_data.b
#         # Use the actual phi range from the data used to fit the line
#         phi_line = sort(final_line_data.phi_vals)
#         mu_line = a .* phi_line .+ b
#         # Plot the line in μ-φ space (or normalise if needed for your x-axis)
#         Plots.plot!(mu_line, phi_line; color=:red, linewidth=2, label="Final MP regression")
#     end

#     if final_poly_data !== nothing
#         p = final_poly_data.poly
#         phi_line = range(minimum(final_poly_data.phi_vals), maximum(final_poly_data.phi_vals), length=200)
#         mu_line = p.(phi_line)
#         Plots.plot!(mu_line, phi_line; color=:red, linewidth=2, label="Final MP poly fit")
#     end

#     Plots.savefig(savepath)
#     # Plots.display(plt)
# end

# # Similarly for plt_plateaus_vs_phi_idop, plt_plateaus_vs_phi_coloured_legend, plt_qled_coloured_gaps_mu_vs_phi
# # Add the same check after filtering in each:

# # In plt_plateaus_vs_phi_idop:
# filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
# if isempty(filtered_df)
#     @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
#     return
# end

# # In plt_plateaus_vs_phi_coloured_legend:
# filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
# if isempty(filtered_df)
#     @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
#     return
# end

# # In plt_qled_coloured_gaps_mu_vs_phi:
# filtered_df = filter(row -> row.Delta ≈ Delta, df)  # Or use fixed_values if updated
# if isempty(filtered_df)
#     @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
#     return
# end


# function plt_idop_vs_mu(df::DataFrame, savepath::String; fixed_values...)
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     plt = Plots.plot(
#         xlabel = "mu",
#         ylabel = "IDOP",
#         legend = :topright,
#         grid = true,
#         size = (800, 600),
#         title = "IDOP vs mu"
#     )

#     for (i, row) in enumerate(eachrow(filtered_df))
#         mu_vals = row.mu_values       # Vector{Float64}
#         idop_vals = row.idop          # Vector{Float64}
#         label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
#         Plots.plot!(plt, mu_vals, idop_vals, label = label_str)
#     end

#     Plots.savefig(savepath)
#     Plots.display(plt)
# end

# function plt_idop_plateaus_check(df::DataFrame, plateaus::Vector{Vector{Float64}}, savepath::String; fixed_values...)
#     """Plot IDOP vs mu for filtered rows and overlay horizontal lines at provided plateau values (per-row)."""
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end
#     if length(plateaus) != nrow(filtered_df)
#         error("Length of plateaus vector must match number of filtered rows.")
#     end

#     plt = Plots.plot(
#         xlabel = "mu",
#         ylabel = "IDOP",
#         legend = :topright,
#         grid = true,
#         size = (900, 600),
#         title = "IDOP vs mu with detected plateaus"
#     )

#     for (i, row) in enumerate(eachrow(filtered_df))
#         mu_vals = row.mu_values
#         idop_vals = row.idop
#         label_str = hasproperty(row, :sequence_name) ? row.sequence_name : "row $i"
#         Plots.plot!(plt, mu_vals, idop_vals, label = label_str)
#         for p in plateaus[i]
#             Plots.hline!(plt, [p], linestyle = :dot, color = :black, label = false)
#         end
#     end

#     Plots.savefig(savepath)
#     # Plots.display(plt)
# end

# function plt_plateaus_vs_phi_idop(df::DataFrame, savepath::String; fixed_values...)
#     """
#     Scatter plot of IDOP plateau values vs phi.
#     - df must contain columns :phi and :plateaus (Vector{Float64} per row).
#     - fixed_values... filters rows before plotting (same pattern as other helpers).
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     sorted_df = sort(filtered_df, :phi)

#     x_points = Float64[]
#     y_points = Float64[]

#     for row in eachrow(sorted_df)
#         phi = row.phi
#         for p in row.plateaus
#             push!(x_points, p)
#             push!(y_points, phi)
#         end
#     end

#     plt = Plots.scatter(
#         x_points, y_points;
#         markersize = 2.5,
#         markerstrokewidth = 0,
#         legend = false,
#         grid = true,
#         xlabel = "IDOP Plateau Value",
#         ylabel = "phi",
#         size = (900, 600),
#         title = "IDOP Plateaus vs phi"
#     )

#     Plots.savefig(savepath)
#     # Plots.display(plt)
# end


# ###########################################
# ########## Gap-labelled plotting ##########
# ###########################################

# """
# Plot IDOP plateaus vs. phi, colored by gap labels (p/q).
# """
# function plt_plateaus_vs_phi_coloured_legend(
#     df::DataFrame,
#     savepath::AbstractString;
#     fixed_values...
# )
#     # Filter by fixed values (flexible, like other plotting functions)
#     _kv_match(row, key, val) = begin
#         rv = getproperty(row, key)
#         (rv isa Number && val isa Number) ? isapprox(rv, val; atol=1e-8, rtol=1e-6) : (rv == val)
#     end
#     _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
#     filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end
    
#     x_points = Float64[]
#     y_points = Float64[]
#     colors = String[]
    
#     for row in eachrow(filtered_df)
#         phi = row.phi
#         for label in row.gap_labels
#             push!(x_points, label.plateau)
#             push!(y_points, phi)
#             push!(colors, "$(label.p)/$(label.q)")  # Color by p/q label
#         end
#     end
    
#     plt = Plots.scatter(
#         x_points, y_points,
#         group=colors,
#         markersize=2.5,
#         markerstrokewidth=0,
#         legend=:outerright,
#         grid=true,
#         xlabel="IDOP Plateau Value",
#         ylabel="phi",
#         size=(900,600),
#         title="IDOP Plateaus vs phi (Gap-Labeled)"
#     )
#     Plots.savefig(savepath)
# end


# ############
# ## unused ##
# ############
# # """
# # Plot colored gaps (plateaus) in energy space vs. phi.
# # Adapts IDOS logic: assumes plateaus represent "gaps" at quantized levels.
# # """
# # function plt_coloured_gaps_energy_vs_phi(
# #     df::DataFrame,
# #     savepath::AbstractString;
# #     cmap=:viridis,
# #     mu::Float64=NaN,
# #     Delta::Float64=NaN,
# #     phason::Float64=NaN,
# #     atol=1e-8,
# #     rtol=1e-6,
# #     verbose=false
# # )
# #     # Simplified: Plot plateau values as "energy" levels vs. phi, colored by p/q
# #     # (Full adaptation would mirror IDOS, but this aligns with IDOP structures)
# #     filtered_df = df
# #     phi_vals = Float64[]
# #     plateau_vals = Float64[]
# #     labels = String[]
    
# #     for row in eachrow(filtered_df)
# #         for label in row.gap_labels
# #             push!(phi_vals, row.phi)
# #             push!(plateau_vals, label.plateau)
# #             push!(labels, "$(label.p)/$(label.q)")
# #         end
# #     end
    
# #     plt = Plots.scatter(
# #         phi_vals, plateau_vals,
# #         group=labels,
# #         colormap=cmap,
# #         markersize=3,
# #         legend=:outerright,
# #         xlabel="phi",
# #         ylabel="Plateau Value (Energy-like)",
# #         title="Gap-Labeled Plateaus vs phi",
# #         size=(900,600)
# #     )
# #     Plots.savefig(savepath)
# # end

# function plt_qled_coloured_gaps_mu_vs_phi(
#     df::DataFrame,
#     savepath::AbstractString;
#     mu_axis::Symbol = :mu_mp_norm,  # Choose :mu_values or :mu_mp_norm for x-axis
#     cmap=:RdBu,
#     atol=1e-8,
#     rtol=1e-6,
#     verbose=false,
#     fixed_values...
# )
#     # Filter by fixed values (flexible, like other plotting functions)
#     _kv_match(row, key, val) = begin
#         rv = getproperty(row, key)
#         (rv isa Number && val isa Number) ? isapprox(rv, val; atol=atol, rtol=rtol) : (rv == val)
#     end
#     _row_matches(row, kvs) = all(_kv_match(row, k, v) for (k,v) in kvs)
    
#     filtered_df = isempty(fixed_values) ? df : filter(row -> _row_matches(row, fixed_values), df)
#     if isempty(filtered_df)
#         if verbose
#             @info "No rows match fixed_values" keys=keys(Dict(fixed_values))
#             for k in keys(Dict(fixed_values))
#                 if k in names(df)
#                     vals = unique(df[:, k])
#                     @info "Unique $k" n=length(vals) sample=first(vals, min(10, length(vals)))
#                 end
#             end
#         end
#         error("No gaps match the specified fixed values.")
#     end

#     # Collect all gaps
#     mu_starts = Float64[]
#     mu_ends = Float64[]
#     phi_vals = Float64[]
#     ps = Int[]
#     qs = Int[]
#     pq_set = Set{Tuple{Int,Int}}()
#     q_to_ps = Dict{Int, Vector{Int}}()

#     for row in eachrow(filtered_df)
#         phi = row.phi
#         gap_labels = row.gap_labels  # Vector{NamedTuple} with :plateau, :p, :q, :err
#         idop_vals = row.idop  # Vector of IDOP values
#         mu_vals = row[mu_axis]  # Vector of μ values (raw or normalized), may contain Missing
        
#         for label in gap_labels
#             plateau_val = label.plateau
#             p = label.p
#             q = label.q
            
#             # Find μ ranges where IDOP ≈ plateau_val (flat regions)
#             ranges = find_idop_ranges(idop_vals, plateau_val, atol)
#             for (start_idx, end_idx) in ranges
#                 mu_start = mu_vals[start_idx]
#                 mu_end = mu_vals[end_idx]
                
#                 # Skip if any μ value is missing
#                 if ismissing(mu_start) || ismissing(mu_end)
#                     continue
#                 end
                
#                 push!(mu_starts, mu_start)
#                 push!(mu_ends, mu_end)
#                 push!(phi_vals, phi)
#                 push!(ps, p)
#                 push!(qs, q)
#                 push!(pq_set, (p, q))
#                 q_to_ps[q] = get(q_to_ps, q, Int[])
#                 if !(p in q_to_ps[q])
#                     push!(q_to_ps[q], p)
#                 end
#             end
#         end
#     end

#     if isempty(mu_starts)
#         @warn "No valid gap labels to plot (possibly due to missing μ values)."
#         return
#     end

#     # Determine unique q ordering (stable/sorted)
#     q_list = sort(collect(keys(q_to_ps)))
#     n_q = length(q_list)

#     # Build base colours for q from cmap
#     base_colors = nothing
#     if isa(cmap, Symbol) || isa(cmap, String)
#         base_colors = cgrad(cmap, n_q; categorical=true).colors
#     elseif isa(cmap, AbstractVector)
#         if length(cmap) >= n_q
#             base_colors = cmap[1:n_q]
#         else
#             base_colors = cgrad(cmap, n_q; categorical=true).colors
#         end
#     else
#         base_colors = cgrad(:viridis, n_q; categorical=true).colors
#     end

#     # Build final color for each (p,q): same hue per q, vary value by index of p within that q
#     pq_to_color = Dict{Tuple{Int,Int}, Colorant}()
#     for (qi, q) in enumerate(q_list)
#         ps_for_q = sort(q_to_ps[q])   # deterministic ordering of p within q
#         n_pq = length(ps_for_q)
#         base_hsv = HSV(base_colors[qi])
#         for (pi, pval) in enumerate(ps_for_q)
#             t = n_pq == 1 ? 0.5 : (pi - 1) / (n_pq - 1)   # 0..1
#             v_new = clamp(0.45 + 0.50 * t, 0.0, 1.0)
#             s_new = clamp(0.6 + 0.35 * (1 - t), 0.0, 1.0)
#             color_pq = RGB(HSV(base_hsv.h, s_new, v_new))
#             pq_to_color[(pval, q)] = color_pq
#         end
#     end

#     # Axes ranges
#     xmin = minimum(mu_starts); xmax = maximum(mu_ends)
#     ymin = minimum(phi_vals); ymax = maximum(phi_vals)
#     x_span = xmax - xmin; x_span = x_span == 0 ? 1.0 : x_span
#     y_span = ymax - ymin; y_span = y_span == 0 ? 1.0 : y_span

#     # Band half-height in phi units
#     h = max(1e-6, 0.5 * min(0.02 * y_span, (0.9 * y_span) / max(1, length(unique(phi_vals)))))

#     plt = plot(xlabel=mu_axis == :mu_mp_norm ? "Normalized mu" : "mu", ylabel="phi", legend=false, grid=true, size=(900,600), title="Q-Led Colored Gaps in $(mu_axis == :mu_mp_norm ? "Normalized " : "")mu vs phi (IDOP)")

#     # Draw each gap rectangle using pq_to_color mapping
#     for i in eachindex(mu_starts)
#         col = pq_to_color[(ps[i], qs[i])]
#         phi = phi_vals[i]
#         xs = [mu_starts[i], mu_ends[i], mu_ends[i], mu_starts[i], mu_starts[i]]
#         ys = [phi - h, phi - h, phi + h, phi + h, phi - h]
#         plot!(plt, xs, ys; seriestype=:shape, fillcolor=col, fillalpha=0.85, linecolor=:transparent, label=false)
#     end

#     # Legend-like annotations (same layout as IDOS)
#     pq_list = sort(collect(pq_set), by = x -> (x[2], x[1]))  # sort by q then p
#     right_margin = 0.40 * x_span
#     plot!(plt, xlim=(xmin, xmax + right_margin))

#     x_col1 = xmax + 0.03 * x_span
#     x_col2 = xmax + 0.20 * x_span

#     # Column 1: pq_list order
#     N1 = length(pq_list)
#     y_start1 = ymax - 0.02 * y_span
#     y_step1 = min(0.03 * y_span, (0.9 * y_span) / max(1, N1))
#     ys1 = [y_start1 - (i-1) * y_step1 for i in 1:N1]
#     for (i, pq) in enumerate(pq_list)
#         col = pq_to_color[pq]
#         scatter!(plt, [x_col1 - 0.01 * x_span], [ys1[i]];
#                  markersize=6, markerstrokewidth=0, color=col, label=false)
#         annotate!(plt, (x_col1, ys1[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
#     end

#     # Column 2: ordered by q
#     pq_by_q = sort(pq_list, by = x -> x[2])
#     N2 = length(pq_by_q)
#     y_start2 = ymax - 0.02 * y_span
#     y_step2 = min(0.03 * y_span, (0.9 * y_span) / max(1, N2))
#     ys2 = [y_start2 - (i-1) * y_step2 for i in 1:N2]
#     for (i, pq) in enumerate(pq_by_q)
#         col = pq_to_color[pq]
#         scatter!(plt, [x_col2 - 0.01 * x_span], [ys2[i]];
#                  markersize=6, markerstrokewidth=0, color=col, label=false)
#         annotate!(plt, (x_col2, ys2[i], Plots.text(" (p=$(pq[1]), q=$(pq[2]))", 8, :left)))
#     end

#     Plots.savefig(savepath)
# end

# # Helper function to find ranges where IDOP ≈ plateau_val
# function find_idop_ranges(idop_vals::Vector{Float64}, plateau_val::Float64, atol::Float64)
#     ranges = Tuple{Int, Int}[]
#     in_range = false
#     start_idx = 0
#     for (i, val) in enumerate(idop_vals)
#         if isapprox(val, plateau_val; atol=atol)
#             if !in_range
#                 start_idx = i
#                 in_range = true
#             end
#         else
#             if in_range
#                 push!(ranges, (start_idx, i-1))
#                 in_range = false
#             end
#         end
#     end
#     if in_range
#         push!(ranges, (start_idx, length(idop_vals)))
#     end
#     return ranges
# end







using Plots
using DataFrames

## Updated to not error on empty sets
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