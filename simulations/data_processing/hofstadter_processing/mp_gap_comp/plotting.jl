module MPGapPlotting

using Plots
using DataFrames
using Printf
using Colors
using BSON
using Statistics

## Colour Colorscheme
function q_anchor_colormap(qs::Vector{Int})
    lerp(c1::RGB, c2::RGB, t) = RGB((1-t)*c1.r + t*c2.r,
                                    (1-t)*c1.g + t*c2.g,
                                    (1-t)*c1.b + t*c2.b)

    qs_sorted = sort(qs)
    cmap = fill(RGB(0.7,0.7,0.7), length(qs_sorted))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs_sorted)
    idx_pos1     = findfirst(==(  1), qs_sorted)
    idx_pos_high = length(qs_sorted)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    if !isnothing(idx_neg_high) && idx_neg_high <= length(cmap); cmap[idx_neg_high] = col_neg_high; end
    if !isnothing(idx_pos_high) && idx_pos_high <= length(cmap); cmap[idx_pos_high] = col_pos_high; end
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end

    # Special case: Set q=0 to grey
    idx_zero = findfirst(==(0), qs_sorted)
    if !isnothing(idx_zero)
        cmap[idx_zero] = RGB(0.5, 0.5, 0.5)  # Grey
    end

    if !isnothing(idx_neg1) && idx_neg1 > idx_neg_high
        span = idx_neg1 - idx_neg_high
        for (j,i) in enumerate(idx_neg_high+1:idx_neg1-1)
            t = j / span
            if t < 1/3
                cmap[i] = lerp(col_neg_high, col_neg_purp, t/(1/3))
            elseif t < 2/3
                cmap[i] = lerp(col_neg_purp, col_neg_dark, (t-1/3)/(1/3))
            else
                cmap[i] = lerp(col_neg_dark, col_neg_light, (t-2/3)/(1/3))
            end
        end
    end

    if !isnothing(idx_pos1) && idx_pos_high > idx_pos1
        span = idx_pos_high - idx_pos1
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high-1)
            t = j / span
            if t < 0.5
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.5)
            else
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.5)/0.5)
            end
        end
    end

    if !isnothing(idx_neg1) && !isnothing(idx_pos1) && idx_pos1 - idx_neg1 > 1
        span = idx_pos1 - idx_neg1
        for (j,i) in enumerate(idx_neg1+1:idx_pos1-1)
            t = j / span
            cmap[i] = lerp(col_neg_light, RGB(1.0,0.0,0.0), t)
        end
    end

    cmap
end

## fully functional for energy and mp gaps labelled + raw_mp plot
function plt_mp_gap_comparison(
    df_slice::DataFrame,
    savepath::AbstractString;
    gap_segments::Union{Nothing,DataFrame}=nothing,
    mp_gap_segments::Union{Nothing,DataFrame}=nothing,
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu,
    mp_col::Symbol=:mp,
    disc_mp_col::Symbol=:disc_mp,
    gap_labels_col::Symbol=:gap_labels,
    full_range::Bool=true,
    q_range::Vector{Int}=collect(-5:5), # Default fixed range
    plot_disc_mp::Bool=true,
    palette = cgrad(:thermal, 11)
)
    # Helper to check column existence safely (Symbol vs String)
    has_col(df, col_sym) = string(col_sym) in names(df)

    # 1. Sort and Extract Basic Data
    df_sorted = sort(df_slice, mu_col)
    mu_vals = df_sorted[!, mu_col]
    eig_sets = df_sorted[!, eig_col]
    num_eigs = length(first(eig_sets))

    # 2. Determine Y-Limits for Eigenvalues
    if full_range
        ylim_eig = nothing
    else
        # mu_zero_idx = argmin(abs.(mu_vals))
        # max_at_zero = maximum(eig_sets[mu_zero_idx])
        # min_at_zero = minimum(eig_sets[mu_zero_idx])
        # ylim_eig = (min_at_zero * 1.2, max_at_zero * 1.2)

        ylim_eig=(-3.0,3.0)
    end

    # 3. Extract Gap Labels
    
    # A) Energy Gaps (Vertical at mu=0)
    energy_gaps = NamedTuple[]
    
    if gap_segments !== nothing && nrow(gap_segments) > 0 && has_col(gap_segments, gap_labels_col)
        # println("Debug: Extracting energy gaps from gap_segments...")
        gs_sorted = sort(gap_segments, mu_col)
        mu_vals_gs = gs_sorted[!, mu_col]
        zero_idx = argmin(abs.(mu_vals_gs))
        gls = getproperty(gs_sorted[zero_idx, :], gap_labels_col)
        if gls !== nothing && !ismissing(gls)
            energy_gaps = gls
        end
    elseif has_col(df_sorted, gap_labels_col)
        # println("Debug: Extracting energy gaps from df_slice...")
        mu_zero_idx = argmin(abs.(mu_vals))
        zero_mu_row = df_sorted[mu_zero_idx, :]
        gls = getproperty(zero_mu_row, gap_labels_col)
        if gls !== nothing && !ismissing(gls)
            energy_gaps = gls
        end
    else
        # println("Debug: No energy gap source found.")
    end
    # println("Debug: Found $(length(energy_gaps)) energy gaps.")

    # B) MP Gaps (Horizontal along mu axis)
    mp_gaps = NamedTuple[]
    if mp_gap_segments !== nothing && nrow(mp_gap_segments) > 0 && has_col(mp_gap_segments, gap_labels_col)
        # println("Debug: Extracting MP gaps...")
        # Access the column directly
        raw_col = mp_gap_segments[!, gap_labels_col]
        if !isempty(raw_col)
            raw_val = raw_col[1] # Get first row's value
            
            if raw_val !== nothing && !ismissing(raw_val)
                if isa(raw_val, Vector) && !isempty(raw_val)
                    if isa(raw_val[1], NamedTuple)
                        mp_gaps = raw_val
                    elseif isa(raw_val[1], Vector)
                        # Handle nested vector case: Vector{Vector{NamedTuple}}
                        println("Debug: Unpacking nested MP gaps...")
                        mp_gaps = raw_val[1]
                    end
                end
            end
        end
    else
        # println("Debug: No MP gap source found or empty.")
    end
    # println("Debug: Found $(length(mp_gaps)) MP gaps.")

    # 4. Prepare Colors
    # # Modified to return nothing if q is 0, effectively ignoring q=0 gaps
    # clean_q(q) = (ismissing(q) || !isfinite(q) || round(Int, q) == 0) ? nothing : Int(round(q))
    # Modified to include q=0 (does not ignore it)
    clean_q(q) = (ismissing(q) || !isfinite(q)) ? nothing : Int(round(q))

    # Use the fixed q_range for the colormap to ensure consistency across plots
    # # Filter out 0 from the range if present, as we don't plot q=0
    # fixed_qs = sort(filter(x -> x != 0, unique(q_range)))
    # Include q=0 now
    fixed_qs = sort(unique(q_range))
    
    cmap = q_anchor_colormap(fixed_qs)
    q_to_color = Dict(q => cmap[i] for (i, q) in enumerate(fixed_qs))

    # Identify which qs are actually present in this specific plot for the legend
    qs_energy = [clean_q(g.q) for g in energy_gaps]
    qs_mp = [clean_q(g.q) for g in mp_gaps]
    present_qs = sort(unique(Int.(filter(!isnothing, vcat(qs_energy, qs_mp)))))

    # 5. Setup Plot Dimensions
    gap_label_width = 0.1
    min_mu = minimum(mu_vals)
    max_mu = maximum(mu_vals)
    xlims_val = (min_mu - 1.5*gap_label_width, max_mu)

    # 6. Create Eigenvalue Plot
    eig_plot = plot(; 
        xlabel="", 
        ylabel="Eigenvalues", 
        legend=:topright, 
        xticks=:auto,
        minorgrid=:x,
        xlims=xlims_val,
        title="Eigen Spectrum vs μ",
        # aspect_ratio=0.5
    )
    
    if ylim_eig !== nothing
        plot!(eig_plot; ylim=ylim_eig)
    end

    # Plot Eigenvalues
    for j in 1:num_eigs
        y_vals = [eig_vec[j] for eig_vec in eig_sets]
        plot!(eig_plot, mu_vals, y_vals; color=:steelblue, alpha=0.45, linewidth=1.4, label=false)
    end

    # Plot Discretized MP Points
    if plot_disc_mp &&has_col(df_sorted, disc_mp_col)
        disc_idxs = findall(x -> x == 1, df_sorted[!, disc_mp_col])
        scatter!(eig_plot, mu_vals[disc_idxs], zeros(length(disc_idxs)); 
            color=:red, markersize=1.5, marker=:circle, alpha=1.0, markerstrokewidth=0, label=false)
    end

    # 7. Plot Energy Gaps (Vertical bars)
    for gap in energy_gaps
        q = clean_q(gap.q)
        if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
            color = q_to_color[q]
            y1, y2 = minmax(gap.E_low, gap.E_high)
            plot!(eig_plot, Shape([min_mu - gap_label_width, min_mu, min_mu, min_mu - gap_label_width], [y1, y1, y2, y2]);
                fillcolor=color, alpha=0.7, strokewidth=0, label=false)
        end
    end

    # 8. Plot MP Gaps (Horizontal bars)
    mp_bar_height = 0.05
    for gap in mp_gaps
        q = clean_q(gap.q)
        if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
            color = q_to_color[q]
            x1, x2 = minmax(gap.E_low, gap.E_high)
            plot!(eig_plot, Shape([x1, x2, x2, x1], [-mp_bar_height, -mp_bar_height, mp_bar_height, mp_bar_height]);
                fillcolor=color, alpha=0.8, strokewidth=0, label=false)
        end
    end

    # 9. Add Legend
    # Always add legend entries for the full q_range to ensure consistency across plots
    if !isempty(fixed_qs)
        for q in fixed_qs
            plot!(eig_plot, [NaN], [NaN]; color=q_to_color[q], label="q = $q")
        end
    end

    # 10. Create MP Plot
    mp_vals = real.(df_sorted[!, mp_col])
    mp_plot = plot(mu_vals, mp_vals; 
        xlabel="μ", 
        ylabel="MP", 
        ylim=(-1.05, 0.05), 
        xlims=xlims_val,
        legend=false, 
        minorgrid=:x,
        color=:forestgreen, 
        linewidth=1.8, 
        title="Majorana Polarization vs μ"
    )

    layout = @layout [a{0.75h}; b]
    combined = plot(eig_plot, mp_plot; layout=layout, size=(900, 850))
    # combined = plot(eig_plot, mp_plot; layout=layout, size=(1200, 900), link=:x)
    savefig(combined, savepath)
    
    return combined
end

# ## v1 of trying to add mu vs mp_tol plot
# function plt_mp_gap_comparison(
#     df_slice::DataFrame,
#     savepath::AbstractString;
#     gap_segments::Union{Nothing,DataFrame}=nothing,
#     mp_gap_segments::Union{Nothing,DataFrame}=nothing,
#     eig_col::Symbol=:eigenvalues,
#     mu_col::Symbol=:mu,
#     mp_col::Symbol=:mp,
#     disc_mp_col::Symbol=:disc_mp,
#     gap_labels_col::Symbol=:gap_labels,
#     full_range::Bool=true,
#     q_range::Vector{Int}=collect(-5:5),
#     plot_disc_mp::Bool=true,
#     mp_sweep_matrix::Union{Nothing, Matrix}=nothing,
#     mp_sweep_tols::Union{Nothing, Vector}=nothing,
#     palette = cgrad(:thermal, 11)
# )
#     # Helper to check column existence safely (Symbol vs String)
#     has_col(df, col_sym) = string(col_sym) in names(df)

#     # 1. Sort and Extract Basic Data
#     df_sorted = sort(df_slice, mu_col)
#     mu_vals = df_sorted[!, mu_col]
#     eig_sets = df_sorted[!, eig_col]
#     num_eigs = length(first(eig_sets))

#     # 2. Determine Y-Limits for Eigenvalues
#     if full_range
#         ylim_eig = nothing
#     else
#         ylim_eig=(-3.0,3.0)
#     end

#     # 3. Extract Gap Labels
    
#     # A) Energy Gaps (Vertical at mu=0)
#     energy_gaps = NamedTuple[]
    
#     if gap_segments !== nothing && nrow(gap_segments) > 0 && has_col(gap_segments, gap_labels_col)
#         gs_sorted = sort(gap_segments, mu_col)
#         mu_vals_gs = gs_sorted[!, mu_col]
#         zero_idx = argmin(abs.(mu_vals_gs))
#         gls = getproperty(gs_sorted[zero_idx, :], gap_labels_col)
#         if gls !== nothing && !ismissing(gls)
#             energy_gaps = gls
#         end
#     elseif has_col(df_sorted, gap_labels_col)
#         mu_zero_idx = argmin(abs.(mu_vals))
#         zero_mu_row = df_sorted[mu_zero_idx, :]
#         gls = getproperty(zero_mu_row, gap_labels_col)
#         if gls !== nothing && !ismissing(gls)
#             energy_gaps = gls
#         end
#     end

#     # B) MP Gaps (Horizontal along mu axis)
#     mp_gaps = NamedTuple[]
#     if mp_gap_segments !== nothing && nrow(mp_gap_segments) > 0 && has_col(mp_gap_segments, gap_labels_col)
#         raw_col = mp_gap_segments[!, gap_labels_col]
#         if !isempty(raw_col)
#             raw_val = raw_col[1] # Get first row's value
            
#             if raw_val !== nothing && !ismissing(raw_val)
#                 if isa(raw_val, Vector) && !isempty(raw_val)
#                     if isa(raw_val[1], NamedTuple)
#                         mp_gaps = raw_val
#                     elseif isa(raw_val[1], Vector)
#                         # Handle nested vector case: Vector{Vector{NamedTuple}}
#                         mp_gaps = raw_val[1]
#                     end
#                 end
#             end
#         end
#     end

#     # 4. Prepare Colors
#     clean_q(q) = (ismissing(q) || !isfinite(q) || round(Int, q) == 0) ? nothing : Int(round(q))
    
#     fixed_qs = sort(filter(x -> x != 0, unique(q_range)))
#     cmap = q_anchor_colormap(fixed_qs)
#     q_to_color = Dict(q => cmap[i] for (i, q) in enumerate(fixed_qs))

#     # 5. Setup Plot Dimensions
#     gap_label_width = 0.1
#     min_mu = minimum(mu_vals)
#     max_mu = maximum(mu_vals)
#     xlims_val = (min_mu - 1.5*gap_label_width, max_mu)

#     # 6. Create Eigenvalue Plot
#     eig_plot = plot(; 
#         xlabel="", 
#         ylabel="Eigenvalues", 
#         legend=:topright, 
#         xticks=:auto,
#         minorgrid=:x,
#         xlims=xlims_val,
#         title="Eigen Spectrum vs μ",
#     )
    
#     if ylim_eig !== nothing
#         plot!(eig_plot; ylim=ylim_eig)
#     end

#     # Plot Eigenvalues
#     for j in 1:num_eigs
#         y_vals = [eig_vec[j] for eig_vec in eig_sets]
#         plot!(eig_plot, mu_vals, y_vals; color=:steelblue, alpha=0.45, linewidth=1.4, label=false)
#     end

#     # Plot Discretized MP Points
#     if plot_disc_mp && has_col(df_sorted, disc_mp_col)
#         disc_idxs = findall(x -> x == 1, df_sorted[!, disc_mp_col])
#         scatter!(eig_plot, mu_vals[disc_idxs], zeros(length(disc_idxs)); 
#             color=:red, markersize=1.5, marker=:circle, alpha=1.0, markerstrokewidth=0, label=false)
#     end

#     # 7. Plot Energy Gaps (Vertical bars)
#     for gap in energy_gaps
#         q = clean_q(gap.q)
#         if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
#             color = q_to_color[q]
#             y1, y2 = minmax(gap.E_low, gap.E_high)
#             plot!(eig_plot, Shape([min_mu - gap_label_width, min_mu, min_mu, min_mu - gap_label_width], [y1, y1, y2, y2]);
#                 fillcolor=color, alpha=0.7, strokewidth=0, label=false)
#         end
#     end

#     # 8. Plot MP Gaps (Horizontal bars)
#     mp_bar_height = 0.05
#     for gap in mp_gaps
#         q = clean_q(gap.q)
#         if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
#             color = q_to_color[q]
#             x1, x2 = minmax(gap.E_low, gap.E_high)
#             plot!(eig_plot, Shape([x1, x2, x2, x1], [-mp_bar_height, -mp_bar_height, mp_bar_height, mp_bar_height]);
#                 fillcolor=color, alpha=0.8, strokewidth=0, label=false)
#         end
#     end

#     # 9. Add Legend
#     if !isempty(fixed_qs)
#         for q in fixed_qs
#             plot!(eig_plot, [NaN], [NaN]; color=q_to_color[q], label="q = $q")
#         end
#     end

#     # 10. Create MP Plot
#     mp_vals = real.(df_sorted[!, mp_col])
#     mp_plot = plot(mu_vals, mp_vals; 
#         xlabel="", # Remove xlabel here as it will be on the bottom plot
#         ylabel="MP", 
#         ylim=(-1.05, 0.05), 
#         xlims=xlims_val,
#         legend=false, 
#         minorgrid=:x,
#         color=:forestgreen, 
#         linewidth=1.8, 
#         title="Majorana Polarization vs μ"
#     )

#     # 11. Create MP Sweep Plot (Heatmap)
#     if mp_sweep_matrix !== nothing && mp_sweep_tols !== nothing
#         # mp_sweep_matrix is (n_mu, n_tols)
#         # We want x-axis = mu, y-axis = log10(tolerance)
        
#         # Transpose for heatmap: rows=y (tol), cols=x (mu)
#         heatmap_data = permutedims(mp_sweep_matrix) # Now (n_tols, n_mu)
        
#         y_vals = log10.(mp_sweep_tols)
        
#         sweep_plot = heatmap(mu_vals, y_vals, heatmap_data;
#             xlabel="μ",
#             ylabel="log10(Tol)",
#             color=cgrad([:white, :red], 2, categorical=true), # 0=white, 1=red
#             legend=false,
#             xlims=xlims_val,
#             title="Discretized MP vs Tolerance"
#         )
        
#         # Combine 3 plots
#         layout = @layout [a{0.6h}; b{0.2h}; c{0.2h}]
#         combined = plot(eig_plot, mp_plot, sweep_plot; layout=layout, size=(900, 1000), link=:x)
#     else
#         # Fallback to 2 plots if no sweep data
#         layout = @layout [a{0.75h}; b]
#         combined = plot(eig_plot, mp_plot; layout=layout, size=(900, 850), link=:x)
#     end

#     savefig(combined, savepath)
    
#     return combined
# end

# ## mu vs mp_tol chekcing plot (plus identifying appropriate mp_tol range)
# function plt_mu_vs_mptol_check(
#     sweep_input::Union{String,Dict};
#     N::Int,
#     mode::Symbol = :full,               # :full | :fixed | :test
#     fixed_tols::Union{Nothing, Vector{Float64}} = nothing,
#     test_index::Union{Nothing, Int} = nothing,
#     rational_tol::Float64 = 1e-6,       # tolerance for rational approximation of phi
#     min_gap_length::Int = 1,            # min consecutive zeros to count as a gap
#     savepath::Union{Nothing, String} = nothing
# )
#     # Load sweep data if a filename provided
#     sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
#     @assert !isnothing(sweep) "sweep data/file must be provided"

#     mu_vals = sweep[:mu]
#     mp_tol_range = sweep[:mp_tol_range]
#     disc_mp_matrix = sweep[:disc_mp_matrix]         # (n_mu, n_tols)
#     params = get(sweep, :params, nothing)
#     phi_val = params !== nothing ? params[:phi] : get(sweep, :phi, nothing)

#     n_mu = length(mu_vals)
#     n_tols = length(mp_tol_range)
#     @assert size(disc_mp_matrix, 1) == n_mu "disc_mp_matrix rows must match mu length"

#     # Determine expected number of gaps: j = denominator of rational approximation of phi
#     j = nothing
#     if !isnothing(phi_val)
#         try
#             r = rationalize(phi_val, rational_tol)
#             j = denominator(r)
#         catch
#             j = nothing
#         end
#     end
#     expected_gaps = isnothing(j) ? (2N - 1) : min(2N - 1, Int(j))

#     # Helper to count zero-runs in a vector (with minimum length)
#     function count_zero_runs(vec::AbstractVector{<:Integer})
#         cnt = 0
#         runlen = 0
#         for v in vec
#             if v == 0
#                 runlen += 1
#             else
#                 if runlen >= min_gap_length
#                     cnt += 1
#                 end
#                 runlen = 0
#             end
#         end
#         # trailing run
#         if runlen >= min_gap_length
#             cnt += 1
#         end
#         return cnt
#     end

#     # Compute number of gaps for each tolerance
#     gap_counts = Vector{Int}(undef, n_tols)
#     for i in 1:n_tols
#         col = disc_mp_matrix[:, i]
#         gap_counts[i] = count_zero_runs(col)
#     end

#     # Decide which tolerances to highlight based on mode
#     highlight_runs = []  # will collect (i_start, i_end) ranges of tol indices to highlight
#     if mode == :full
#         # find contiguous ranges where gap_counts == expected_gaps
#         i = 1
#         while i <= n_tols
#             if gap_counts[i] == expected_gaps
#                 starti = i
#                 while i <= n_tols && gap_counts[i] == expected_gaps
#                     i += 1
#                 end
#                 push!(highlight_runs, (starti, i - 1))
#             else
#                 i += 1
#             end
#         end
#     elseif mode == :fixed && !isnothing(fixed_tols)
#         # highlight each fixed tol if it yields expected_gaps
#         for (i, tol) in enumerate(mp_tol_range)
#             if tol in fixed_tols && gap_counts[i] == expected_gaps
#                 push!(highlight_runs, (i, i))
#             end
#         end
#     elseif mode == :test && !isnothing(test_index)
#         ti = test_index
#         if ti >= 1 && ti <= n_tols && gap_counts[ti] == expected_gaps
#             push!(highlight_runs, (ti, ti))
#         end
#     else
#         # unsupported combination => no highlights
#     end

#     # Prepare heatmap data: Plots.heatmap expects z dims = (length(y), length(x))
#     heatmap_z = permutedims(disc_mp_matrix)  # (n_tols, n_mu)
#     y_vals = log10.(mp_tol_range)            # plot tolerances on log scale

#     # Base heatmap: use extrema of viridis (dark blue for 1 (positive discretisation, so MP~-1.0 if mp_targ=-1.0), yellow for 0)
#     cmap = cgrad(:viridis, 2, categorical=true, rev=true)

#     p = heatmap(mu_vals, y_vals, heatmap_z;
#         xlabel="μ",
#         ylabel="log10(mp_tol)",
#         color=cmap,
#         xlims=(minimum(mu_vals), maximum(mu_vals)),
#         yticks = :auto,
#         title = "μ vs mp_tol (disc MP)",
#         aspect_ratio = :auto,
#         legend=true
#     )

#     # Overlay highlighted tol windows as translucent horizontal bands
#     xmin = minimum(mu_vals)
#     xmax = maximum(mu_vals)
#     for (i1, i2) in highlight_runs
#         tol_low = mp_tol_range[i1]
#         tol_high = mp_tol_range[i2]
#         y1 = log10(tol_low)
#         y2 = log10(tol_high)
#         # draw translucent rectangle
#         plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
#             fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.25, strokewidth=0, label=false)
#         # add text label in center of band
#         midy = (y1 + y2) / 2
#         midx = (xmin + xmax) / 2
#         ann = @sprintf("tol idx %d-%d\ntol %.3g-%.3g\ngaps=%d", i1, i2, tol_low, tol_high, expected_gaps)
#         annotate!(p, midx, midy, text(ann, :black, 8, :center))
#     end

#     # If test_index provided, also add side-panel showing mp_disc vs mu for that tol
#     if mode == :test && !isnothing(test_index)
#         ti = test_index
#         if 1 <= ti <= n_tols
#             vall = disc_mp_matrix[:, ti]
#             # create small plot of mp_disc along mu
#             p2 = plot(mu_vals, vall;
#                 xlabel="μ",
#                 ylabel="disc_mp",
#                 ylim = (-0.1, 1.1),
#                 title = @sprintf("disc_mp @ tol idx %d (tol=%.3g), gaps=%d", ti, mp_tol_range[ti], gap_counts[ti]),
#                 legend=false,
#                 color=:black
#             )
#             layout = @layout [a; b{0.25h}]
#             combined = plot(p, p2; layout=layout, size=(1000,800), link=:x)
#             if !isnothing(savepath)
#                 savefig(combined, savepath)
#             end
#             return combined
#         end
#     end

#     # annotate overall info: expected gaps, j (if found)
#     info_str = isnothing(j) ? @sprintf("expected_gaps = min(2N-1, j) with j unknown; using 2N-1=%d", 2N-1) :
#                               @sprintf("expected_gaps = min(2N-1, j=%d) = %d", j, expected_gaps)
#     annotate!(p, minimum(mu_vals), maximum(y_vals), text(info_str, :left, 8, :black))

#     if !isnothing(savepath)
#         savefig(p, savepath)
#     end

#     return p
# end


## mu vs mp_tol chekcing plot (plus annotating appropriate mp_tol range -- processing done before plotting)
 ## bad, but functional, legend layout
function plt_mu_vs_mptol_check(
    sweep_input::Union{String,Dict};
    N::Int,
    mode::Symbol = :full,               # :full | :fixed | :test
    fixed_tols::Union{Nothing, Vector{Float64}} = nothing,
    test_index::Union{Nothing, Int} = nothing,
    savepath::Union{Nothing, String} = nothing
)
    # Load sweep data if a filename provided
    sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
    @assert !isnothing(sweep) "sweep data/file must be provided"

    mu_vals = sweep[:mu]
    mp_tol_range = sweep[:mp_tol_range]
    disc_mp_matrix = sweep[:disc_mp_matrix]         # (n_mu, n_tols)
    params = get(sweep, :params, nothing)
    mp_targ = params !== nothing ? params[:mp_targ] : get(sweep, :mp_targ, nothing)
    
    # Extract Analysis Data
    analysis = get(sweep, :analysis, nothing)
    
    if !isnothing(analysis)
        j = analysis[:j]
        expected_gaps = analysis[:expected_gaps]
        gap_counts = analysis[:gap_counts]
        valid_indices = analysis[:valid_indices]
        mu_c_vec = get(analysis, :mu_c_vec, nothing) # Extract mu_c vector
    else
        @warn "Analysis data missing in BSON. Plotting without gap info."
        j = nothing
        expected_gaps = -1
        gap_counts = zeros(Int, length(mp_tol_range))
        valid_indices = Int[]
        mu_c_vec = nothing
    end

    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)

    # 1. Calculate representative scalar mu_c (median of valid indices)
    mu_c_scalar = nothing
    if !isnothing(mu_c_vec)
        # Use valid indices if available, otherwise use all
        relevant_mu_c = !isempty(valid_indices) ? mu_c_vec[valid_indices] : mu_c_vec
        if !isempty(relevant_mu_c)
            mu_c_scalar = median(relevant_mu_c)
        end
    end

    # Decide which tolerances to highlight based on mode
    highlight_runs = []
    
    if mode == :full
        if !isempty(valid_indices)
            sort!(valid_indices)
            start_i = valid_indices[1]
            prev_i = valid_indices[1]
            for i in valid_indices[2:end]
                if i == prev_i + 1
                    prev_i = i
                else
                    push!(highlight_runs, (start_i, prev_i))
                    start_i = i
                    prev_i = i
                end
            end
            push!(highlight_runs, (start_i, prev_i))
        end
    elseif mode == :fixed && !isnothing(fixed_tols)
        for (i, tol) in enumerate(mp_tol_range)
            if tol in fixed_tols && (i in valid_indices)
                push!(highlight_runs, (i, i))
            end
        end
    elseif mode == :test && !isnothing(test_index)
        ti = test_index
        if ti >= 1 && ti <= n_tols && (ti in valid_indices)
            push!(highlight_runs, (ti, ti))
        end
    end

    # Prepare heatmap data
    heatmap_z = permutedims(disc_mp_matrix)  # (n_tols, n_mu)
    y_vals = log10.(mp_tol_range)            # plot tolerances on log scale

    # Base heatmap colors
    # rev=true: 0 -> Yellow (Gap), 1 -> Dark (MBS)
    cmap = cgrad(:viridis, 2, categorical=true, rev=true)
    col_mbs = RGB(0.267, 0.004, 0.329) # Viridis start (Dark)
    col_gap = RGB(0.993, 0.906, 0.144) # Viridis end (Yellow)

    p = heatmap(mu_vals, y_vals, heatmap_z;
        xlabel="μ",
        ylabel="log10(mp_tol)",
        color=cmap,
        ylims=(minimum(y_vals), maximum(y_vals)),
        xlims=(minimum(mu_vals), maximum(mu_vals)),
        yticks = :auto,
        title = "μ vs mp_tol (disc MP)",
        aspect_ratio = :auto,
        colorbar=false,         # Disable colorbar
        legend=:outerbottom      # Legend outside
    )

    # --- Custom Legend Entries ---
    # 1. mu_c line
    if !isnothing(mu_c_scalar)
        vline!(p, [mu_c_scalar], color=:red, linestyle=:dot, linewidth=2, label="μ_c")
    end
    # 2. MBS Color
    mbs_label = isnothing(mp_targ) ? "MBS" : @sprintf("MBS (MP ≈ %.1f)", mp_targ)
    scatter!(p, [NaN], [NaN], color=col_mbs, label=mbs_label, shape=:square, markersize=5, legend_font_pointsize=8)
    # 3. Gap Color
    scatter!(p, [NaN], [NaN], color=col_gap, label="Gap", shape=:square, markersize=5)

    # Overlay highlighted tol windows (visual boxes only)
    xmin = minimum(mu_vals)
    xmax = maximum(mu_vals)
    for (i1, i2) in highlight_runs
        tol_low = mp_tol_range[i1]
        tol_high = mp_tol_range[i2]
        y1 = log10(tol_low)
        y2 = log10(tol_high)
        
        plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
            fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.25, strokewidth=0, label=false)
    end

    # --- Add Annotation Directly to Main Plot ---
    # Position at bottom left (adjust x/y as needed for exact placement)
    lines = String[]
    push!(lines, @sprintf("Parameters: N=%d, j=%s", N, isnothing(j) ? "?" : string(j)))
    push!(lines, @sprintf("Expected Gaps: %d", expected_gaps))
    if !isnothing(mu_c_scalar)
        push!(lines, @sprintf("μ_c ≈ %.4f", mu_c_scalar))
    end
    
    if !isempty(highlight_runs)
        push!(lines, "Valid Tolerance Ranges:")
        for (i1, i2) in highlight_runs
            t1 = mp_tol_range[i1]
            t2 = mp_tol_range[i2]
            push!(lines, @sprintf("  Idx %d-%d: %.1e to %.1e", i1, i2, t1, t2))
        end
    else
        push!(lines, "No tolerance ranges found matching expected gaps.")
    end
    
    info_text = join(lines, "\n")
    
    # Annotate at bottom left of plot area (x=minimum(mu_vals), y=minimum(y_vals) - offset for below)
    offset_y = 0.2 * (maximum(y_vals) - minimum(y_vals))  # Adjust offset to place below plot
    annotate!(p, minimum(mu_vals), minimum(y_vals) - offset_y, text(info_text, :left, 8, :black))

    # Return the single plot (no layout needed)
    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end

end  # module MPGapPlotting











