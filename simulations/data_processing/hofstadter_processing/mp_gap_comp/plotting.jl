module MPGapPlotting

using Plots
using DataFrames
using Printf
using Colors
using BSON
using Statistics
using LaTeXStrings
using Measures

## Colour Colorscheme
function q_anchor_colormap(qs::Vector{Int})
    lerp(c1::RGB, c2::RGB, t) = RGB((1-t)*c1.r + t*c2.r,
                                    (1-t)*c1.g + t*c2.g,
                                    (1-t)*c1.b + t*c2.b)

    qs_sorted = sort(qs)
    
    # --- SENTINEL HANDLING ---
    SC_GAP_Q = 9999
    has_sentinel = SC_GAP_Q in qs_sorted
    
    # Filter out sentinel for gradient calculation
    qs_phys = filter(!=(SC_GAP_Q), qs_sorted)
    # -------------------------

    cmap = fill(RGB(0.7,0.7,0.7), length(qs_phys))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs_phys)
    idx_pos1     = findfirst(==(  1), qs_phys)
    idx_pos_high = length(qs_phys)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    cmap[idx_neg_high] = col_neg_high
    
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end

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
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high)
            t = j / span
            if t < 0.33
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.33)
            elseif t < 0.66
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.33)/0.33)
            else
                cmap[i] = lerp(col_pos_yell, col_pos_high, (t-0.66)/0.34)
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

    # --- RECONSTRUCT MAP WITH SENTINEL ---
    # Map q -> color
    color_dict = Dict(q => c for (q,c) in zip(qs_phys, cmap))
    
    if has_sentinel
        color_dict[SC_GAP_Q] = RGB(0.85, 0.85, 0.85) # Light Grey for SC gaps
    end
    
    # Return colors in the order of the original sorted input
    return [color_dict[q] for q in qs_sorted]
end

## fully functional for energy and mp gaps labelled + raw_mp plot
function plt_mp_gap_comparison(
    df_slice::DataFrame,
    savepath::AbstractString;
    gap_segments::Union{Nothing,DataFrame}=nothing,
    mp_gap_segments::Union{Nothing,DataFrame}=nothing,
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu,
    _mp_col::Symbol=:mp,
    disc_mp_col::Symbol=:disc_mp,
    gap_labels_col::Symbol=:gap_labels,
    full_range::Bool=true,
    mu_range_lim::Union{Nothing, Tuple{Float64,Float64}}=nothing,
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

    if mu_range_lim !== nothing
        mu_min, mu_max = mu_range_lim
        valid_idxs = findall(x -> (x >= mu_min) && (x <= mu_max), mu_vals)
        df_sorted = df_sorted[valid_idxs, :]
        mu_vals = df_sorted[!, mu_col]
        eig_sets = df_sorted[!, eig_col]
        num_eigs = length(first(eig_sets))
    end

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
                        # println("Debug: Unpacking nested MP gaps...")
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
        yguidefontrotation=0,
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
        yguidefontrotation=0,
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

# Including mp_tol plot as third pannel
function plt_mp_gap_comparison_full(
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
    mu_range_lim::Union{Nothing, Tuple{Float64,Float64}}=nothing,
    q_range::Vector{Int}=collect(-5:5), # Default fixed range
    plot_disc_mp::Bool=true,
    mp_sweep_matrix::Union{Nothing, Matrix}=nothing,
    mp_sweep_tols::Union{Nothing, Vector}=nothing,
    axis_label_fontsize::Int=20,
    axis_number_fontsize::Int=14,
    raw_mp_tol::Float64=5e-2,
    raw_mp_target::Float64= -1.0,
    mp_gaps_top::Bool=true
)
    # Helper to check column existence safely (Symbol vs String)
    has_col(df, col_sym) = string(col_sym) in names(df)

    # 1. Sort and Extract Basic Data
    perm = sortperm(df_slice[!, mu_col])
    df_sorted = df_slice[perm, :]
    
    sweep_mat_sorted = nothing
    if mp_sweep_matrix !== nothing
        if size(mp_sweep_matrix, 1) == nrow(df_slice)
            sweep_mat_sorted = mp_sweep_matrix[perm, :]
        else
            @warn "mp_sweep_matrix rows $(size(mp_sweep_matrix, 1)) != df rows $(nrow(df_slice)). Ignoring sweep data."
        end
    end

    mu_vals = df_sorted[!, mu_col]
    eig_sets = df_sorted[!, eig_col]
    num_eigs = length(first(eig_sets))

    # Calculate mu_critical (last point where MP is within tolerance of target)
    mu_critical = nothing
    if has_col(df_sorted, mp_col)
        mp_calc_vals = real.(df_sorted[!, mp_col])
        in_range_mask = abs.(mp_calc_vals .- raw_mp_target) .< raw_mp_tol
        if any(in_range_mask)
            # The critical mu is the maximum mu where the condition holds
            mu_critical = maximum(mu_vals[in_range_mask])
        end
    end

    if mu_range_lim !== nothing
        mu_min, mu_max = mu_range_lim
        valid_idxs = findall(x -> (x >= mu_min) && (x <= mu_max), mu_vals)
        df_sorted = df_sorted[valid_idxs, :]
        mu_vals = df_sorted[!, mu_col]
        eig_sets = df_sorted[!, eig_col]
        num_eigs = length(first(eig_sets))
        if sweep_mat_sorted !== nothing; sweep_mat_sorted = sweep_mat_sorted[valid_idxs, :]; end
    end

    # 2. Determine Y-Limits for Eigenvalues
    if full_range
        ylim_eig = nothing
    else
        ylim_eig=(-3.0,3.0)
    end

    # 3. Extract Gap Labels
    energy_gaps = NamedTuple[]
    if gap_segments !== nothing && nrow(gap_segments) > 0 && has_col(gap_segments, gap_labels_col)
        gs_sorted = sort(gap_segments, mu_col)
        mu_vals_gs = gs_sorted[!, mu_col]
        zero_idx = argmin(abs.(mu_vals_gs))
        gls = getproperty(gs_sorted[zero_idx, :], gap_labels_col)
        if gls !== nothing && !ismissing(gls); energy_gaps = gls; end
    elseif has_col(df_sorted, gap_labels_col)
        mu_zero_idx = argmin(abs.(mu_vals))
        zero_mu_row = df_sorted[mu_zero_idx, :]
        gls = getproperty(zero_mu_row, gap_labels_col)
        if gls !== nothing && !ismissing(gls); energy_gaps = gls; end
    end

    mp_gaps = NamedTuple[]
    if mp_gap_segments !== nothing && nrow(mp_gap_segments) > 0 && has_col(mp_gap_segments, gap_labels_col)
        raw_col = mp_gap_segments[!, gap_labels_col]
        if !isempty(raw_col)
            raw_val = raw_col[1]
            if raw_val !== nothing && !ismissing(raw_val)
                if isa(raw_val, Vector) && !isempty(raw_val)
                    if isa(raw_val[1], NamedTuple); mp_gaps = raw_val;
                    elseif isa(raw_val[1], Vector); mp_gaps = raw_val[1]; end
                end
            end
        end
    end

    # 4. Prepare Colors
    clean_q(q) = (ismissing(q) || !isfinite(q)) ? nothing : Int(round(q))
    # Helper for latex formatting axis numbers
    latex_num_fmt(x) = (abs(x - round(x)) < 1e-5) ? L"%$(Int(round(x)))" : L"%$(round(x, digits=2))"

    fixed_qs = sort(unique(q_range))
    cmap = q_anchor_colormap(fixed_qs)
    q_to_color = Dict(q => cmap[i] for (i, q) in enumerate(fixed_qs))

    # 5. Setup Plot Dimensions
    gap_label_width = 0.1
    min_mu = minimum(mu_vals)
    max_mu = maximum(mu_vals)
    xlims_val = (min_mu - 1.5*gap_label_width, max_mu)

    # 6. Create Eigenvalue Plot
    eig_plot = plot(; 
        xlabel="", ylabel=L"E", yguidefontrotation=0, legend=false, xticks=:auto, minorgrid=:x,
        xlims=xlims_val, bottom_margin = -1.0mm,
        xformatter = _ -> "",
        yformatter = latex_num_fmt,
        framestyle=:box,
        guidefontsize=axis_label_fontsize,
        tickfontsize=axis_number_fontsize
    )
    if ylim_eig !== nothing; plot!(eig_plot; ylim=ylim_eig); end

    for j in 1:num_eigs
        y_vals = [eig_vec[j] for eig_vec in eig_sets]
        plot!(
            eig_plot, 
            mu_vals, 
            y_vals; 
            color=:steelblue, 
            alpha=0.45, 
            linewidth=1.4, 
            label=false
        )
    end

    if mu_critical !== nothing
        vline!(eig_plot, [mu_critical], color=:red, linestyle=:dash, label="")
    end

    if plot_disc_mp && has_col(df_sorted, disc_mp_col)
        disc_idxs = findall(x -> x == 1, df_sorted[!, disc_mp_col])
        scatter!(
            eig_plot, 
            mu_vals[disc_idxs], 
            zeros(length(disc_idxs)); 
            color=:red,
            markersize=1.5, 
            marker=:circle, 
            alpha=1.0, 
            markerstrokewidth=0,
            label=false
        )
    end

    # 7 & 8. Plot Gaps
    for gap in energy_gaps
        q = clean_q(gap.q)
        if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
            y1, y2 = minmax(gap.E_low, gap.E_high)
            plot!(eig_plot, Shape([min_mu - gap_label_width, min_mu, min_mu, min_mu - gap_label_width], [y1, y1, y2, y2]);
                fillcolor=q_to_color[q], alpha=1.0, strokewidth=0, label=false)
        end
    end

    if mp_gaps_top
        if ylim_eig !== nothing
            y_top = ylim_eig[2]
        else
            y_top = maximum(maximum(e) for e in eig_sets)
        end
        bar_h = 0.2
        y_bottom_edge = y_top - bar_h
        y_top_edge = y_top
    else
        bar_h = 0.1
        y_bottom_edge = -bar_h
        y_top_edge = bar_h
    end

    for gap in mp_gaps
        q = clean_q(gap.q)
        if q !== nothing && haskey(q_to_color, q) && isfinite(gap.E_low) && isfinite(gap.E_high)
            x1, x2 = minmax(gap.E_low, gap.E_high)
            plot!(eig_plot, Shape([x1, x2, x2, x1], [y_bottom_edge, y_bottom_edge, y_top_edge, y_top_edge]);
                fillcolor=q_to_color[q], alpha=1.0, strokewidth=0, label=false)
        end
    end

    # 9. (Legend/Colorbar handled via inset in step 12)

    # 10. Create MP Plot
    # Formatter to only show 0 and -1 labels, hiding others but keeping ticks
    mp_y_formatter = y -> (abs(y) < 0.05 || abs(y + 1.0) < 0.05) ? L"%$(round(Int, y))" : ""

    mp_vals = real.(df_sorted[!, mp_col])
    mp_plot = plot(mu_vals, mp_vals; 
        xlabel="", ylabel=L"\mathcal{M}", yguidefontrotation=0, ylim=(-1.05, 0.05), xlims=xlims_val,
        legend=false, minorgrid=:x, 
        color=:black, #:forestgreen, 
        linewidth=1.8, 
        bottom_margin = -1.0mm,
        yformatter = mp_y_formatter, # Apply custom formatter
        xformatter = _ -> "",
        framestyle=:box,
        guidefontsize=axis_label_fontsize,
        tickfontsize=axis_number_fontsize
    )
    hline!(mp_plot, [raw_mp_target+raw_mp_tol], color=:red, linestyle=:dot, label="")
    if mu_critical !== nothing
        vline!(mp_plot, [mu_critical], color=:red, linestyle=:dash, label="")
    end

    # 11. Create Sweep Plot & Layout
    println("sweep_mat_sorted: ", sweep_mat_sorted === nothing ? "nothing" : "present")
    if sweep_mat_sorted !== nothing && mp_sweep_tols !== nothing
        println("Creating MP Tolerance Sweep Plot...")
        y_vals_sweep = log10.(mp_sweep_tols)
        z_vals_sweep = permutedims(sweep_mat_sorted)
        cmap_sweep = cgrad(:viridis, 2, categorical=true, rev=true, alpha=0.7)
        
        # Formatter: show 0 and -10
        sweep_y_fmt = y -> (abs(y) < 0.1 || abs(y + 10.0) < 0.1) ? L"%$(Int(round(y)))" : ""

        sweep_plot = heatmap(mu_vals, y_vals_sweep, z_vals_sweep;
            xlabel=L"\mu'", ylabel=L"\log{\epsilon}", yguidefontrotation=0, color=cmap_sweep,
            xlims=xlims_val,
            guidefontsize=axis_label_fontsize,
            tickfontsize=axis_number_fontsize, 
            ylims=(-3.0, maximum(y_vals_sweep)),
            legend=false, colorbar=false,
            yformatter=sweep_y_fmt,
            xformatter=latex_num_fmt,
            framestyle=:box,
        )
        if mu_critical !== nothing
            vline!(sweep_plot, [mu_critical], color=:red, linestyle=:dash, label="")
        end

        data_layout = @layout [a{0.7h}; b{0.15h}; c{0.15h}]
        combined = plot(eig_plot, mp_plot, sweep_plot; layout=data_layout, link=:x, size=(900, 1200))
    else
        # Formatting for mp_plot x-axis as it is the bottom plot
        plot!(mp_plot, xformatter=latex_num_fmt, xlabel=L"\mu", bottom_margin=3Plots.mm)

        data_layout = @layout [a{0.75h}; b]
        combined = plot(eig_plot, mp_plot; layout=data_layout, link=:x, size=(900, 850))
    end

    # 12. Add Inset Colorbar (Overlay on Eig Plot)
    # Exclude 0 from colorbar as it is meaningless
    cb_qs = filter(!=(0), fixed_qs)
    
    if !isempty(cb_qs)
        n_qs = length(cb_qs)
        c_vals = [q_to_color[q] for q in cb_qs]
        cmap_q = cgrad(c_vals, n_qs, categorical=true)
        
        # Calculate specific ticks: min, -1, 1, max
        target_qs = unique(sort([minimum(cb_qs), -1, 1, maximum(cb_qs)]))
        # Keep only targets that actually exist in the data
        show_qs = filter(x -> x in cb_qs, target_qs)
        tick_idxs = [findfirst(==(q), cb_qs) for q in show_qs]
        tick_labs = string.(show_qs)

        # Adjustable Inset Parameters (relative to Subplot 1)
        # To change position: adjust cb_x (horizontal, 0-1) and cb_y (vertical, 0-1)
        cb_w = 0.025  # Slightly wider to include labels in the box
        cb_h = 0.45
        cb_x = 0.95
        cb_y = 0.025
        
        # Add inset to subplot 1 (eig_plot)
        # We add it to the 'combined' plot.
        # Subplot index for the inset will be length(combined.subplots) + 1
        next_sub = length(combined.subplots) + 1
        
        # Create inset with white background for the WHOLE subplot area
        plot!(combined, inset=(1, bbox(cb_x, cb_y, cb_w, cb_h)), subplot=next_sub,
              bg_subplot=:white,    # White background for the entire inset box (including labels)
              bg_inside=:white,     # White background for the plot area
              framestyle=:box,      # Black border around the heatmap axis
              grid=false, 
              xticks=false, yticks=false, minorticks=false
        )
        
        # Plot Heatmap in Inset
        z_mat = reshape(1:n_qs, n_qs, 1)
        
        heatmap!(combined, [1], 1:n_qs, z_mat, subplot=next_sub,
                 color=cmap_q, clims=(1, n_qs+1),
                 xticks=false,
                 yticks=(tick_idxs, tick_labs),
                 colorbar=false,
                #  title=nothing, #L"q", 
                 titlefontsize=10,
                 tickfontsize=8,
                 right_margin=3Plots.mm, # Padding inside the box for labels
                 left_margin=1Plots.mm,
                 bottom_margin=1Plots.mm,
                 top_margin=1Plots.mm
        )
    end

    savefig(combined, savepath)
    return combined
end


function plt_eigs_and_mptol_with_sc_gap_condition(
    df_slice::DataFrame,
    sweep_dict::Dict,
    savepath::AbstractString;
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu,
    mp_col::Symbol=:mp,
    axis_label_fontsize::Int=18,
    axis_number_fontsize::Int=12,
    plot_full_y_range::Bool=false,
    energy_gap_ranges_override=nothing
)
    @assert nrow(df_slice) > 0 "df_slice must not be empty"
    df_sorted = sort(df_slice, mu_col)
    mu_vals = df_sorted[!, mu_col]
    eig_sets = df_sorted[!, eig_col]
    n_mu = length(mu_vals)
    num_eigs = length(first(eig_sets))
    params = get(sweep_dict, :params, nothing)
    @assert params !== nothing "sweep_dict missing :params"
    delta_val = get(params, :Delta, nothing)
    @assert delta_val !== nothing ":Delta not found in params"
    mp_targ_val = get(params, :mp_targ, nothing)
    mp_tol_range = sweep_dict[:mp_tol_range]
    disc_mp_matrix = sweep_dict[:disc_mp_matrix]
    heatmap_mu = sweep_dict[:mu]
    analysis = get(sweep_dict, :analysis, Dict{Symbol,Any}())
    gap_counts = haskey(analysis, :gap_counts) ? Vector{Int}(analysis[:gap_counts]) : fill(0, length(mp_tol_range))
    expected_gaps = get(analysis, :expected_gaps, -1)
    valid_indices = haskey(analysis, :valid_indices) ? Vector{Int}(analysis[:valid_indices]) : Int[]
    mu_c_vec = get(analysis, :mu_c_vec, nothing)
    mu_c_scalar = get(analysis, :mu_c_value, nothing)
    energy_gap_counts = get(analysis, :energy_gap_counts, nothing)
    energy_gap_ranges = energy_gap_ranges_override === nothing ? get(analysis, :energy_gap_ranges, nothing) : energy_gap_ranges_override
    if mu_c_scalar === nothing && mu_c_vec !== nothing
        relevant_mu_c = !isempty(valid_indices) ? mu_c_vec[valid_indices] : mu_c_vec
        if !isempty(relevant_mu_c)
            mu_c_scalar = median(relevant_mu_c)
        end
    end

    zero_idx = argmin(abs.(mu_vals))
    zero_eigs = zero_idx <= n_mu ? collect(real.(df_sorted[zero_idx, eig_col])) : Float64[]

    overlap_tol = 1e-6
    zero_filter_tol = 1e-8

    # Normalise gap ranges depending on whether we plot the full symmetric window or positive branch only
    function normalise_gap_range(low::Float64, high::Float64)
        l, h = minmax(low, high)
        if !plot_full_y_range
            if h < -zero_filter_tol
                return nothing
            end
            l = max(l, 0.0)
            h = max(h, l)
            if h - l < zero_filter_tol
                return nothing
            end
        end
        return (l, h)
    end

    normalised_energy_ranges = Vector{Vector{Tuple{Float64, Float64}}}(undef, n_mu)
    for idx in 1:n_mu
        normalised_energy_ranges[idx] = Tuple{Float64, Float64}[]
        ranges_i = (energy_gap_ranges isa AbstractVector && idx <= length(energy_gap_ranges)) ? energy_gap_ranges[idx] : nothing
        if ranges_i === nothing || ismissing(ranges_i)
            continue
        end
        for range_i in ranges_i
            if range_i === nothing
                continue
            end
            low_i, high_i = range_i
            norm_range = normalise_gap_range(low_i, high_i)
            if norm_range !== nothing
                push!(normalised_energy_ranges[idx], norm_range)
            end
        end
    end

    zero_gap_ranges = zero_idx <= length(normalised_energy_ranges) ? copy(normalised_energy_ranges[zero_idx]) : Tuple{Float64, Float64}[]
    sort!(zero_gap_ranges, by=first)

    raw_mu_segments = Tuple{Float64, Float64}[]
    if !isempty(zero_gap_ranges)
        @assert length(normalised_energy_ranges) == n_mu "energy_gap_ranges length must match μ grid"
        @assert length(heatmap_mu) == n_mu "μ grid mismatch between analysis and plotting data"

        for (gap_low, gap_high) in zero_gap_ranges
            start_idx = nothing
            for idx in 1:n_mu
                ranges_i = normalised_energy_ranges[idx]
                gap_present = any(range_i -> min(range_i[2], gap_high) >= max(range_i[1], gap_low) - overlap_tol, ranges_i)

                if gap_present
                    if start_idx === nothing
                        start_idx = idx
                    end
                elseif start_idx !== nothing
                    mu_start = heatmap_mu[start_idx]
                    mu_end = heatmap_mu[idx - 1]
                    push!(raw_mu_segments, (min(mu_start, mu_end), max(mu_start, mu_end)))
                    start_idx = nothing
                end
            end
            if start_idx !== nothing
                mu_start = heatmap_mu[start_idx]
                mu_end = heatmap_mu[end]
                push!(raw_mu_segments, (min(mu_start, mu_end), max(mu_start, mu_end)))
            end
        end
    end

    sort!(raw_mu_segments, by=first)
    mu_gap_segments = Tuple{Float64, Float64}[]
    for seg in raw_mu_segments
        if isempty(mu_gap_segments)
            push!(mu_gap_segments, seg)
        else
            last_seg = mu_gap_segments[end]
            if seg[1] <= last_seg[2]
                mu_gap_segments[end] = (last_seg[1], max(last_seg[2], seg[2]))
            else
                push!(mu_gap_segments, seg)
            end
        end
    end

    gap_color = RGBA(0.85, 0.33, 0.1, 0.6)
    gap_color_solid = RGB(0.85, 0.33, 0.1)
    gap_label_width = 0.12
    min_mu = minimum(mu_vals)
    max_mu = maximum(mu_vals)
    full_eig_min = minimum(minimum(real.(e)) for e in eig_sets)
    full_eig_max = maximum(maximum(real.(e)) for e in eig_sets)
    if !plot_full_y_range && !isempty(zero_eigs)
        plot_eig_min = minimum(zero_eigs)
        plot_eig_max = maximum(zero_eigs)
    else
        plot_eig_min = full_eig_min
        plot_eig_max = full_eig_max
    end
    eig_range = max(plot_eig_max - plot_eig_min, 1e-6)
    bar_h = 0.12 * eig_range
    y_bottom_edge = plot_eig_max
    y_top_edge = plot_eig_max + bar_h
    xlims_combined = (min_mu - 1.5 * gap_label_width, max_mu)
    ylims_combined = (plot_eig_min - 0.05 * eig_range, plot_eig_max + 1.2 * bar_h)

    eig_plot = plot(
        xlabel="",
        ylabel=L"E",
        yguidefontrotation=0,
        legend=false,
        xticks=:auto,
        minorgrid=:x,
        xlims=xlims_combined,
        ylim=ylims_combined,
        framestyle=:box,
        guidefontsize=axis_label_fontsize,
        tickfontsize=axis_number_fontsize
    )

    for j in 1:num_eigs
        y_vals = [real(eig_sets[i][j]) for i in 1:n_mu]
        plot!(eig_plot, mu_vals, y_vals; color=:steelblue, alpha=0.45, linewidth=1.4, label=false)
    end

    for (_, (E_low, E_high)) in enumerate(zero_gap_ranges)
        plot!(eig_plot,
            Shape([min_mu - gap_label_width, min_mu, min_mu, min_mu - gap_label_width],
                  [E_low, E_low, E_high, E_high]);
            fillcolor=gap_color_solid,
            alpha=0.9,
            strokewidth=0,
            label=false
        )
    end

    for (mu_start, mu_end) in mu_gap_segments
        plot!(eig_plot,
            Shape([mu_start, mu_end, mu_end, mu_start],
                  [y_bottom_edge, y_bottom_edge, y_top_edge, y_top_edge]);
            fillcolor=gap_color,
            strokewidth=0,
            label=false
        )
    end

    annotate!(eig_plot, xlims_combined[1], y_top_edge, text(@sprintf("Δ = %.3g", delta_val), :left, 8, :black))

    highlight_runs = Tuple{Int, Int}[]
    if !isempty(valid_indices)
        sort!(valid_indices)
        start_i = valid_indices[1]
        prev_i = start_i
        for idx in valid_indices[2:end]
            if idx == prev_i + 1
                prev_i = idx
            else
                push!(highlight_runs, (start_i, prev_i))
                start_i = idx
                prev_i = idx
            end
        end
        push!(highlight_runs, (start_i, prev_i))
    end

    col_mbs = RGB(0.267, 0.004, 0.329)
    col_gap = RGB(0.993, 0.906, 0.144)
    y_vals_heat = log10.(mp_tol_range)
    heatmap_z = permutedims(disc_mp_matrix)
    cmap = cgrad(:viridis, 2, categorical=true, rev=true)
    heatmap_plot = heatmap(heatmap_mu, y_vals_heat, heatmap_z;
        xlabel=L"\mu",
        ylabel=L"\log \epsilon",
        yguidefontrotation=0,
        color=cmap,
        ylims=(minimum(y_vals_heat), maximum(y_vals_heat)),
        xlims=xlims_combined,
        yticks=:auto,
        legend=:outerbottom,
        framestyle=:box,
        guidefontsize=axis_label_fontsize,
        tickfontsize=axis_number_fontsize,
        colorbar=false
    )
    mbs_label = isnothing(mp_targ_val) ? "MBS" : @sprintf("MBS (MP ≈ %.1f)", mp_targ_val)
    scatter!(heatmap_plot, [NaN], [NaN]; color=col_mbs, label=mbs_label, shape=:square, markersize=5)
    scatter!(heatmap_plot, [NaN], [NaN]; color=col_gap, label="Gap", shape=:square, markersize=5)

    xmin = minimum(heatmap_mu)
    xmax = maximum(heatmap_mu)
    for (i1, i2) in highlight_runs
        tol_low = mp_tol_range[i1]
        tol_high = mp_tol_range[i2]
        y1 = log10(tol_low)
        y2 = log10(tol_high)
        plot!(heatmap_plot,
            Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
            fillcolor=gap_color,
            alpha=0.25,
            strokewidth=0,
            label=false
        )
    end

    if mu_c_scalar === nothing
        mu_c_scalar = full_eig_max
    end

    if mu_c_scalar !== nothing
        vline!(heatmap_plot, [mu_c_scalar]; color=:red, linestyle=:dot, linewidth=2, label="μ_c")
    end

    info_lines = String[]
    if expected_gaps >= 0
        push!(info_lines, @sprintf("|ΔE|>Δ gaps: %d", expected_gaps))
    else
        push!(info_lines, "|ΔE|>Δ gaps: ?")
    end
    if energy_gap_counts !== nothing && !isempty(energy_gap_counts)
        eg_min = minimum(energy_gap_counts)
        eg_max = maximum(energy_gap_counts)
        if eg_min == eg_max
            push!(info_lines, @sprintf("per-μ gaps: %d", eg_min))
        else
            push!(info_lines, @sprintf("per-μ gaps range: %d-%d", eg_min, eg_max))
        end
    end
    if mu_c_scalar !== nothing
        push!(info_lines, @sprintf("μ_c ≈ %.4f", mu_c_scalar))
    end
    annotate!(heatmap_plot, xlims_combined[1], maximum(y_vals_heat), text(join(info_lines, "\n"), :left, 8, :black))

    layout = @layout [a{0.65h}; b]
    combined = plot(eig_plot, heatmap_plot; layout=layout, size=(900, 950), link=:x)
    savefig(combined, savepath)
    # return combined
end


## mu vs mp_tol chekcing plot (plus identifying appropriate mp_tol range)
function plt_mu_vs_mptol_check_and_identify(
    sweep_input::Union{String,Dict};
    N::Int,
    mode::Symbol = :full,               # :full | :fixed | :test
    fixed_tols::Union{Nothing, Vector{Float64}} = nothing,
    test_index::Union{Nothing, Int} = nothing,
    rational_tol::Float64 = 1e-6,       # tolerance for rational approximation of phi
    min_gap_length::Int = 1,            # min consecutive zeros to count as a gap
    savepath::Union{Nothing, String} = nothing
)
    # Load sweep data if a filename provided
    sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
    @assert !isnothing(sweep) "sweep data/file must be provided"

    mu_vals = sweep[:mu]
    mp_tol_range = sweep[:mp_tol_range]
    disc_mp_matrix = sweep[:disc_mp_matrix]         # (n_mu, n_tols)
    params = get(sweep, :params, nothing)
    phi_val = params !== nothing ? params[:phi] : get(sweep, :phi, nothing)
    analysis = get(sweep, :analysis, nothing)

    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)
    @assert size(disc_mp_matrix, 1) == n_mu "disc_mp_matrix rows must match mu length"

    # Determine expected gaps from analysis if available; fallback to rational approximation
    j = nothing
    expected_gaps = nothing
    energy_gap_counts = nothing
    if analysis !== nothing
        if haskey(analysis, :expected_gaps)
            expected_gaps = analysis[:expected_gaps]
        end
        if haskey(analysis, :energy_gap_counts)
            energy_gap_counts = analysis[:energy_gap_counts]
        end
        if haskey(analysis, :j)
            j = analysis[:j]
        end
    end
    if expected_gaps === nothing
        if !isnothing(phi_val)
            try
                r = rationalize(phi_val, rational_tol)
                j = denominator(r)
            catch
                j = nothing
            end
        end
        expected_gaps = isnothing(j) ? (2 * N - 1) : min(2 * N - 1, Int(j))
    end

    # Helper to count zero-runs in a vector (with minimum length)
    function count_zero_runs(vec::AbstractVector{<:Integer})
        cnt = 0
        runlen = 0
        for v in vec
            if v == 0
                runlen += 1
            else
                if runlen >= min_gap_length
                    cnt += 1
                end
                runlen = 0
            end
        end
        # trailing run
        if runlen >= min_gap_length
            cnt += 1
        end
        return cnt
    end

    # Compute number of gaps for each tolerance
    gap_counts = Vector{Int}(undef, n_tols)
    for i in 1:n_tols
        col = disc_mp_matrix[:, i]
        gap_counts[i] = count_zero_runs(col)
    end

    # Decide which tolerances to highlight based on mode
    highlight_runs = []  # will collect (i_start, i_end) ranges of tol indices to highlight
    if mode == :full
        # find contiguous ranges where gap_counts == expected_gaps
        i = 1
        while i <= n_tols
            if gap_counts[i] == expected_gaps
                starti = i
                while i <= n_tols && gap_counts[i] == expected_gaps
                    i += 1
                end
                push!(highlight_runs, (starti, i - 1))
            else
                i += 1
            end
        end
    elseif mode == :fixed && !isnothing(fixed_tols)
        # highlight each fixed tol if it yields expected_gaps
        for (i, tol) in enumerate(mp_tol_range)
            if tol in fixed_tols && gap_counts[i] == expected_gaps
                push!(highlight_runs, (i, i))
            end
        end
    elseif mode == :test && !isnothing(test_index)
        ti = test_index
        if ti >= 1 && ti <= n_tols && gap_counts[ti] == expected_gaps
            push!(highlight_runs, (ti, ti))
        end
    else
        # unsupported combination => no highlights
    end

    # Prepare heatmap data: Plots.heatmap expects z dims = (length(y), length(x))
    heatmap_z = permutedims(disc_mp_matrix)  # (n_tols, n_mu)
    y_vals = log10.(mp_tol_range)            # plot tolerances on log scale

    # Base heatmap: use extrema of viridis (dark blue for 1 (positive discretisation, so MP~-1.0 if mp_targ=-1.0), yellow for 0)
    cmap = cgrad(:viridis, 2, categorical=true, rev=true)

    p = heatmap(mu_vals, y_vals, heatmap_z;
        xlabel="μ",
        ylabel="log10(mp_tol)",
        color=cmap,
        xlims=(minimum(mu_vals), maximum(mu_vals)),
        yticks = :auto,
        title = "μ vs mp_tol (disc MP)",
        aspect_ratio = :auto,
        legend=true
    )

    # Overlay highlighted tol windows as translucent horizontal bands
    xmin = minimum(mu_vals)
    xmax = maximum(mu_vals)
    for (i1, i2) in highlight_runs
        tol_low = mp_tol_range[i1]
        tol_high = mp_tol_range[i2]
        y1 = log10(tol_low)
        y2 = log10(tol_high)
        # draw translucent rectangle
        plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
            fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.25, strokewidth=0, label=false)
        # add text label in center of band
        midy = (y1 + y2) / 2
        midx = (xmin + xmax) / 2
        ann = @sprintf("tol idx %d-%d\ntol %.3g-%.3g\ngaps=%d", i1, i2, tol_low, tol_high, expected_gaps)
        annotate!(p, midx, midy, text(ann, :black, 8, :center))
    end

    # If test_index provided, also add side-panel showing mp_disc vs mu for that tol
    if mode == :test && !isnothing(test_index)
        ti = test_index
        if 1 <= ti <= n_tols
            vall = disc_mp_matrix[:, ti]
            # create small plot of mp_disc along mu
            p2 = plot(mu_vals, vall;
                xlabel="μ",
                ylabel="disc_mp",
                ylim = (-0.1, 1.1),
                title = @sprintf("disc_mp @ tol idx %d (tol=%.3g), gaps=%d", ti, mp_tol_range[ti], gap_counts[ti]),
                legend=false,
                color=:black
            )
            layout = @layout [a; b{0.25h}]
            combined = plot(p, p2; layout=layout, size=(1000,800), link=:x)
            if !isnothing(savepath)
                savefig(combined, savepath)
            end
            return combined
        end
    end

    # annotate overall info using energy-gap-derived expectation when available
    info_lines = String[]
    push!(info_lines, @sprintf("|ΔE|>Δ gaps = %d", expected_gaps))
    if energy_gap_counts !== nothing && !isempty(energy_gap_counts)
        local_min = minimum(energy_gap_counts)
        local_max = maximum(energy_gap_counts)
        if local_min == local_max
            push!(info_lines, @sprintf("per-μ gaps = %d", local_min))
        else
            push!(info_lines, @sprintf("per-μ gaps range = %d-%d", local_min, local_max))
        end
    end
    if !isnothing(j)
        push!(info_lines, @sprintf("rational j = %d", j))
    end
    info_str = join(info_lines, "\n")
    annotate!(p, minimum(mu_vals), maximum(y_vals), text(info_str, :left, 8, :black))

    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end


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

    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)

    j = nothing
    expected_gaps = -1
    gap_counts = zeros(Int, n_tols)
    valid_indices = Int[]
    mu_c_vec = nothing
    mu_c_scalar = nothing
    energy_gap_counts = nothing

    if !isnothing(analysis)
        j = get(analysis, :j, j)
        expected_gaps = get(analysis, :expected_gaps, expected_gaps)
        gap_counts = get(analysis, :gap_counts, gap_counts)
        valid_indices = get(analysis, :valid_indices, valid_indices)
        mu_c_scalar = get(analysis, :mu_c_value, nothing)
        mu_c_vec = get(analysis, :mu_c_vec, mu_c_vec) # Extract mu_c vector
        energy_gap_counts = get(analysis, :energy_gap_counts, nothing)
    else
        @warn "Analysis data missing in BSON. Plotting without gap info."
    end

    # 1. Calculate representative scalar mu_c (median of valid indices if value missing)
    if isnothing(mu_c_scalar) && !isnothing(mu_c_vec)
        # Use valid indices if available, otherwise use all
        relevant_mu_c = !isempty(valid_indices) ? mu_c_vec[valid_indices] : mu_c_vec
        if !isempty(relevant_mu_c)
            mu_c_scalar = median(relevant_mu_c)
        end
    end

    if isnothing(mu_c_scalar)
        mu_c_scalar = maximum(mu_vals)
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
        yguidefontrotation=0,
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
    push!(lines, @sprintf("Parameters: N=%d", N))
    if !isnothing(j)
        push!(lines, @sprintf("Rational j: %d", j))
    end
    if expected_gaps >= 0
        push!(lines, @sprintf("|ΔE|>Δ gaps: %d", expected_gaps))
    else
        push!(lines, "|ΔE|>Δ gaps: ?")
    end
    if energy_gap_counts !== nothing && !isempty(energy_gap_counts)
        eg_min = minimum(energy_gap_counts)
        eg_max = maximum(energy_gap_counts)
        if eg_min == eg_max
            push!(lines, @sprintf("per-μ gaps: %d", eg_min))
        else
            push!(lines, @sprintf("per-μ gaps range: %d-%d", eg_min, eg_max))
        end
    end
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

function plt_mu_vs_mptol_check_normalized(
    sweep_input::Union{String,Dict};
    N::Int,
    mode::Symbol = :full,               # :full | :fixed | :test
    fixed_tols::Union{Nothing, Vector{Float64}} = nothing,
    test_index::Union{Nothing, Int} = nothing,
    savepath::Union{Nothing, String} = nothing,
    noise_threshold_override::Union{Nothing, Float64} = nothing  # Optional: override noise threshold
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

    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)

    j = nothing
    expected_gaps = -1
    gap_counts = zeros(Int, n_tols)
    valid_indices = Int[]
    mu_c_vec = nothing
    mu_c_scalar = nothing
    energy_gap_counts = nothing

    if !isnothing(analysis)
        j = get(analysis, :j, j)
        expected_gaps = get(analysis, :expected_gaps, expected_gaps)
        gap_counts = get(analysis, :gap_counts, gap_counts)
        valid_indices = get(analysis, :valid_indices, valid_indices)
        mu_c_scalar = get(analysis, :mu_c_value, nothing)
        mu_c_vec = get(analysis, :mu_c_vec, mu_c_vec) # Extract mu_c vector
        energy_gap_counts = get(analysis, :energy_gap_counts, nothing)
    else
        @warn "Analysis data missing in BSON. Plotting without gap info."
    end

    # 1. Calculate representative scalar mu_c (median of valid indices if value missing)
    if isnothing(mu_c_scalar) && !isnothing(mu_c_vec)
        # Use valid indices if available, otherwise use all
        relevant_mu_c = !isempty(valid_indices) ? mu_c_vec[valid_indices] : mu_c_vec
        if !isempty(relevant_mu_c)
            mu_c_scalar = median(relevant_mu_c)
        end
    end

    if isnothing(mu_c_scalar)
        mu_c_scalar = maximum(mu_vals)
    end

    # Compute noise threshold from boundary curve (as in plt_boundary_curve_with_peak_heights)
    boundary_tol = Vector{Union{Float64, Missing}}(missing, n_mu)
    for mu_idx in 1:n_mu
        for tol_idx in 1:n_tols
            if disc_mp_matrix[mu_idx, tol_idx] == 1
                boundary_tol[mu_idx] = mp_tol_range[tol_idx]
                break
            end
        end
        if ismissing(boundary_tol[mu_idx])
            boundary_tol[mu_idx] = maximum(mp_tol_range)
        end
    end
    log_boundary = log10.(boundary_tol)
    noise_threshold = isnothing(noise_threshold_override) ? minimum(skipmissing(log_boundary)) : noise_threshold_override

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

    # Prepare heatmap data with normalized y-axis
    heatmap_z = permutedims(disc_mp_matrix)  # (n_tols, n_mu)
    y_vals = log10.(mp_tol_range)
    # y_vals = y_vals_raw .- noise_threshold  # Normalize to noise threshold

    # Base heatmap colors
    # rev=true: 0 -> Yellow (Gap), 1 -> Dark (MBS)
    cmap = cgrad(:viridis, 2, categorical=true, rev=true)
    col_mbs = RGB(0.267, 0.004, 0.329) # Viridis start (Dark)
    col_gap = RGB(0.993, 0.906, 0.144) # Viridis end (Yellow)

    p = heatmap(mu_vals, y_vals, heatmap_z;
        xlabel="μ",
        ylabel="log10(mp_tol) - noise_threshold",
        yguidefontrotation=0,
        color=cmap,
        ylims=(-1.25, maximum(y_vals)), #minimum(noise_threshold), maximum(y_vals)),
        xlims=(minimum(mu_vals), maximum(mu_vals)),
        yticks = :auto,
        title = "μ vs Normalized mp_tol (disc MP)",
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
        y1 = log10(tol_low) - noise_threshold
        y2 = log10(tol_high) - noise_threshold
        
        plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
            fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.25, strokewidth=0, label=false)
    end

    # --- Add Annotation Directly to Main Plot ---
    # Position at bottom left (adjust x/y as needed for exact placement)
    lines = String[]
    push!(lines, @sprintf("Parameters: N=%d", N))
    if !isnothing(j)
        push!(lines, @sprintf("Rational j: %d", j))
    end
    if expected_gaps >= 0
        push!(lines, @sprintf("|ΔE|>Δ gaps: %d", expected_gaps))
    else
        push!(lines, "|ΔE|>Δ gaps: ?")
    end
    if energy_gap_counts !== nothing && !isempty(energy_gap_counts)
        eg_min = minimum(energy_gap_counts)
        eg_max = maximum(energy_gap_counts)
        if eg_min == eg_max
            push!(lines, @sprintf("per-μ gaps: %d", eg_min))
        else
            push!(lines, @sprintf("per-μ gaps range: %d-%d", eg_min, eg_max))
        end
    end
    if !isnothing(mu_c_scalar)
        push!(lines, @sprintf("μ_c ≈ %.4f", mu_c_scalar))
    end
    push!(lines, @sprintf("Noise Threshold: %.3f", noise_threshold))
    
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

function plt_boundary_curve_with_peak_heights(
    sweep_input::Union{String,Dict};
    N::Int,
    savepath::Union{Nothing, String} = nothing,
    noise_threshold_override::Union{Nothing, Float64} = nothing  # Optional: override noise threshold
)
    # Load sweep data if a filename provided
    sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
    @assert !isnothing(sweep) "sweep data/file must be provided"

    mu_vals = sweep[:mu]
    mp_tol_range = sweep[:mp_tol_range]
    disc_mp_matrix = sweep[:disc_mp_matrix]  # (n_mu, n_tols)
    params = get(sweep, :params, nothing)
    mp_targ = params !== nothing ? params[:mp_targ] : get(sweep, :mp_targ, nothing)
    
    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)

    # Extract boundary: for each mu, smallest tol where disc_mp == 1
    boundary_tol = Vector{Union{Float64, Missing}}(missing, n_mu)
    for mu_idx in 1:n_mu
        for tol_idx in 1:n_tols
            if disc_mp_matrix[mu_idx, tol_idx] == 1
                boundary_tol[mu_idx] = mp_tol_range[tol_idx]
                break
            end
        end
        # If never 1, set to max tol (or leave missing)
        if ismissing(boundary_tol[mu_idx])
            boundary_tol[mu_idx] = maximum(mp_tol_range)
        end
    end

    # Log-transform boundary
    log_boundary = log10.(boundary_tol)

    # Noise threshold: global min, or override
    noise_threshold = isnothing(noise_threshold_override) ? minimum(skipmissing(log_boundary)) : noise_threshold_override

    # Heights: y - noise_threshold
    heights = log_boundary .- noise_threshold

    # Identify peaks: local maxima (skip edges)
    peak_indices = Int[]
    for i in 2:(n_mu-1)
        if !ismissing(log_boundary[i]) && log_boundary[i] > log_boundary[i-1] && log_boundary[i] > log_boundary[i+1]
            push!(peak_indices, i)
        end
    end

    # Plot the boundary curve with gradient coloring
    p = plot(mu_vals, log_boundary;
        xlabel = "μ",
        ylabel = "log10(mp_tol)",
        title = "Boundary Curve with Peak Heights",
        line_z = heights,  # Color gradient based on heights
        color = :viridis,  # Gradient colormap
        linewidth = 2,
        colorbar = true,
        colorbar_title = "Height",
        legend = false
    )

    # Highlight peaks with markers colored by height
    if !isempty(peak_indices)
        peak_mu = mu_vals[peak_indices]
        peak_y = log_boundary[peak_indices]
        peak_heights = heights[peak_indices]
        scatter!(p, peak_mu, peak_y;
            marker_z = peak_heights,
            color = :viridis,
            markersize = 6,
            markerstrokewidth = 1,
            label = "Peaks"
        )
    end

    # Optional: Add mu_c if available (from analysis)
    analysis = get(sweep, :analysis, nothing)
    mu_c_scalar = nothing
    if !isnothing(analysis)
        mu_c_scalar = get(analysis, :mu_c_value, nothing)
        if isnothing(mu_c_scalar) && haskey(analysis, :mu_c_vec)
            mu_c_vec = analysis[:mu_c_vec]
            if !isempty(mu_c_vec)
                mu_c_scalar = median(mu_c_vec)
            end
        end
    end
    if !isnothing(mu_c_scalar)
        vline!(p, [mu_c_scalar], color=:red, linestyle=:dash, label="μ_c")
    end

    # Annotations
    info_lines = String[]
    push!(info_lines, @sprintf("N = %d", N))
    push!(info_lines, @sprintf("Noise Threshold: %.3f", noise_threshold))
    push!(info_lines, @sprintf("Num Peaks: %d", length(peak_indices)))
    if !isnothing(mp_targ)
        push!(info_lines, @sprintf("MP Target: %.1f", mp_targ))
    end
    info_text = join(info_lines, "\n")
    annotate!(p, minimum(mu_vals), minimum(log_boundary) - 0.1 * (maximum(log_boundary) - minimum(log_boundary)), 
              text(info_text, :left, 8, :black))

    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end

function plt_mu_vs_mptol_check_with_precomputed_ranges(
    sweep_input::Union{String,Dict},
    tol_range::Tuple{Float64,Float64},
    mu_count_range::Tuple{Float64,Float64},
    savepath::String
)
    # Load sweep data if a filename provided
    sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
    @assert !isnothing(sweep) "sweep data/file must be provided"

    mu_vals = sweep[:mu]
    mp_tol_range = sweep[:mp_tol_range]
    disc_mp_matrix = sweep[:disc_mp_matrix]         # (n_mu, n_tols)

    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)

    highlight_runs = []
    tol_low, tol_high = tol_range
    for (i, tol) in enumerate(mp_tol_range)
        if tol >= tol_low && tol <= tol_high
            push!(highlight_runs, (i, i))
        end
    end

    # Assume mu_vals and mu_count_range are defined
    if mu_count_range[1] == 0.0
        xlims = (0.0, maximum(mu_vals))
    else
        xlims = (minimum(mu_vals), maximum(mu_vals))
    end

    # Prepare heatmap data
    heatmap_z = permutedims(disc_mp_matrix)  # (n_tols, n_mu)
    y_vals = log10.(mp_tol_range)            # plot tolerances on log scale

    # Base heatmap colors
    cmap = cgrad(:viridis, 2, categorical=true, rev=true)

    p = heatmap(mu_vals, y_vals, heatmap_z;
        xlabel="μ",
        ylabel="log10(mp_tol)",
        yguidefontrotation=0,
        color=cmap,
        ylims=(minimum(y_vals), maximum(y_vals)),
        xlims=xlims,
        yticks = :auto,
        title = "μ vs mp_tol (disc MP)",
        aspect_ratio = :auto,
        colorbar=false
    )

    vline!([mu_count_range[2]]; color=:red, linestyle=:dash, linewidth=2, label="µ_c")

    # Overlay highlighted tol windows (visual boxes only)
    xmin = minimum(mu_vals)
    xmax = maximum(mu_vals)
    for (i1, i2) in highlight_runs
        tol_low = mp_tol_range[i1]
        tol_high = mp_tol_range[i2]
        y1 = log10(tol_low)
        y2 = log10(tol_high)
        
        plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
            fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.5, strokewidth=0, label=false)
    end

    if !isnothing(savepath)
        savefig(p, savepath)
    end

    return p
end


# ## mu vs mp_tol plot is supplied the expected number of gaps to identify valid tol range for
# function plt_mu_vs_mptol_check_ext_gaps(
#     sweep_input::Union{String,Dict},
#     external_expected_gaps::Int;
#     N::Int,
#     mode::Symbol = :full,               # :full | :fixed | :test
#     fixed_tols::Union{Nothing, Vector{Float64}} = nothing,
#     test_index::Union{Nothing, Int} = nothing,
#     savepath::Union{Nothing, String} = nothing
# )
#     # Load sweep data if a filename provided
#     sweep = isa(sweep_input, String) ? BSON.load(sweep_input) : sweep_input
#     @assert !isnothing(sweep) "sweep data/file must be provided"

#     mu_vals = sweep[:mu]
#     mp_tol_range = sweep[:mp_tol_range]
#     disc_mp_matrix = sweep[:disc_mp_matrix]         # (n_mu, n_tols)
#     params = get(sweep, :params, nothing)
#     mp_targ = params !== nothing ? params[:mp_targ] : get(sweep, :mp_targ, nothing)
    
#     # Extract Analysis Data
#     analysis = get(sweep, :analysis, nothing)

#     n_mu = length(mu_vals)
#     n_tols = length(mp_tol_range)

#     gap_counts = zeros(Int, n_tols)
#     valid_indices = Int[]
#     mu_c_vec = nothing
#     mu_c_scalar = nothing
#     energy_gap_counts = nothing
#     j_from_analysis = nothing

#     if !isnothing(analysis)
#         j_from_analysis = get(analysis, :j, nothing)
#         gap_counts = get(analysis, :gap_counts, gap_counts)
#         # We re-calculate valid_indices based on external_expected_gaps
#         # instead of trusting the one from analysis which might be based on a different expectation
#         mu_c_scalar = get(analysis, :mu_c_value, nothing)
#         mu_c_vec = get(analysis, :mu_c_vec, mu_c_vec) 
#         energy_gap_counts = get(analysis, :energy_gap_counts, nothing)
#     end

#     # Re-calculate valid indices matching external requirement
#     # We assume 'gap_counts' in analysis is correct (count of gaps per tolerance row)
#     # If not present, we would need to recompute it from disc_mp_matrix
#     if all(gap_counts .== 0)
#         # Helper to count zero-runs in a vector (re-implementation as it's not global)
#         function count_zero_runs_local(vec::AbstractVector{<:Integer})
#             cnt = 0
#             runlen = 0
#             # Assuming min_gap_length=1 as default if not specified
#             min_gap_length = 1 
#             for v in vec
#                 if v == 0
#                     runlen += 1
#                 else
#                     if runlen >= min_gap_length
#                         cnt += 1
#                     end
#                     runlen = 0
#                 end
#             end
#             if runlen >= min_gap_length
#                 cnt += 1
#             end
#             return cnt
#         end
        
#         for i in 1:n_tols
#             col = disc_mp_matrix[:, i]
#             gap_counts[i] = count_zero_runs_local(col)
#         end
#     end

#     valid_indices = findall(x -> x == external_expected_gaps, gap_counts)

#     # 1. Calculate representative scalar mu_c (median of valid indices if value missing)
#     if isnothing(mu_c_scalar) && !isnothing(mu_c_vec)
#         relevant_mu_c = !isempty(valid_indices) ? mu_c_vec[valid_indices] : mu_c_vec
#         if !isempty(relevant_mu_c)
#             mu_c_scalar = median(relevant_mu_c)
#         end
#     end

#     if isnothing(mu_c_scalar)
#         mu_c_scalar = maximum(mu_vals)
#     end

#     # Decide which tolerances to highlight based on mode
#     highlight_runs = []
    
#     if mode == :full
#         if !isempty(valid_indices)
#             sort!(valid_indices)
#             start_i = valid_indices[1]
#             prev_i = valid_indices[1]
#             for i in valid_indices[2:end]
#                 if i == prev_i + 1
#                     prev_i = i
#                 else
#                     push!(highlight_runs, (start_i, prev_i))
#                     start_i = i
#                     prev_i = i
#                 end
#             end
#             push!(highlight_runs, (start_i, prev_i))
#         end
#     elseif mode == :fixed && !isnothing(fixed_tols)
#         # Not fully supported for "largest range" logic, but keeping for compatibility
#         for (i, tol) in enumerate(mp_tol_range)
#             if tol in fixed_tols && (i in valid_indices)
#                 push!(highlight_runs, (i, i))
#             end
#         end
#     elseif mode == :test && !isnothing(test_index)
#          # Not fully supported for "largest range" logic, but keeping for compatibility
#         ti = test_index
#         if ti >= 1 && ti <= n_tols && (ti in valid_indices)
#             push!(highlight_runs, (ti, ti))
#         end
#     end

#     # Select the range with the largest mp_tol values if multiple ranges exist
#     # Assuming mp_tol_range is sorted? Or at least we look at indices.
#     # Typically sweep logic sorts tolerances ascending or descending.
#     # Let's assume we want the range that contains the largest tolerance values.
    
#     selected_range_tol = Float64[]
    
#     if !isempty(highlight_runs)
#         if length(highlight_runs) > 1
#              # We want the run whose MAXIMUM tolerance index maps to the LARGEST tolerance value
#              # Let's check max value in each run
#              run_max_vals = [maximum(mp_tol_range[r[1]:r[2]]) for r in highlight_runs]
#              best_run_idx = argmax(run_max_vals)
#              highlight_runs = [highlight_runs[best_run_idx]]
#         end
        
#         # Now we calculate the [min, max] tolerance for the selected single range
#         r = highlight_runs[1]
#         tols_in_run = mp_tol_range[r[1]:r[2]]
#         selected_range_tol = [minimum(tols_in_run), maximum(tols_in_run)]
#     end


#     # Prepare heatmap data
#     heatmap_z = permutedims(disc_mp_matrix)
#     y_vals = log10.(mp_tol_range)

#     # Base heatmap colors
#     cmap = cgrad(:viridis, 2, categorical=true, rev=true)
#     col_mbs = RGB(0.267, 0.004, 0.329) # Viridis start (Dark)
#     col_gap = RGB(0.993, 0.906, 0.144) # Viridis end (Yellow)

#     p = heatmap(mu_vals, y_vals, heatmap_z;
#         xlabel="μ",
#         ylabel="log10(mp_tol)",
#         yguidefontrotation=0,
#         color=cmap,
#         ylims=(minimum(y_vals), maximum(y_vals)),
#         xlims=(minimum(mu_vals), maximum(mu_vals)),
#         yticks = :auto,
#         title = "μ vs mp_tol (disc MP) [External Gap Target]",
#         aspect_ratio = :auto,
#         colorbar=false,
#         legend=:outerbottom
#     )

#     # --- Legend Entries ---
#     if !isnothing(mu_c_scalar)
#         vline!(p, [mu_c_scalar], color=:red, linestyle=:dot, linewidth=2, label="μ_c")
#     end
#     mbs_label = isnothing(mp_targ) ? "MBS" : @sprintf("MBS (MP ≈ %.1f)", mp_targ)
#     scatter!(p, [NaN], [NaN], color=col_mbs, label=mbs_label, shape=:square, markersize=5, legend_font_pointsize=8)
#     scatter!(p, [NaN], [NaN], color=col_gap, label="Gap", shape=:square, markersize=5)

#     # Overlay highlighted tol windows
#     xmin = minimum(mu_vals)
#     xmax = maximum(mu_vals)
#     for (i1, i2) in highlight_runs
#         tol_low = mp_tol_range[i1]
#         tol_high = mp_tol_range[i2]
#         y1 = log10(tol_low)
#         y2 = log10(tol_high)
        
#         plot!(p, Shape([xmin, xmax, xmax, xmin], [y1, y1, y2, y2]);
#             fillcolor = RGB(1.0, 0.85, 0.3), alpha = 0.25, strokewidth=0, label=false)
#     end

#     # --- Add Annotation ---
#     lines = String[]
#     push!(lines, @sprintf("Parameters: N=%d", N))
#     push!(lines, @sprintf("EXTERNAL TARGET GAPS: %d", external_expected_gaps))
    
#     if !isnothing(j_from_analysis)
#         push!(lines, @sprintf("Rational j: %d", j_from_analysis))
#     end

#     if energy_gap_counts !== nothing && !isempty(energy_gap_counts)
#         eg_min = minimum(energy_gap_counts)
#         eg_max = maximum(energy_gap_counts)
#         if eg_min == eg_max
#             push!(lines, @sprintf("per-μ energy gaps: %d", eg_min))
#         else
#             push!(lines, @sprintf("per-μ energy gaps range: %d-%d", eg_min, eg_max))
#         end
#     end
#     if !isnothing(mu_c_scalar)
#         push!(lines, @sprintf("μ_c ≈ %.4f", mu_c_scalar))
#     end
    
#     if !isempty(highlight_runs)
#         push!(lines, "Selected Max Tolerance Range:")
#         for (i1, i2) in highlight_runs
#             t1 = mp_tol_range[i1]
#             t2 = mp_tol_range[i2]
#             push!(lines, @sprintf("  Idx %d-%d: %.1e to %.1e", i1, i2, t1, t2))
#         end
#     else
#         push!(lines, "No tolerance ranges found matching target.")
#     end
    
#     info_text = join(lines, "\n")
#     offset_y = 0.2 * (maximum(y_vals) - minimum(y_vals))
#     annotate!(p, minimum(mu_vals), minimum(y_vals) - offset_y, text(info_text, :left, 8, :black))

#     if !isnothing(savepath)
#         savefig(p, savepath)
#     end

#     return selected_range_tol
# end




#############################
### QC vs Bulk Gap criterion plots
#############################
function plt_bulk_vs_qc_gaps(
    mu_vals::AbstractVector{Float64},
    bulk_gaps::AbstractVector{Float64},
    E_vals::AbstractVector{Float64},
    qc_gaps::AbstractVector{Float64},
    savepath::AbstractString;
    plot_title::String="",
    xlims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    logyscale::Bool=false,
    size::Tuple{Int,Int}=(1000,200)
)
    if logyscale
        bulk_gaps = replace(bulk_gaps, 0.0=>1e-10)
        qc_gaps = replace(qc_gaps, 0.0=>1e-10)
        p= plot(
            [mu_vals, E_vals], 
            [bulk_gaps, qc_gaps],
            xlabel=L"\mu \quad / \quad E",
            ylabel="Energy Gap Size",
            yscale=:log10,
            # title=plot_title,
            label=["Bulk Gaps (vs μ)" "QC Gaps (vs E)"],
            linewidth=1.5,
            alpha=0.8,
            legend=:bottomright,
            minorgrid=true,
            xlims=xlims,
            size=size,
            framestyle=:box
        )
    else
        p = plot(
            [mu_vals, E_vals], 
            [bulk_gaps, qc_gaps],
            xlabel=L"\mu \quad / \quad E",
            ylabel="Energy Gap Size",
            # title=plot_title,
            label=["Bulk Gaps (vs μ)" "QC Gaps (vs E)"],
            linewidth=1.5,
            alpha=0.8,
            legend=:topright,
            minorgrid=true,
            # xlims=xlims,
            ylims=(-0.01, 0.6),
            xlims=(-0.1, 3.0), #maximum(mu_vals)*1.01),
            size=size,
            framestyle=:box
        )
    end

    if !isempty(savepath)
        mkpath(dirname(savepath))
        savefig(p, savepath)
    end
    
    return p
end


function plt_all_gap_components(
    mu_vals::AbstractVector{Float64},
    bulk_gaps::AbstractVector{Float64},
    E_vals::AbstractVector{Float64},
    qc_gaps::AbstractVector{Float64},
    savepath::AbstractString;
    plot_title::String="",
    xlims::Union{Nothing, Tuple{Float64, Float64}}=nothing,
    logyscale::Bool=false,
    mid_gaps::Union{Nothing, AbstractVector{Float64}}=nothing,
    outer_gaps::Union{Nothing, AbstractVector{Float64}}=nothing
)
    # Prepare data for plotting
    x_data = [mu_vals, E_vals]
    y_data = [bulk_gaps, qc_gaps]
    labels_vec = String["Bulk Gaps (vs μ)", "QC Gaps (vs E)"]
    alphas_vec = Float64[1.0, 1.0]
    
    if mid_gaps !== nothing
        push!(x_data, mu_vals)
        push!(y_data, mid_gaps)
        push!(labels_vec, "Mid Gaps (vs μ)")
        push!(alphas_vec, 0.3)
    end
    
    if outer_gaps !== nothing
        push!(x_data, mu_vals)
        push!(y_data, outer_gaps)
        push!(labels_vec, "Outer Gaps (vs μ)")
        push!(alphas_vec, 0.3)
    end

    labels_row = reshape(labels_vec, 1, :)
    alphas_row = reshape(alphas_vec, 1, :)
    
    if logyscale
        y_data = [replace(y, 0.0=>1e-10) for y in y_data]
        p = plot(
            x_data, 
            y_data,
            xlabel=L"\mu \quad / \quad E",
            ylabel="Energy Gap Size",
            yscale=:log10,
            title=plot_title,
            label=labels_row,
            linewidth=1.5,
            alpha=alphas_row,
            legend=:bottomright,
            minorgrid=true,
            xlims=xlims
        )
    else
        p = plot(
            x_data, 
            y_data,
            xlabel=L"\mu \quad / \quad E",
            ylabel="Energy Gap Size",
            title=plot_title,
            label=labels_row,
            linewidth=1.5,
            alpha=alphas_row,
            legend=:topright,
            minorgrid=true,
            xlims=xlims
        )
    end

    if !isempty(savepath)
        mkpath(dirname(savepath))
        savefig(p, savepath)
    end
    
    return p
end

function plt_bulk_vs_qc_diff(
    mu_vals::AbstractVector{Float64},
    bulk_gaps::AbstractVector{Float64},
    E_qc::AbstractVector{Float64},
    qc_gaps::AbstractVector{Float64},
    savepath::AbstractString;
    plot_title::String="",
    xlims::Union{Nothing, Tuple{Float64, Float64}}=nothing
)

    @assert length(mu_vals) == length(bulk_gaps) "mu_vals and bulk_gaps must have same length"
    @assert length(E_qc) == length(qc_gaps) "E_qc and qc_gaps must have same length"

    # Interpolate QC gaps onto the mu grid
    # QC gaps are defined on E_qc. We want to compare them at values of mu.
    # We only compute the difference where mu is within the range of E_qc.
    
    qc_interp = fill(NaN, length(mu_vals))
    E_min, E_max = minimum(E_qc), maximum(E_qc)
    
    for (i, mu) in enumerate(mu_vals)
        if mu >= E_min && mu <= E_max
            # Find surrounding indices in E_qc
            # E_qc is sorted
            k = searchsortedlast(E_qc, mu)
            
            val = 0.0
            if k == 0
                val = qc_gaps[1] # Should not happen given bounds check
            elseif k >= length(E_qc)
                val = qc_gaps[end]
            else
                x0 = E_qc[k]
                x1 = E_qc[k+1]
                y0 = qc_gaps[k]
                y1 = qc_gaps[k+1]
                
                # Linear interpolation
                if x1 == x0
                    val = y0
                else
                    val = y0 + (mu - x0) * (y1 - y0) / (x1 - x0)
                end
            end
            qc_interp[i] = val
        end
    end

    # Difference = Bulk - QC
    # This will contain NaNs where mu is outside E_range
    diff_vals = bulk_gaps .- qc_interp

    # 1. Identify Negative Regions for Counting & Marking
    regions = Float64[]
    region_count = 0
    in_region = false
    start_mu = 0.0

    # Define counting bounds based on plotting xlims
    count_min, count_max = isnothing(xlims) ? (-Inf, Inf) : xlims

    # Ensure mu_vals are sorted for region identification logic
    for i in eachindex(diff_vals)
        cur_mu = mu_vals[i]
        
        # Skip NaNs (outside valid interpolation range)
        if isnan(diff_vals[i])
            if in_region
                in_region = false
                push!(regions, start_mu)
                push!(regions, mu_vals[i-1])
            end
            continue
        end

        # Skip values outside the specified plotting range
        if cur_mu < count_min || cur_mu > count_max
            # If we were strictly inside a region and stepped out, close it at the boundary/last point
            if in_region
                in_region = false
                push!(regions, start_mu)
                push!(regions, mu_vals[i-1])
            end
            continue
        end

        if diff_vals[i] < 0
            if !in_region
                in_region = true
                start_mu = cur_mu
                region_count += 1
            end
        else
            if in_region
                in_region = false
                push!(regions, start_mu)
                push!(regions, mu_vals[i-1])
            end
        end
    end
    # Close any open region at the end
    if in_region
        push!(regions, start_mu)
        # Find the last valid index within range to close properly
        valid_idxs = findall(i -> !isnan(diff_vals[i]) && count_min <= mu_vals[i] <= count_max, eachindex(mu_vals))
        if !isempty(valid_idxs)
             push!(regions, mu_vals[valid_idxs[end]])
        else
             push!(regions, start_mu)
        end
    end
    
    # Construct label with info
    x_lbl = L"\mu" * "\n" * "Number of QC>Bulk gap regions: $(region_count)"

    p = plot(
        mu_vals, 
        diff_vals,
        xlabel=x_lbl,
        ylabel=L"\Delta_{Bulk} - \Delta_{QC}",
        title=plot_title,
        label="Difference",
        linewidth=1.5,
        color=:purple,
        legend=:topright,
        minorgrid=true,
        xlims=xlims,
        bottom_margin=10mm
    )
    
    # Shade negative regions yellow (50% transparent)
    # We filter diff_vals for shading to handle NaNs (min(NaN, 0) is NaN)
    shade_vals = [isnan(x) ? NaN : min(x, 0.0) for x in diff_vals]
    
    plot!(p, mu_vals, shade_vals, 
        fillrange=0.0, 
        fillcolor=:yellow, 
        fillalpha=0.5, 
        linewidth=0, 
        label=""
    )

    # Mark ranges on x-axis (y=min_val)
    # Use global minimum for the mark y-position, ignoring NaNs
    valid_diffs = filter(!isnan, diff_vals)
    y_mark = isempty(valid_diffs) ? 0.0 : minimum(valid_diffs)
    
    for k in 1:2:length(regions)
        r_start = regions[k]
        r_end = regions[k+1]
        # Draw a segment on the axis
        plot!(p, [r_start, r_end], [y_mark, y_mark], linewidth=3.0, color=:orange, label="")
    end

    hline!(p, [0.0], color=:black, linestyle=:dash, label="")

    # Save plot
    if !isempty(savepath)
        mkpath(dirname(savepath))
        savefig(p, savepath)
    end
    
    return p
end



end  # module MPGapPlotting











