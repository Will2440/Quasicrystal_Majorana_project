using BSON
using DataFrames
using Colors
using Plots

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "hof_style_slopes_N500_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-3p0-5001)_Delta(0p05-0p05-1)"
run_set_name = "specific_slopes"
unpack_folder = "mptarg-1.0_mptol0.5_IDOPminsamp2_N_rangenothing_Delta_rangenothing_tn_rangenothing_phason_rangenothing_mpduplicated-true"

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

gap_label_file = joinpath(base_results_dir, "idop", unpack_folder, "all_phason_idop_gap_label_records.bson")


######################################
############## Load data #############
######################################

function load_gap_label_records(path::AbstractString)
    isfile(path) || error("Gap label file not found: $(path)")
    payload = BSON.load(path)

    if haskey(payload, :label_dataset_df)
        df = payload[:label_dataset_df]
    elseif haskey(payload, :label_records)
        df = DataFrame(payload[:label_records])
    else
        error("Unexpected structure in gap label BSON; available keys: $(collect(keys(payload)))")
    end

    return df isa DataFrame ? df : DataFrame(df)
end

data = load_gap_label_records(gap_label_file)



#######################################
########## Plot Butterfly #############
#######################################

## Colour Colorscheme
function q_anchor_colormap(
    qs::Vector{Int}
)

    lerp(c1::RGB, c2::RGB, t) = RGB(
        (1-t)*c1.r + t*c2.r,
        (1-t)*c1.g + t*c2.g,
        (1-t)*c1.b + t*c2.b
    )
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

    
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end

    if !isnothing(idx_neg1) && idx_neg1 > idx_neg_high
        cmap[idx_neg_high] = col_neg_high
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

function plt_gap_labels_vs_phi(
    df::DataFrame;
    label_col::Symbol = :idop_labels,
    savepath::Union{Nothing, AbstractString} = nothing,
    linewidth::Real = 2,
    title::Union{Nothing, AbstractString} = nothing,
    q_max::Union{Int, Nothing} = nothing,
    fixed_values...
)
    if !(label_col in propertynames(df))
        error("DataFrame does not contain column $(label_col). Available columns: $(propertynames(df))")
    end

    # Filter logic considering column existence and float tolerance
    local filtered_df = df
    if !isempty(fixed_values)
        try
             filtered_df = filter(row -> begin
                all(begin
                    if !(k in propertynames(row))
                        # If a requested filter column isn't in the df, we can't filter by it.
                        # Decide: Fail (as it seems the user expects it) or Skip?
                        # The error suggests 'getproperty' fails if column missing.
                        # Let's check existence first.
                        false # Row fails if it doesn't have the column
                    else
                        val = row[k]
                        if val isa Number && v isa Number
                            abs(val - v) < 1e-6 # Float tolerance
                        else
                            val == v
                        end
                    end
                end for (k, v) in fixed_values)
            end, df)
        catch e
            @warn "Filtering failed" error=e keys(fixed_values) propertynames(df)
            return nothing
        end
    end
    
    if isempty(filtered_df)
        @warn "No data matches fixed_values for plotting: $fixed_values. Skipping plot."
        return nothing
    end

    gap_records = NamedTuple{(:mu_low, :mu_high, :phi, :q), Tuple{Float64, Float64, Float64, Int}}[]
    qs_encountered = Int[]

    for row in eachrow(filtered_df)
        phi_val = :phi in propertynames(row) ? row[:phi] : missing
        if ismissing(phi_val)
            continue
        end

        labels = row[label_col]
        if labels === nothing || isempty(labels)
            continue
        end

        for label in labels
            mu_low, mu_high, q = label
            push!(gap_records, (mu_low=Float64(mu_low), mu_high=Float64(mu_high), phi=Float64(phi_val), q=Int(q)))
            push!(qs_encountered, Int(q))
        end
    end

    if isempty(gap_records)
        @warn "No gap labels found in column $(label_col); skipping plot."
        return nothing
    end

    unique_qs = sort(unique(qs_encountered))
    
    # Determine the set of qs used for the colour scale
    if q_max !== nothing
        # Create a full dense range from -q_max to q_max
        scale_qs = collect(-q_max:q_max)
        
        # If the data contains the sentinel 9999, make sure it's included in our scale/map
        if 9999 in unique_qs && !(9999 in scale_qs)
            push!(scale_qs, 9999)
        end
        
        # Use strictly this scale, ensuring it is sorted for the color function
        cbar_qs = sort(unique(scale_qs))
    else
        cbar_qs = unique_qs
    end

    color_array = q_anchor_colormap(cbar_qs)
    color_dict = Dict(cbar_qs[i] => color_array[i] for i in eachindex(cbar_qs))

    plt_title = isnothing(title) ? "Gap labels ($(label_col))" : String(title)
    main_plt = plot(; xlabel="μ", ylabel="φ", legend=false, grid=true, title=plt_title)

    for rec in gap_records
        q_target = rec.q
        
        # If clamping is active, clamp values exceeding q_max 
        # (but preserve sentinel 9999)
        if q_max !== nothing && rec.q != 9999
            if q_target > q_max
                q_target = q_max
            elseif q_target < -q_max
                q_target = -q_max
            end
        end

        c = get(color_dict, q_target, RGB(0.85, 0.85, 0.85))
        plot!(main_plt, [rec.mu_low, rec.mu_high], [rec.phi, rec.phi]; color=c, linewidth=linewidth, label=false)
    end

    # Create colorbar subplot
    n_q = length(cbar_qs)
    colorbar_data = reshape(1:n_q, n_q, 1)  # Use indices to ensure uniform spacing

    # Handle ticks: if too many, sparse them
    if n_q > 25
        step_idx = ceil(Int, n_q / 20)
        tick_indices = 1:step_idx:n_q
        
        # Ensure the last one is included if it's the sentinel or just the max
        if !(n_q in tick_indices)
            # If the last element is sentinel, we definitely want it. 
            # Or just append the last index.
            tick_indices = vcat(collect(tick_indices), n_q)
            unique!(tick_indices)
        end
        
        cbar_yticks = (tick_indices, string.(cbar_qs[tick_indices]))
    else
        cbar_yticks = (1:n_q, string.(cbar_qs))
    end

    colorbar_plt = heatmap(
        colorbar_data,
        c=cgrad(color_array, categorical=true),
        xticks=false,
        yticks=cbar_yticks,
        ylabel="q",
        title="",
        colorbar=false,
        clims=(1, n_q),
        framestyle=:box,
        grid=false
    )

    # Combine plots
    # Use a custom layout to reduce the width of the colorbar subplot
    l = @layout [a b{0.1w}]
    combined_plt = plot(main_plt, colorbar_plt, layout=l, size=(1000, 600), margin=5Plots.mm)

    if savepath !== nothing
        mkpath(dirname(savepath))
        savefig(combined_plt, savepath)
    end

    return combined_plt
end



q_max = 20

plot_outdir = joinpath(base_results_dir, "idop", unpack_folder, "label_plots_NEW")
isdir(plot_outdir) || mkpath(plot_outdir)

phason_range = collect(0.0:0.01:1.0)
for phason in phason_range
    try
        idop_plot_path = joinpath(plot_outdir, "gap_labels_idop_vs_phi_phason$(round(phason,digits=2))_qmax$(q_max).png")
        plt_gap_labels_vs_phi(
            data;
            label_col=:idop_labels,
            savepath=idop_plot_path,
            title="IDOP Gap Labels vs φ (phason=$(round(phason,digits=2)))",
            q_max=q_max,
            phason=phason
        )
        println("Saved IDOP gap label plot to: ", idop_plot_path)
    catch e
        @warn "Failed to plot IDOP gap labels for phason=$(phason)" error=e
        showerror(stderr, e, catch_backtrace())
    end
end



# try
#     winding_plot_path = joinpath(plot_outdir, "gap_labels_winding_vs_phi_qmax$(q_max).png")
#     plt_gap_labels_vs_phi(
#         data;
#         label_col=:winding_labels,
#         savepath=winding_plot_path,
#         title="Winding Gap Labels vs φ",
#         q_max=q_max
#     )
#     println("Saved winding gap label plot to: ", winding_plot_path)
# catch e
#     @warn "Failed to plot winding gap labels" error=e
#     showerror(stderr, e, catch_backtrace())
# end
