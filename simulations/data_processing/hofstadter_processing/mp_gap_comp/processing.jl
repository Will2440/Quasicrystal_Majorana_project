module MPGapProcessing

using DataFrames
using StatsBase
using DelimitedFiles

include("../idos/processing.jl")
using .IDOSProcessing

################################
### Range changing functions ###
################################
const MU_RANGE_TOL = 1e-9

function mask_negative_mu_df(df::DataFrame, mu_col::Symbol; tol::Float64=MU_RANGE_TOL)
    mask = df[!, mu_col] .>= -tol
    return copy(df[mask, :])
end

function mirror_mu_df(df::DataFrame, mu_col::Symbol; tol::Float64=MU_RANGE_TOL)
    df_sorted = sort(df, mu_col)
    mu_vals = df_sorted[!, mu_col]
    pos_mask = mu_vals .> tol
    if !any(pos_mask)
        return df_sorted
    end
    df_pos = df_sorted[pos_mask, :]
    df_pos_rev = df_pos[end:-1:1, :]
    df_neg = deepcopy(df_pos_rev)
    df_neg[!, mu_col] .= -df_neg[!, mu_col]
    return vcat(df_neg, df_sorted; cols=:union)
end

function mask_negative_mu_sweep!(sweep_dict::Dict; tol::Float64=MU_RANGE_TOL)
    mu_vals = sweep_dict[:mu]
    mask = mu_vals .>= -tol
    sweep_dict[:mu] = mu_vals[mask]
    if haskey(sweep_dict, :mp_raw)
        sweep_dict[:mp_raw] = sweep_dict[:mp_raw][mask]
    end
    sweep_dict[:disc_mp_matrix] = sweep_dict[:disc_mp_matrix][mask, :]
    if haskey(sweep_dict, :analysis)
        analysis = sweep_dict[:analysis]
        if analysis !== nothing && haskey(analysis, :energy_gap_counts)
            analysis[:energy_gap_counts] = analysis[:energy_gap_counts][mask]
        end
        if analysis !== nothing && haskey(analysis, :energy_gap_ranges)
            analysis[:energy_gap_ranges] = analysis[:energy_gap_ranges][mask]
        end
    end
    return sweep_dict
end

function mirror_mu_sweep!(sweep_dict::Dict; tol::Float64=MU_RANGE_TOL)
    mu_vals = sweep_dict[:mu]
    pos_mask = mu_vals .> tol
    if !any(pos_mask)
        return sweep_dict
    end
    mu_pos = mu_vals[pos_mask]
    mu_neg = -reverse(mu_pos)
    sweep_dict[:mu] = vcat(mu_neg, mu_vals)
    if haskey(sweep_dict, :mp_raw)
        mp_raw = sweep_dict[:mp_raw]
        sweep_dict[:mp_raw] = vcat(reverse(mp_raw[pos_mask]), mp_raw)
    end
    disc = sweep_dict[:disc_mp_matrix]
    disc_neg = disc[pos_mask, :]
    disc_neg = disc_neg[end:-1:1, :]
    sweep_dict[:disc_mp_matrix] = vcat(disc_neg, disc)
    if haskey(sweep_dict, :analysis)
        analysis = sweep_dict[:analysis]
        if analysis !== nothing && haskey(analysis, :energy_gap_counts)
            counts = analysis[:energy_gap_counts]
            analysis[:energy_gap_counts] = vcat(reverse(counts[pos_mask]), counts)
        end
        if analysis !== nothing && haskey(analysis, :energy_gap_ranges)
            ranges = analysis[:energy_gap_ranges]
            ranges_pos = ranges[pos_mask]
            if isempty(ranges_pos)
                mirrored_ranges = Vector{Vector{Tuple{Float64, Float64}}}()
            else
                mirrored_ranges = [copy(r) for r in ranges_pos[end:-1:1]]
            end
            analysis[:energy_gap_ranges] = vcat(mirrored_ranges, ranges)
        end
    end
    return sweep_dict
end

function adjust_mu_range(df_slice::DataFrame, sweep_dict::Dict, data_full::Bool, plot_full::Bool; mu_col::Symbol=:mu)
    df_adj = df_slice
    if data_full && !plot_full
        df_adj = mask_negative_mu_df(df_adj, mu_col)
        mask_negative_mu_sweep!(sweep_dict)
    elseif (!data_full) && plot_full
        df_adj = mirror_mu_df(df_adj, mu_col)
        mirror_mu_sweep!(sweep_dict)
    end
    return df_adj, sweep_dict
end


#####################################
### Energy gap analysis functions ###
#####################################
function extract_gap_ranges_from_plateaus(
    plateaus,
    min_width::Float64,
    zero_tol::Float64,
    bridge_fraction::Float64
)
    ranges = Tuple{Float64, Float64}[]
    if plateaus === nothing || ismissing(plateaus) || isempty(plateaus)
        return ranges
    end

    sorted_gaps = sort(collect(plateaus); by = g -> g.N_gap)
    current_cluster = [sorted_gaps[1]]
    tol = max(bridge_fraction, 0.0) + 1e-12

    function flush_cluster!(cluster)
        isempty(cluster) && return
        low = cluster[1].E_low
        high = cluster[end].E_high
        width = abs(high - low)
        if width >= min_width && !(low <= zero_tol && high >= -zero_tol)
            push!(ranges, (low, high))
        end
    end

    prev_gap = sorted_gaps[1]
    for gap in Base.Iterators.drop(sorted_gaps, 1)
        idos_step = abs(gap.N_gap - prev_gap.N_gap)
        if idos_step <= tol
            push!(current_cluster, gap)
        else
            flush_cluster!(current_cluster)
            empty!(current_cluster)
            push!(current_cluster, gap)
        end
        prev_gap = gap
    end

    flush_cluster!(current_cluster)
    return ranges
end

function dominant_gap_count(counts::Vector{Int})
    isempty(counts) && return 0
    freq = Dict{Int, Int}()
    best_val = counts[1]
    best_freq = 0
    for c in counts
        new_freq = get!(freq, c, 0) + 1
        freq[c] = new_freq
        if new_freq > best_freq || (new_freq == best_freq && c > best_val)
            best_freq = new_freq
            best_val = c
        end
    end
    return best_val
end

# Helper to count gaps (runs of zeros) in discretized MP data
function count_zero_runs(vec::AbstractVector{<:Integer})
    cnt = 0
    runlen = 0
    min_gap_length = 1
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
    if runlen >= min_gap_length
        cnt += 1
    end
    return cnt
end

function process_eigenvalue_analysis!(
    df_slice::DataFrame,
    sweep_dict::Dict;
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu,
    plot_mu_range_full::Bool=true,
    zero_tol::Float64=1e-9,
    idos_bridge_fraction::Union{Nothing, Float64}=nothing
)
    """
        process_eigenvalue_analysis!(df_slice::DataFrame, sweep_dict::Dict;
                         eig_col=:eigenvalues, mu_col=:mu,
                         plot_mu_range_full=true, zero_tol=1e-9,
                         idos_bridge_fraction=nothing)

        Analyze eigenvalue data from df_slice and add analysis results to sweep_dict[:analysis].
        This includes:
        - Energy gap counts and ranges for each μ
        - Expected gap count (dominant across μ)
        - mu_c (maximum eigenvalue at μ≈0)
        - MP gap counts respecting mu_c
        - Valid tolerance indices

    The counting window is adjusted based on `plot_mu_range_full`:
    - true  ⇒ use μ ∈ [-μ_c, μ_c]
    - false ⇒ use μ ∈ [0, μ_c]

    IDOS-derived gaps are clustered so that up to `floor(frac * N)` intermediate
    states can be ignored when measuring a single energy gap, where `frac` is
    `idos_bridge_fraction` (defaulting to `1/N`).

        Modifies sweep_dict in place by adding or updating sweep_dict[:analysis].
    """
    # Extract parameters
    params = sweep_dict[:params]
    N_val = params.N
    Delta_val = params.Delta
    phi_val = params.phi
    
    # Get sweep data
    mu_vals = sweep_dict[:mu]
    mp_tol_range = sweep_dict[:mp_tol_range]
    disc_mp_matrix = sweep_dict[:disc_mp_matrix]
    
    n_mu = length(mu_vals)
    n_tols = length(mp_tol_range)
    
    # Sort df_slice by mu to align with sweep_dict
    df_sorted = sort(df_slice, mu_col)
    
    # Rationalize phi to find j
    r = rationalize(phi_val, tol=1e-10)
    j = denominator(r)

    bridge_fraction = isnothing(idos_bridge_fraction) ? (1.0 / max(N_val, 1)) : idos_bridge_fraction
    
    # Analyze energy gaps for each μ
    energy_gap_counts = Vector{Int}(undef, n_mu)
    energy_gap_ranges = Vector{Vector{Tuple{Float64, Float64}}}(undef, n_mu)
    
    df_sorted = IDOSProcessing.compute_idos_df!(df_sorted)
    df_sorted = IDOSProcessing.compute_plateaus_on_idos_df_extra_thresh!(
        df_sorted;
        threshold=Delta_val,
        # idos_threshold=nothing,
        eigs_col=eig_col
    )

    plateaus_col = :plateaus_idxd in names(df_sorted) ? df_sorted[!, :plateaus_idxd] : fill(nothing, n_mu)

    for (idx, plateaus) in enumerate(plateaus_col)
        gap_ranges = extract_gap_ranges_from_plateaus(
            plateaus,
            Delta_val,
            zero_tol,
            bridge_fraction
        )
        energy_gap_counts[idx] = length(gap_ranges)
        energy_gap_ranges[idx] = gap_ranges
    end
    
    # Calculate mu_c from eigenvalue at μ≈0
    zero_mu_idx = argmin(abs.(mu_vals))
    zero_eigs = collect(real.(df_sorted[zero_mu_idx, eig_col]))
    mu_c_value = isempty(zero_eigs) ? mu_vals[end] : maximum(zero_eigs)
    pos_mu_idx = searchsortedlast(mu_vals, mu_c_value + zero_tol)
    neg_mu_idx = searchsortedfirst(mu_vals, -mu_c_value - zero_tol)
    pos_mu_idx = clamp(pos_mu_idx, 1, n_mu)
    neg_mu_idx = clamp(neg_mu_idx, 1, n_mu)
    nonneg_idx = searchsortedfirst(mu_vals, -zero_tol)
    nonneg_idx = clamp(nonneg_idx, 1, n_mu)
    
    # Determine expected gaps
    zero_gap_ranges = energy_gap_ranges[zero_mu_idx]
    full_energy_gap_count = (isnothing(zero_gap_ranges) || ismissing(zero_gap_ranges)) ? 0 : length(zero_gap_ranges)
    expected_gaps_full = full_energy_gap_count
    expected_gaps = plot_mu_range_full ? expected_gaps_full : div(expected_gaps_full, 2)
    max_possible_gaps = max(0, 2 * N_val - 1)
    expected_gaps = clamp(expected_gaps, 0, max_possible_gaps)
    
    # Count MP gaps for each tolerance, respecting mu_c
    gap_counts = Vector{Int}(undef, n_tols)
    mu_c_vec = Vector{Float64}(undef, n_tols)
    
    for i in 1:n_tols
        col = disc_mp_matrix[:, i]
        last_one_idx = findlast(x -> x == 1, col)
        
        mu_c_vec[i] = mu_c_value

        if pos_mu_idx <= 0 || isnothing(last_one_idx)
            gap_counts[i] = 0
            continue
        end
        
        limit_idx = min(last_one_idx, pos_mu_idx)
        if limit_idx <= 0
            gap_counts[i] = 0
        else
            start_idx = plot_mu_range_full ? neg_mu_idx : nonneg_idx
            start_idx = min(max(start_idx, 1), limit_idx)
            if limit_idx < start_idx
                gap_counts[i] = 0
            else
                gap_counts[i] = count_zero_runs(col[start_idx:limit_idx])
            end
        end
    end
    
    # Identify valid tolerance indices
    valid_indices = findall(x -> x == expected_gaps, gap_counts)
    
    # Store analysis results
    sweep_dict[:analysis] = Dict(
        :j => j,
        :expected_gaps => expected_gaps,
        :expected_gaps_full => expected_gaps_full,
        :gap_counts => gap_counts,
        :valid_indices => valid_indices,
        :mu_c_vec => mu_c_vec,
        :mu_c_value => mu_c_value,
        :energy_gap_counts => energy_gap_counts,
        :energy_gap_ranges => energy_gap_ranges
    )
    
    return sweep_dict
end



#####################################
### Bulk/QC energy gap comparison ###
#####################################

function compute_bulk_gaps!(
    df::DataFrame;
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu,
    tol::Float64=1e-9
)
    # Ensure dataframe is sorted by mu to prevent plotting artifacts (connecting lines)
    if mu_col in propertynames(df)
        sort!(df, mu_col)
    end

    bulk_gaps = Float64[]
    sizehint!(bulk_gaps, nrow(df))

    # Iterate over each row's eigenvalues
    for evals in df[!, eig_col]
        # Assuming evals is sorted length is 2*N, so middle is at n and n+1
        n = length(evals) ÷ 2
        
        # Verify we have enough eigenvalues for the outer indices (n-1 and n+2)
        if n < 2
            push!(bulk_gaps, NaN)
            continue
        end

        # Calculate the gap between the middle two states (candidate for zero energy difference)
        mid_gap = abs(evals[n+1] - evals[n])

        if mid_gap < tol
            # If the middle gap is within tolerance (degenerate zero modes),
            # define bulk gap as difference between the next outer pair of states
            outer_gap = abs(evals[n+2] - evals[n-1])
            push!(bulk_gaps, outer_gap)
        else
            # Otherwise, use the middle gap
            push!(bulk_gaps, mid_gap)
        end
    end

    # println("Computed bulk gaps for ", length(bulk_gaps), " mu values.")

    df[!, :bulk_energy_gaps] = bulk_gaps
    return df, bulk_gaps
end

function compute_bulk_gaps_robust(
    df::DataFrame;
    eig_col::Symbol=:eigenvalues,
    mu_col::Symbol=:mu
)
    # Ensure dataframe is sorted by mu to prevent plotting artifacts (connecting lines)
    if mu_col in propertynames(df)
        sort!(df, mu_col)
    end

    n_rows = nrow(df)
    bulk_gaps = Vector{Float64}(undef, n_rows)
    mid_gaps = Vector{Float64}(undef, n_rows)
    outer_gaps = Vector{Float64}(undef, n_rows)

    # Iterate over each row's eigenvalues
    for (i, evals) in enumerate(df[!, eig_col])
        # Assuming evals is sorted length is 2*N, so middle is at n and n+1
        n = length(evals) ÷ 2

        # Calculate the gap between the middle two states (candidate for zero energy difference)
        mid_gap = abs(evals[n+1] - evals[n])
        outer_gap = abs(evals[n+2] - evals[n-1])

        bulk_gap = max(mid_gap, outer_gap)

        mid_gaps[i] = mid_gap
        outer_gaps[i] = outer_gap
        bulk_gaps[i] = bulk_gap
    end

    writedlm("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/access_results/bulk_gaps.txt", bulk_gaps)
    writedlm("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/access_results/mid_gaps.txt", mid_gaps)
    writedlm("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/access_results/outer_gaps.txt", outer_gaps)
    writedlm("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/access_results/mu_vals.txt", df[!, mu_col])

    return mid_gaps, outer_gaps, bulk_gaps
end

function compute_qc_gaps_direct(
    df::DataFrame;
    mu_val::Float64=0.0,
    mu_col::Symbol=:mu,
    eig_col::Symbol=:eigenvalues,
    eig_tol::Float64=1e-9,
    N_val=nothing
)
    # Find index for closest mu in the dataframe
    mu_vals = df[!, mu_col]
    _, idx = findmin(abs.(mu_vals .- mu_val))

    # Get sorted eigenvalues at that mu
    eigs_at_mu = sort(df[idx, eig_col])
    writedlm("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/access_results/eigs_at_mu_mu$(mu_val)_N$(N_val).txt", eigs_at_mu)

    n_eigs = length(eigs_at_mu)
    
    # Define E_range based on spectrum spread
    # Length matches the number of mu points (nrow(df)) to maintain consistent x-axis resolution for plotting
    # x_length = length(df.mu)
    x_length = 15000
    # println("x_length: ", x_length)
    # x_max = maximum(df.mu)
    # println("x_max: ", x_max)
    E_min = minimum(eigs_at_mu)
    E_max = maximum(eigs_at_mu)
    E_range = collect(range(E_min, E_max, length=x_length))
    # println("E_range: ", E_range)
    
    qc_gaps = Vector{Float64}(undef, x_length)

    for (i, E) in enumerate(E_range)
        # Find index of the first eigenvalue >= E
        # searchsortedfirst is O(log N)
        k = searchsortedfirst(eigs_at_mu, E)
        
        # Check if we hit an eigenvalue exactly (gap = 0)
        is_on_state = false
        if k <= n_eigs && abs(eigs_at_mu[k] - E) < eig_tol
            is_on_state = true
        elseif k > 1 && abs(eigs_at_mu[k-1] - E) < eig_tol
            # Check k-1 in case E is slightly above eigenval due to float precision
            is_on_state = true
        end

        if is_on_state
            qc_gaps[i] = 0.0
        elseif k == 1 || k > n_eigs
            # E is outside the spectral range (below lowest or above highest)
            qc_gaps[i] = 0.0
        else
            # E is strictly inside a gap between eigs[k-1] and eigs[k]
            qc_gaps[i] = eigs_at_mu[k] - eigs_at_mu[k-1]
        end
    end
    
    return E_range, qc_gaps
end

function compute_qc_gaps_idos(
    df::DataFrame;
    mu_val::Float64=0.0,
    mu_col::Symbol=:mu,
    eig_col::Symbol=:eigenvalues,
    plateaus_col::Symbol=:plateaus_idxd
)
    # Find index for closest mu
    mu_vals = df[!, mu_col]
    _, idx = findmin(abs.(mu_vals .- mu_val))

    # Retrieve plateaus for the selected mu
    # Try the provided column name, or fallback to :plateaus_idxd if specific one not found
    plats = nothing
    if plateaus_col in propertynames(df)
        plats = df[idx, plateaus_col]
    elseif :plateaus_idxd in propertynames(df)
        plats = df[idx, :plateaus_idxd]
    elseif :gap_labels in propertynames(df)
        plats = df[idx, :gap_labels]
    end

    # Determine E_range from eigenvalues to match direct method's range
    eigs_at_mu = df[idx, eig_col]
    
    if isempty(eigs_at_mu)
        return Float64[], Float64[]
    end

    E_min = minimum(eigs_at_mu)
    E_max = maximum(eigs_at_mu)
    x_length = nrow(df)
    E_range = collect(range(E_min, E_max, length=x_length))
    
    # Initialize gaps with zeros (representing 'on-band' regions)
    qc_gaps = zeros(Float64, x_length)
    
    if plats === nothing || ismissing(plats) || isempty(plats)
        return E_range, qc_gaps
    end

    # Iterate over IDOS plateaus and mark gap sizes
    for p in plats
        # Handle different structures
        # 1. NamedTuple with E_low, E_high (from compute_plateaus_on_idos_df_extra_thresh! or gap_labels)
        # 2. Tuple (j, IDOS) - Old format, likely needs eigenvalues to look up E_low/High
        
        low, high = 0.0, 0.0
        
        if hasproperty(p, :E_low) && hasproperty(p, :E_high)
            low = p.E_low
            high = p.E_high
        elseif p isa Tuple && length(p) == 2 && p[1] isa Int
             # Fallback for old (j, IDOS) tuple format if eigs available
             # p[1] is the index j below the gap
             j = p[1]
             if j >= 1 && j < length(eigs_at_mu)
                 low = eigs_at_mu[j]
                 high = eigs_at_mu[j+1]
             else
                 continue
             end
        else
            continue
        end

        gap_size = high - low
        
        # Optimisation: Find index range in E_range that could overlap with this plateau
        i_start = searchsortedfirst(E_range, low)
        i_end = searchsortedlast(E_range, high)
        
        if i_start > i_end
            continue
        end

        # Assign gap size only if E is strictly inside the plateau
        for i in i_start:i_end
            val = E_range[i]
            if val > low && val < high
                qc_gaps[i] = gap_size
            end
        end
    end

    return E_range, qc_gaps
end

## Merges (just one following) gap if next gap is larger than qc_gap_tol
function compute_qc_gaps_direct_robust(
    df::DataFrame;
    mu_val::Float64=0.0,
    mu_col::Symbol=:mu,
    eig_col::Symbol=:eigenvalues,
    hit_eig_tol::Float64=1e-9,
    qc_gap_tol::Float64=1e-8
)
    # Find index for closest mu in the dataframe
    mu_vals = df[!, mu_col]
    _, idx = findmin(abs.(mu_vals .- mu_val))

    # Get sorted eigenvalues at that mu
    eigs_at_mu = sort(df[idx, eig_col])
    n_eigs = length(eigs_at_mu)
    
    # Define E_range based on spectrum spread
    x_length = length(df.mu)
    
    E_min = minimum(eigs_at_mu)
    E_max = maximum(eigs_at_mu)
    E_range = collect(range(E_min, E_max, length=x_length))

    # 1. Identify Gap Regions and Merge Small Intervening "Islands"
    # We store regions as (start_E, end_E, gap_size)
    gap_regions = Tuple{Float64, Float64, Float64}[]
    merged_count = 0

    if n_eigs > 1
        i = 1
        while i < n_eigs
            # Current gap bounds
            g1_start = eigs_at_mu[i]
            g1_end = eigs_at_mu[i+1]
            g1_size = g1_end - g1_start

            merged = false
            if i + 1 < n_eigs
                g2_start = eigs_at_mu[i+1]
                g2_end = eigs_at_mu[i+2]
                g2_size = g2_end - g2_start
                
                # Check merge condition
                if g2_size > qc_gap_tol
                    # Check if there's a 3rd gap to merge as well
                    if i + 2 < n_eigs
                        # We merge g1, g2, AND g3
                        g3_start = eigs_at_mu[i+2]
                        g3_end = eigs_at_mu[i+3]
                        # g3_size = g3_end - g3_start 
                        
                        # Merge g1, g2, g3 => start of g1 to end of g3
                        total_gap_size = g3_end - g1_start
                        push!(gap_regions, (g1_start, g3_end, total_gap_size))
                        
                        i += 3 # Skip g2 and g3 processing
                        merged_count += 1 # Considers this block 1 merge event (or 2?)
                        merged = true
                    else
                        # Only g2 available (at end of list), so just merge g1 and g2
                        total_gap_size = g2_end - g1_start
                        push!(gap_regions, (g1_start, g2_end, total_gap_size))
                        i += 2 # Skip g2 processing
                        merged_count += 1
                        merged = true
                    end
                end
            end
            
            if !merged
                push!(gap_regions, (g1_start, g1_end, g1_size))
                i += 1
            end
        end
    end
    
    # 2. Map regions to E_range
    qc_gaps = zeros(Float64, x_length)
    
    for (g_start, g_end, g_size) in gap_regions
        # Find indices in E_range strictly inside (g_start, g_end)
        i_start = searchsortedfirst(E_range, g_start + hit_eig_tol)
        i_end_cand = searchsortedfirst(E_range, g_end - hit_eig_tol)
        
        # Handle searchsorted returning index >= val
        # We want index where E <= g_end - tol, so if E[i] >= val, we take i-1
        if i_end_cand <= x_length && E_range[i_end_cand] >= g_end - hit_eig_tol
             i_end = i_end_cand - 1
        else
             i_end = i_end_cand
             if i_end > x_length; i_end = x_length; end
        end
        
        if i_start <= i_end
            qc_gaps[i_start:i_end] .= g_size
        end
    end
    
    return E_range, qc_gaps, merged_count
end


function count_qc_beats_sc_regions(
    mu_vals::AbstractVector{Float64},
    bulk_gaps::AbstractVector{Float64},
    E_qc::AbstractVector{Float64},
    qc_gaps::AbstractVector{Float64};
    xlims::Union{Nothing, Tuple{Float64, Float64}}=nothing
)
    @assert length(mu_vals) == length(bulk_gaps) "mu_vals and bulk_gaps must have same length"
    @assert length(E_qc) == length(qc_gaps) "E_qc and qc_gaps must have same length"

    # Interpolate QC gaps onto the mu grid
    qc_interp = fill(NaN, length(mu_vals))
    E_min, E_max = minimum(E_qc), maximum(E_qc)
    
    for (i, mu) in enumerate(mu_vals)
        if mu >= E_min && mu <= E_max
            k = searchsortedlast(E_qc, mu)
            
            val = 0.0
            if k == 0
                val = qc_gaps[1]
            elseif k >= length(E_qc)
                val = qc_gaps[end]
            else
                x0 = E_qc[k]
                x1 = E_qc[k+1]
                y0 = qc_gaps[k]
                y1 = qc_gaps[k+1]
                
                # Linear Interpolation
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
    diff_vals = bulk_gaps .- qc_interp

    # Count disconnected negative regions
    region_count = 0.0
    in_region = false
    count_min, count_max = isnothing(xlims) ? (-Inf, Inf) : xlims

    for i in eachindex(diff_vals)
        cur_mu = mu_vals[i]
        val = diff_vals[i]
        
        # Check if we are inside valid range
        is_in_range = (cur_mu >= count_min) && (cur_mu <= count_max)
        
        # Determine if this point represents a negative gap condition
        # Must be valid (not NaN), negative, and within xlims
        is_negative = !isnan(val) && val < 0 && is_in_range

        if is_negative
            if !in_region
                in_region = true
                region_count += 1.0
            end
        else
            # Reset region if we hit non-negative, NaN, or out-of-bounds
            in_region = false
        end
    end

    return region_count
end


##############################
#### mp_tol range finding ####
##############################

# function find_mptol_ranges_for_expected_gaps(
#     sweep_data::Dict,
#     mu_count_range::Tuple{Float64, Float64},
#     expected_gaps::Int
# )
#     mu_vals = sweep_data[:mu]
#     mp_tol_range = sweep_data[:mp_tol_range]
#     disc_mp_matrix = sweep_data[:disc_mp_matrix]

#     # Filter mu indices within the range
#     min_mu, max_mu = mu_count_range
#     if min_mu > max_mu; min_mu, max_mu = max_mu, min_mu; end

#     mask = (mu_vals .>= min_mu) .& (mu_vals .<= max_mu)
#     valid_idxs = findall(mask)
    
#     if isempty(valid_idxs)
#         return Tuple{Float64, Float64}[]
#     end

#     matched_indices = Int[]
#     n_tols = length(mp_tol_range)

#     # Iterate through all tolerances and count gaps in the valid mu window
#     for i in 1:n_tols
#         col_slice = view(disc_mp_matrix, valid_idxs, i)
        
#         # Use existing helper to count gaps (runs of zeros)
#         # Note: count_zero_runs requires AbstractVector{<:Integer}
#         # If disc_mp_matrix is Float, this relies on it being Integer-like or compatible
#         gaps = count_zero_runs(col_slice)
        
#         if gaps == expected_gaps
#             push!(matched_indices, i)
#         end
#     end

#     # Group consecutive indices into ranges
#     ranges = Tuple{Float64, Float64}[]
#     if isempty(matched_indices)
#         return ranges
#     end

#     start_idx = matched_indices[1]
#     prev_idx = start_idx

#     for i in 2:length(matched_indices)
#         curr_idx = matched_indices[i]
#         # Check if indices are consecutive to form a range
#         if curr_idx == prev_idx + 1
#             prev_idx = curr_idx
#         else
#             # End of a contiguous block
#             push!(ranges, (mp_tol_range[start_idx], mp_tol_range[prev_idx]))
#             start_idx = curr_idx
#             prev_idx = curr_idx
#         end
#     end
#     # Push the final range
#     push!(ranges, (mp_tol_range[start_idx], mp_tol_range[prev_idx]))

#     # Return only highest mp_tol valued range if multiple found
#     if isempty(ranges)
#         return (0.0, 0.0)
#     else
#         highest_range = ranges[argmax(last.(ranges))]
#         return highest_range
#     end
# end

function find_mptol_ranges_for_expected_gaps(
    sweep_data::Dict,
    mu_count_range::Tuple{Float64, Float64},
    expected_gaps::Int
)
    mu_vals = sweep_data[:mu]
    mp_tol_range = sweep_data[:mp_tol_range]
    disc_mp_matrix = sweep_data[:disc_mp_matrix]

    # Filter mu indices within the range
    min_mu, max_mu = mu_count_range
    if min_mu > max_mu; min_mu, max_mu = max_mu, min_mu; end

    mask = (mu_vals .>= min_mu) .& (mu_vals .<= max_mu)
    valid_idxs = findall(mask)
    
    if isempty(valid_idxs)
        return nothing  # Changed from Tuple{Float64, Float64}[] to nothing
    end

    matched_indices = Int[]
    n_tols = length(mp_tol_range)

    # Iterate through all tolerances and count gaps in the valid mu window
    for i in 1:n_tols
        col_slice = view(disc_mp_matrix, valid_idxs, i)
        
        # Use existing helper to count gaps (runs of zeros)
        # Note: count_zero_runs requires AbstractVector{<:Integer}
        # If disc_mp_matrix is Float, this relies on it being Integer-like or compatible
        gaps = count_zero_runs(col_slice)
        
        if gaps == expected_gaps
            push!(matched_indices, i)
        end
    end

    # Group consecutive indices into ranges
    ranges = Tuple{Float64, Float64}[]
    if isempty(matched_indices)
        return (0.0,0.0)  # Or (0.0, 0.0) if you prefer a default tuple
    end

    start_idx = matched_indices[1]
    prev_idx = start_idx

    for i in 2:length(matched_indices)
        curr_idx = matched_indices[i]
        # Check if indices are consecutive to form a range
        if curr_idx == prev_idx + 1
            prev_idx = curr_idx
        else
            # End of a contiguous block
            push!(ranges, (mp_tol_range[start_idx], mp_tol_range[prev_idx]))
            start_idx = curr_idx
            prev_idx = curr_idx
        end
    end
    # Push the final range
    push!(ranges, (mp_tol_range[start_idx], mp_tol_range[prev_idx]))

    # Return only highest mp_tol valued range if multiple found
    if isempty(ranges)
        return nothing  # Or (0.0, 0.0)
    else
        highest_range = ranges[argmax(last.(ranges))]
        return highest_range
    end
end


end  # module MPGapProcessing