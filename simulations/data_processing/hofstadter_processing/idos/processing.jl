module IDOSProcessing

using DataFrames
using ProgressMeter

############################################
######## Eigenvalue Normalisations #########
############################################

function compute_eigs_bandwidth_norm!(
    df::DataFrame; 
    eigs_col::Symbol=:eigenvalues, 
    out_col::Symbol=:eigs_global_norm, 
    tol::Float64=1e-6,
    min_out_col::Union{Symbol, Nothing}=nothing,
    max_out_col::Union{Symbol, Nothing}=nothing
)
    df[!, out_col] = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df))
    
    # Initialize new columns if requested
    if !isnothing(min_out_col)
        df[!, min_out_col] = Vector{Union{Float64, Missing}}(undef, nrow(df))
    end
    if !isnothing(max_out_col)
        df[!, max_out_col] = Vector{Union{Float64, Missing}}(undef, nrow(df))
    end
    
    for (i, row) in enumerate(eachrow(df))
        es = row[eigs_col]
        if es === missing || isempty(es)
            df[i, out_col] = missing
            if !isnothing(min_out_col)
                df[i, min_out_col] = missing
            end
            if !isnothing(max_out_col)
                df[i, max_out_col] = missing
            end
            continue
        end

        vals = sort!(collect(real(es)))
        # Drop energies in 0 ± tol
        vals = [x for x in vals if abs(x) > tol]
        if isempty(vals)
            df[i, out_col] = Float64[]
            if !isnothing(min_out_col)
                df[i, min_out_col] = missing
            end
            if !isnothing(max_out_col)
                df[i, max_out_col] = missing
            end
            continue
        end

        # Find global min and max
        global_min = minimum(vals)
        global_max = maximum(vals)

        # Store min and max if requested
        if !isnothing(min_out_col)
            df[i, min_out_col] = global_min
        end
        if !isnothing(max_out_col)
            df[i, max_out_col] = global_max
        end

        # Linear mapping: min -> -1, max -> 1
        if global_max == global_min
            # All values are the same, map to 0
            normed = zeros(Float64, length(vals))
        else
            normed = [2 * (x - global_min) / (global_max - global_min) - 1 for x in vals]
        end

        df[i, out_col] = normed
    end
    return df
end

function compute_eigs_norm!(df::DataFrame; eigs_col::Symbol=:eigenvalues, out_col::Symbol=:eigs_norm, tol::Float64=1e-6)
    df[!, out_col] = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        es = row[eigs_col]
        if es === missing || isempty(es)
            df[i, out_col] = missing
            continue
        end

        vals = sort!(collect(real(es)))
        # drop energies in 0 ± tol
        vals = [x for x in vals if abs(x) > tol]
        if isempty(vals)
            df[i, out_col] = Float64[]
            continue
        end

        negs = [x for x in vals if x < 0.0]
        poss = [x for x in vals if x > 0.0]

        # precompute per-side scales so that:
        # - min |neg| -> 0, max |neg| -> -1
        # - min pos   -> 0, max pos   -> 1
        amin, amax = if isempty(negs)
            (0.0, 1.0)
        else
            a = abs.(negs)
            (minimum(a), maximum(a))
        end
        pmin, pmax = if isempty(poss)
            (0.0, 1.0)
        else
            (minimum(poss), maximum(poss))
        end

        normed = [ x < 0 ?
                   (amax == amin ? -1.0 : -((abs(x) - amin) / (amax - amin))) :
                   (pmax == pmin ?  1.0 :  ((x      - pmin) / (pmax - pmin)))
                   for x in vals ]

        df[i, out_col] = normed
    end
    return df
end

function compute_eigs_inner_norm!(df::DataFrame; eigs_col::Symbol=:eigenvalues, out_col::Symbol=:eigs_inner_norm, tol::Float64=1e-6)
    df[!, out_col] = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        es = row[eigs_col]
        if es === missing || isempty(es)
            df[i, out_col] = missing
            continue
        end

        vals = sort!(collect(real(es)))
        # drop energies in 0 ± tol
        vals = [x for x in vals if abs(x) > tol]
        if isempty(vals)
            df[i, out_col] = Float64[]
            continue
        end

        negs = [x for x in vals if x < 0.0]  # in [-amax, -amin]
        poss = [x for x in vals if x > 0.0]  # in [ pmin,  pmax]

        # inner/outer magnitudes
        amin, amax = isempty(negs) ? (0.0, 0.0) : (minimum(abs.(negs)), maximum(abs.(negs)))
        pmin, pmax = isempty(poss) ? (0.0, 0.0) : (minimum(poss),          maximum(poss))

        # linear maps that set inner edge -> 0 and keep outer edge unchanged
        # negatives: f(x) = a*x + b, with a = amax/(amax-amin), b = a*amin
        # positives: g(x) = c*x + d, with c = pmax/(pmax-pmin), d = -c*pmin
        a = (amax > amin) ? (amax / (amax - amin)) : 1.0
        b = (amax > amin) ? (a * amin)             : 0.0
        c = (pmax > pmin) ? (pmax / (pmax - pmin)) : 1.0
        d = (pmax > pmin) ? (-c * pmin)            : 0.0

        mapped = similar(vals)
        @inbounds for k in eachindex(vals)
            x = vals[k]
            mapped[k] = x < 0 ? (a*x + b) : (c*x + d)
        end

        df[i, out_col] = mapped
    end
    return df
end


"""
    compute_interpolated_normalization!(df; mode=:both, ref_indices=(1,1), out_col=:eigs_interp_norm)

Normalizes eigenvalues with interpolated scaling (Outer) and bulk-edge zeroing (Inner).

# Arguments
- `mode`: :inner, :outer, or :both.
- `ref_indices`: Tuple (start_offset, end_offset). 
   (1, 1) uses the very first and very last slope values.
   (1, 5) uses the first and the 5th-from-last.
- `interp_func`: Optional function (phi_start, m_start, phi_end, m_end) -> (phi -> multiplier).
   Defaults to linear interpolation of the multiplier.
"""
function compute_interpolated_normalisation!(
    df::DataFrame; 
    out_col::Symbol=:eigs_interp_norm,
    mode::Symbol=:both, # :inner, :outer, :both
    ref_indices::Tuple{Int, Int}=(1, 1), 
    interp_func=nothing
)
    
    # 1. Get unique phis and sort
    phis = sort(unique(df.phi))
    n_phis = length(phis)
    
    # 2. Resolve Reference Indices
    i_start = clamp(ref_indices[1], 1, n_phis)
    i_end   = clamp(n_phis - ref_indices[2] + 1, 1, n_phis)
    phi_start = phis[i_start]
    phi_end   = phis[i_end]

    # 3. Pre-calculate Bulk Edges (for Inner Norm) if needed
    # We need to know the bulk edge BEFORE scaling to do Inner first.
    # Or, we can do it row-by-row.
    
    # We need to determine the SCALING factors based on the "Inner-Normalized" data
    # if we want the final result to be exactly [-1, 1].
    
    # Helper to get bulk edge
    get_bulk_edge = (eigs) -> begin
        abs_e = abs.(eigs)
        candidates = filter(x -> x > 1e-4, abs_e) # Threshold for MBS
        return isempty(candidates) ? 0.0 : minimum(candidates)
    end

    # Helper to apply inner norm
    apply_inner = (eigs, edge) -> [sign(e) * max(0.0, abs(e) - edge) for e in eigs]

    # Calculate Reference Maxima *after* hypothetical inner norm
    row_start = first(filter(r -> r.phi == phi_start, df))
    row_end   = first(filter(r -> r.phi == phi_end, df))
    
    eigs_start = real(row_start.eigenvalues)
    eigs_end   = real(row_end.eigenvalues)
    
    if mode == :inner || mode == :both
        edge_start = get_bulk_edge(eigs_start)
        edge_end   = get_bulk_edge(eigs_end)
        eigs_start = apply_inner(eigs_start, edge_start)
        eigs_end   = apply_inner(eigs_end, edge_end)
    end
    
    max_E_start = maximum(abs.(eigs_start))
    max_E_end   = maximum(abs.(eigs_end))
    
    m_start = 1.0 / (max_E_start > 1e-9 ? max_E_start : 1.0)
    m_end   = 1.0 / (max_E_end > 1e-9 ? max_E_end : 1.0)
    
    # 4. Interpolation Function
    if isnothing(interp_func)
        get_mult = (phi) -> begin
            if phi_end == phi_start; return m_start; end
            t = (phi - phi_start) / (phi_end - phi_start)
            return m_start + t * (m_end - m_start)
        end
    else
        get_mult = interp_func(phi_start, m_start, phi_end, m_end)
    end

    # 5. Apply to all rows
    new_col = Vector{Vector{Float64}}(undef, nrow(df))
    
    for (i, row) in enumerate(eachrow(df))
        eigs = real(row.eigenvalues)
        phi = row.phi
        
        # A. INNER NORMALIZATION (Shift) - FIRST
        if mode == :inner || mode == :both
            edge = get_bulk_edge(eigs)
            eigs = apply_inner(eigs, edge)
        end
        
        # B. OUTER NORMALIZATION (Scale) - SECOND
        if mode == :outer || mode == :both
            mult = get_mult(phi)
            eigs = eigs .* mult
        end
        
        new_col[i] = eigs
    end
    
    df[!, out_col] = new_col
    return df
end



"""
    compute_smooth_normalisation!(df; out_col=:eigs_smooth_norm, eigs_col=:eigenvalues, poly_deg=4)

Fits a polynomial envelope to the maximum eigenvalues across phi and normalizes by this smooth envelope.
This avoids shearing artifacts caused by local cusps at rational phi values.
"""
function compute_smooth_normalisation!(df::DataFrame; out_col::Symbol=:eigs_smooth_norm, eigs_col::Symbol=:eigenvalues, poly_deg::Int=4, coeffs_out_file::Union{String, Nothing}=nothing)
    # 1. Extract max energy for each phi
    data_xy = []
    for row in eachrow(df)
        es = row[eigs_col]
        if !ismissing(es) && !isempty(es)
            push!(data_xy, (row.phi, maximum(abs.(real(es)))))
        end
    end
    
    if isempty(data_xy)
        df[!, out_col] = [Float64[] for _ in 1:nrow(df)]
        return df
    end

    phis = [x[1] for x in data_xy]
    ymax = [x[2] for x in data_xy]

    # 2. Fit Polynomial using Least Squares: y ≈ c0 + c1*phi + ...
    A = zeros(length(phis), poly_deg + 1)
    for i in 1:length(phis)
        for p in 0:poly_deg
            A[i, p+1] = phis[i]^p
        end
    end
    
    # Solve for coeffs c
    c = A \ ymax

    # Save coefficients to file if requested
    if !isnothing(coeffs_out_file)
        # Write c (coefficients) to file
        open(coeffs_out_file, "w") do io
            println(io, "degree,coeff_index,value")
            for (idx, val) in enumerate(c)
                println(io, "$poly_deg,$(idx-1),$val")
            end
        end
    end

    # 3. Apply Normalization
    new_col = Vector{Vector{Float64}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        es = row[eigs_col]
        phi = row.phi
        if ismissing(es) || isempty(es)
            new_col[i] = Float64[]
            continue
        end

        # Evaluate polynomial envelope at this phi
        norm_val = 0.0
        for p in 0:poly_deg
            norm_val += c[p+1] * phi^p
        end

        # Normalize
        factor = norm_val > 1e-9 ? 1.0 / norm_val : 1.0
        new_col[i] = real(es) .* factor
    end
    
    df[!, out_col] = new_col
    return df
end


############################################
####### Phason-independent Spectrum ########
############################################
function compute_phason_independent_spectrum_gridTol(
    df::DataFrame; 
    eigs_col::Symbol=:eigenvalues, 
    tolerance::Float64=0.01, 
    grid_spacing::Float64=0.005,
    verbose::Bool=false
)

    """
        compute_phason_independent_spectrum(df; eigs_col=:eigenvalues, tolerance=1e-2, grid_spacing=1e-3, verbose=false)

        Calculates the set of energies that are persistently present (within `tolerance`) across all phason configurations in the DataFrame.
        This corresponds to the intersection of all spectra, effectively retaining the "bulk" Cantor dust bands while filtering out 
        subgap states that move as a function of phason.

        # Arguments
        - `df`: DataFrame containing eigenvalues (one row per phason, all rows should have same slope).
        - `eigs_col`: Column symbol for eigenvalues.
        - `tolerance`: Maximum allowable distance to the nearest eigenvalue. If an energy requires a jump > tolerance to find a match in ANY phason slice, it is discarded.
        - `grid_spacing`: Resolution of the output spectrum.

        # Returns
        - `Vector{Float64}`: A dense collection of energy points representing the stable bands.
    """
    # Filter for rows that actually have data
    valid_rows = filter(r -> !ismissing(r[eigs_col]) && !isempty(r[eigs_col]), df)
    
    if nrow(valid_rows) == 0
        return Float64[]
    end

    # 1. Determine global energy range efficiently
    local global_min, global_max
    
    first_eigs = real(valid_rows[1, eigs_col])
    global_min, global_max = extrema(first_eigs)
    
    # scan remaining rows to expand range
    for row in eachrow(valid_rows)
        es = real(row[eigs_col])
        mi, ma = extrema(es)
        if mi < global_min; global_min = mi; end
        if ma > global_max; global_max = ma; end
    end
    
    # 2. Define search grid
    # Expand slightly to process edges correctly
    grid = (global_min - tolerance):grid_spacing:(global_max + tolerance)
    n_points = length(grid)
    
    if n_points == 0
        return Float64[]
    end

    # max_dists[i] will store the Maximum distance found so far for grid point i across the rows processed
    # We define persistent spectrum as energies where max_dists[i] <= tolerance
    max_dists = zeros(Float64, n_points) 
    
    # Helper: efficiently find min distance from x to a sorted vector
    function dist_to_sorted_set(x::Float64, sorted_set::AbstractVector{<:Real})
        # searchsortedfirst returns the index of the first element >= x
        idx = searchsortedfirst(sorted_set, x)
        
        if idx == 1
            return abs(x - sorted_set[1])
        elseif idx > length(sorted_set)
            return abs(x - sorted_set[end])
        else
            # x is between sorted_set[idx-1] and sorted_set[idx]
            d_left = abs(x - sorted_set[idx-1])
            d_right = abs(x - sorted_set[idx])
            return min(d_left, d_right)
        end
    end

    if verbose
        p = Progress(nrow(valid_rows), desc="Computing independent spectrum: ")
    end

    # 3. Iterate over every phason configuration (row)
    for row in eachrow(valid_rows)
        # Ensure eigenvalues are real and sorted (usually they are, but safety first)
        eigs = sort(real(row[eigs_col]))
        
        # Check every grid point against this spectrum
        for i in 1:n_points
            
            # Optimization: If this grid point is already disqualified by a previous row, skip it
            # (Note: This works because we strictly require ALL rows to effectively match. 
            # However, since we want the actual max_dist for debugging/tuning, we might want to keep calculating.
            # Here we skip for speed if we only care about the binary result.)
            if max_dists[i] > tolerance
                continue
            end
            
            d = dist_to_sorted_set(grid[i], eigs)
            
            # We want the worst-case match across all phasons
            if d > max_dists[i]
                max_dists[i] = d
            end
        end
        
        if verbose; next!(p); end
    end
    
    # 4. Filter the grid to return only the "robust" energies
    robust_indices = findall(d -> d <= tolerance, max_dists)
    
    return collect(grid[robust_indices])
end


function compute_phason_independent_spectrum_specDens(
    df::DataFrame; 
    eigs_col::Symbol=:eigenvalues, 
    threshold::Float64=0.9, 
    n_bins::Int=1000,
    verbose::Bool=false,
)

    """
        compute_phason_independent_spectrum(df; eigs_col=:eigenvalues, n_bins=1000, threshold=0.9, verbose=false)

        Calculates the "stable" spectrum by computing the occupancy probability of energy bins across all phason configurations.
        1. The energy axis is divided into `n_bins` (or bins of `grid_spacing` if provided).
        2. For each phason slice (row), we mark which bins contain at least one eigenvalue.
        3. We sum these marks to get an 'occupancy count' for each bin.
        4. Bins with `occupancy / n_slices >= threshold` are retained.

        This method is robust to small energetic jitters and occasional outliers.

        # Arguments
        - `threshold`: Fraction of phason slices that must contain a state in a bin for it to be kept (default 0.9).
        - `n_bins`: Number of bins for the histogram. Higher N -> finer resolution but requires more stability.
    """
    # Filter for valid rows
    valid_rows = filter(r -> !ismissing(r[eigs_col]) && !isempty(r[eigs_col]), df)
    n_slices = nrow(valid_rows)
    
    if n_slices == 0
        return Float64[]
    end

    if verbose
        println("Computing independent spectrum across $n_slices phason slices...")
    end

    # 1. Determine global range
    all_eigs_iter = Iterators.flatten(valid_rows[!, eigs_col])
    global_min, global_max = extrema(real, all_eigs_iter)
    
    # Add small margin
    span = global_max - global_min
    margin = span * 0.01
    emin = global_min - margin
    emax = global_max + margin
    
    # Define bins
    edges = range(emin, emax, length=n_bins+1)
    bin_width = step(edges)
    
    # 2. Count Occupancy
    # occupancy[i] counts how many rows have at least one eigenvalue in bin i
    occupancy = zeros(Int, n_bins)
    
    for row in eachrow(valid_rows)
        es = real(row[eigs_col])
        
        # Find which bins are hit by this row's eigenvalues
        # Using Set to ensure we count a row only once per bin
        hit_bins = Set{Int}()
        
        for e in es
            # Map e to bin index
            # idx = 1 + floor((e - emin) / bin_width)
            idx = floor(Int, (e - emin) / bin_width) + 1
            
            # clamp to be safe with float precision at edges
            idx = clamp(idx, 1, n_bins)
            push!(hit_bins, idx)
        end
        
        for b in hit_bins
            occupancy[b] += 1
        end
    end
    
    # 3. Filter
    min_count = ceil(Int, threshold * n_slices)
    valid_indices = findall(c -> c >= min_count, occupancy)
    
    valid_centers = [(edges[i] + edges[i+1])/2 for i in valid_indices]
    
    if verbose
        println("  -> Retained $(length(valid_centers)) / $n_bins bins (threshold $(threshold*100)%)")
    end

    return valid_centers
end


############################################
########### IDOS Computation ###############
############################################

function compute_idos_df!(df::DataFrame)
    idos_col = Vector{Vector{Float64}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        eigenvalues = sort(real(row.eigenvalues))
        N = length(eigenvalues)
        idos = collect(1:N) ./ N
        idos_col[i] = idos
    end
    df.idos = idos_col
    return df
end

function find_idos_plateaus(df::DataFrame; threshold=0.05, fixed_values...)
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plateau_idos = Vector{Vector{Float64}}(undef, nrow(filtered_df))
    for (i, row) in enumerate(eachrow(filtered_df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        gaps = Int[]
        for j in 1:length(energies)-1
            if energies[j+1] - energies[j] > threshold
                push!(gaps, j)
            end
        end
        # IDOS value at the gap is between j and j+1
        plateau_idos[i] = [(j + 0.5) / length(idos) for j in gaps]
    end
    return plateau_idos
end

function compute_plateaus_df!(df::DataFrame; threshold=0.01)
    plateaus_col = Vector{Vector{Float64}}(undef, nrow(df))
    @showprogress for (i, row) in enumerate(eachrow(df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        gaps = Int[]
        for j in 1:length(energies)-1
            if energies[j+1] - energies[j] > threshold
                push!(gaps, j)
            end
        end
        plateaus_col[i] = [(j + 0.5) / length(idos) for j in gaps]
    end
    df.plateaus = plateaus_col
    return df
end


#######################################
######## Gap Labelling Theorem ########
#######################################

function compute_plateaus_on_idos_df!(df::DataFrame; threshold=0.01)
    plateaus_idxd = Vector{Vector{Tuple{Int,Float64}}}(undef, nrow(df))
    plateaus = Vector{Vector{Float64}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos

        gap_list = Tuple{Int,Float64}[]
        plateaus_idos = Float64[]
        for j in 1:length(energies)-1
            if energies[j+1] - energies[j] > threshold
                # store (index j, IDOS[j])
                push!(gap_list, (j, idos[j]))
                push!(plateaus_idos, idos[j])
            end
        end
        plateaus_idxd[i] = gap_list
        plateaus[i] = plateaus_idos
    end

    df.plateaus = plateaus
    df.plateaus_idxd = plateaus_idxd
    return df
end

##########################################################################################
##### New plateaus function with second threshold on IDOS proximity for merging gaps #####
##########################################################################################

function compute_plateaus_on_idos_df_extra_thresh!(df::DataFrame; threshold=0.01, eigs_col::Symbol=:eigenvalues)
    # Store explicit bounds: (E_low, E_high, N_gap)
    plateaus_idxd = Vector{Vector{NamedTuple{(:E_low, :E_high, :N_gap), Tuple{Float64, Float64, Float64}}}}(undef, nrow(df))
    plateaus = Vector{Vector{Float64}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        # Use the specified column (e.g. :eigs_interp_norm)
        energies = sort(real(row[eigs_col]))
        N = length(energies)
        
        raw_gaps = NamedTuple{(:E_low, :E_high, :N_gap), Tuple{Float64, Float64, Float64}}[]
        for j in 1:N-1
            dE = energies[j+1] - energies[j]
            if dE > threshold
                push!(raw_gaps, (E_low=energies[j], E_high=energies[j+1], N_gap=j/N))
            end
        end

        plateaus_idxd[i] = raw_gaps
        plateaus[i] = [g.N_gap for g in raw_gaps]
    end

    df.plateaus = plateaus
    df.plateaus_idxd = plateaus_idxd
    return df
end


#####################################################################################

function compute_gap_labels!(df::DataFrame; p_range::Vector{Int}=[0,1])
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        energies = sort(real(row.eigenvalues))
        plateaus = row.plateaus_idxd     # assumed as Vector of (j, N_gap)
        phi = row.phi                    # your geometric parameter (interpretation by mode)

        label_list = NamedTuple[]
        for (j, N_gap) in plateaus
            E_low  = energies[j]
            E_high = energies[j+1]

            best_err = Inf
            best_p = 0
            best_q = 0

            factor = phi

            # typically p in {0,1} for normalized IDOS in [0,1]
            for p_try in p_range
                q_candidate = round(Int, (N_gap - p_try) / factor)
                val = p_try + q_candidate * factor

                err = abs(N_gap - val)
                if err < best_err
                    best_err = err
                    best_p = p_try
                    best_q = q_candidate
                end
            end

            push!(label_list, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=best_p, q=best_q, err=best_err))
        end
        gap_labels[i] = label_list
    end

    df.gap_labels = gap_labels
    return df
end

function compute_gap_labels_qlim!(df::DataFrame; p_range::Vector{Int}=[0,1], q_max::Int=5)
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        energies = sort(real(row.eigenvalues))
        plateaus = row.plateaus_idxd     # Vector of (j, N_gap)
        phi = row.phi                    # slope parameter

        label_list = NamedTuple[]
        for (j, N_gap) in plateaus
            E_low  = energies[j]
            E_high = energies[j+1]

            best_err = Inf
            best_p = 0
            best_q = 0

            factor = phi

            # Guard against phi ≈ 0 to avoid division by zero and Inf errors
            if !isfinite(phi) || abs(phi) < eps(Float64)
                push!(label_list, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=missing, q=missing, err=missing))
                continue
            end

            for p_try in p_range
                q_candidate = round(Int, (N_gap - p_try) / factor)

                # enforce q_max cutoff
                if abs(q_candidate) > q_max
                    continue
                end

                val = p_try + q_candidate * factor
                err = abs(N_gap - val)

                if err < best_err
                    best_err = err
                    best_p = p_try
                    best_q = q_candidate
                end
            end

            push!(label_list, (E_low=E_low, E_high=E_high,
                               N_gap=N_gap, p=best_p, q=best_q,
                               err=best_err))
        end

        gap_labels[i] = label_list
    end

    df.gap_labels = gap_labels
    return df
end

function compute_gap_labels_parsimonious!(
    df::DataFrame; 
    p_range::Vector{Int}=[0,1], 
    q_max::Int=5, 
    tol::Float64=1e-9, 
    merge_threshold::Union{Nothing, Float64}=nothing,
    check_mp::Bool=false,
    mp_col::Symbol=:mp,
    mp_threshold::Float64=0.1,
    always_mask_central::Bool=false # <--- New Argument
)
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))
    q_search_limit = 1000 
    
    # Sentinel value for SC gaps to be coloured Grey
    SC_GAP_Q = 9999

    for (i, row) in enumerate(eachrow(df))
        plateaus = row.plateaus_idxd     
        phi = row.phi                    

        # --- Determine System Size N ---
        if hasproperty(row, :N)
            N = row.N
        elseif hasproperty(row, :eigenvalues)
            N = length(row.eigenvalues)
        else
            error("Cannot determine system size N. Ensure 'N' or 'eigenvalues' column exists.")
        end

        # --- Check for MBS (Optional if always_mask_central is true) ---
        has_mbs = false
        if check_mp && hasproperty(row, mp_col) && !ismissing(row[mp_col])
            if abs(row[mp_col] + 1.0) < mp_threshold
                has_mbs = true
            end
        end
        
        # Decide whether to mask central gaps for this row
        should_mask_central = always_mask_central || has_mbs

        # --- Step A: Calculate Labels for ALL raw segments ---
        raw_labels = NamedTuple[]
        for gap in plateaus
            E_low, E_high, N_gap = gap.E_low, gap.E_high, gap.N_gap

            # --- SC GAP FILTER (INDEX BASED) ---
            if should_mask_central
                # Recover integer index j from IDOS N_gap
                # N_gap = j / N  =>  j = round(N_gap * N)
                j_gap = round(Int, N_gap * N)
                
                # Central indices for N states:
                # Gaps of interest: N/2 - 1, N/2, N/2 + 1
                center_idx = div(N, 2)
                
                if j_gap >= (center_idx - 1) && j_gap <= (center_idx + 1)
                    push!(raw_labels, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=0, q=SC_GAP_Q, err=0.0))
                    continue
                end
            end
            # ---------------------

            best_err = Inf
            best_p = missing
            best_q = missing

            factor = phi

            if isfinite(phi) && abs(phi) > eps(Float64)
                for p_try in p_range
                    q_ideal = (N_gap - p_try) / factor
                    q_candidate = round(Int, q_ideal)

                    if abs(q_candidate) > q_search_limit; continue; end

                    val = p_try + q_candidate * factor
                    err = abs(N_gap - val)

                    if err < best_err - tol
                        best_err = err
                        best_p = p_try
                        best_q = q_candidate
                    elseif err <= best_err + tol
                        if abs(q_candidate) < abs(best_q)
                            best_err = err
                            best_p = p_try
                            best_q = q_candidate
                        elseif abs(q_candidate) == abs(best_q) && abs(p_try) < abs(best_p)
                            best_err = err
                            best_p = p_try
                            best_q = q_candidate
                        end
                    end
                end
            end

            # Clamp q
            final_q = best_q
            if !ismissing(best_q) && abs(best_q) > q_max
                final_q = Int(sign(best_q) * q_max)
            end

            push!(raw_labels, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=best_p, q=final_q, err=best_err))
        end

        # --- Step B: Merge adjacent segments with SAME Label ---
        if merge_threshold !== nothing && !isempty(raw_labels)
            merged_labels = NamedTuple[]
            current_cluster = [raw_labels[1]]
            
            for k in 2:length(raw_labels)
                lbl = raw_labels[k]
                prev = current_cluster[end]
                
                diff_idos = lbl.N_gap - prev.N_gap
                
                same_label = !ismissing(lbl.p) && !ismissing(prev.p) && 
                             (lbl.p == prev.p) && (lbl.q == prev.q)

                if diff_idos <= merge_threshold && same_label
                    push!(current_cluster, lbl)
                else
                    push!(merged_labels, merge_labeled_cluster(current_cluster))
                    current_cluster = [lbl]
                end
            end
            push!(merged_labels, merge_labeled_cluster(current_cluster))
            gap_labels[i] = merged_labels
        else
            gap_labels[i] = raw_labels
        end
    end

    df.gap_labels = gap_labels
    return df
end


function merge_labeled_cluster(cluster)
    # E_low from first, E_high from last
    E_low = cluster[1].E_low
    E_high = cluster[end].E_high
    
    # Use properties of the "best" fit (smallest error) or widest gap?
    # Usually widest gap has best N_gap estimate.
    # But since labels are identical, we just keep the label.
    # We can pick the N_gap from the segment with lowest error.
    
    # argmin returns the ELEMENT that minimizes the function
    rep = argmin(x -> ismissing(x.err) ? Inf : x.err, cluster)
    
    return (E_low=E_low, E_high=E_high, N_gap=rep.N_gap, p=rep.p, q=rep.q, err=rep.err)
end

end # module IDOSProcessing