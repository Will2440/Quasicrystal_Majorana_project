using Statistics
using DataFrames
using Polynomials



############################################
######## DataFrame Preparation #############
############################################

function prep_df_for_IDOP(df::DataFrame; fixed_values...)
    """
    Collect mu-range vectors and corresponding discretised indicators per unique (phi, N, t_n, Delta).
    Assumes :disc_mp exists (scalar 0/1 per row) and :disc_eigenvalues exists (Vector{Int} per row).
    Returns DataFrame with columns:
      :phi, :N, :t_n, :Delta, :mu_values, :disc_mp, :disc_evals_first
    """
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    rows = Vector{Dict{Symbol,Any}}()

    for g in DataFrames.groupby(filtered_df, [:phi, :N, :t_n, :Delta, :phason])
        g_sorted = sort(g, :mu)

        mu_values = [r.mu for r in eachrow(g_sorted)]
        # assume disc_mp is scalar 0/1 (may be Float64 but integer-like)
        disc_mp_vals = [Int(r.disc_mp) for r in eachrow(g_sorted)]
        # assume disc_eigenvalues is a vector and we want its first element
        disc_eval_first_vals = [Int(r.disc_eigenvalues[1]) for r in eachrow(g_sorted)]

        first_row = first(g_sorted)

        push!(rows, Dict(
            :phi => first_row.phi,
            :N => first_row.N,
            :t_n => first_row.t_n,
            :Delta => first_row.Delta,
            :phason => first_row.phason,
            :mu_values => mu_values,
            :disc_mp => disc_mp_vals,
            :disc_evals => disc_eval_first_vals
        ))
    end

    return DataFrame(rows)
end

## original, based on indexing disc_mp to mu_values and stretching mu_values ot fit
# function compute_phase_norm!(df::DataFrame; mu_col::Symbol=:mu_values, disc_col::Symbol=:disc_mp, mu_norm_col::Symbol=:mu_values_norm, disc_norm_col::Symbol=:mp_disc_norm)
#     """
#     For each row, normalize :mu_values so that the mu at the first disc_mp=1 maps to 0,
#     and the mu at the last disc_mp=1 maps to 1. Add :mu_values_norm and :mp_disc_norm (unchanged disc_mp).
#     Assumes df has :mu_values (Vector) and :disc_mp (Vector) per row.
#     """
#     # if !(mu_col in names(df)) || !(disc_col in names(df))
#     #     error("Required columns $mu_col and $disc_col not found in DataFrame")
#     # end
    
#     mu_norm_vec = Vector{Vector{Float64}}(undef, nrow(df))
#     disc_norm_vec = Vector{Vector{Int}}(undef, nrow(df))
    
#     for (i, row) in enumerate(eachrow(df))
#         mu_vals = row[mu_col]
#         disc_vals = row[disc_col]
        
#         # Find indices where disc_mp == 1
#         idxs = findall(x -> x == 1, disc_vals)
#         if isempty(idxs)
#             # No 1s, perhaps set norm to original or handle
#             mu_norm_vec[i] = mu_vals  # or zeros, but let's assume there are
#             disc_norm_vec[i] = disc_vals
#             continue
#         end
        
#         idx_first = first(idxs)
#         idx_last = last(idxs)
        
#         mu_first = mu_vals[idx_first]
#         mu_last = mu_vals[idx_last]
        
#         if mu_last == mu_first
#             # Degenerate, set to 0.5 or something
#             mu_norm_vec[i] = fill(0.5, length(mu_vals))
#         else
#             mu_norm_vec[i] = (mu_vals .- mu_first) ./ (mu_last - mu_first)
#         end
        
#         disc_norm_vec[i] = disc_vals  # unchanged
#     end
    
#     df[!, mu_norm_col] = mu_norm_vec
#     df[!, disc_norm_col] = disc_norm_vec
    
#     return df
# end

## revised version creating :mu_mp and :mu_mp_norm, a list of mu values where disc_mp == 1 (else missing), and normalized version
function compute_phase_norm!(df::DataFrame; mu_col::Symbol=:mu_values, disc_col::Symbol=:disc_mp, mu_mp_col::Symbol=:mu_mp, mu_mp_norm_col::Symbol=:mu_mp_norm, disc_norm_col::Symbol=:mp_disc_norm)
    """
    For each row, create :mu_mp (mu values where disc_mp == 1, else missing) and :mu_mp_norm (normalized to [0,1] per row based on min/max of mu_mp).
    Also keeps :mp_disc_norm as the original disc_mp vector.
    Assumes df has :mu_values (Vector) and :disc_mp (Vector) per row.
    """
    # if !(mu_col in names(df)) || !(disc_col in names(df))
    #     error("Required columns $mu_col and $disc_col not found in DataFrame")
    # end
    
    mu_mp_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    mu_mp_norm_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    disc_norm_vec = Vector{Vector{Int}}(undef, nrow(df))
    
    for (i, row) in enumerate(eachrow(df))
        mu_vals = row[mu_col]
        disc_vals = row[disc_col]
        
        # Create mu_mp: mu where disc == 1, else missing
        mu_mp = [disc_vals[j] == 1 ? mu_vals[j] : missing for j in 1:length(mu_vals)]
        mu_mp_vec[i] = mu_mp
        
        # Normalize mu_mp to [0,1] based on its own min/max (non-missing values)
        non_missing_mu = filter(!ismissing, mu_mp)
        if !isempty(non_missing_mu)
            min_mu = minimum(non_missing_mu)
            max_mu = maximum(non_missing_mu)
            if max_mu > min_mu
                mu_mp_norm = [ismissing(m) ? missing : (m - min_mu) / (max_mu - min_mu) for m in mu_mp]
            else
                # Degenerate case: all mu_mp the same, set to 0.5
                mu_mp_norm = [ismissing(m) ? missing : 0.5 for m in mu_mp]
            end
        else
            # No matches: keep as missing
            mu_mp_norm = mu_mp
        end
        mu_mp_norm_vec[i] = mu_mp_norm
        
        # Keep original disc_mp
        disc_norm_vec[i] = disc_vals
    end
    
    df[!, mu_mp_col] = mu_mp_vec
    df[!, mu_mp_norm_col] = mu_mp_norm_vec
    df[!, disc_norm_col] = disc_norm_vec
    
    return df
end

## revised version normalising based on a line of best fit to the endge of mp_disc==1
function fit_final_mp_line(df::DataFrame; mu_col=:mu_values, disc_col=:disc_mp, phi_col=:phi, ignore_n_largest_phi::Int=0)
    phi_vals = Float64[]
    mu_final_mp = Float64[]
    for row in eachrow(df)
        mu_vec = row[mu_col]
        disc_vec = row[disc_col]
        idxs = findall(x -> x == 1, disc_vec)
        if !isempty(idxs)
            push!(phi_vals, row[phi_col])
            push!(mu_final_mp, mu_vec[last(idxs)])
        end
    end
    # Ignore the largest n phi values if requested
    if ignore_n_largest_phi > 0 && length(phi_vals) > ignore_n_largest_phi
        inds = sortperm(phi_vals)
        keep_inds = inds[1:end-ignore_n_largest_phi]
        phi_vals = phi_vals[keep_inds]
        mu_final_mp = mu_final_mp[keep_inds]
    end
    # Linear regression: μ = a*φ + b
    X = hcat(ones(length(phi_vals)), phi_vals)
    coeffs = X \ mu_final_mp
    b, a = coeffs[1], coeffs[2]
    return (a=a, b=b, phi_vals=phi_vals, mu_final_mp=mu_final_mp)
end

function compute_phase_norm_with_line!(
    df::DataFrame;
    mu_col::Symbol=:mu_values,
    disc_col::Symbol=:disc_mp,
    phi_col::Symbol=:phi,
    mu_mp_col::Symbol=:mu_mp,
    mu_mp_norm_col::Symbol=:mu_mp_norm,
    disc_norm_col::Symbol=:mp_disc_norm,
    line_coeffs=nothing
)
    if line_coeffs === nothing
        error("Must provide regression line coefficients as line_coeffs=(a=..., b=...)")
    end
    a, b = line_coeffs.a, line_coeffs.b

    mu_mp_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    mu_mp_norm_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    disc_norm_vec = Vector{Vector{Int}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        mu_vals = row[mu_col]
        disc_vals = row[disc_col]
        phi = row[phi_col]
        mu_mp = [disc_vals[j] == 1 ? mu_vals[j] : missing for j in 1:length(mu_vals)]
        nonmissing_mu = filter(!ismissing, mu_mp)
        mu_start = isempty(nonmissing_mu) ? missing : first(nonmissing_mu)
        mu_end = a * phi + b
        if ismissing(mu_start) || mu_end == mu_start
            mu_mp_norm = [ismissing(m) ? missing : 0.5 for m in mu_mp]
        else
            mu_mp_norm = [ismissing(m) ? missing : (m - mu_start) / (mu_end - mu_start) for m in mu_mp]
        end
        mu_mp_vec[i] = mu_mp
        mu_mp_norm_vec[i] = mu_mp_norm
        disc_norm_vec[i] = disc_vals
    end

    df[!, mu_mp_col] = mu_mp_vec
    df[!, mu_mp_norm_col] = mu_mp_norm_vec
    df[!, disc_norm_col] = disc_norm_vec
    return df
end

## rervised version based on polynomial of best fitfunction fit_final_mp_poly(df::DataFrame; mu_col=:mu_values, disc_col=:disc_mp, phi_col=:phi, order::Int=1)
function fit_final_mp_poly(df::DataFrame; mu_col=:mu_values, disc_col=:disc_mp, phi_col=:phi, order::Int=1, ignore_n_largest_phi::Int=0)
    phi_vals = Float64[]
    mu_final_mp = Float64[]
    for row in eachrow(df)
        mu_vec = row[mu_col]
        disc_vec = row[disc_col]
        idxs = findall(x -> x == 1, disc_vec)
        if !isempty(idxs)
            push!(phi_vals, row[phi_col])
            push!(mu_final_mp, mu_vec[last(idxs)])
        end
    end
    # Ignore the largest n phi values if requested
    if ignore_n_largest_phi > 0 && length(phi_vals) > ignore_n_largest_phi
        inds = sortperm(phi_vals)
        keep_inds = inds[1:end-ignore_n_largest_phi]
        phi_vals = phi_vals[keep_inds]
        mu_final_mp = mu_final_mp[keep_inds]
    end
    # Polynomial fit: μ = p(φ)
    p = fit(phi_vals, mu_final_mp, order)
    return (poly=p, phi_vals=phi_vals, mu_final_mp=mu_final_mp)
end

function compute_phase_norm_with_poly!(
    df::DataFrame;
    mu_col::Symbol=:mu_values,
    disc_col::Symbol=:disc_mp,
    phi_col::Symbol=:phi,
    mu_mp_col::Symbol=:mu_mp,
    mu_mp_norm_col::Symbol=:mu_mp_norm,
    disc_norm_col::Symbol=:mp_disc_norm,
    poly=nothing
)
    if poly === nothing
        error("Must provide polynomial as poly=...")
    end

    mu_mp_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    mu_mp_norm_vec = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df))
    disc_norm_vec = Vector{Vector{Int}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        mu_vals = row[mu_col]
        disc_vals = row[disc_col]
        phi = row[phi_col]
        mu_mp = [disc_vals[j] == 1 ? mu_vals[j] : missing for j in 1:length(mu_vals)]
        nonmissing_mu = filter(!ismissing, mu_mp)
        mu_start = isempty(nonmissing_mu) ? missing : first(nonmissing_mu)
        mu_end = poly(phi)
        if ismissing(mu_start) || mu_end == mu_start
            mu_mp_norm = [ismissing(m) ? missing : 0.5 for m in mu_mp]
        else
            mu_mp_norm = [ismissing(m) ? missing : (m - mu_start) / (mu_end - mu_start) for m in mu_mp]
        end
        mu_mp_vec[i] = mu_mp
        mu_mp_norm_vec[i] = mu_mp_norm
        disc_norm_vec[i] = disc_vals
    end

    df[!, mu_mp_col] = mu_mp_vec
    df[!, mu_mp_norm_col] = mu_mp_norm_vec
    df[!, disc_norm_col] = disc_norm_vec
    return df
end

############################################
########### IDOS Computation ###############
############################################

function compute_idop_df!(df::DataFrame; disc_variable::Symbol = :disc_mp)
    """
    Compute Integrated Density Of Phase-gaps (IDOP) from a discretised column.
    Normalises each cumulative count by the total number of positive entries so the final value = 1
    when there is at least one match. If there are zero matches the IDOP is a zero-vector.
    """
    if !(disc_variable in names(df) || string(disc_variable) in names(df))
        error("discretised column $(disc_variable) not found in DataFrame")
    end

    n = nrow(df)
    idop_col = Vector{Vector{Float64}}(undef, n)

    for (i, row) in enumerate(eachrow(df))
        vec = row[disc_variable]                  # expected AbstractVector{Int} or 0/1-like
        v = Int.(vec)                             # coerce elementwise

        total = sum(v)
        if total == 0
            # no matches -> IDOP stays zero (length preserved)
            idop_col[i] = zeros(Float64, length(v))
        else
            idop_col[i] = cumsum(v) ./ float(total)
        end
    end

    df.idop = idop_col
    return df
end

function find_idop_plateaus_from_disc(
    df::DataFrame;
    mu_gap_thresh::Float64 = 1e-2,
    disc_symbol::Symbol = :disc_mp,
    fixed_values...
)
    """
    Find IDOP plateau positions from a discretised indicator vector (e.g. :disc_mp).

    For each grouped row (after filtering by fixed_values), we:
      - take mu_values (Vector{Float64}) and disc_vec (Vector{Int} / 0-1)
      - find indices idxs where disc_vec == 1
      - if consecutive matched indices idxs[k], idxs[k+1] are separated by
        mu[idx_next] - mu[idx_cur] > mu_gap_thresh we treat that as a gap
        in mu and record the IDOP plateau at (idx_cur + 0.5) / L where L = length(mu_values).

    Returns Vector{Vector{Float64}} aligned with filtered rows.
    """
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plateau_idop = Vector{Vector{Float64}}(undef, nrow(filtered_df))
    for (i, row) in enumerate(eachrow(filtered_df))
        mu_vals = get(row, :mu_values, nothing)
        disc_vec = get(row, disc_symbol, nothing)

        if mu_vals === nothing || disc_vec === nothing || length(mu_vals) != length(disc_vec)
            plateau_idop[i] = Float64[]    # length mismatch / missing -> no plateaus
            continue
        end

        L = length(mu_vals)
        idxs = findall(x -> Int(round(x)) == 1, disc_vec)
        if isempty(idxs)
            plateau_idop[i] = Float64[]
            continue
        end

        plateaus = Float64[]
        # look for large mu-gaps between consecutive matched indices
        for k in 1:length(idxs)-1
            idx_cur = idxs[k]
            idx_next = idxs[k+1]
            if mu_vals[idx_next] - mu_vals[idx_cur] > mu_gap_thresh
                push!(plateaus, (idx_cur + 0.5) / L)
            end
        end

        # (optional) check gaps at edges if desired:
        # if mu_vals[first_matched] - mu_vals[1] > mu_gap_thresh -> push!((0 + 0.5)/L)
        # if mu_vals[end] - mu_vals[last_matched]  > mu_gap_thresh -> push!((last_matched + 0.5)/L)

        plateau_idop[i] = plateaus
    end

    return plateau_idop
end

function find_idop_plateaus_from_idop(
    df::DataFrame;
    delta_idop_thresh::Float64 = 1e-4,   # max allowed |ΔIDOP| between consecutive samples inside a flat region
    min_mu_span::Float64 = 0.0,          # minimal mu span (absolute) for a plateau; 0.0 disables mu-span check
    min_samples::Int = 3,                # minimal number of IDOP points in plateau (>=2 diffs -> points = diffs+1)
    fixed_values...
)
    """
    Detect IDOP plateaus by scanning idop (Vector{Float64}) vs mu_values (Vector{Float64}).
    Returns Vector{Vector{Float64}} aligned with filtered rows where each inner vector
    contains plateau y-values (mean IDOP over the flat region).

    Plateaus are contiguous runs where abs(diff(idop)) <= delta_idop_thresh.
    A run is reported if either:
      - its mu span >= min_mu_span (and min_mu_span > 0), OR
      - number of IDOP samples in the run >= min_samples.
    """
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plateau_idop = Vector{Vector{Float64}}(undef, nrow(filtered_df))

    for (i, row) in enumerate(eachrow(filtered_df))
        mu_vals = get(row, :mu_values, nothing)
        idop = get(row, :idop, nothing)

        if mu_vals === nothing || idop === nothing || length(mu_vals) != length(idop) || length(idop) < 2
            plateau_idop[i] = Float64[]   # not enough data or missing
            continue
        end

        L = length(idop)
        diffs = abs.(diff(idop))  # length L-1
        # boolean mask where consecutive points are "flat"
        flat_mask = diffs .<= delta_idop_thresh

        plateaus = Float64[]
        k = 1
        while k <= length(flat_mask)
            if !flat_mask[k]
                k += 1
                continue
            end
            # start of run at diff index `k`
            run_start = k
            run_end = k
            while run_end + 1 <= length(flat_mask) && flat_mask[run_end + 1]
                run_end += 1
            end
            # run covers diffs indices run_start:run_end -> corresponds to idop indices run_start:(run_end+1)
            idx_a = run_start
            idx_b = run_end + 1  # inclusive idop index
            n_points = idx_b - idx_a + 1
            mu_span = mu_vals[idx_b] - mu_vals[idx_a]   # span between first and last mu in plateau

            if (min_mu_span > 0.0 && mu_span >= min_mu_span) || (n_points >= min_samples)
                # report plateau y-value as mean (robust to small residual noise)
                push!(plateaus, mean(idop[idx_a:idx_b]))
            end

            k = run_end + 1
        end

        plateau_idop[i] = plateaus
    end

    return plateau_idop
end

function compute_idop_plateaus_all!(
    df::DataFrame;
    delta_idop_thresh::Float64 = 1e-4,
    min_mu_span::Float64 = 0.0,
    min_samples::Int = 3
)::DataFrame
    plateaus = find_idop_plateaus_from_idop(
        df;
        delta_idop_thresh = delta_idop_thresh,
        min_mu_span = min_mu_span,
        min_samples = min_samples
    )
    @assert length(plateaus) == nrow(df) "find_idop_plateaus_from_idop returned mismatched length"
    df.plateaus = plateaus
    return df
end

function filter_idop_plateaus_at_1!(df::DataFrame; atol::Float64=1e-6)
    """
    Filter out plateaus in :plateaus column that are approximately equal to 1.0.
    Assumes :plateaus is a Vector{Vector{Float64}} where each inner vector contains plateau IDOP values.
    Modifies df in-place.
    """
    # if !(:plateaus in names(df))
    #     error(":plateaus column not found in DataFrame")
    # end
    
    for row in eachrow(df)
        plateaus = row.plateaus
        # Filter out values close to 1.0
        filtered_plateaus = filter(p -> !isapprox(p, 1.0; atol=atol), plateaus)
        row.plateaus = filtered_plateaus
    end
    
    return df
end

############################################
############# Gap Labelling ################
############################################

function compute_gap_labels_qlim!(df::DataFrame; p_range::Vector{Int}=collect(-5:5), q_max::Int=5)
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))
    
    for (i, row) in enumerate(eachrow(df))
        plateaus = row.plateaus  # Vector of Float64 plateau values from compute_idop_plateaus_all!
        phi = row.phi            # Float64 parameter for rational approximation
        
        label_list = NamedTuple[]
        for plateau in plateaus
            adjusted_plateau = plateau + 0.5
            best_err = Inf
            best_p = 0
            best_q = 0
            
            # Search for best p/q such that plateau ≈ p + q * phi
            for p_try in p_range
                for q_try in -q_max:q_max
                    if q_try == 0
                        val = p_try  # Avoid division by zero
                    else
                        val = p_try + q_try * phi
                    end
                    err = abs(adjusted_plateau - val)
                    if err < best_err
                        best_err = err
                        best_p = p_try
                        best_q = q_try
                    end
                end
            end
            
            push!(label_list, (plateau=adjusted_plateau, p=best_p, q=best_q, err=best_err))
        end
        gap_labels[i] = label_list
    end
    
    df[!, :gap_labels] .= gap_labels
    return df
end