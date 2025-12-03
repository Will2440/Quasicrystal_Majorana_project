module IDOPProcessing

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

function prep_df_for_IDOP_single(df_slice::DataFrame, phi_val::Float64, phason_val::Float64, N_val::Int, Delta_val::Float64, t_n::Vector{Float64})
    """
    Prepare a single-row DataFrame for IDOP computation from a slice DataFrame.
    Assumes df_slice has :mu, :disc_mp, :eigenvalues columns.
    Returns DataFrame with columns: :phi, :N, :t_n, :Delta, :phason, :mu_values, :disc_mp, :disc_evals
    """
    mu_values = df_slice.mu
    disc_mp_vals = Int.(df_slice.disc_mp)  # Ensure Int

    row = Dict(
        :phi => phi_val,
        :N => N_val,
        :t_n => t_n,
        :Delta => Delta_val,
        :phason => phason_val,
        :mu_values => mu_values,
        :disc_mp => disc_mp_vals,
    )

    return DataFrame([row])
end

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

function compute_phase_norm_with_poly_inplace!(
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

    for row in eachrow(df)
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
        row[mu_mp_col] = mu_mp
        row[mu_mp_norm_col] = mu_mp_norm
        row[disc_norm_col] = disc_vals
    end
    return df
end

############################################
########### IDOP Computation ###############
############################################

# function compute_idop_df!(df::DataFrame; disc_variable::Symbol = :disc_mp)
#     """
#     Compute Integrated Density Of Phase-gaps (IDOP) from a discretised column.
#     Normalises each cumulative count by the total number of positive entries so the final value = 1
#     when there is at least one match. If there are zero matches the IDOP is a zero-vector.
#     """
#     if !(disc_variable in names(df) || string(disc_variable) in names(df))
#         error("discretised column $(disc_variable) not found in DataFrame")
#     end

#     n = nrow(df)
#     idop_col = Vector{Vector{Float64}}(undef, n)

#     for (i, row) in enumerate(eachrow(df))
#         vec = row[disc_variable]                  # expected AbstractVector{Int} or 0/1-like
#         v = Int.(vec)                             # coerce elementwise

#         total = sum(v)
#         if total == 0
#             # no matches -> IDOP stays zero (length preserved)
#             idop_col[i] = zeros(Float64, length(v))
#         else
#             idop_col[i] = cumsum(v) ./ float(total)
#         end

#         idop_col[i] = 0.5 .+ 0.5 .* idop_col[i]   # rescale to [0.5, 1.0]
#     end

#     df.idop = idop_col
#     return df
# end


function compute_idop_df!(
    df::DataFrame;
    disc_variable::Symbol = :disc_mp,
    rescale_idop::Bool = true
)
    """
    Compute Integrated Density Of Phase-gaps (IDOP) from a discretised column.
    For binary (Int/Bool) disc_variable: Normalises each cumulative count by the total number of positive entries so the final value = 1
    when there is at least one match. If there are zero matches the IDOP is a zero-vector. Gaps are identified as 0s.
    For Float64/Missing disc_variable (e.g. mu_mp): IDOP is a cumulative sum over non-missing values, holding constant (plateau) over missing values.
    Normalizes by the total number of non-missing entries, then rescales to [0.5, 1.0] if rescale_idop=true.
    """
    if !(disc_variable in names(df) || string(disc_variable) in names(df))
        error("discretised column $(disc_variable) not found in DataFrame")
    end

    n = nrow(df)
    idop_col = Vector{Vector{Float64}}(undef, n)

    for (i, row) in enumerate(eachrow(df))
        vec = row[disc_variable]

        if eltype(vec) <: Union{Int, Bool}
            # Binary case: gaps as 0s, accumulate 1s
            v = Int.(vec)
            total = sum(v)
            if total == 0
                idop = zeros(Float64, length(v))
            else
                idop = cumsum(v) ./ float(total)
            end
            if rescale_idop
                idop_col[i] = 0.5 .+ 0.5 .* idop   # rescale to [0.5, 1.0]
            else
                idop_col[i] = idop
            end
        elseif eltype(vec) <: Union{Float64, Missing} || eltype(vec) == Any
            # mu_mp case: cumulative sum over non-missing, plateau over missing
            is_valid = .!ismissing.(vec)
            total = count(is_valid)
            idop = zeros(Float64, length(vec))
            acc = 0
            if total == 0
                # No non-missing values, IDOP stays zero
                idop .= 0.0
            else
                for j in 1:length(vec)
                    if is_valid[j]
                        acc += 1
                    end
                    idop[j] = acc / total
                end
            end
            if rescale_idop
                idop_col[i] = 0.5 .+ 0.5 .* idop
            else
                idop_col[i] = idop
            end
        else
            error("Unsupported element type for disc_variable: $(eltype(vec))")
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


#############################################
# function find_idop_plateaus_from_idop(
#     df::DataFrame;
#     delta_idop_thresh::Float64 = 1e-4,   # max allowed |ΔIDOP| between consecutive samples inside a flat region
#     min_mu_span::Float64 = 0.0,          # minimal mu span (absolute) for a plateau; 0.0 disables mu-span check
#     min_samples::Int = 3,                # minimal number of IDOP points in plateau (>=2 diffs -> points = diffs+1)
#     fixed_values...
# )
#     """
#     Detect IDOP plateaus by scanning idop (Vector{Float64}) vs mu_values (Vector{Float64}).
#     Returns Vector{Vector{Float64}} aligned with filtered rows where each inner vector
#     contains plateau y-values (mean IDOP over the flat region).

#     Plateaus are contiguous runs where abs(diff(idop)) <= delta_idop_thresh.
#     A run is reported if either:
#       - its mu span >= min_mu_span (and min_mu_span > 0), OR
#       - number of IDOP samples in the run >= min_samples.
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     plateau_idop = Vector{Vector{Float64}}(undef, nrow(filtered_df))

#     for (i, row) in enumerate(eachrow(filtered_df))
#         mu_vals = get(row, :mu_values, nothing)
#         idop = get(row, :idop, nothing)

#         if mu_vals === nothing || idop === nothing || length(mu_vals) != length(idop) || length(idop) < 2
#             plateau_idop[i] = Float64[]   # not enough data or missing
#             continue
#         end

#         L = length(idop)
#         diffs = abs.(diff(idop))  # length L-1
#         # boolean mask where consecutive points are "flat"
#         flat_mask = diffs .<= delta_idop_thresh

#         plateaus = Float64[]
#         k = 1
#         while k <= length(flat_mask)
#             if !flat_mask[k]
#                 k += 1
#                 continue
#             end
#             # start of run at diff index `k`
#             run_start = k
#             run_end = k
#             while run_end + 1 <= length(flat_mask) && flat_mask[run_end + 1]
#                 run_end += 1
#             end
#             # run covers diffs indices run_start:run_end -> corresponds to idop indices run_start:(run_end+1)
#             idx_a = run_start
#             idx_b = run_end + 1  # inclusive idop index
#             n_points = idx_b - idx_a + 1
#             mu_span = mu_vals[idx_b] - mu_vals[idx_a]   # span between first and last mu in plateau

#             if (min_mu_span > 0.0 && mu_span >= min_mu_span) || (n_points >= min_samples)
#                 # report plateau y-value as mean (robust to small residual noise)
#                 push!(plateaus, mean(idop[idx_a:idx_b]))
#             end

#             k = run_end + 1
#         end

#         plateau_idop[i] = plateaus
#     end

#     return plateau_idop
# end

# function compute_idop_plateaus_all!(
#     df::DataFrame;
#     delta_idop_thresh::Float64 = 1e-4,
#     min_mu_span::Float64 = 0.0,
#     min_samples::Int = 3
# )::DataFrame
#     plateaus = find_idop_plateaus_from_idop(
#         df;
#         delta_idop_thresh = delta_idop_thresh,
#         min_mu_span = min_mu_span,
#         min_samples = min_samples
#     )
#     @assert length(plateaus) == nrow(df) "find_idop_plateaus_from_idop returned mismatched length"
#     df.plateaus = plateaus
#     return df
# end

# function filter_idop_plateaus_at_1!(df::DataFrame; atol::Float64=1e-6)
#     """
#     Filter out plateaus in :plateaus column that are approximately equal to 1.0.
#     Assumes :plateaus is a Vector{Vector{Float64}} where each inner vector contains plateau IDOP values.
#     Modifies df in-place.
#     """
#     # if !(:plateaus in names(df))
#     #     error(":plateaus column not found in DataFrame")
#     # end
    
#     for row in eachrow(df)
#         plateaus = row.plateaus
#         # Filter out values close to 1.0
#         filtered_plateaus = filter(p -> !isapprox(p, 1.0; atol=atol), plateaus)
#         row.plateaus = filtered_plateaus
#     end
    
#     return df
# end
##############################################
function find_idop_plateaus_from_idop(
    df::DataFrame;
    delta_idop_thresh::Float64 = 1e-4,   # max allowed |ΔIDOP| between consecutive samples inside a flat region
    min_mu_span::Float64 = 0.0,          # minimal mu span (absolute) for a plateau; 0.0 disables mu-span check
    min_samples::Int = 3,                # minimal number of IDOP points in plateau (>=2 diffs -> points = diffs+1)
    fixed_values...
)
    """
    Detect IDOP plateaus by scanning idop (Vector{Float64}) vs mu_values (Vector{Float64}).
    Returns Vector{Vector{NamedTuple}} aligned with filtered rows where each inner vector
    contains plateau info: (plateau=mean IDOP, mu_low=min mu, mu_high=max mu) over the flat region.

    Plateaus are contiguous runs where abs(diff(idop)) <= delta_idop_thresh.
    A run is reported if either:
      - its mu span >= min_mu_span (and min_mu_span > 0), OR
      - number of IDOP samples in the run >= min_samples.
    """
    filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plateau_idop = Vector{Vector{NamedTuple}}(undef, nrow(filtered_df))

    for (i, row) in enumerate(eachrow(filtered_df))
        mu_vals = get(row, :mu_values, nothing)
        idop = get(row, :idop, nothing)

        if mu_vals === nothing || idop === nothing || length(mu_vals) != length(idop) || length(idop) < 2
            plateau_idop[i] = NamedTuple[]   # not enough data or missing
            continue
        end

        L = length(idop)
        diffs = abs.(diff(idop))  # length L-1
        # boolean mask where consecutive points are "flat"
        flat_mask = diffs .<= delta_idop_thresh

        plateaus = NamedTuple[]
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
                # report plateau info
                push!(plateaus, (plateau=mean(idop[idx_a:idx_b]), mu_low=mu_vals[idx_a], mu_high=mu_vals[idx_b]))
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
    for row in eachrow(df)
        plateaus = row.plateaus
        # Filter out values close to 1.0
        filtered_plateaus = filter(p -> !isapprox(p.plateau, 1.0; atol=atol), plateaus)
        row.plateaus = filtered_plateaus
    end
    return df
end

function filter_idop_plateaus_at_0!(df::DataFrame; atol::Float64=1e-6)
    for row in eachrow(df)
        plateaus = row.plateaus
        # Filter out values close to 0.0
        filtered_plateaus = filter(p -> !isapprox(p.plateau, 0.0; atol=atol), plateaus)
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
                    err = abs(plateau - val)
                    if err < best_err
                        best_err = err
                        best_p = p_try
                        best_q = q_try
                    end
                end
            end
            
            push!(label_list, (plateau=plateau, p=best_p, q=best_q, err=best_err))
        end
        gap_labels[i] = label_list
    end
    
    df[!, :gap_labels] .= gap_labels
    return df
end

function compute_gap_labels_qlim_new!(df::DataFrame; p_range::Vector{Int}=[0,1], q_max::Int=5)
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        plateaus = row.plateaus  # Vector of NamedTuple (plateau, mu_low, mu_high) from compute_idop_plateaus_all!
        phi = row.phi            # Float64 parameter for rational approximation

        label_list = NamedTuple[]
        for gap in plateaus
            plateau = gap.plateau
            mu_low = gap.mu_low
            mu_high = gap.mu_high
            best_err = Inf
            best_p = 0
            best_q = 0

            factor = phi

            # Guard against phi ≈ 0 to avoid division by zero and Inf errors
            if !isfinite(phi) || abs(phi) < eps(Float64)
                push!(label_list, (E_low=mu_low, E_high=mu_high, N_gap=plateau, p=missing, q=missing, err=missing))
                continue
            end

            for p_try in p_range
                q_candidate = round(Int, (plateau - p_try) / factor)

                # enforce q_max cutoff
                if abs(q_candidate) > q_max
                    continue
                end

                val = p_try + q_candidate * factor
                err = abs(plateau - val)

                if err < best_err
                    best_err = err
                    best_p = p_try
                    best_q = q_candidate
                end
            end

            push!(label_list, (E_low=mu_low, E_high=mu_high,
                               N_gap=plateau, p=best_p, q=best_q,
                               err=best_err))
        end

        gap_labels[i] = label_list
    end

    df.gap_labels = gap_labels
    return df
end





############################################
########### Phason Processing ##############
############################################

## currently not working as expected
function process_gap_eigenvalues(df_slice::DataFrame)
    """
    Process df_slice to compute gap eigenvalues (closest to zero in gaps where disc_mp == 0).
    Adds :gap_pos_eigs and :gap_neg_eigs columns.
    Assumes df_slice has :min_pos_eigs, :min_neg_eigs, and :mu_mp (to detect gaps).
    """
    # if !(:min_pos_eigs in names(df_slice)) || !(:min_neg_eigs in names(df_slice)) || !(:mu_mp in names(df_slice))
    #     @warn "df_slice missing required columns. Skipping gap processing."
    #     return df_slice
    # end

    println("DataFrame keynames: ", names(df_slice))

    gap_pos_eigs_all = []
    gap_neg_eigs_all = []

    for row in eachrow(df_slice)
        min_pos_eigs = getproperty(row, :min_pos_eigs)
        min_neg_eigs = getproperty(row, :min_neg_eigs)
        mu_mp_vec = getproperty(row, :mu_mp)
        mu_range = getproperty(row, :mu_range)

        gap_pos_eigs = []
        gap_neg_eigs = []

        for (i, mu) in enumerate(mu_range)
            is_gap = !(mu in mu_mp_vec)
            if is_gap
                # Closest to zero is the one with smaller abs between min_pos and min_neg
                pos_abs = abs(min_pos_eigs[i])
                neg_abs = abs(min_neg_eigs[i])
                if pos_abs < neg_abs
                    push!(gap_pos_eigs, min_pos_eigs[i])
                    push!(gap_neg_eigs, missing)
                else
                    push!(gap_pos_eigs, missing)
                    push!(gap_neg_eigs, min_neg_eigs[i])
                end
            else
                push!(gap_pos_eigs, missing)
                push!(gap_neg_eigs, missing)
            end
        end

        push!(gap_pos_eigs_all, gap_pos_eigs)
        push!(gap_neg_eigs_all, gap_neg_eigs)
    end

    df_slice[!, :gap_pos_eigs] = gap_pos_eigs_all
    df_slice[!, :gap_neg_eigs] = gap_neg_eigs_all

    return df_slice
end


end # module IDOPProcessing