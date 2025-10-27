using Statistics

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

    for g in DataFrames.groupby(filtered_df, [:phi, :N, :t_n, :Delta])
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
            :mu_values => mu_values,
            :disc_mp => disc_mp_vals,
            :disc_evals => disc_eval_first_vals
        ))
    end

    return DataFrame(rows)
end


# function compute_idop_df!(df::DataFrame; disc_variable::Symbol = :disc_mp)
#     """
#     Compute Integrated Density Of Phase-gaps (IDOP) from a discretised column.
#     Assumes `disc_variable` is present in `df` and contains Vector{Int} for each row.
#     """
#     if !(disc_variable in names(df) || string(disc_variable) in names(df))
#         error("discretised column $(disc_variable) not found in DataFrame")
#     end

#     n = nrow(df)
#     idop_col = Vector{Vector{Float64}}(undef, n)

#     for (i, row) in enumerate(eachrow(df))
#         vec = row[disc_variable]            # expect AbstractVector{Int}
#         v = Int.(vec)                       # coerce elementwise to Int
#         idop_col[i] = cumsum(v) ./ length(v)
#     end

#     df.idop = idop_col
#     return df
# end

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

# function compute_idop_plateaus_all!(
#     df::DataFrame;
#     delta_idop_thresh::Float64 = 1e-4,
#     min_mu_span::Float64 = 0.0,
#     min_samples::Int = 3
# )::DataFrame
#     """
#     Compute IDOP plateaus for every row in `df` and store result in column :plateaus.
#     Uses the same plateau-detection logic as `find_idop_plateaus_from_idop` but operates
#     on the full dataframe (no filtering). Returns the modified dataframe.
#     """
#     n = nrow(df)
#     plateaus_col = Vector{Vector{Float64}}(undef, n)

#     for (i, row) in enumerate(eachrow(df))
#         mu_vals = get(row, :mu_values, nothing)
#         idop = get(row, :idop, nothing)

#         if mu_vals === nothing || idop === nothing || length(mu_vals) != length(idop) || length(idop) < 2
#             plateaus_col[i] = Float64[]
#             continue
#         end

#         diffs = abs.(diff(idop))
#         flat_mask = diffs .<= delta_idop_thresh

#         plateaus = Float64[]
#         k = 1
#         while k <= length(flat_mask)
#             if !flat_mask[k]
#                 k += 1
#                 continue
#             end
#             run_start = k
#             run_end = k
#             while run_end + 1 <= length(flat_mask) && flat_mask[run_end + 1]
#                 run_end += 1
#             end

#             # idop indices covered: run_start : (run_end+1)
#             idx_a = run_start
#             idx_b = run_end + 1
#             n_points = idx_b - idx_a + 1
#             mu_span = mu_vals[idx_b] - mu_vals[idx_a]

#             if (min_mu_span > 0.0 && mu_span >= min_mu_span) || (n_points >= min_samples)
#                 push!(plateaus, mean(idop[idx_a:idx_b]))
#             end

#             k = run_end + 1
#         end

#         plateaus_col[i] = plateaus
#     end

#     df.plateaus = plateaus_col
#     return df
# end


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