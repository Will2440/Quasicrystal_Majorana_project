using DataFrames
using ProgressMeter

# function prep_df_for_IDOS(df::DataFrame; fixed_values...)
#     """
#     Collect mu-range vectors and corresponding discretised indicators per unique (phi, N, t_n, Delta).
#     Assumes :disc_mp exists (scalar 0/1 per row) and :disc_eigenvalues exists (Vector{Int} per row).
#     Returns DataFrame with columns:
#       :phi, :N, :t_n, :Delta, :mu_values, :disc_mp, :disc_evals_first
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     rows = Vector{Dict{Symbol,Any}}()

#     for g in DataFrames.groupby(filtered_df, [:phi, :N, :t_n, :Delta])
#         g_sorted = sort(g, :mu)

#         mu_values = [r.mu for r in eachrow(g_sorted)]
#         # assume disc_mp is scalar 0/1 (may be Float64 but integer-like)
#         disc_mp_vals = [Int(r.disc_mp) for r in eachrow(g_sorted)]
#         # assume disc_eigenvalues is a vector and we want its first element
#         disc_eval_first_vals = [Int(r.disc_eigenvalues[1]) for r in eachrow(g_sorted)]

#         eigenvalues_per_mu = [r.eigenvalues for r in eachrow(g_sorted)]

#         first_row = first(g_sorted)

#         push!(rows, Dict(
#             :phi => first_row.phi,
#             :N => first_row.N,
#             :t_n => first_row.t_n,
#             :Delta => first_row.Delta,
#             :mu_values => mu_values,
#             :disc_mp => disc_mp_vals,
#             :disc_evals => disc_eval_first_vals,
#             :eigenvalues => eigenvalues_per_mu
#         ))
#     end

#     return DataFrame(rows)
# end

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

## for computing the areas in energy space coresponding to the idos plateaus (gap labelled).
# function compute_gap_intervals_df(df::DataFrame; fixed_values...)
#     """
#     Build a DataFrame with one row per identified gap (from df.gap_labels).
#     Returned DataFrame columns:
#       :phi, :N, :t_n, :Delta, :mu, :phason,
#       :E_low, :E_high, :N_gap, :p, :q

#     Use this output to draw coloured gap bands in energy space.
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     rows = Vector{Dict{Symbol,Any}}()
#     for row in eachrow(filtered_df)
#         # Expect row.gap_labels to be a Vector{NamedTuple} with fields E_low, E_high, N_gap, p, q
#         if !haskey(row, :gap_labels) || isempty(row.gap_labels)
#             continue
#         end
#         for gap in row.gap_labels
#             push!(rows, Dict(
#                 :phi   => row.phi,
#                 :N     => row.N,
#                 :t_n   => row.t_n,
#                 :Delta => row.Delta,
#                 :mu    => row.mu,
#                 :phason=> row.phason,
#                 :E_low => gap.E_low,
#                 :E_high=> gap.E_high,
#                 :N_gap => gap.N_gap,
#                 :p     => gap.p,
#                 :q     => gap.q
#             ))
#         end
#     end

#     return DataFrame(rows)
# end

function flatten_gap_intervals(df::DataFrame)
    # expects df has :gap_labels and a :phi column (and optionally :mu, :Delta, :phason)
    rows = NamedTuple[]
    has_mu = :mu in names(df); has_Delta = :Delta in names(df); has_phason = :phason in names(df)
    for r in eachrow(df)
        gls = r.gap_labels
        (gls === nothing || isempty(gls)) && continue
        for g in gls
            push!(rows, (
                E_low   = Float64(g.E_low),
                E_high  = Float64(g.E_high),
                phi     = Float64(r.phi),
                p       = Int(g.p),
                q       = Int(g.q),
                N_gap   = Float64(g.N_gap),
                mu      = has_mu ? Float64(r.mu) : NaN,
                Delta   = has_Delta ? Float64(r.Delta) : NaN,
                phason  = has_phason ? Float64(r.phason) : NaN,
            ))
        end
    end
    return DataFrame(rows)
end