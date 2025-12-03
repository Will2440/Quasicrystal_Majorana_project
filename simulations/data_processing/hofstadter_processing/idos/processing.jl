module IDOSProcessing

using DataFrames
using ProgressMeter

############################################
######## Eigenvalue Normalisations #########
############################################

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

end # module IDOSProcessing