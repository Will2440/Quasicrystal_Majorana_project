module IDOPProcessing

using Statistics
using DataFrames
using Polynomials
using LinearAlgebra
using FFTW
using CSV
using BSON



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

function normalise_scatter!(
    df::DataFrame,
    unormalised_column::Symbol,
    normalised_column::Symbol;
    norm_range::Tuple{Float64, Float64}=(-1.0, 1.0)
)
    """
    Takes in df and normalises the data in unormalised_column PER ROW.
    The normalised data is stored in normalised_column -- a new column created in the df.
    The normalisation rescales the data of THAT ROW to fit within norm_range (default: (-1, 1)).
    """
    range_min, range_max = norm_range
    range_span = range_max - range_min

    # Apply normalization row-wise
    df[!, normalised_column] = map(df[!, unormalised_column]) do vec
        if ismissing(vec)
            return missing
        end
        
        # Find local min/max for this row, ignoring missings
        valid_vals = collect(skipmissing(vec))
        
        if isempty(valid_vals)
             return vec 
        end

        min_val = minimum(valid_vals)
        max_val = maximum(valid_vals)
        data_span = max_val - min_val

        # Normalize elements
        map(vec) do x
            if ismissing(x)
                return missing
            end
            
            if data_span == 0.0
                 return range_min + range_span / 2.0
            end
            
            return range_min + ((x - min_val) / data_span) * range_span
        end
    end

    return df
end

function apply_smooth_normalisation!(
    df::DataFrame, 
    coeffs_file::String; 
    target_col::Symbol=:mu_mp, 
    out_col::Symbol=:mu_mp_norm, 
    phi_col::Symbol=:phi
)
    
    if !isfile(coeffs_file)
        error("Normalization coefficients file not found: $coeffs_file")
    end

    # 1. Read Coefficients
    coeffs = Dict{Int, Float64}() # coeff_index => value
    poly_deg = 0
    
    open(coeffs_file, "r") do io
        header = readline(io) # skip header
        for line in eachline(io)
            parts = split(line, ",")
            if length(parts) >= 3
                deg = parse(Int, parts[1])
                idx = parse(Int, parts[2])
                val = parse(Float64, parts[3])
                coeffs[idx] = val
                poly_deg = max(poly_deg, deg)
            end
        end
    end
    
    if isempty(coeffs)
        error("No coefficients found in file: $coeffs_file")
    end

    # 2. Apply Normalization
    new_col = Vector{Vector{Union{Float64, Missing}}}(undef, nrow(df))
    
    for (i, row) in enumerate(eachrow(df))
        vals = row[target_col]
        # Handle cases where phi might be a fixed parameter not in the row
        # We assume for now it is in the row or we fail.
        if !hasproperty(row, phi_col)
             error("Column $phi_col not found in DataFrame for normalization.")
        end
        phi = row[phi_col] 
        
        if ismissing(vals) || isempty(vals)
            new_col[i] = Float64[]
            continue
        end

        # Evaluate polynomial
        norm_val = 0.0
        for p in 0:poly_deg
            c = get(coeffs, p, 0.0)
            norm_val += c * phi^p
        end

        # Normalize (divide by envelope)
        factor = norm_val > 1e-9 ? 1.0 / norm_val : 1.0
        # new_col[i] = vals .* factor
        new_col[i] = [ismissing(v) ? missing : Float64(v) * factor for v in vals]
    end
    
    df[!, out_col] = new_col
    return df
end

function apply_outer_normalisation!(
    df::DataFrame, 
    norm_filename::String; 
    target_col::Symbol=:mu_mp, 
    out_col::Symbol=:mu_mp_norm
)
    """
        apply_outer_normalisation!(df::DataFrame, norm_filename::String; target_col::Symbol=:mu_mp, out_col::Symbol=:mu_mp_norm)

        Applies the stored global min/max normalization (from eigenvalue bandwidth normalization) to a target column in the DataFrame.

        # Arguments
        - `df`: DataFrame to modify in-place.
        - `norm_filename`: Path to the CSV file containing normalization parameters (:global_min, :global_max, and matching parameter columns like :phi, :phason, etc.).
        - `target_col`: Column in `df` containing vectors of values to normalize (default :mu_mp).
        - `out_col`: New column in `df` to store the normalized vectors (default :mu_mp_norm).

        # Behavior
        - Loads the normalization DataFrame from `norm_filename`.
        - For each row in `df`, matches based on parameter columns (:phi, :phason, :mu, :Delta, :t_n, :N).
        - Applies linear mapping: 2 * (x - global_min) / (global_max - global_min) - 1 to each value in `target_col`.
        - If global_max == global_min, maps all values to 0.0.
        - Skips rows with no match or missing/invalid data, logging warnings.
    """
    # Load the normalization parameters from CSV
    norm_df = CSV.read(norm_filename, DataFrame)
    
    # Define the parameter columns used for matching (adjust if your DataFrames have different columns)
    param_cols = [:phi, :phason]
    
    # Ensure required columns exist
    required_cols = [param_cols..., :global_min, :global_max]
    for col in required_cols
        if !(col in propertynames(norm_df))
            error("Normalization CSV missing required column: $col")
        end
    end
    
    # Initialize the output column (allow Missing in vector elements)
    df[!, out_col] = Vector{Union{Vector{Union{Missing, Float64}}, Missing}}(undef, nrow(df))

    for (i, row) in enumerate(eachrow(df))
        # Extract parameter values for matching
        param_values = NamedTuple(param_cols .=> getproperty.(Ref(row), param_cols))
        
        # Find matching row in norm_df
        matching_rows = filter(r -> all(getproperty(r, col) == param_values[col] for col in param_cols), norm_df)
        
        if nrow(matching_rows) == 0
            @warn "No matching normalization parameters found for row $i with params $param_values. Skipping normalization."
            df[i, out_col] = missing
            continue
        elseif nrow(matching_rows) > 1
            @warn "Multiple matching normalization parameters found for row $i with params $param_values. Using the first match."
        end
        
        # Get min/max from the first match
        global_min = matching_rows[1, :global_min]
        global_max = matching_rows[1, :global_max]
        
        # Check for missing or invalid min/max
        if ismissing(global_min) || ismissing(global_max) || !isfinite(global_min) || !isfinite(global_max)
            @warn "Invalid global_min ($global_min) or global_max ($global_max) for row $i. Skipping normalization."
            df[i, out_col] = missing
            continue
        end
        
        # Get the target values (vector, may contain Missing)
        target_vals = row[target_col]
        if ismissing(target_vals) || isempty(target_vals)
            df[i, out_col] = missing
            continue
        end

        # Apply normalization (handle Missing in elements)
        if global_max == global_min
            # All values map to 0, but preserve Missing
            normed = [ismissing(x) ? missing : 0.0 for x in target_vals]
        else
            normed = [ismissing(x) ? missing : 2 * (x - global_min) / (global_max - global_min) - 1 for x in target_vals]
        end
        
        df[i, out_col] = normed
    end
    
    return df
end

############################################
########### IDOP Computation ###############
############################################

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

function compute_gap_labels_parsimonious!(
    df::DataFrame; 
    p_range::Vector{Int}=[0,1], 
    q_max::Int=5, 
    tol::Float64=1e-9, 
    merge_threshold::Union{Nothing, Float64}=nothing,
    check_mp::Bool=false,
    mp_col::Symbol=:mp,
    mp_threshold::Float64=0.1,
    always_mask_central::Bool=false
)
    gap_labels = Vector{Vector{NamedTuple}}(undef, nrow(df))
    q_search_limit = 1000 
    
    # Sentinel value for SC gaps to be coloured Grey
    SC_GAP_Q = 9999

    for (i, row) in enumerate(eachrow(df))
        plateaus = row.plateaus  # Vector of NamedTuple (plateau, mu_low, mu_high)
        phi = row.phi

        # --- Determine System Size N ---
        N = hasproperty(row, :N) ? row.N : nothing
        
        # --- Check for MBS (Optional) ---
        has_mbs = false
        if check_mp && hasproperty(row, mp_col) && !ismissing(row[mp_col])
            if abs(row[mp_col] + 1.0) < mp_threshold
                has_mbs = true
            end
        end
        
        should_mask_central = always_mask_central || has_mbs

        raw_labels = NamedTuple[]
        for gap in plateaus
            # IDOP specific mapping
            N_gap = gap.plateau
            E_low = gap.mu_low
            E_high = gap.mu_high

            # --- Central Gap Filtering ---
            if should_mask_central && !isnothing(N)
                j_gap = round(Int, N_gap * N)
                center_idx = div(N, 2)
                if j_gap >= (center_idx - 1) && j_gap <= (center_idx + 1)
                     push!(raw_labels, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=0, q=SC_GAP_Q, err=0.0))
                     continue
                end
            end

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

                    # Parsimonious selection
                    if err < best_err - tol
                        best_err = err
                        best_p = p_try
                        best_q = q_candidate
                    elseif err <= best_err + tol
                        # Prefer smaller |q|
                        if ismissing(best_q) || abs(q_candidate) < abs(best_q)
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

            final_q = best_q
            if !ismissing(best_q) && abs(best_q) > q_max
                final_q = Int(sign(best_q) * q_max)
            end

            push!(raw_labels, (E_low=E_low, E_high=E_high, N_gap=N_gap, p=best_p, q=final_q, err=best_err))
        end

        # --- Merge adjacent ---
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
    E_low = cluster[1].E_low
    E_high = cluster[end].E_high
    rep = argmin(x -> ismissing(x.err) ? Inf : x.err, cluster)
    return (E_low=E_low, E_high=E_high, N_gap=rep.N_gap, p=rep.p, q=rep.q, err=rep.err)
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


function calculate_winding_gaps(df, mu_values, eig_col::Symbol)
    n_mu = length(mu_values)
    is_majorana = falses(n_mu)
    
    mu_min, mu_max = minimum(mu_values), maximum(mu_values)
    mu_step = (mu_max - mu_min) / (n_mu - 1)
    
    # Check consistency of eigenvalue vectors
    if isempty(df); return []; end
    
    # Verify all rows have same length for the target column
    lengths = [length(row[eig_col]) for row in eachrow(df) if !ismissing(row[eig_col])]
    if !isempty(lengths) && !all(l -> l == n_mu, lengths)
        @warn "Data inconsistency: Eigenvalue vector lengths $(unique(lengths)) do not match mu_values length $n_mu"
        # We proceed, but this suggests the mu_values reconstruction might be slightly off 
        # or the data has varying sizes.
    end

    # Mark Majorana regions (Union of all phason slices)
    for row in eachrow(df)
        if !ismissing(row.mu_mp)
            # Handle sparse mu_mp (list of mu values)
            for val in row.mu_mp
                if !ismissing(val)
                    # Map value to index
                    idx = round(Int, (val - mu_min) / mu_step) + 1
                    if 1 <= idx <= n_mu
                        is_majorana[idx] = true
                    end
                end
            end
        end
    end
    
    # Find gaps (runs of false in is_majorana)
    gaps = []
    in_gap = false
    start_idx = 1
    
    for i in 1:n_mu
        if !is_majorana[i]
            if !in_gap
                in_gap = true
                start_idx = i
            end
        else
            if in_gap
                push!(gaps, (start_idx, i-1))
                in_gap = false
            end
        end
    end
    if in_gap
        push!(gaps, (start_idx, n_mu))
    end
    
    # Calculate Winding for each gap
    gap_windings = []
    
    # Sort df by phason to ensure correct tracing
    df_sorted = sort(df, :phason)
    phasons = df_sorted.phason
    
    for (idx_start, idx_end) in gaps
        idx_mid = div(idx_start + idx_end, 2)
        
        # Extract trace of eigenvalues at the gap midpoint
        # Use abs() because we look for zero-crossings (dips to 0)
        trace = [abs(coalesce(row[eig_col][idx_mid], 1.0)) for row in eachrow(df_sorted)]
        
        # Count dips (zero crossings approximation)
        thresh = 1e-2 
        n_dips = 0
        is_low = trace[1] < thresh
        
        for i in 2:length(trace)
            curr_low = trace[i] < thresh
            if curr_low && !is_low
                n_dips += 1
            end
            is_low = curr_low
        end
        
        # Check wrap around
        if length(phasons) > 1 && (phasons[end] - phasons[1]) >= 0.9 
             if (trace[1] < thresh) && !(trace[end] < thresh)
                 n_dips += 1
             end
        end
        
        q = n_dips
        push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], q))
    end
    
    return gap_windings
end

# ...existing code...
function calculate_winding_gaps_combined(
    df, 
    mu_values, 
    col_pos::Symbol, 
    col_neg::Symbol; 
    log_process::Bool=false
)
    # Sort by phason
    df_sorted = sort(df, :phason)
    n_phason = nrow(df_sorted)
    n_mu = length(mu_values)
    
    # 1. Build Combined Energy Matrix Z
    Z = zeros(n_phason, n_mu)
    for (j, row) in enumerate(eachrow(df_sorted))
        p_vec = coalesce.(row[col_pos], 1.0)
        n_vec = coalesce.(row[col_neg], 1.0)
        for k in 1:n_mu
            val_p = ismissing(p_vec[k]) ? 1.0 : abs(p_vec[k])
            val_n = ismissing(n_vec[k]) ? 1.0 : abs(n_vec[k])
            Z[j, k] = min(val_p, val_n)
        end
    end

    # 2. Identify Gaps (Robust Index-Based)
    is_majorana = falses(n_mu)
    
    # Check alignment
    if nrow(df_sorted) > 0
        len_mp = length(df_sorted.mu_mp[1])
        if len_mp != n_mu
            n_mu = min(len_mp, n_mu)
        end
    end

    for row in eachrow(df_sorted)
        if hasproperty(row, :mu_mp) && !ismissing(row.mu_mp)
            for k in 1:n_mu
                if !ismissing(row.mu_mp[k])
                    is_majorana[k] = true
                end
            end
        end
    end
    
    gap_regions = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if !is_majorana[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gap_regions, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gap_regions, (start_idx, n_mu)); end

    gap_windings = []
    
    # 3. Calculate Winding
    for (idx_start, idx_end) in gap_regions
        width = idx_end - idx_start + 1
        
        # --- FIX: Trim Margins ---
        # We exclude the outer 10% (or at least 1 pixel) to avoid bulk band edges
        # oscillating into the detection window.
        margin = max(1, div(width, 10))
        
        # Ensure we don't trim the whole gap if it's tiny
        if width <= 2
            scan_start = idx_start
            scan_end = idx_end
        else
            scan_start = idx_start + margin
            scan_end = idx_end - margin
        end
        
        # Scan for minimum energy in the "Deep Gap"
        trace = vec(minimum(Z[:, scan_start:scan_end], dims=2))
        
        # Count zero crossings (dips)
        # Adjust threshold based on log_process flag
        # If log_process is true, we expect to resolve much smaller energies, 
        # so we tighten the threshold to avoid counting "soft" gap closings.
        thresh = log_process ? 1e-3 : 0.02
        
        n_dips = 0
        
        in_dip = false
        for val in trace
            if val < thresh
                if !in_dip
                    n_dips += 1
                    in_dip = true
                end
            else
                in_dip = false
            end
        end
        
        # Handle wrap-around
        if trace[end] < thresh && trace[1] < thresh
             nothing
        end
        
        push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], n_dips))
    end
    
    return gap_windings
end

function calculate_winding_gaps_refined(
    df, 
    mu_values, 
    col_pos::Symbol, 
    col_neg::Symbol; 
    log_process::Bool=false,
    threshold::Float64=0.0
)
    # Sort by phason
    df_sorted = sort(df, :phason)
    n_phason = nrow(df_sorted)
    n_mu = length(mu_values)
    
    # 1. Build Combined Energy Matrix Z
    Z = zeros(n_phason, n_mu)
    for (j, row) in enumerate(eachrow(df_sorted))
        p_vec = coalesce.(row[col_pos], 1.0)
        n_vec = coalesce.(row[col_neg], 1.0)
        for k in 1:n_mu
            val_p = ismissing(p_vec[k]) ? 1.0 : abs(p_vec[k])
            val_n = ismissing(n_vec[k]) ? 1.0 : abs(n_vec[k])
            Z[j, k] = min(val_p, val_n)
        end
    end

    # 2. Identify Gaps (Robust Index-Based)
    is_majorana = falses(n_mu)
    
    # Check alignment
    if nrow(df_sorted) > 0
        len_mp = length(df_sorted.mu_mp[1])
        if len_mp != n_mu
            n_mu = min(len_mp, n_mu)
        end
    end

    for row in eachrow(df_sorted)
        if hasproperty(row, :mu_mp) && !ismissing(row.mu_mp)
            for k in 1:n_mu
                if !ismissing(row.mu_mp[k])
                    is_majorana[k] = true
                end
            end
        end
    end
    
    gap_regions = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if !is_majorana[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gap_regions, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gap_regions, (start_idx, n_mu)); end

    gap_windings = []
    
    # 3. Calculate Winding using Slice Counting with Interval Logic
    # We count the number of separated regions below threshold at each phason slice.
    # This correctly handles flat-bottomed minima and avoids merging crossing modes 
    # except at the exact crossing point (which median filters out).
    for (idx_start, idx_end) in gap_regions
        # Extract sub-matrix for this gap
        sub_Z = Z[:, idx_start:idx_end]
        width = size(sub_Z, 2)
        
        # Determine threshold logic
        if threshold > 0.0
            thresh = threshold
        else
            thresh = log_process ? 5e-5 : 5e-3
        end
        
        # Smart Edge Filtering:
        # Check if the edges are "polluted" by bulk bands.
        # If the first or last column is below threshold for > 90% of phason values,
        # we treat it as a bulk band artifact rather than a winding mode.
        left_polluted = false
        right_polluted = false
        if width > 0
            left_polluted = count(sub_Z[:, 1] .< thresh) > 0.9 * n_phason
            right_polluted = count(sub_Z[:, end] .< thresh) > 0.9 * n_phason
        end
        
        counts = Int[]
        
        for r in 1:n_phason
            row_vals = sub_Z[r, :]
            
            # Find indices below threshold
            below_thresh = findall(x -> x < thresh, row_vals)
            
            if isempty(below_thresh)
                push!(counts, 0)
                continue
            end
            
            # Identify connected intervals and filter them
            n_intervals = 0
            
            # Helper to process an interval [min_idx, max_idx]
            # Indices are 1-based relative to sub_Z width
            function check_interval(min_idx, max_idx)
                is_left = (min_idx == 1)
                is_right = (max_idx == width)
                
                if (is_left && left_polluted) || (is_right && right_polluted)
                    return 0
                end
                return 1
            end
            
            if !isempty(below_thresh)
                curr_start = below_thresh[1]
                curr_prev = below_thresh[1]
                
                for k in 2:length(below_thresh)
                    idx = below_thresh[k]
                    if idx == curr_prev + 1
                        # Continue interval
                        curr_prev = idx
                    else
                        # Scan finished interval
                        n_intervals += check_interval(curr_start, curr_prev)
                        # Start new
                        curr_start = idx
                        curr_prev = idx
                    end
                end
                # Check last interval
                n_intervals += check_interval(curr_start, curr_prev)
            end
            
            push!(counts, n_intervals)
        end
        
        # Use median to filter out transient events like crossings or noise
        q = round(Int, median(counts))
        consistency = isempty(counts) ? 0.0 : count(x -> x == q, counts) / length(counts)
        
        push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], q, consistency))
    end
    
    # Post-processing: Filter to ensure specific labels (winding numbers) appear only once per mu-range.
    # We treat positive and negative mu ranges separately due to symmetry.
    final_gaps = []
    
    # Split gaps into positive and negative ranges (based on center)
    gaps_pos = filter(g -> (g[1] + g[2])/2 >= 0, gap_windings)
    gaps_neg = filter(g -> (g[1] + g[2])/2 < 0, gap_windings)
    
    for subgroup in [gaps_pos, gaps_neg]
        # Always keep q=0 (trivial gaps)
        kept = filter(g -> g[3] == 0, subgroup)
        
        # For q >= 1, enforce uniqueness by picking the gap with highest consistency
        candidates = filter(g -> g[3] > 0, subgroup)
        best_by_q = Dict{Int, Any}()
        
        for gap in candidates
            q_val = gap[3]
            # If we haven't seen this q, or if this gap is more consistent, keep it
            if !haskey(best_by_q, q_val) || gap[4] > best_by_q[q_val][4]
                best_by_q[q_val] = gap
            end
        end
        
        append!(kept, values(best_by_q))
        append!(final_gaps, kept)
    end
    
    # Sort final list by mu start
    sort!(final_gaps, by = x -> x[1])
    
    return [(g[1], g[2], g[3]) for g in final_gaps]
end

function calculate_winding_gaps_sine_fit(
    df, 
    mu_values, 
    col_pos::Symbol, 
    col_neg::Symbol; 
    q_max::Int=10,
    unique_q::Bool=false,
    filter_edges::Bool=false,
    threshold::Float64=5e-3,
    log_process::Bool=false,
    fit_strategy::Symbol=:grid, # :grid or :fft
    fft_periods::Int=3, # of repetitions for FFT
    fft_samples::Int=256 # number of samples for FFT
)
    # Sort by phason
    df_sorted = sort(df, :phason)
    n_phason = nrow(df_sorted)
    n_mu = length(mu_values)
    phason_vals = df_sorted.phason

    fit_strategy ∈ (:grid, :fft) || error("fit_strategy must be :grid or :fft")
    if fit_strategy == :fft
        fft_periods >= 1 || error("fft_periods must be ≥ 1")
        fft_samples >= 8 || error("fft_samples must be ≥ 8")
    end

    wrap01(x) = x - floor(x)

    function prepare_phason_parameter(phi_vals::Vector{Float64})
        n = length(phi_vals)
        if n <= 1
            return wrap01.(phi_vals), 1.0
        end

        phi_mod = wrap01.(phi_vals)
        sorted_phi = sort(phi_mod)
        diffs = diff(sorted_phi)
        wrap_gap = wrap01(sorted_phi[1] - sorted_phi[end])
        gap_list = vcat(diffs, wrap_gap)
        max_gap, idx = findmax(gap_list)
        coverage = clamp(1.0 - max_gap, 1e-6, 1.0)

        if coverage >= 0.95
            return phi_mod, coverage
        end

        start_idx = idx == length(sorted_phi) ? 1 : idx + 1
        start_phi = sorted_phi[start_idx]
        shifted = wrap01.(phi_mod .- start_phi)
        scaled = shifted ./ coverage
        return clamp.(scaled, 0.0, 1.0), coverage
    end

    function dedup_sorted_phi(phi_vec::Vector{Float64}, y_vec::Vector{Float64})
        phi_unique = Float64[]
        y_unique = Float64[]
        i = 1
        n = length(phi_vec)
        while i <= n
            j = i
            acc = y_vec[i]
            count = 1
            while j < n && abs(phi_vec[j+1] - phi_vec[j]) < 1e-8
                j += 1
                acc += y_vec[j]
                count += 1
            end
            push!(phi_unique, phi_vec[i])
            push!(y_unique, acc / count)
            i = j + 1
        end
        return phi_unique, y_unique
    end

    function resample_uniform_series(phi_param::Vector{Float64}, y_vals::Vector{Float64}; n_samples::Int=256)
        if length(phi_param) <= 1
            val = isempty(y_vals) ? 0.0 : y_vals[1]
            return fill(val, n_samples)
        end

        idx = sortperm(phi_param)
        phi_sorted = phi_param[idx]
        y_sorted = y_vals[idx]
        phi_unique, y_unique = dedup_sorted_phi(phi_sorted, y_sorted)

        sample_points = collect(range(0.0, 1.0, length=n_samples))
        y_interp = similar(sample_points)

        for (idx_sp, sp) in enumerate(sample_points)
            if sp <= phi_unique[1]
                y_interp[idx_sp] = y_unique[1]
            elseif sp >= phi_unique[end]
                y_interp[idx_sp] = y_unique[end]
            else
                hi = searchsortedfirst(phi_unique, sp)
                lo = max(hi - 1, 1)
                if hi > length(phi_unique)
                    hi = length(phi_unique)
                    lo = hi - 1
                end
                x0 = phi_unique[lo]
                x1 = phi_unique[hi]
                y0 = y_unique[lo]
                y1 = y_unique[hi]
                t = isapprox(x1, x0; atol=1e-12) ? 0.0 : (sp - x0) / (x1 - x0)
                y_interp[idx_sp] = y0 + t * (y1 - y0)
            end
        end

        return y_interp
    end

    function compute_fft_scores_for_gap(g, q_max, periods::Int, n_samples::Int)
        y_uniform = resample_uniform_series(g.phi_param, g.y; n_samples=n_samples)
        y_centered = y_uniform .- mean(y_uniform)
        y_ext = repeat(y_centered, periods)
        N = length(y_ext)
        N <= 2 && return Dict{Int, Float64}()

        window = sin.(range(0, π, length=N))
        y_windowed = y_ext .* window

        spectrum = fft(y_windowed)
        half = fld(N, 2)
        scores = Dict{Int, Float64}()

        for k in 1:(half-1)
            freq = k / periods
            q_est = round(Int, 4 * freq)
            if 1 <= q_est <= q_max
                mag = abs(spectrum[k+1])
                if mag > get(scores, q_est, 0.0)
                    scores[q_est] = mag
                end
            end
        end

        return scores
    end

    function fit_sine_parameters(g, q)
        omega = q * pi / 2
        theta = omega .* g.phi_param
        s = sin.(theta)
        c = cos.(theta)
        M = hcat(s, c)
        coeffs = M \ g.y
        A = coeffs[1]
        B = coeffs[2]
        delta = atan(B, A)
        return g.amp, delta, theta
    end
    
    # 1. Build Combined Energy Matrix Z
    Z = zeros(n_phason, n_mu)
    for (j, row) in enumerate(eachrow(df_sorted))
        p_vec = coalesce.(row[col_pos], 1.0)
        n_vec = coalesce.(row[col_neg], 1.0)
        for k in 1:n_mu
            val_p = ismissing(p_vec[k]) ? 1.0 : abs(p_vec[k])
            val_n = ismissing(n_vec[k]) ? 1.0 : abs(n_vec[k])
            val = min(val_p, val_n)
            
            if log_process
                # Transform to log scale, avoiding log(0)
                Z[j, k] = log10(val + 1e-20)
            else
                Z[j, k] = val
            end
        end
    end

    # 2. Identify Gaps (Robust Index-Based)
    is_majorana = falses(n_mu)
    if nrow(df_sorted) > 0
        len_mp = length(df_sorted.mu_mp[1])
        if len_mp != n_mu; n_mu = min(len_mp, n_mu); end
    end
    for row in eachrow(df_sorted)
        if hasproperty(row, :mu_mp) && !ismissing(row.mu_mp)
            for k in 1:n_mu
                if !ismissing(row.mu_mp[k]); is_majorana[k] = true; end
            end
        end
    end
    
    gap_regions = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if !is_majorana[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gap_regions, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gap_regions, (start_idx, n_mu)); end

    # 3. Extract Traces and Fit
    # Structure to hold gap info
    gap_data = [] 
    
    for (idx, (idx_start, idx_end)) in enumerate(gap_regions)
        sub_Z = Z[:, idx_start:idx_end]
        width = size(sub_Z, 2)
        if width < 3; continue; end
        
        mu_start = mu_values[idx_start]
        mu_end = mu_values[idx_end]
        mu_center = (mu_start + mu_end) / 2.0
        gap_amp = (mu_end - mu_start) / 2.0 # Half-width
        
        # Edge Pollution Check
        left_polluted = false
        right_polluted = false
        if filter_edges
            # Adjust threshold for log scale if needed
            eff_thresh = log_process ? log10(threshold + 1e-20) : threshold
            
            # Check if edges are consistently low energy (>90% of phason slices)
            left_polluted = count(sub_Z[:, 1] .< eff_thresh) > 0.9 * n_phason
            right_polluted = count(sub_Z[:, end] .< eff_thresh) > 0.9 * n_phason
        end

        trace_phi = Float64[]
        trace_y = Float64[]
        
        for r in 1:n_phason
            local_vals = sub_Z[r, :]
            min_val, min_k = findmin(local_vals)
            
            # If filtering is ON and the minimum is at a polluted edge, ignore this point
            if filter_edges
                if (min_k == 1 && left_polluted) || (min_k == width && right_polluted)
                    continue
                end
            end
            
            # Helper to get mu value relative to center
            mu_val = mu_values[idx_start + min_k - 1]
            
            push!(trace_phi, phason_vals[r])
            push!(trace_y, mu_val - mu_center)
        end
        
        # Only add if we have enough points
        if length(trace_phi) > 0.5 * n_phason
            phi_param, coverage = prepare_phason_parameter(trace_phi)
            push!(gap_data, (
                id = idx,
                mu_start = mu_start,
                mu_end = mu_end,
                amp = gap_amp,
                phi = trace_phi,
                phi_param = phi_param,
                coverage = coverage,
                y = trace_y
            ))
        end
    end
    
    assignments = [] # (mu_start, mu_end, q)
    fit_results = [] # (phi, mu, q)
    
    # Process Positive and Negative ranges separately if required (assumption: yes)
    groups = [
        filter(g -> (g.mu_start + g.mu_end)/2 >= 0, gap_data),
        filter(g -> (g.mu_start + g.mu_end)/2 < 0, gap_data)
    ]
    
    for group in groups
        if isempty(group); continue; end
        
        available_gaps = Set([g.id for g in group])
        gap_map = Dict(g.id => g for g in group)
        
        # Pre-calculate errors for ALL (gap, q) pairs to handle unique vs non-unique logic
        # Store as (mse, gap_id) for each q
        
        # Function to compute MSE for a gap and a given q
        # RETURNS: (mse, best_delta, best_model_func)
        function compute_fit_loss_detailed(g, q)
            if q == 0; return (Inf, 0.0, nothing); end
            
            # Frequency parameter with adaptive phason scaling
            omega = q * pi / 2
            theta = omega .* g.phi_param
            
            # Optimize phase delta over grid
            best_mse = Inf
            best_delta = 0.0
            
            phases = range(0, stop=2*pi, length=40) # Broader search or analytical?
            
            for delta in phases
                # Model based on user request: SINE function
                model_y = g.amp .* sin.(theta .+ delta)
                mse = mean((g.y .- model_y).^2)
                
                if mse < best_mse
                    best_mse = mse
                    best_delta = delta
                end
            end
            
            return (best_mse / (g.amp^2), best_delta, theta)
        end

        # Calculate losses / scores for candidate q values
        q_fits = Dict{Int, Vector{Tuple{Int, Float64}}}()
        fit_details_cache = Dict{Tuple{Int, Int}, Tuple{Float64, Vector{Float64}, Float64}}()

        if fit_strategy == :grid
            for q in 1:q_max
                for gid in available_gaps
                    loss, delta, theta = compute_fit_loss_detailed(gap_map[gid], q)
                    push!(get!(q_fits, q, []), (gid, loss))
                    fit_details_cache[(gid, q)] = (delta, theta, gap_map[gid].amp)
                end
            end
        else
            for gid in available_gaps
                scores = compute_fft_scores_for_gap(gap_map[gid], q_max, fft_periods, fft_samples)
                isempty(scores) && continue
                for (q, score) in scores
                    amp_fit, delta, theta = fit_sine_parameters(gap_map[gid], q)
                    loss = -score
                    push!(get!(q_fits, q, []), (gid, loss))
                    fit_details_cache[(gid, q)] = (delta, theta, amp_fit)
                end
            end
        end
        
        # Assignment Logic
        if unique_q
            # Greedy hierarchy: q=1, then q=2...
            for q in 1:q_max
                if isempty(available_gaps); break; end
                
                candidates = get(q_fits, q, nothing)
                candidates === nothing && continue
                # Filter to available
                candidates = filter(c -> c[1] in available_gaps, candidates)
                if isempty(candidates); continue; end
                
                # Find best fit for this q
                best_g, best_loss = candidates[argmin(c -> c[2], candidates)]
                
                # Assign
                g = gap_map[best_g]
                push!(assignments, (g.mu_start, g.mu_end, q))
                
                # Reconstruct best fit for logging/plotting
                delta, theta, amp_fit = fit_details_cache[(best_g, q)]
                center = (g.mu_start + g.mu_end) / 2.0
                fitted_y = amp_fit .* sin.(theta .+ delta)
                fitted_mu = fitted_y .+ center
                push!(fit_results, (
                    phi = copy(g.phi),
                    phi_param = copy(g.phi_param),
                    theta = copy(theta),
                    mu = fitted_mu,
                    q = q,
                    amp = amp_fit,
                    center = center,
                    delta = delta,
                    delta_shift = mod2pi(delta + pi),
                    coverage = g.coverage
                ))

                delete!(available_gaps, best_g)
            end
        else
            # Independent assignment: each gap picks its best q
            for gid in available_gaps
                best_q = 0
                min_loss = Inf
                
                # Check q=1..q_max for this specific gap
                for q in 1:q_max
                    loss_list = get(q_fits, q, nothing)
                    loss_list === nothing && continue
                    idx = findfirst(x -> x[1] == gid, loss_list)
                    idx === nothing && continue
                    loss = loss_list[idx][2]
                    
                    if loss < min_loss
                        min_loss = loss
                        best_q = q
                    end
                end
                
                # Assign
                g = gap_map[gid]
                push!(assignments, (g.mu_start, g.mu_end, best_q))
                
                if best_q > 0
                    delta, theta, amp_fit = fit_details_cache[(gid, best_q)]
                    center = (g.mu_start + g.mu_end) / 2.0
                    fitted_y = amp_fit .* sin.(theta .+ delta)
                    fitted_mu = fitted_y .+ center
                    push!(fit_results, (
                        phi = copy(g.phi),
                        phi_param = copy(g.phi_param),
                        theta = copy(theta),
                        mu = fitted_mu,
                        q = best_q,
                        amp = amp_fit,
                        center = center,
                        delta = delta,
                        delta_shift = mod2pi(delta + pi),
                        coverage = g.coverage
                    ))
                end
            end
            empty!(available_gaps)
        end
        
        # Fill remaining with 0
        for gid in available_gaps
            g = gap_map[gid]
            push!(assignments, (g.mu_start, g.mu_end, 0))
        end
    end
    sort!(assignments, by = x -> x[1])
    
    return assignments, fit_results
end


# --- Custom Connected Components (BFS) ---
function connected_components(binary_mask::AbstractMatrix{Bool})
    rows, cols = size(binary_mask)
    labels = zeros(Int, rows, cols)
    current_label = 0
    queue = Tuple{Int, Int}[]
    sizehint!(queue, rows * cols)

    for c in 1:cols, r in 1:rows
        if binary_mask[r, c] && labels[r, c] == 0
            current_label += 1
            labels[r, c] = current_label
            push!(queue, (r, c))
            
            while !isempty(queue)
                curr_r, curr_c = popfirst!(queue)
                
                # Check 8-connected neighbors
                for dr in -1:1, dc in -1:1
                    (dr == 0 && dc == 0) && continue
                    nr, nc = curr_r + dr, curr_c + dc
                    
                    if 1 <= nr <= rows && 1 <= nc <= cols
                        if binary_mask[nr, nc] && labels[nr, nc] == 0
                            labels[nr, nc] = current_label
                            push!(queue, (nr, nc))
                        end
                    end
                end
            end
        end
    end
    return labels, current_label
end

function calculate_winding_gaps_nodal_lines(df, mu_values, eig_col::Symbol)
    # 1. Construct the 2D Matrix Z[phason, mu]
    # Ensure df is sorted by phason
    df_sorted = sort(df, :phason)
    phasons = df_sorted.phason
    n_phason = length(phasons)
    n_mu = length(mu_values)
    
    Z = zeros(n_phason, n_mu)
    for (j, row) in enumerate(eachrow(df_sorted))
        val_vec = row[eig_col]
        if length(val_vec) == n_mu
            # Use abs() to treat pos/neg symmetrically if needed, 
            # though min_pos_eigs is already positive.
            Z[j, :] .= replace(abs.(val_vec), missing => 1.0)
        else
            Z[j, :] .= 1.0 # Fill with high energy if missing
        end
    end

    # 2. Identify Gaps (Columns where MP != -1 generally)
    # Heuristic: If median energy of a column is high, it's a gap (no Majorana)
    # If median energy is ~0, it's a Majorana region.
    is_gap_col = [median(Z[:, i]) > 1e-3 for i in 1:n_mu]
    
    # Group into contiguous gap regions
    gap_regions = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if is_gap_col[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gap_regions, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gap_regions, (start_idx, n_mu)); end

    gap_windings = []

    # 3. Process each gap region
    for (idx_start, idx_end) in gap_regions
        # Extract sub-matrix for this gap
        sub_Z = Z[:, idx_start:idx_end]
        
        # Adaptive threshold: look for "nodal lines" (near zero)
        # We want to count distinct lines that traverse the phason axis.
        local_min = minimum(sub_Z)
        # Threshold: slightly above minimum, but capped to avoid noise
        thresh = max(local_min * 5, 1e-2) 
        
        binary_img = sub_Z .< thresh
        
        # Label connected components using our custom function
        labels, n_components = label_components_simple(binary_img)
        
        q_count = 0
        for k in 1:n_components
            # Find all pixels belonging to this component
            indices = findall(x -> x == k, labels)
            if isempty(indices); continue; end
            
            rows_involved = [idx[1] for idx in indices]
            min_r, max_r = extrema(rows_involved)
            
            # Check if it spans the phason range (allowing for small edge effects)
            # If it covers > 90% of the phason axis, we count it as a winding mode.
            if (max_r - min_r) >= 0.9 * n_phason
                q_count += 1
            end
        end
        
        push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], q_count))
    end

    return gap_windings
end

function calculate_winding_gaps_slice_counting(df, mu_values, eig_col::Symbol)
    # Sort by phason to ensure rows correspond to sequential phason values
    df_sorted = sort(df, :phason)
    n_phason = nrow(df_sorted)
    n_mu = length(mu_values)
    
    # Build Matrix Z[phason, mu]
    Z = zeros(n_phason, n_mu)
    for (j, row) in enumerate(eachrow(df_sorted))
        val_vec = row[eig_col]
        if length(val_vec) == n_mu
            Z[j, :] .= replace(abs.(val_vec), missing => 1.0)
        else
            Z[j, :] .= 1.0
        end
    end

    # Identify Gaps: Columns where median energy is high (not a bulk Majorana region)
    # Threshold 1e-3 separates "pinned at zero" regions from "gapped but dipping" regions
    is_gap_col = [median(Z[:, i]) > 1e-3 for i in 1:n_mu]
    
    gap_regions = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if is_gap_col[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gap_regions, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gap_regions, (start_idx, n_mu)); end

    gap_windings = []
    
    # Threshold for a "zero mode" dip
    # Must be significantly smaller than the gap, but allow for finite-size effects
    zero_mode_thresh = 0.02 

    for (idx_start, idx_end) in gap_regions
        sub_Z = Z[:, idx_start:idx_end]
        width = size(sub_Z, 2)
        if width < 3; continue; end # Too narrow to find peaks reliably

        counts = Int[]
        for r in 1:n_phason
            row = sub_Z[r, :]
            n_min = 0
            # Find local minima in mu direction
            for k in 2:(width-1)
                val = row[k]
                # Check if local minimum
                if val < row[k-1] && val < row[k+1]
                    # Check if it's a "zero mode"
                    if val < zero_mode_thresh
                        n_min += 1
                    end
                end
            end
            push!(counts, n_min)
        end
        
        # The winding number is the typical number of zero modes observed
        # Median is robust against noise or transient crossings
        q = round(Int, median(counts))
        
        if q > 0
            push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], q))
        end
    end
    
    return gap_windings
end

function calculate_winding_wilson_loop(df, mu_values, evec_col::Symbol)
    """
    Calculates the winding number of the subgap states using the discretized Wilson loop method.
    Requires 'evec_col' to contain the eigenvectors (e.g. N x 2 matrix) for the relevant states.
    """
    # Sort by phason to ensure correct order for integration
    df_sorted = sort(df, :phason)
    phasons = df_sorted.phason
    n_phason = length(phasons)
    n_mu = length(mu_values)
    
    # 1. Identify Gaps (regions where we want to calculate winding)
    # We reuse the logic: "Gap" = region where system is NOT in the bulk Majorana phase (MP != -1)
    # We assume 'mu_mp' column exists to identify these regions.
    
    is_majorana = falses(n_mu)
    mu_min, mu_max = minimum(mu_values), maximum(mu_values)
    mu_step = (mu_max - mu_min) / (n_mu - 1)
    
    for row in eachrow(df)
        if hasproperty(row, :mu_mp) && !ismissing(row.mu_mp)
            for val in row.mu_mp
                if !ismissing(val)
                    idx = round(Int, (val - mu_min) / mu_step) + 1
                    if 1 <= idx <= n_mu; is_majorana[idx] = true; end
                end
            end
        end
    end
    
    gaps = []
    in_gap = false; start_idx = 1
    for i in 1:n_mu
        if !is_majorana[i]
            if !in_gap; in_gap = true; start_idx = i; end
        else
            if in_gap; push!(gaps, (start_idx, i-1)); in_gap = false; end
        end
    end
    if in_gap; push!(gaps, (start_idx, n_mu)); end
    
    gap_windings = []
    
    # 2. Calculate Wilson Loop for each gap
    for (idx_start, idx_end) in gaps
        # We pick the middle of the gap to trace
        idx_mid = div(idx_start + idx_end, 2)
        
        total_phase = 0.0
        valid_path = true
        
        for i in 1:n_phason
            # Wrap around index for next step
            next_i = (i % n_phason) + 1
            
            # Extract eigenvector matrices for the chosen mu index
            # Assumes row[evec_col] is a Vector of Matrices (one per mu)
            # or a Matrix of Matrices.
            
            # Check if data exists
            if ismissing(df_sorted[i, evec_col]) || ismissing(df_sorted[next_i, evec_col])
                valid_path = false; break
            end
            
            # Get the specific eigenvectors for the gap midpoint
            # We assume the column contains a Vector{Any} of length n_mu
            evecs_curr = df_sorted[i, evec_col][idx_mid]
            evecs_next = df_sorted[next_i, evec_col][idx_mid]
            
            if ismissing(evecs_curr) || ismissing(evecs_next)
                valid_path = false; break
            end
            
            # Calculate Overlap Matrix M = U_i^dagger * U_{i+1}
            # evecs should be N_sites x k matrices (e.g. k=2)
            M = evecs_curr' * evecs_next
            
            # Accumulate phase: Im(ln(det(M)))
            # This tracks the rotation of the subspace
            det_M = det(M)
            if abs(det_M) < 1e-6
                # Singularity or orthogonality (should not happen in adiabatic evolution)
                valid_path = false; break
            end
            
            phase_step = imag(log(det_M))
            total_phase += phase_step
        end
        
        if valid_path
            # Winding number = Total Phase / 2pi
            q = round(Int, total_phase / (2π))
            # Only record non-trivial windings if desired, or all
            if q != 0
                push!(gap_windings, (mu_values[idx_start], mu_values[idx_end], q))
            end
        end
    end
    
    return gap_windings
end



############################################
###### Process slope-dependent MP ##########
############################################


# function extract_mp_tol(filepath::String)
#     # Load data
#     df = CSV.read(filepath, DataFrame)

#     # Check required columns
#     required_cols = [:N, :Delta, :t_n, :phi, :phason, :qc_region_count, :mp_tol_ranges]
#     for col in required_cols
#         hasproperty(df, col) || error("DataFrame must contain column: $col")
#     end

#     df.phi = Float64.(df.phi)
#     df.phason = Float64.(df.phason)

#     # Parse mp_tol_ranges and compute midpoint
#     mp_tols = Float64[]
#     for r in df.mp_tol_ranges
#         s = string(r)
#         # Remove parentheses and spaces
#         s_clean = replace(s, "(" => "", ")" => "", " " => "")
#         parts = split(s_clean, ',')
#         if length(parts) == 2
#             v1 = parse(Float64, parts[1])
#             v2 = parse(Float64, parts[2])
#             push!(mp_tols, (v1 + v2) / 2)
#         else
#             error("Invalid mp_tol_ranges format: $s")
#         end
#     end
    
#     # Create output df with :phi, :phason, :mp_tol
#     out_df = select(df, [:phi, :phason])
#     out_df.mp_tol = mp_tols

#     return out_df
# end

function extract_mp_tol(filepath::String)
    # Load data from BSON
    data = BSON.load(filepath)
    found_mp_tol_ranges = data[:found_mp_tol_ranges]  # Vector{MPTolResult}

    # Extract phi, phason, and compute mp_tol
    phis = Float64[]
    phasons = Float64[]
    mp_tols = Float64[]

    for result in found_mp_tol_ranges
        push!(phis, result.phi)
        push!(phasons, result.phason)
        if result.mp_tol_ranges === nothing
            push!(mp_tols, NaN)  # Use NaN for missing ranges
        else
            v1, v2 = result.mp_tol_ranges
            # Clamp the range to [1e-10, 1e-2]
            v1 = max(v1, 1e-10)
            v2 = min(v2, 1e-2)
            push!(mp_tols, (v1 + v2) / 2)
        end
    end

    # Create output df with :phi, :phason, :mp_tol
    out_df = DataFrame(phi=phis, phason=phasons, mp_tol=mp_tols)

    return out_df
end

function add_mp_disc_column!(
    df_slice::DataFrame, 
    mp_tol_df::Union{DataFrame, Float64}; 
    mp_targ::Float64=-1.0, 
    raw_col::Symbol=:mp_raw
)
    disc_mp_all = Vector{Vector{Int}}(undef, nrow(df_slice))
    for i in 1:nrow(df_slice)
        row = df_slice[i, :]
        if mp_tol_df isa Float64
            mp_tol = mp_tol_df
        else
            phi_val = row.phi
            phason_val = row.phason
            # find matching row in mp_tol_df
            matching_rows = filter(r -> r.phi ≈ phi_val && r.phason ≈ phason_val, mp_tol_df)
            if isempty(matching_rows)
                error("No matching mp_tol found for phi=$phi_val, phason=$phason_val")
            end
            mp_tol = first(matching_rows).mp_tol
        end
        # Discretize mp_raw: 1 if mp_raw <= mp_targ + mp_tol, else 0
        disc_mp = [mp <= mp_targ + mp_tol ? 1 : 0 for mp in row[raw_col]]
        disc_mp_all[i] = disc_mp
    end
    df_slice[!, :disc_mp] = disc_mp_all
    return df_slice
end

function add_mu_mp_column!(df_slice::DataFrame)
    mu_mp_all = Vector{Vector{Union{Missing, Float64}}}(undef, nrow(df_slice))
    for i in 1:nrow(df_slice)
        row = df_slice[i, :]
        disc_mp = row.disc_mp
        mu_mp = [disc_mp[j] == 1 ? row.mu_range[j] : missing for j in 1:length(disc_mp)]
        mu_mp_all[i] = mu_mp
    end
    df_slice[!, :mu_mp] = mu_mp_all
    return df_slice
end

end # module IDOPProcessing