"""
    file name:   processing.jl
    created:     2026-09-22

    overview:
        Processing module for spectral localizer butterfly data.
        Extracts spectral localizer invariant values from raw JLD2 chunk files,
        aggregates them by mu and phi (slope), and normalizes using eigenvalue
        bandwidth data from the IDOP pipeline.

    normalization strategy:
        Instead of computing bounds from specloc data alone, this uses pre-computed
        eigenvalue bandwidth data from the continuous MB IDOP calculation. This
        provides a more physical normalization based on the topological region
        boundaries (eigenvalue extrema) rather than data-driven bounds.

    structure:
        - load_chunk_files: Load all chunk files for a given data folder
        - extract_specloc_data: Extract specloc values and organize by mu and phi/slope
        - aggregate_specloc_data: Combine and optionally average over x0
        - downsample_data: Reduce grid resolution for fast plotting
        - mirror_mu_data: Create symmetric data about mu=0
        - load_eigs_bandwidth_norms: Load pre-computed eigenvalue bandwidth from IDOP CSV
        - apply_eigs_bandwidth_normalisation!: Apply eigenvalue-based normalization
"""

module SpecLocProcessing

using JLD2
using DataFrames
using ProgressMeter
using Statistics
using CSV

# Include bson_unpacker at module level to avoid world age issues
include("../bson_unpacker.jl")

export load_chunk_files, extract_specloc_data, aggregate_specloc_data, downsample_data, mirror_mu_data, compute_mu_bounds_specloc, apply_mu_normalisation!, compute_specloc_boundaries, apply_boundary_normalisation!, load_eigs_bandwidth_norms, apply_eigs_bandwidth_normalisation!, extract_eigs_spectrum_at_mu0, compute_eigs_normalization_bounds_from_spectrum, extract_and_compute_eigs_bounds


"""
    load_chunk_files(folder_path::String)

Load all JLD2 chunk files from a folder. Returns a vector of DataFrames.
"""
function load_chunk_files(folder_path::String)
    # Get all .jld2 files
    files = filter(f -> endswith(f, ".jld2"), readdir(folder_path))
    sort!(files)  # Sort to ensure consistent ordering
    
    dataframes::Vector{DataFrame} = DataFrame[]
    
    @showprogress "Loading chunks..." for filename in files
        filepath = joinpath(folder_path, filename)
        try
            data = JLD2.load(filepath)
            # The JLD2 file should contain a DataFrame or dictionary with a DataFrame
            if haskey(data, "result")
                df = data["result"]
            elseif length(data) == 1
                # If there's only one key, use that
                df = first(values(data))
            else
                df = data
            end
            
            # Make sure it's a DataFrame
            if isa(df, DataFrame)
                push!(dataframes, df)
            end
        catch e
            @warn "Failed to load $filename: $e"
        end
    end
    
    return dataframes
end


"""
    extract_specloc_data(dataframes::Vector{DataFrame})

Extract spectral localizer invariant values from DataFrames.
Returns organized data structure with columns: N, t_n, mu, phi, specloc_gap, specloc_winding
"""
function extract_specloc_data(dataframes::Vector{DataFrame})
    # Combine all DataFrames
    combined_df = vcat(dataframes...)
    
    # Extract relevant columns - look for specloc-related columns
    # Common names: specloc_gap, specloc_winding, winding, localiser_gap
    
    extracted_data = DataFrame(
        N = Int64[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        phi = Float64[],  # slope from sequence_id[1]
        phason = Float64[],  # phason from sequence_id[2]
        x0 = Float64[],
        Delta = Float64[],
        specloc_invariant = Float64[],
        specloc_gap = Float64[]
    )

    # Diagnose x0 column
    x0_col = if "x0" in names(combined_df)
        :x0
    elseif "x_0" in names(combined_df)
        :x_0
    else
        nothing
    end
    
    if isnothing(x0_col)
        @warn "x0 column not found in combined DataFrame. Available columns: $(names(combined_df))"
    else
        sample_x0 = combined_df[1:min(10, nrow(combined_df)), x0_col]
        @info "x0 column detected: $x0_col. Sample values: $sample_x0"
    end
    
    for row in eachrow(combined_df)
        # Extract basic parameters
        N = get(row, :N, missing)
        t_n = get(row, :t_n, missing)
        mu = get(row, :mu, missing)
        Delta = get(row, :Delta, missing)
        sequence_id = get(row, :sequence_id, missing)
        x0 = missing
        if !isnothing(x0_col)
            x0 = get(row, x0_col, missing)
        end
        
        # Extract phi (slope) and phason from sequence_id tuple
        phi = missing
        phason = missing
        if !ismissing(sequence_id) && isa(sequence_id, Tuple) && length(sequence_id) >= 2
            phi = convert(Float64, sequence_id[1])
            phason = convert(Float64, sequence_id[2])
        end
        
        # Extract specloc values
        specloc_invariant = missing
        specloc_gap = missing
        
        # Try to get specloc_invariant first (main quantity)
        if hasproperty(row, :specloc_invariant)
            specloc_invariant = convert(Float64, get(row, :specloc_invariant, missing))
        end
        
        # Also get specloc_gap if available
        if hasproperty(row, :specloc_gap)
            specloc_gap = convert(Float64, get(row, :specloc_gap, missing))
        end
        
        # Only add if we have the essential data
        if !ismissing(phi) && !ismissing(phason) && !ismissing(specloc_invariant)
            push!(extracted_data, (
                N = ismissing(N) ? 500 : N,
                t_n = ismissing(t_n) ? [1.0, 2.0] : t_n,
                mu = ismissing(mu) ? 0.0 : mu,
                phi = phi,
                phason = phason,
                x0 = ismissing(x0) ? NaN : convert(Float64, x0),
                Delta = ismissing(Delta) ? 0.05 : Delta,
                specloc_invariant = specloc_invariant,
                specloc_gap = ismissing(specloc_gap) ? specloc_invariant : specloc_gap
            ))
        end
    end
    
    # Diagnostic: check extracted x0 values
    n_valid_x0 = count(!isnan, extracted_data.x0)
    @info "Extraction complete: $(nrow(extracted_data)) rows extracted, $n_valid_x0 rows with valid x0 (non-NaN)"
    if n_valid_x0 > 0
        @info "x0 range in extracted data: $(minimum(extracted_data.x0[.!isnan.(extracted_data.x0)])) to $(maximum(extracted_data.x0[.!isnan.(extracted_data.x0)]))"
    end
    
    return extracted_data
end


"""
    aggregate_specloc_data(df::DataFrame; by_mean::Bool=true, average_over_x0::Bool=false)

Aggregate specloc data by (mu, phi, Delta, x0) groups for efficient plotting.
Takes statistics (mean) if multiple entries per grid point, reducing data points significantly.

Returns aggregated DataFrame with columns: mu, phi, phason, specloc_invariant, specloc_gap, N, Delta
"""
function aggregate_specloc_data(df::DataFrame; by_mean::Bool=true, average_over_x0::Bool=false)
    if by_mean
        # Group either per-x0 or averaged-over-x0.
        # Per-x0 keeps one 1201x1000 map per x0 value.
        grouping_cols = average_over_x0 ? [:mu, :phi, :Delta] : [:mu, :phi, :Delta, :x0]
        grouped = combine(
            groupby(df, grouping_cols),
            :specloc_invariant => mean => :specloc_invariant,
            :specloc_gap => mean => :specloc_gap,
            :phason => first => :phason,  # phason should be constant for a group
            :N => first => :N,  # N should be constant for a group
            :t_n => (values -> Ref(first(values))) => :t_n,
            :x0 => first => :x0  # Explicitly preserve x0
        )
        if average_over_x0
            grouped.x0 = fill(NaN, nrow(grouped))
        end
        return grouped
    else
        # No aggregation, return as-is
        return df
    end
end


"""
    downsample_data(df::DataFrame, factor::Int=10)

Downsample data by taking every Nth unique (mu, phi) point.
Useful for quick visualization and testing.

Returns downsampled DataFrame.
"""
function downsample_data(df::AbstractDataFrame, factor::Int=10)
    unique_mus = sort(unique(df.mu))
    unique_phis = sort(unique(df.phi))
    
    # Subsample by taking every factor-th point
    downsampled_mus = unique_mus[1:factor:end]
    downsampled_phis = unique_phis[1:factor:end]
    downsampled_mu_set = Set(downsampled_mus)
    downsampled_phi_set = Set(downsampled_phis)
    
    # Filter to keep only downsampled points
    mask = [row.mu in downsampled_mu_set && row.phi in downsampled_phi_set for row in eachrow(df)]
    downsampled_df = df[mask, :]
    
    n_original = nrow(df)
    n_downsampled = nrow(downsampled_df)
    @info "Downsampled from $n_original to $n_downsampled points (factor $factor)"
    
    return downsampled_df
end


"""
    mirror_mu_data(df::AbstractDataFrame)

Mirror the data about mu=0. If all mu values are >= 0, create mirror copies with -mu.
Useful for creating symmetric butterfly plots around the origin.

Returns the original DataFrame if already symmetric, or vcat of original + mirrored.
"""
function mirror_mu_data(df::AbstractDataFrame)
    min_mu = minimum(df.mu)
    max_mu = maximum(df.mu)
    
    # Check if already symmetric (has both positive and negative)
    has_neg = any(df.mu .< 0)
    has_pos = any(df.mu .> 0)
    
    if has_neg && has_pos
        @info "Data already symmetric around mu=0 (range: $min_mu to $max_mu). No mirroring needed."
        return df
    end
    
    if min_mu >= 0
        @info "Mirroring data about mu=0 (original range: $min_mu to $max_mu)"
        df_mirror = copy(df)
        df_mirror.mu .= .-df_mirror.mu  # Negate mu values
        return vcat(df, df_mirror)
    else
        @warn "Unexpected mu range: $min_mu to $max_mu. No mirroring applied."
        return df
    end
end


"""
    compute_mu_bounds_specloc(df::AbstractDataFrame, by_cols::Vector{Symbol})

Compute global_min_mu and global_max_mu bounds per parameter group.
Faithful to IDOP's approach: groups data and extracts mu range for each combination.

Args:
    df: DataFrame with :mu, :phi, :Delta, and optionally :x0
    by_cols: Parameter columns to group by, e.g. [:phi, :Delta, :x0] for per-x0
             or [:phi, :Delta] for x0-averaged

Returns:
    DataFrame with parameter columns + :global_min_mu + :global_max_mu
    (Analogous to IDOP's norm CSV structure)
"""
function compute_mu_bounds_specloc(df::AbstractDataFrame, by_cols::Vector{Symbol})
    # Ensure all grouping columns exist
    for col in by_cols
        if !(string(col) in names(df))
            error("Column $col not found in DataFrame. Available: $(names(df))")
        end
    end
    
    # Create mu_specloc_invariant column (the x-axis values to be normalized)
    if !("mu_specloc_invariant" in names(df))
        df = copy(df)
        df[!, :mu_specloc_invariant] = df.mu
    end
    
    # Group and compute bounds per parameter combination
    bounds_df = combine(
        groupby(df, by_cols),
        :mu_specloc_invariant => minimum => :global_min_mu,
        :mu_specloc_invariant => maximum => :global_max_mu,
        :N => first => :N,
        :Delta => first => :Delta
    )
    
    @info "Computed mu bounds for $(nrow(bounds_df)) parameter groups"
    
    return bounds_df
end


"""
    apply_mu_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame; 
                             by_cols::Vector{Symbol}=[:phi, :Delta, :x0])

Apply mu normalization using pre-computed bounds (IDOP-faithful approach).
Normalizes mu values to [-1, 1] range using: 2 * (x - global_min_mu) / (global_max_mu - global_min_mu) - 1

Creates column :mu_specloc_invariant_norm with normalized mu values.

Args:
    df: Data DataFrame with :mu column
    bounds_df: Bounds DataFrame from compute_mu_bounds_specloc()
    by_cols: Parameter columns to match on (must match those used in bounds computation)

Returns:
    df with added :mu_specloc_invariant_norm column
"""
function apply_mu_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame; 
                                  by_cols::Vector{Symbol}=[:phi, :Delta, :x0])
    # Ensure mu_specloc_invariant exists
    if !("mu_specloc_invariant" in names(df))
        df = copy(df)
        df[!, :mu_specloc_invariant] = df.mu
    end
    
    # Create lookup dict for O(1) bound retrieval (only ~100-200 entries)
    bounds_dict = Dict()
    
    @inbounds for row in eachrow(bounds_df)
        key = Tuple(row[col] for col in by_cols)
        bounds_dict[key] = (row.global_min_mu, row.global_max_mu)
    end
    
    # Initialize output vector (avoids combine() row reordering issue)
    norm_values = Vector{Union{Float64, Missing}}(undef, nrow(df))
    n_groups_processed = 0
    n_unmatched = 0
    
    # Use groupby to process each group but preserve row indices
    # parentindices(group)[1] gives original row positions without copying
    for group in groupby(df, by_cols)
        key = Tuple(group[1, col] for col in by_cols)
        group_indices = parentindices(group)[1]  # Get original row indices for this group
        
        if haskey(bounds_dict, key)
            vmin, vmax = bounds_dict[key]
            if vmax ≈ vmin
                norm_values[group_indices] .= 0.0
            else
                # Vectorized normalization for entire group at once
                @. norm_values[group_indices] = 2 * (group.mu_specloc_invariant - vmin) / (vmax - vmin) - 1
            end
            n_groups_processed += 1
        else
            # No bounds found for this group
            norm_values[group_indices] .= missing
            n_unmatched += length(group_indices)
        end
    end
    
    # Quantize normalized values to 4 decimal places to eliminate floating-point artifacts
    # from eigenvalue data (which has much higher mu resolution than specloc data).
    # This reduces thousands of unique floating-point values to ~200-400 discrete values
    # while preserving physical accuracy in the [-1, 1] range.
    @inbounds for i in eachindex(norm_values)
        if !ismissing(norm_values[i])
            norm_values[i] = round(norm_values[i], digits=4)
        end
    end
    
    # Assign pre-computed vector as column (once, not in loop)
    df[!, :mu_specloc_invariant_norm] = norm_values
    
    if n_unmatched > 0
        @warn "Could not match $n_unmatched rows to bounds (missing parameter combination)"
    end
    
    @info "Applied mu normalization to $(nrow(df)) rows via $n_groups_processed parameter groups (quantized to 4 decimals)"
    
    return df
end


"""
    compute_specloc_boundaries(df::AbstractDataFrame, by_cols::Vector{Symbol})

Extract critical mu values where specloc_invariant transitions from 1 to 0.
Faithful to IDOP's mu_mp extraction: for each parameter group, collect all mu values
where specloc_invariant == 1. These form a vector representing the "safe" region boundary.

Args:
    df: DataFrame with :mu, :specloc_invariant, and grouping columns
    by_cols: Parameter columns to group by, e.g. [:phi, :Delta, :x0]

Returns:
    Dictionary mapping parameter tuples to vectors of boundary mu values
    Format: (phi, Delta, x0) => [mu_values_where_specloc==1]
    
Similar structure to IDOP's mu_mp extraction where mu_mp contains mu values where disc_mp == 1.
"""
function compute_specloc_boundaries(df::AbstractDataFrame, by_cols::Vector{Symbol})
    # Ensure all required columns exist
    for col in by_cols
        if !(string(col) in names(df))
            error("Column $col not found in DataFrame. Available: $(names(df))")
        end
    end
    
    if !("specloc_invariant" in names(df))
        error("Column specloc_invariant not found in DataFrame")
    end
    
    # Dictionary to store boundary vectors: key => [mu_values_where_invariant==1]
    boundary_dict = Dict()
    
    # Group by parameter columns and extract boundaries per group (vectorized)
    for group in groupby(df, by_cols)
        key = Tuple(group[1, col] for col in by_cols)
        
        # Vectorized filter: get all mu where invariant == 1.0 (faster than indexing)
        mask = @. group.specloc_invariant == 1.0
        boundary_mus = group.mu[mask]
        
        # Store as vector (even if empty)
        boundary_dict[key] = copy(boundary_mus)  # copy ensures it's a standalone vector
    end
    
    @info "Extracted specloc boundaries for $(length(boundary_dict)) parameter groups"
    
    return boundary_dict
end


"""
    apply_boundary_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame, 
                                   boundary_dict::Dict; by_cols::Vector{Symbol}=[:phi, :Delta, :x0])

Normalize extracted boundary mu values using the same min/max bounds as the mu axis.
Faithful to IDOP's apply_outer_normalisation! applied to mu_mp.

Args:
    df: Data DataFrame
    bounds_df: Bounds DataFrame from compute_mu_bounds_specloc()
    boundary_dict: Boundary dictionary from compute_specloc_boundaries()
    by_cols: Parameter columns used for grouping

Creates columns:
    :mu_specloc_boundary - Original boundary mu values where invariant == 1
    :mu_specloc_boundary_norm - Normalized boundary values using axis bounds
    
The normalized boundary marks where |mu_norm| > |mu_boundary_norm| should be masked/zeroed.
"""
function apply_boundary_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame, 
                                        boundary_dict::Dict; by_cols::Vector{Symbol}=[:phi, :Delta, :x0])
    
    # PRE-COMPUTE all normalizations once per group (not per row)
    # This avoids 1.2M dict lookups and instead does ~100-200 lookups
    combined_lookup = Dict()
    
    @inbounds for row in eachrow(bounds_df)
        key = Tuple(row[col] for col in by_cols)
        vmin, vmax = row.global_min_mu, row.global_max_mu
        boundary_mus = get(boundary_dict, key, Float64[])
        
        # Pre-normalize boundaries once per group, not per row (vectorized @. macro)
        if !isempty(boundary_mus)
            if vmax ≈ vmin
                boundary_norm = zeros(length(boundary_mus))
            else
                boundary_norm = @. 2 * (boundary_mus - vmin) / (vmax - vmin) - 1
            end
        else
            boundary_norm = Float64[]
        end
        
        combined_lookup[key] = (boundary_mus, boundary_norm)
    end
    
    # Initialize output vectors (assign once after loop, not per row)
    boundary_col = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df))
    boundary_norm_col = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df))
    
    n_unmatched = 0
    
    # KEY OPTIMIZATION: Use groupby to process groups, not individual rows
    # Reduces from 1.2M iterations to ~100-200 iterations
    for group in groupby(df, by_cols)
        key = Tuple(group[1, col] for col in by_cols)
        group_indices = axes(group, 1)  # Get row indices for this group
        
        if haskey(combined_lookup, key)
            boundary_mus, boundary_norm = combined_lookup[key]
            # Assign same values to all rows in group (vectorized operation)
            for idx in group_indices
                boundary_col[idx] = boundary_mus
                boundary_norm_col[idx] = boundary_norm
            end
        else
            # No match: mark all rows in group as missing
            for idx in group_indices
                boundary_col[idx] = missing
                boundary_norm_col[idx] = missing
            end
            n_unmatched += length(group_indices)
        end
    end
    
    # Assign pre-computed columns ONCE (not in loop)
    df[!, :mu_specloc_boundary] = boundary_col
    df[!, :mu_specloc_boundary_norm] = boundary_norm_col
    
    if n_unmatched > 0
        @warn "Could not match $n_unmatched rows to boundary data (missing parameter combination)"
    end
    
    @info "Applied boundary normalization to $(nrow(df)) rows via $(length(combined_lookup)) groups"
    
    return df
end



"""
    extract_eigs_spectrum_at_mu0(butterfly_data_path::String; mu_tol::Float64=0.01)

Extract eigenvalue spectra from raw butterfly_data at mu ≈ 0.
The eigenspectrum at mu=0 provides the reference for normalization:
for each phi (slope), the min/max of eigenvalues defines the normalization scale.

This approach:
- Avoids CSV resolution mismatch (uses only needed mu points)
- Provides per-phi independent normalization
- Directly matches eigenspectrum at mu=0 to normalize mu axis

Args:
    butterfly_data_path: Path to butterfly_data folder containing BSON files
    mu_tol: Tolerance for selecting mu ≈ 0 (default 0.01)

Returns:
    DataFrame with columns: phi, Delta, N, t_n, phason, eigenvalues
    Each row contains the full eigenspectrum for (phi, Delta, N, t_n, phason) at mu≈0
"""
function extract_eigs_spectrum_at_mu0(butterfly_data_path::String; mu_tol::Float64=0.01)
    @info "Extracting eigenvalue spectra at mu≈0 from butterfly_data..."
    
    df_raw = process_bson_files(butterfly_data_path)
    
    if isempty(df_raw)
        error("No data loaded from butterfly_data folder: $butterfly_data_path")
    end
    
    @info "  Loaded $(nrow(df_raw)) raw records"
    
    # Extract phi and phason from sequence_id FIRST (like IDOP does)
    # phi = sequence_id[1], phason = sequence_id[2]
    df_raw[!, :phi] = getindex.(df_raw.sequence_id, 1)
    df_raw[!, :phason] = getindex.(df_raw.sequence_id, 2)
    
    @info "  Extracted phi and phason from sequence_id"
    
    # Filter to mu ≈ 0 (within tolerance)
    df_mu0 = filter(row -> abs(row.mu) <= mu_tol, df_raw)
    
    if nrow(df_mu0) == 0
        @warn "No records found near mu=0 (within tolerance $mu_tol). Trying broader tolerance..."
        mu_tol_broad = 0.1
        df_mu0 = filter(row -> abs(row.mu) <= mu_tol_broad, df_raw)
        if nrow(df_mu0) == 0
            error("Could not find any records near mu=0")
        end
        @warn "  Using broader tolerance $mu_tol_broad instead"
    end
    
    @info "  Filtered to $(nrow(df_mu0)) records with |mu| ≤ $(mu_tol)"
    
    eigs_spectra = DataFrame()
    
    for group in groupby(df_mu0, [:phi, :Delta, :N, :t_n])
        if nrow(group) > 0
            first_row = first(group)
            
            # phason already extracted above, just get it from the group
            phason_val = first_row.phason
            
            # Collect all eigenvalues from this group
            all_eigs = Float64[]
            for row in eachrow(group)
                if !ismissing(row.eigenvalues) && !isempty(row.eigenvalues)
                    append!(all_eigs, collect(real(row.eigenvalues)))
                end
            end
            
            if !isempty(all_eigs)
                push!(eigs_spectra, (
                    phi = first_row.phi,
                    Delta = first_row.Delta,
                    N = first_row.N,
                    t_n = first_row.t_n,
                    phason = phason_val,
                    eigenvalues = unique(sort!(all_eigs))
                ))
            end
        end
    end
    
    @info "Extracted eigenvalue spectra for $(nrow(eigs_spectra)) (phi, Delta, N, t_n) combinations"
    if nrow(eigs_spectra) > 0
        @info "  phi range: $(extrema(eigs_spectra.phi))"
        @info "  Delta range: $(extrema(eigs_spectra.Delta))"
    end
    
    return eigs_spectra
end


"""
    compute_eigs_normalization_bounds_from_spectrum(eigs_spectra_df::DataFrame; eig_tol::Float64=1e-6)

Compute normalization bounds from eigenvalue spectra using IDOP approach.
For each (phi, Delta) pair, extract min/max of eigenvalues (filtering near-zero).
The maximum eigenvalue bandwidth (abs(max)) scales the mu axis multiplicatively.

Normalization: norm_mu = mu / abs(max_eigenvalue)
This uniformly scales all mu values where specloc_invariant=1 by a constant factor,
matching the mp_disc normalization approach without distorting data point positions.

Args:
    eigs_spectra_df: DataFrame from extract_eigs_spectrum_at_mu0()
    eig_tol: Tolerance for ignoring near-zero eigenvalues

Returns:
    DataFrame with columns: phi, Delta, N, t_n, phason, global_min_phi, global_max_phi
    (Compatible with apply_eigs_bandwidth_normalisation!)
"""
function compute_eigs_normalization_bounds_from_spectrum(eigs_spectra_df::DataFrame; eig_tol::Float64=1e-6)
    @info "Computing eigenvalue normalization bounds from spectra..."
    
    bounds_list = []
    
    for row in eachrow(eigs_spectra_df)
        eigs = row.eigenvalues
        
        if ismissing(eigs) || isempty(eigs)
            continue
        end
        
        # Filter to non-zero eigenvalues
        eigs_nonzero = filter(x -> abs(x) > eig_tol, eigs)
        
        if isempty(eigs_nonzero)
            @warn "No non-zero eigenvalues for phi=$(row.phi), Delta=$(row.Delta)"
            continue
        end
        
        # Compute min and max
        eig_min = minimum(eigs_nonzero)
        eig_max = maximum(eigs_nonzero)
        
        push!(bounds_list, (
            phi = row.phi,
            Delta = row.Delta,
            N = row.N,
            t_n = row.t_n,
            phason = row.phason,
            global_min_phi = eig_min,
            global_max_phi = eig_max
        ))
    end
    
    if isempty(bounds_list)
        error("Could not compute bounds from eigenvalue spectra")
    end
    
    bounds_df = DataFrame(bounds_list)
    
    @info "Computed bounds for $(nrow(bounds_df)) (phi, Delta) combinations"
    @info "  global_min_phi range: $(extrema(bounds_df.global_min_phi))"
    @info "  global_max_phi range: $(extrema(bounds_df.global_max_phi))"
    
    return bounds_df
end


"""
    extract_and_compute_eigs_bounds(butterfly_data_path::String; mu_tol::Float64=0.01)

Convenience function: extract eigenvalue spectra and compute bounds in one step.

Args:
    butterfly_data_path: Path to butterfly_data folder with BSON files
    mu_tol: Tolerance for mu≈0 selection

Returns:
    DataFrame with normalization bounds (compatible with apply_eigs_bandwidth_normalisation!)
"""
function extract_and_compute_eigs_bounds(butterfly_data_path::String; mu_tol::Float64=0.01)
    eigs_spectra = extract_eigs_spectrum_at_mu0(butterfly_data_path; mu_tol=mu_tol)
    bounds_df = compute_eigs_normalization_bounds_from_spectrum(eigs_spectra)
    return bounds_df
end


"""
    load_eigs_bandwidth_norms(csv_path::String)

Load the IDOP eigenvalue bandwidth cache and retain the mu=0 bounds per parameter set.
These are the bounds used by IDOP's outer normalization when it matches the cache
without using mu as a key.

Args:
    csv_path: Path to eigs_bandwidth_norms_generated.csv from IDOP results

Returns:
    DataFrame with columns: phi, Delta, phason, N, t_n, global_min_phi, global_max_phi

Example:
    bounds_df = load_eigs_bandwidth_norms("eigs_bandwidth_norms_generated.csv")
"""
function load_eigs_bandwidth_norms(csv_path::String; mu0_tol::Float64=1e-10)
    if !isfile(csv_path)
        error("Eigenvalue bandwidth CSV not found at: $csv_path")
    end
    
    @info "Loading eigenvalue bandwidth normalization data from: $csv_path"
    
    # Read CSV
    df_csv = CSV.read(csv_path, DataFrame)
    
    required_cols = [:phi, :phason, :mu, :Delta, :t_n, :N, :global_min, :global_max]
    missing_cols = filter(col -> !(col in propertynames(df_csv)), required_cols)
    isempty(missing_cols) || error("Normalization CSV missing columns: $missing_cols")

    # IDOP's lookup excludes mu from its matching key; the cache is ordered with
    # the mu=0 row first for each parameter set. Select that reference explicitly.
    df_mu0 = filter(row -> abs(row.mu) <= mu0_tol, df_csv)
    isempty(df_mu0) && error("Normalization CSV contains no mu=0 rows: $csv_path")

    bounds_df = select(
        df_mu0,
        :phi, :phason, :Delta, :N, :t_n,
        :global_min => :global_min_phi,
        :global_max => :global_max_phi
    )
    key_cols = [:phi, :phason, :Delta, :N, :t_n]
    if nrow(unique(bounds_df, key_cols)) != nrow(bounds_df)
        error("Normalization CSV has duplicate mu=0 rows for a parameter key")
    end
    
    @info "Loaded mu=0 eigenvalue bounds for $(nrow(bounds_df)) parameter combinations"
    @info "  phi range: $(extrema(bounds_df.phi))"
    @info "  global_min range: $(extrema(bounds_df.global_min_phi))"
    @info "  global_max range: $(extrema(bounds_df.global_max_phi))"
    
    return bounds_df
end


"""
        apply_eigs_bandwidth_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame)

Apply eigenvalue-bandwidth-based normalization to specloc data.
Uses the exact affine mapping in IDOP's `apply_outer_normalisation!` on the actual mu
coordinate for each parameter set. The invariant remains the plotted value; it must
not be multiplied into mu, since that collapses all zero-invariant samples to mu=0.

Args:
    df: Specloc data DataFrame with :mu, :phi, :Delta, and optionally :x0
    bounds_df: mu=0 eigenvalue bounds from load_eigs_bandwidth_norms()

Creates columns:
    :mu_specloc_invariant_norm - the IDOP-normalized mu coordinate

Returns:
    df with added :mu_specloc_invariant_norm column
"""
function apply_eigs_bandwidth_normalisation!(df::AbstractDataFrame, bounds_df::AbstractDataFrame)
    match_cols = [col for col in [:phi, :phason, :Delta, :N, :t_n]
                  if col in propertynames(df) && col in propertynames(bounds_df)]
    isempty(match_cols) && error("No shared parameter columns found for eigenvalue normalization")

    _canon(v) = if ismissing(v)
        "__MISSING__"
    elseif v isa Number
        string(round(Float64(v); digits=8))
    elseif v isa AbstractVector || v isa Tuple
        string(collect(v))
    else
        string(v)
    end
    _key(row, cols) = Tuple(_canon(getproperty(row, col)) for col in cols)

    bounds_dict = Dict{Tuple, Tuple{Float64, Float64}}()
    for row in eachrow(bounds_df)
        bounds_dict[_key(row, match_cols)] = (row.global_min_phi, row.global_max_phi)
    end

    norm_values = Vector{Union{Float64, Missing}}(undef, nrow(df))
    n_groups_processed = 0
    n_unmatched = 0

    for group in groupby(df, match_cols)
        key = _key(first(group), match_cols)
        indices = parentindices(group)[1]
        if haskey(bounds_dict, key)
            vmin, vmax = bounds_dict[key]
            if vmax == vmin
                norm_values[indices] .= 0.0
            else
                @. norm_values[indices] = 2 * (group.mu - vmin) / (vmax - vmin) - 1
            end
            n_groups_processed += 1
        else
            norm_values[indices] .= missing
            n_unmatched += length(indices)
        end
    end

    df[!, :mu_specloc_invariant_norm] = norm_values

    if n_unmatched > 0
        @warn "Could not match $n_unmatched rows to eigenvalue bounds" match_cols=match_cols
    end

    @info "Applied IDOP mu normalization to $(nrow(df)) rows via $n_groups_processed parameter groups"
    return df
end


end  # end module SpecLocProcessing
