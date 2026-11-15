"""
    diagnose_bounds.jl
    
Diagnostic script to verify eigenvalue bounds extraction and matching against specloc data.
"""

using JLD2
using DataFrames
using ProgressMeter
using Logging
using Statistics

level = Logging.Info
global_logger(ConsoleLogger(stderr, level))

include("processing.jl")
include("plotting.jl")

using .SpecLocProcessing
using .SpecLocPlotting


# Load specloc data
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p0-1201)_Delta(0p05-0p05-1)_x0(218p75-281p25-5)_kappa(0p0001-0p0001-1)"
folder_path_specloc = normpath(joinpath(@__DIR__, "..", "..", "raw_data", "specloc_butterfly", data_folder))

@info "Loading specloc data..."
dataframes = SpecLocProcessing.load_chunk_files(folder_path_specloc)
specloc_df = SpecLocProcessing.extract_specloc_data(dataframes)
agg_by_x0_df = SpecLocProcessing.aggregate_specloc_data(specloc_df; by_mean=true, average_over_x0=false)
agg_by_x0_df = SpecLocProcessing.mirror_mu_data(agg_by_x0_df)

@info "Specloc data parameters:"
@info "  Unique phi values: $(length(unique(agg_by_x0_df.phi)))"
@info "  Unique Delta values: $(unique(agg_by_x0_df.Delta))"
@info "  Unique N values: $(unique(agg_by_x0_df.N))"
@info "  Unique x0 values: $(sort(unique(agg_by_x0_df.x0)))"
@info "  Unique mu values: $(length(unique(agg_by_x0_df.mu)))"
@info "  mu range: $(extrema(agg_by_x0_df.mu))"
@info "  phi range: $(extrema(agg_by_x0_df.phi))"

# Extract eigenvalue bounds
butterfly_data_path = normpath(joinpath(@__DIR__, "..", "..", "raw_data", "butterfly_data", 
    "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)_FINALDATASHARED"))

@info "\nExtracting eigenvalue bounds..."
eigs_bounds_df = SpecLocProcessing.extract_and_compute_eigs_bounds(butterfly_data_path; mu_tol=0.01)

@info "Eigenvalue bounds parameters:"
@info "  Number of (phi, Delta) combinations: $(nrow(eigs_bounds_df))"
@info "  Unique phi values: $(length(unique(eigs_bounds_df.phi)))"
@info "  Unique Delta values: $(unique(eigs_bounds_df.Delta))"
@info "  Unique N values: $(unique(eigs_bounds_df.N))"
@info "  phi range in bounds: $(extrema(eigs_bounds_df.phi))"
@info "  Delta range in bounds: $(extrema(eigs_bounds_df.Delta))"
@info "  global_min_phi range: $(extrema(eigs_bounds_df.global_min_phi))"
@info "  global_max_phi range: $(extrema(eigs_bounds_df.global_max_phi))"

# Check for matching
@info "\nMatching analysis:"
specloc_phi_delta = Set([(row.phi, row.Delta) for row in eachrow(agg_by_x0_df)])
bounds_phi_delta = Set([(row.phi, row.Delta) for row in eachrow(eigs_bounds_df)])

matched = intersect(specloc_phi_delta, bounds_phi_delta)
only_in_specloc = setdiff(specloc_phi_delta, bounds_phi_delta)
only_in_bounds = setdiff(bounds_phi_delta, specloc_phi_delta)

@info "  Matched (phi, Delta) pairs: $(length(matched))"
@info "  Only in specloc: $(length(only_in_specloc))"
@info "  Only in bounds: $(length(only_in_bounds))"

if length(only_in_specloc) > 0
    @warn "Missing bounds for some specloc (phi, Delta) pairs:"
    for (phi, Delta) in collect(only_in_specloc)[1:min(5, length(only_in_specloc))]
        @info "    phi=$phi, Delta=$Delta"
    end
end

# Sample some bounds to check they look reasonable
@info "\nSample eigenvalue bounds (first 5):"
for i in 1:min(5, nrow(eigs_bounds_df))
    row = eigs_bounds_df[i, :]
    @info "  phi=$(row.phi), Delta=$(row.Delta): min=$(row.global_min_phi), max=$(row.global_max_phi), span=$(row.global_max_phi - row.global_min_phi)"
end

# Check scaling factors
@info "\nScaling factors (abs(max_eigenvalue)):"
scale_factors = abs.(eigs_bounds_df.global_max_phi)
@info "  Range: $(extrema(scale_factors))"
@info "  Mean: $(mean(scale_factors))"
@info "  Std: $(std(scale_factors))"

# Check what mu values would look like after normalization (sample)
@info "\nSample normalized mu values (for matched pairs):"
sample_indices = findall(row -> (row.phi, row.Delta) in matched, eachrow(agg_by_x0_df))[1:min(20, sum(1 for row in eachrow(agg_by_x0_df) if (row.phi, row.Delta) in matched))]
for idx in sample_indices
    row = agg_by_x0_df[idx, :]
    # Find corresponding bounds
    bounds_rows = filter(r -> r.phi ≈ row.phi && r.Delta ≈ row.Delta, eigs_bounds_df)
    if nrow(bounds_rows) > 0
        bounds_row = bounds_rows[1, :]
        mu_eff = row.mu * row.specloc_invariant
        mu_norm = mu_eff / abs(bounds_row.global_max_phi)
        @info "  phi=$(round(row.phi, digits=3)), mu=$(round(row.mu, digits=3)), specloc_inv=$(row.specloc_invariant), mu_eff=$(round(mu_eff, digits=3)), mu_norm=$(round(mu_norm, digits=3))"
    end
end

@info "\n✓ Diagnostic complete"
