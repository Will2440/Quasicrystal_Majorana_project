"""
    check_specloc_distribution.jl
    
Check where specloc_invariant=1 and how mu values cluster.
"""

using JLD2
using DataFrames
using ProgressMeter
using Logging

level = Logging.Info
global_logger(ConsoleLogger(stderr, level))

include("processing.jl")

using .SpecLocProcessing

# Load specloc data
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p0-1201)_Delta(0p05-0p05-1)_x0(218p75-281p25-5)_kappa(0p0001-0p0001-1)"
folder_path_specloc = normpath(joinpath(@__DIR__, "..", "..", "raw_data", "specloc_butterfly", data_folder))

@info "Loading specloc data..."
dataframes = SpecLocProcessing.load_chunk_files(folder_path_specloc)
specloc_df = SpecLocProcessing.extract_specloc_data(dataframes)
agg_by_x0_df = SpecLocProcessing.aggregate_specloc_data(specloc_df; by_mean=true, average_over_x0=false)
agg_by_x0_df = SpecLocProcessing.mirror_mu_data(agg_by_x0_df)

# Check distribution of specloc_invariant
@info "Specloc invariant distribution:"
inv_0_count = sum(agg_by_x0_df.specloc_invariant .== 0.0)
inv_1_count = sum(agg_by_x0_df.specloc_invariant .== 1.0)
@info "  Count with invariant=0: $inv_0_count"
@info "  Count with invariant=1: $inv_1_count"
@info "  Ratio: $(round(inv_1_count / nrow(agg_by_x0_df) * 100, digits=2))% have invariant=1"

# For a few different phi values, show where invariant=1
@info "\nFor selected phi values, show mu distribution where invariant=1:"
unique_phis = sort(unique(agg_by_x0_df.phi))
for phi_idx in [1, div(length(unique_phis), 4), div(length(unique_phis), 2), div(3*length(unique_phis), 4), length(unique_phis)]
    phi = unique_phis[phi_idx]
    phi_data = filter(row -> row.phi ≈ phi, agg_by_x0_df)
    invariant_1_data = filter(row -> row.specloc_invariant == 1.0, phi_data)
    
    if nrow(invariant_1_data) > 0
        mu_vals = sort(unique(invariant_1_data.mu))
        @info "  phi=$(round(phi, digits=3)): $(nrow(invariant_1_data)) points with invariant=1"
        @info "    mu range: $(extrema(mu_vals))"
        @info "    mu values: $(mu_vals[1:min(10, length(mu_vals))])" 
        if length(mu_vals) > 10
            @info "    ... ($(length(mu_vals)-10) more)"
        end
    else
        @info "  phi=$(round(phi, digits=3)): NO points with invariant=1"
    end
end

@info "\n✓ Distribution check complete"
