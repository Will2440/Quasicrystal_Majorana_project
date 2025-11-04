"""
    file name:   hpc_compatible/serial_batch/param_prep.jl
    created:     04/11/2025
    last edited: 04/11/2025

    Prepares parameter .dat files for serial batch HPC simulations.
    Loads sequences from BSON files produced by cut_and_project.jl,
    and combines them with user-defined parameter ranges into a
    tab-separated values (TSV) file where each row corresponds to
    one simulation entry.

    Each row contains:
      - Ns: Vector{Int} of system sizes
      - mus: Vector{Float64} of chemical potentials
      - Deltas: Vector{Float64} of pairing strengths
      - t1s, t2s, t3s: Vector{Float64} of hopping amplitudes
      - sequence_name: String identifier of the sequence
      - sequence_id: Tuple{Float64,Float64} (phi, phason)
      - sequence: Vector{Int} representing the quasicrystal sequence

    The mus vector can be split into chunks per row to distribute
    chemical potential values across multiple jobs.

    Usage:
      - Edit the user config section below to set input/output paths
        and parameter ranges.
      - Run this script in Julia to generate the .dat file.
"""

using BSON
using DataFrames
using DelimitedFiles
using Printf

project_root = @__DIR__

# =====================================================================
# User (EDIT THESE)
# =====================================================================

# Input sequences (produced by cut_and_project.jl).
# Place the file in: serial_batch/batch_params/sequences/
sequences_dir = joinpath(project_root, "batch_params", "sequences")
sequence_bson_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0.bson"
indir = normpath(joinpath(sequences_dir, sequence_bson_file))

# Parameter ranges to embed per-row (raw vectors)
Ns      = [500]              # Vector{Int}
mus_all = collect(range(0.0, 3.0; length=301))  # Vector{Float64}
Deltas  = collect(range(0.0, 1.0, 3))                         # Vector{Float64}
t1s     = collect(range(1.0, 1.0, 1))                         # Vector{Float64}
t2s     = collect(range(1.5, 1.5, 1))                         # Vector{Float64}
t3s     = Float64[]                     # Vector{Float64} (empty means “not used”)

# Decide whether to split mus_all into chunks per row
n_mus_per_row = 10  # set to 1 to disable splitting
split_mode = :contiguous  # :contiguous or :roundrobin





# =====================================================================
# Helpers
# =====================================================================

# mu policy: define how mus_all can be distributed accross rows
function mu_window_policy(mus_full::Vector{Float64}, row_index::Int, nrows::Int;
                          n_chunks::Int=10, mode::Symbol=:contiguous)
    # split mus_full into n_chunks contiguous parts (last gets remainder)
    chunk_sizes = fill(length(mus_full) ÷ n_chunks, n_chunks)
    rem = length(mus_full) % n_chunks
    for i in 1:rem
        chunk_sizes[i] += 1
    end
    # compute chunk boundaries
    bounds = cumsum([0; chunk_sizes])
    if mode == :contiguous
        # map rows to chunks evenly: chunk_idx in 1..n_chunks
        chunk_idx = clamp(ceil(Int, row_index * n_chunks / max(1, nrows)), 1, n_chunks)
    elseif mode == :roundrobin
        chunk_idx = ((row_index - 1) % n_chunks) + 1
    else
        error("unknown mode: $mode")
    end
    a = bounds[chunk_idx] + 1
    b = bounds[chunk_idx+1]
    return mus_full[a:b]
end

# Output .dat path
ns_info = isempty(Ns) ? "none" : "$(Ns[1])-$(Ns[end])-$(length(Ns))"
mus_info = isempty(mus_all) ? "none" : string(replace(@sprintf("%.3g", mus_all[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.3g", mus_all[end]), "."=>"p"),
                                                "-", length(mus_all))
deltas_info = isempty(Deltas) ? "none" : "$(Deltas[1])-$(Deltas[end])-$(length(Deltas))"
t1_info = isempty(t1s) ? "none" : string(replace(@sprintf("%.3g", t1s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.3g", t1s[end]), "."=>"p"),
                                                "-", length(t1s))
t2_info = isempty(t2s) ? "none" : string(replace(@sprintf("%.3g", t2s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.3g", t2s[end]), "."=>"p"),
                                                "-", length(t2s))
t3_info = isempty(t3s) ? "none" : string(replace(@sprintf("%.3g", t3s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.3g", t3s[end]), "."=>"p"),
                                                "-", length(t3s))
out_folder = joinpath(project_root, "batch_params/param_sets")
isdir(out_folder) || mkpath(out_folder)
sequence_name = splitext(sequence_bson_file)[1]
out_filename = "params_$(sequence_name)_Ns_$(ns_info)_mus_$(mus_info)_Deltas_$(deltas_info)_t1_$(t1_info)_t2_$(t2_info)_t3_$(t3_info)_nmus_$(n_mus_per_row).dat"
outdir = joinpath(out_folder, "params_$(sequence_name)_ranges.dat")

# Load sequences from cut_and_project bson
function load_sturmian_seq_bson(path::AbstractString)
    @assert isfile(path) "BSON file not found: $path"
    raw = BSON.load(path)
    haskey(raw, :seqs)  || error("BSON missing :seqs")
    haskey(raw, :phis)  || error("BSON missing :phis")
    # :phasons optional in older files
    seqs    = Vector{Vector{Int}}(raw[:seqs])
    phis    = Float64.(raw[:phis])
    phasons = haskey(raw, :phasons) ? Float64.(raw[:phasons]) : fill(0.0, length(phis))
    length(seqs) == length(phis) == length(phasons) || error("Mismatched lengths in BSON arrays")
    return DataFrame(seq = seqs, phi = phis, phason = phasons)
end

# =====================================================================
# Generate param .dat
# =====================================================================

println("Loading sequences from: $indir")
sturm_df = load_sturmian_seq_bson(indir)
println("Loaded $(nrow(sturm_df)) sequences")

open(outdir, "w") do io
    # Header
    header = join((
        "Ns", "mus", "Deltas", "t1s", "t2s", "t3s",
        "sequence_name", "sequence_id", "sequence"
    ), '\t')
    write(io, header * "\n")

    nrows = 0
    total_rows = nrow(sturm_df)
    for (row_idx, row) in enumerate(eachrow(sturm_df))
        phi    = Float64(row.phi)
        phason = Float64(row.phason)
        seq    = Vector{Int}(row.seq)

        # Per-entry mus (subset of global mus_all, if desired)
        mus_entry = mu_window_policy(mus_all, row_idx, total_rows; n_chunks=n_mus_per_row, mode=split_mode)

        row_vals = (
            literal(Ns),
            literal(mus_entry),
            literal(Deltas),
            literal(t1s),
            literal(t2s),
            literal(t3s),
            sequence_name,
            literal((phi, phason)),
            literal(seq)
        )
        write(io, join(row_vals, '\t') * "\n")
        nrows += 1
    end
    println("Wrote $nrows rows to $outdir")
end