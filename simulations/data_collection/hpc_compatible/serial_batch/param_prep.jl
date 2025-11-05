"""
Param prep for serial_batch runs:
- Builds a .dat with one row per parameter combination group (grouped sequences × chunked params).
- Each row gets sub-ranges for mus, Deltas, t1s, t2s, t3s via policies; groups sequences.
- Columns: Ns, mus, Deltas, t1s, t2s, t3s, sequence_names, sequence_ids, sequences
- sequence_names: Vector{String}, sequence_ids: Vector{Tuple{Float64,Float64}}, sequences: Vector{Vector{Int}}
"""

using BSON
using DataFrames
using DelimitedFiles
using Printf

project_root = @__DIR__

# =====================================================================
# User config (EDIT THESE)
# =====================================================================

sequences_dir = joinpath(project_root, "batch_params", "sequences")
sequence_bson_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0.bson"

# Parameter ranges (full vectors; will be chunked per row)
Ns      = [500]                       # Vector{Int} (not chunked)
mus_all = collect(range(0.0, 3.0; length=601))  # Vector{Float64}
Deltas  = [0.0] #[0.0, 0.01, 0.05, 0.1, 0.2]             # Vector{Float64}
t1s     = [1.0]                  # Vector{Float64}
t2s     = [1.5]#[1.5, 2.0, 2.5]             # Vector{Float64}
t3s     = Float64[]                   # Vector{Float64} (empty = no t3)

# Chunking parameters for each
mu_vals_per_row = 10
mu_split_mode = :contiguous  # :contiguous or :roundrobin
delta_vals_per_row = length(Deltas)
delta_split_mode = :contiguous
t1_vals_per_row = 1
t1_split_mode = :contiguous
t2_vals_per_row = length(t2s)
t2_split_mode = :contiguous
t3_vals_per_row = 1  # if t3s not empty
t3_split_mode = :contiguous
n_sequences_per_row = 1
seq_split_mode = :contiguous

# =====================================================================
# Helpers
# =====================================================================

# Generic chunking policy (now takes vals_per_chunk: number of values per chunk/row)
function chunk_policy(full_vec::Vector, row_index::Int, nrows::Int;
                      vals_per_chunk::Int, mode::Symbol=:contiguous)
    isempty(full_vec) && return full_vec
    n_chunks = ceil(Int, length(full_vec) / vals_per_chunk)
    chunk_sizes = fill(vals_per_chunk, n_chunks)
    # Adjust last chunk if remainder
    total_assigned = (n_chunks - 1) * vals_per_chunk
    chunk_sizes[end] = length(full_vec) - total_assigned
    bounds = cumsum([0; chunk_sizes])
    if mode == :contiguous
        chunk_idx = clamp(ceil(Int, row_index * n_chunks / max(1, nrows)), 1, n_chunks)
    elseif mode == :roundrobin
        chunk_idx = ((row_index - 1) % n_chunks) + 1
    else
        error("unknown mode: $mode")
    end
    a = bounds[chunk_idx] + 1
    b = bounds[chunk_idx+1]
    return full_vec[a:b]
end

# Sequence grouping policy (unchanged)
function sequence_group_policy(seqs::Vector{Vector{Int}}, phis::Vector{Float64}, phasons::Vector{Float64}, row_index::Int, nrows::Int;
                               n_groups::Int, mode::Symbol=:contiguous)
    n_seqs = length(seqs)
    group_sizes = fill(n_seqs ÷ n_groups, n_groups)
    rem = n_seqs % n_groups
    for i in 1:rem
        group_sizes[i] += 1
    end
    bounds = cumsum([0; group_sizes])
    if mode == :contiguous
        group_idx = clamp(ceil(Int, row_index * n_groups / max(1, nrows)), 1, n_groups)
    elseif mode == :roundrobin
        group_idx = ((row_index - 1) % n_groups) + 1
    else
        error("unknown mode: $mode")
    end
    a = bounds[group_idx] + 1
    b = bounds[group_idx+1]
    return seqs[a:b], phis[a:b], phasons[a:b]
end

# Output
outdir = joinpath(project_root, "batch_params", "param_sets")
isdir(outdir) || mkpath(outdir)
sequence_name = splitext(sequence_bson_file)[1]

# Build filename with range summaries
ns_info     = isempty(Ns)      ? "none" : "$(Ns[1])-$(Ns[end])-$(length(Ns))"
mus_info    = isempty(mus_all) ? "none" : string(replace(@sprintf("%.6g", mus_all[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.6g", mus_all[end]), "."=>"p"),
                                                "-", length(mus_all))
deltas_info = isempty(Deltas)  ? "none" : "$(Deltas[1])-$(Deltas[end])-$(length(Deltas))"
t1_info     = isempty(t1s)     ? "none" : string(replace(@sprintf("%.6g", t1s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.6g", t1s[end]), "."=>"p"),
                                                "-", length(t1s))
t2_info     = isempty(t2s)     ? "none" : string(replace(@sprintf("%.6g", t2s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.6g", t2s[end]), "."=>"p"),
                                                "-", length(t2s))
t3_info     = isempty(t3s)     ? "none" : string(replace(@sprintf("%.6g", t3s[1]), "."=>"p"),
                                                "-", replace(@sprintf("%.6g", t3s[end]), "."=>"p"),
                                                "-", length(t3s))

out_file = joinpath(outdir,
    "params_$(sequence_name)_Ns_$(ns_info)_mus_$(mus_info)_Deltas_$(deltas_info)_t1_$(t1_info)_t2_$(t2_info)_t3_$(t3_info)_nseq_$(n_sequences_per_row).dat")

# =====================================================================
# Load sequences
# =====================================================================

sequence_bson_path = normpath(joinpath(sequences_dir, sequence_bson_file))
@assert isfile(sequence_bson_path) "BSON file not found: $sequence_bson_path"

raw = BSON.load(sequence_bson_path)
@assert haskey(raw, :seqs) "BSON missing :seqs"
@assert haskey(raw, :phis) "BSON missing :phis"
seqs = Vector{Vector{Int}}(raw[:seqs])
phis = Float64.(raw[:phis])
phasons = haskey(raw, :phasons) ? Float64.(raw[:phasons]) : fill(0.0, length(phis))
length(seqs) == length(phis) == length(phasons) || error("Mismatched lengths in BSON arrays")

# =====================================================================
# Generate sequence groups
# =====================================================================

n_seq_groups = ceil(Int, length(seqs) / n_sequences_per_row)
seq_groups = []
for g in 1:n_seq_groups
    group_seqs, group_phis, group_phasons = sequence_group_policy(seqs, phis, phasons, g, n_seq_groups; n_groups=n_seq_groups, mode=seq_split_mode)
    push!(seq_groups, (seqs=group_seqs, phis=group_phis, phasons=group_phasons))
end

# =====================================================================
# Helpers
# =====================================================================

literal(x) = repr(x)

# =====================================================================
# Generate rows and write .dat
# =====================================================================

open(out_file, "w") do io
    header = join(("Ns", "mus", "Deltas", "t1s", "t2s", "t3s", "sequence_names", "sequence_ids", "sequences"), '\t')
    write(io, header * "\n")
    nrows = 0
    total_rows = length(seq_groups)

    for (row_idx, group) in enumerate(seq_groups)
        group_seqs = group.seqs
        group_phis = group.phis
        group_phasons = group.phasons

        # Truncate sequences to max(Ns) for this row
        max_N = maximum(Ns)
        truncated_sequences = Vector{Vector{Int}}()
        for seq in group_seqs
            if length(seq) < max_N
                error("Sequence length $(length(seq)) < max_N=$max_N for row $row_idx")
            end
            push!(truncated_sequences, seq[1:max_N])
        end

        # Collect sequence data for this group
        sequence_names = fill(sequence_name, length(group_seqs))  # Vector{String}
        sequence_ids = [(group_phis[i], group_phasons[i]) for i in 1:length(group_seqs)]  # Vector{Tuple{Float64,Float64}}
        sequences = truncated_sequences  # Vector{Vector{Int}} (truncated to max_N)

        # Chunk params for this row
        mus_entry = chunk_policy(mus_all, row_idx, total_rows; vals_per_chunk=mu_vals_per_row, mode=mu_split_mode)
        deltas_entry = chunk_policy(Deltas, row_idx, total_rows; vals_per_chunk=delta_vals_per_row, mode=delta_split_mode)
        t1s_entry = chunk_policy(t1s, row_idx, total_rows; vals_per_chunk=t1_vals_per_row, mode=t1_split_mode)
        t2s_entry = chunk_policy(t2s, row_idx, total_rows; vals_per_chunk=t2_vals_per_row, mode=t2_split_mode)
        t3s_entry = chunk_policy(t3s, row_idx, total_rows; vals_per_chunk=t3_vals_per_row, mode=t3_split_mode)

        row_vals = (
            literal(Ns),
            literal(mus_entry),
            literal(deltas_entry),
            literal(t1s_entry),
            literal(t2s_entry),
            literal(t3s_entry),
            literal(sequence_names),
            literal(sequence_ids),
            literal(sequences)
        )
        write(io, join(row_vals, '\t') * "\n")
        nrows += 1
    end
    # println("Wrote $nrows rows to $out_file")
end

println("Done. Param file at: $out_file")
println("Ns (start, end, length): $(Ns[1]), $(Ns[end]), $(length(Ns))")
println("mus (start, end, length): $(mus_all[1]), $(mus_all[end]), $(length(mus_all))")
println("Deltas (start, end, length): $(Deltas[1]), $(Deltas[end]), $(length(Deltas))")
println("t1s (start, end, length): $(t1s[1]), $(t1s[end]), $(length(t1s))")
println("t2s (start, end, length): $(t2s[1]), $(t2s[end]), $(length(t2s))")
println("t3s (start, end, length): $(isempty(t3s) ? "none" : string(t3s[1], ", ", t3s[end], ", ", length(t3s)))")
println("Number of sequences: $(length(seq_groups))")
println("______________________________________________")
println("Number of rows: $(countlines(out_file) - 1)")
println("Combinations per row: $(mu_vals_per_row*delta_vals_per_row*t1_vals_per_row*t2_vals_per_row*max(1, t3_vals_per_row)*n_sequences_per_row)")