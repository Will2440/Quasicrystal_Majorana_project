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
using Dates

project_root = @__DIR__

# =====================================================================
# User config (EDIT THESE)
# =====================================================================

sequences_dir = joinpath(project_root, "batch_params", "sequences")
sequence_bson_file ="hof_style_slopes_N500_phason_0.0-101-1.0_nbins100_npb1.bson"# "sturmian_slopes_K100_L4_balanced_bins5000_mpb2_r50_comp-false_tailoption-0_N500_phason_0.0-1-0.0_const_mapping.bson"

# Sequence selection
sequence_selection = :target_phason  # :all, :user_range, :target_slope, :target_phason, :slope_range
sequence_indices = 100:100  # only used if :user_range
target_slope = 2/(1 + sqrt(5))  # only used if :target_slope
slope_tolerance = 0.001  # only used if :target_slope
target_slope_range = (0.0, 0.5)  # only used if :slope_range
target_phason = 0.0  # only used if :target_phason
phason_tolerance = 0.0001  # only used if :target_phason

# Parameter ranges (full vectors; will be chunked per row)
Ns      = [500]                       # Vector{Int} (not chunked)
mus_all = collect(range(0.0, 3.0, 601)) #vcat(collect(range(0.6, 1.4, 401)), collect(range(1.8, 2.2, 201)))  # Vector{Float64}
Deltas  = [0.05] #, 0.1, 0.2, 0.3, 0.4, 0.5]             # Vector{Float64}
t1s     = [1.0]                  # Vector{Float64}
t2s     = [1.5] #, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]             # Vector{Float64}
t3s     = Float64[]                   # Vector{Float64} (empty = no t3)

# Chunking parameters for each
mu_vals_per_row = 301 #ceil(Int, length(mus_all)/30)
mu_split_mode = :contiguous  # :contiguous or :roundrobin
delta_vals_per_row = length(Deltas)
delta_split_mode = :contiguous
t1_vals_per_row = 1
t1_split_mode = :contiguous
t2_vals_per_row = length(t2s)
t2_split_mode = :contiguous
t3_vals_per_row = 1  # if t3s not empty
t3_split_mode = :contiguous
slope_vals_per_row = 1
phason_vals_per_row = 101
phason_split_mode = :contiguous  # :contiguous or :roundrobin
sequences_per_row = slope_vals_per_row * phason_vals_per_row

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

# Build filename with range summaries (shortened to avoid filename length limits)
ns_info     = isempty(Ns)      ? "none" : "$(Ns[1])-$(Ns[end])-$(length(Ns))"
mus_info    = isempty(mus_all) ? "none" : string(mus_all[1], "-", mus_all[end], "-", length(mus_all))
deltas_info = isempty(Deltas)  ? "none" : "$(Deltas[1])-$(Deltas[end])-$(length(Deltas))"
t1_info     = isempty(t1s)     ? "none" : string(t1s[1], "-", t1s[end], "-", length(t1s))
t2_info     = isempty(t2s)     ? "none" : string(t2s[1], "-", t2s[end], "-", length(t2s))
t3_info     = isempty(t3s)     ? "none" : "none"

# Sequence selection info for filename (shortened)
if sequence_selection == :all
    seq_sel_info = "_all"
elseif sequence_selection == :user_range
    seq_sel_info = "_r$(sequence_indices[1])-$(sequence_indices[end])"
elseif sequence_selection == :target_slope
    seq_sel_info = "_s$(target_slope)_t$(slope_tolerance)"
elseif sequence_selection == :slope_range
    seq_sel_info = "_srange$(target_slope_range[1])-$(target_slope_range[2])"
elseif sequence_selection == :target_phason
    seq_sel_info = "_p$(target_phason)_t$(phason_tolerance)"
else
    seq_sel_info = "_unk"
end

# Shortened chunk info
chunk_info = "_mu$(mu_vals_per_row)_d$(delta_vals_per_row)_t1$(t1_vals_per_row)_t2$(t2_vals_per_row)_slope$(slope_vals_per_row)_p$(phason_vals_per_row)"

# Create a timestamp-based unique identifier
timestamp = replace(string(now()), r"[^\d]" => "")[1:14]  # YYYYMMDDHHMMSS format

out_file = joinpath(outdir,
    "params_$(sequence_name)_$(timestamp)$(seq_sel_info)$(chunk_info)_Ns$(ns_info)_mus$(mus_info)_D$(deltas_info)_t1$(t1_info)_t2$(t2_info).dat")

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

total_sequences = length(seqs)  # Store total before filtering

# =====================================================================
# Sequence selection
# =====================================================================

if sequence_selection == :user_range
    max_idx = length(seqs)
    original_indices = sequence_indices
    sequence_indices = clamp.(sequence_indices, 1, max_idx)
    if original_indices != sequence_indices
        println("Warning: Some sequence indices were out of range. Valid range is 1:$(max_idx). Clamped to: $(sequence_indices)")
    end
    seqs = seqs[sequence_indices]
    phis = phis[sequence_indices]
    phasons = phasons[sequence_indices]
elseif sequence_selection == :target_slope
    mask = abs.(phis .- target_slope) .<= slope_tolerance
    seqs = seqs[mask]
    phis = phis[mask]
    phasons = phasons[mask]
elseif sequence_selection == :slope_range
    mask = (phis .>= target_slope_range[1]) .& (phis .<= target_slope_range[2])
    seqs = seqs[mask]
    phis = phis[mask]
    phasons = phasons[mask]
elseif sequence_selection == :target_phason
    mask = abs.(phasons .- target_phason) .<= phason_tolerance
    seqs = seqs[mask]
    phis = phis[mask]
    phasons = phasons[mask]
elseif sequence_selection != :all
    error("Unknown sequence_selection: $sequence_selection. Use :all, :user_range, :target_slope, :slope_range, or :target_phason")
end

# =====================================================================
# Generate sequence groups
# =====================================================================

# Sort sequences by slope (phi) then phason
sorted_indices = sortperm(collect(zip(phis, phasons)))
seqs = seqs[sorted_indices]
phis = phis[sorted_indices]
phasons = phasons[sorted_indices]

# Group into chunks of sequences_per_row
seq_groups = []
n_groups = ceil(Int, length(seqs) / sequences_per_row)
for g in 1:n_groups
    start_idx = (g-1)*sequences_per_row + 1
    end_idx = min(g*sequences_per_row, length(seqs))
    group_seqs = seqs[start_idx:end_idx]
    group_phis = phis[start_idx:end_idx]
    group_phasons = phasons[start_idx:end_idx]
    push!(seq_groups, (seqs=group_seqs, phis=group_phis, phasons=group_phasons))
end

# # Estimate output size and warn if too large
# estimated_rows = length(mu_chunks) * length(seq_groups)
# println("Estimated number of parameter rows: $estimated_rows")
# if estimated_rows > 1000
#     error("Estimated $(estimated_rows) rows would generate a very large file (disk quota exceeded). Please reduce parameter ranges (e.g., smaller mu range, fewer sequences, or larger chunks).")
# end

# =====================================================================
# Compute chunks for chunked variables
# =====================================================================

# Mus chunks
mu_chunks = []
if mu_vals_per_row < length(mus_all)
    n_mu_chunks = ceil(Int, length(mus_all) / mu_vals_per_row)
    for i in 1:n_mu_chunks
        start_idx = (i-1)*mu_vals_per_row + 1
        end_idx = min(i*mu_vals_per_row, length(mus_all))
        push!(mu_chunks, mus_all[start_idx:end_idx])
    end
else
    push!(mu_chunks, mus_all)
end

# Deltas chunks (not chunked)
delta_chunks = [Deltas]

# t1s chunks (not chunked)
t1_chunks = [t1s]

# t2s chunks (not chunked)
t2_chunks = [t2s]

# t3s chunks (not chunked)
t3_chunks = [t3s]

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

    for mu_chunk in mu_chunks
        for group in seq_groups
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

            # Use chunks for this row
            mus_entry = literal(mu_chunk)
            deltas_entry = literal(Deltas)
            t1s_entry = literal(t1s)
            t2s_entry = literal(t2s)
            t3s_entry = literal(t3s)

            row_vals = (
                literal(Ns),
                mus_entry,
                deltas_entry,
                t1s_entry,
                t2s_entry,
                t3s_entry,
                literal(sequence_names),
                literal(sequence_ids),
                literal(sequences)
            )
            write(io, join(row_vals, '\t') * "\n")
            nrows += 1
        end
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
println("Number of sequences: $(length(seqs)) out of $(total_sequences)")
if sequence_selection == :user_range
    println("Using filtering: $(sequence_selection) with indices: $(sequence_indices)")
elseif sequence_selection == :target_slope
    println("Using filtering: $(sequence_selection) with target slope: $(target_slope) ± $(slope_tolerance)")
elseif sequence_selection == :slope_range
    println("Using filtering: $(sequence_selection) with range: [$(target_slope_range[1]), $(target_slope_range[2])]")
elseif sequence_selection == :target_phason
    println("Using filtering: $(sequence_selection) with target phason: $(target_phason) ± $(phason_tolerance)")
elseif sequence_selection == :all
    println("Using filtering: $(sequence_selection)")
end
println("______________________________________________")
println("Chunking info:")
println("  mus: $mu_vals_per_row per row, mode: $mu_split_mode")
println("  Deltas: $delta_vals_per_row per row, mode: $delta_split_mode")
println("  t1s: $t1_vals_per_row per row, mode: $t1_split_mode")
println("  t2s: $t2_vals_per_row per row, mode: $t2_split_mode")
println("  t3s: $t3_vals_per_row per row, mode: $t3_split_mode")
println("  slopes: $slope_vals_per_row per row")
println("  phasons: $phason_vals_per_row per row, mode: $phason_split_mode")
println("  sequences: $sequences_per_row per row")
println("______________________________________________")
println("Number of rows: $(countlines(out_file) - 1)")
println("Combinations per row: $(mu_vals_per_row * length(Deltas) * length(t1s) * length(t2s) * (isempty(t3s) ? 1 : length(t3s)) * sequences_per_row)")
