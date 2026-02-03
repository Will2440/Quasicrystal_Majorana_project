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

function log_range(base::Float64, start_power::Float64, end_power::Float64, length::Int)
    powers = range(start_power, stop=end_power, length=length)
    values = base .^ powers
    return values
end

# =====================================================================
# User config (EDIT THESE)
# =====================================================================

sequences_dir = joinpath(project_root, "batch_params", "sequences")
sequence_bson_file ="hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-501-1.0_nbins1000_npb1.bson"

# Sequence selection
sequence_selection = :all  # :all, :user_range, :target_slope, :target_phason
sequence_indices = 100:100  # only used if :user_range
target_slope =  2/(1 + sqrt(5))  # only used if :target_slope
slope_tolerance = 0.001  # only used if :target_slope
target_phason = 0.0  # only used if :target_phason
phason_tolerance = 0.0001  # only used if :target_phason

# Parameter ranges (full vectors; will be chunked per row)
Ns      = [1000] #[5,10,50,100,500,1000,5000,10000]                       # Vector{Int} (not chunked)
mus_all = [0.0] #collect(range(0.0, 3.0, 1201)) #vcat(collect(range(0.6, 1.4, 401)), collect(range(1.8, 2.2, 201)))  # Vector{Float64}
Deltas  = [0.05] #[0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]             # Vector{Float64}
t1s     = [1.0]                  # Vector{Float64}
t2s     = [1.5] #[1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0] #[1.5] #, 3.0, 3.5, 4.0, 4.5, 5.0]             # Vector{Float64}
t3s     = Float64[]                   # Vector{Float64} (empty = no t3)

# Auto-mu configuration
mu_range_mode = :user_defined       # :user_defined or :auto_mu
auto_mu_full_range = false           # true: -mu_max to mu_max, false: 0 to mu_max
auto_mu_preserve_resolution = true  # true: fixed mu resolution (variable number of mu points), false: fixed length(mus) (variable mu step size)
auto_mu_resolution_param = 0.001    # Target step size (resolution)

# Smart grouping configuration
smart_sequence_grouping = false
target_number_of_rows = 200 # Target number of rows (jobs) to generate
###########################################################################################################################
## N.B. when using :auto_mu mode, ensure slope_vals_per_row is small, as one mu range will be taken for all these slopes ##
###########################################################################################################################

# Chunking parameters for each
N_vals_per_row = length(Ns)
N_split_mode = :contiguous
mu_vals_per_row = 1 #(length(mus_all)÷2) + 1 # This can safely exceed the number of mu values per row when using :auto_mu mode
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
phason_vals_per_row = 51
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
if mu_range_mode == :auto_mu
    mus_info = "auto_$(auto_mu_full_range ? "full" : "half")_$(auto_mu_preserve_resolution ? "presRes" : "fixedN")_res$(auto_mu_resolution_param)"
else
    mus_info = isempty(mus_all) ? "none" : string(mus_all[1], "-", mus_all[end], "-", length(mus_all))
end
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
chunk_info = "_Ns$(N_vals_per_row)_mu$(mu_vals_per_row)_d$(delta_vals_per_row)_t1$(t1_vals_per_row)_t2$(t2_vals_per_row)_slope$(slope_vals_per_row)_p$(phason_vals_per_row)"

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

# Group into chunks
seq_groups = []

if smart_sequence_grouping && mu_range_mode == :auto_mu
    
    max_t2_global = isempty(t2s) ? 1.0 : maximum(t2s)
    
    # Pre-calculate global N if needed for fallback
    max_phi_global = maximum(phis)
    mu_max_global = max_t2_global * (1.5 + max_phi_global)
    mu_min_global = auto_mu_full_range ? -mu_max_global : 0.0
    range_width_global = mu_max_global - mu_min_global
    global_N_fallback = ceil(Int, range_width_global / auto_mu_resolution_param) + 1

    # Calculate total combinations first to determine target per row
    total_combinations_global = 0
    for i in 1:length(seqs)
        local_phi = phis[i]
        local_mu_max = max_t2_global * (1.5 + local_phi)
        local_mu_min = auto_mu_full_range ? -local_mu_max : 0.0
        
        if auto_mu_preserve_resolution
             local_n_points = ceil(Int, (local_mu_max - local_mu_min) / auto_mu_resolution_param) + 1
        else
             local_n_points = global_N_fallback
        end
        global total_combinations_global += local_n_points
    end

    target_combinations_per_row = ceil(Int, total_combinations_global / target_number_of_rows)
    println("Using smart sequence grouping.")
    println("  Total combinations across all sequences: $total_combinations_global")
    println("  Target rows: $target_number_of_rows")
    println("  Calculated target combinations per row: $target_combinations_per_row")

    current_seqs = Vector{Vector{Int}}()
    current_phis = Vector{Float64}()
    current_phasons = Vector{Float64}()

    for i in 1:length(seqs)
        # Add candidate
        push!(current_seqs, seqs[i])
        push!(current_phis, phis[i])
        push!(current_phasons, phasons[i])
        
        # Calculate expected mu points for this group (driven by max slope in group)
        # Since sorted by phi, the current one (or last one) has the max phi
        local_max_phi = maximum(current_phis) 
        
        local_mu_max = max_t2_global * (1.5 + local_max_phi)
        local_mu_min = auto_mu_full_range ? -local_mu_max : 0.0
        
        if auto_mu_preserve_resolution
             local_n_points = ceil(Int, (local_mu_max - local_mu_min) / auto_mu_resolution_param) + 1
        else
             local_n_points = global_N_fallback
        end
        
        current_combinations = length(current_seqs) * local_n_points
        
        if current_combinations >= target_combinations_per_row
            # Finalize group
            push!(seq_groups, (seqs=copy(current_seqs), phis=copy(current_phis), phasons=copy(current_phasons)))
            empty!(current_seqs)
            empty!(current_phis)
            empty!(current_phasons)
        end
    end
    # Push remaining
    if !isempty(current_seqs)
        push!(seq_groups, (seqs=current_seqs, phis=current_phis, phasons=current_phasons))
    end

else
    # Standard fixed chunking
    n_groups = ceil(Int, length(seqs) / sequences_per_row)
    for g in 1:n_groups
        start_idx = (g-1)*sequences_per_row + 1
        end_idx = min(g*sequences_per_row, length(seqs))
        group_seqs = seqs[start_idx:end_idx]
        group_phis = phis[start_idx:end_idx]
        group_phasons = phasons[start_idx:end_idx]
        push!(seq_groups, (seqs=group_seqs, phis=group_phis, phasons=group_phasons))
    end
end

println("DEBUG: Generated $(length(seq_groups)) sequence groups.")

# =====================================================================
# Compute chunks for chunked variables
# =====================================================================

# Deltas chunks
delta_chunks = []
if delta_vals_per_row < length(Deltas)
    n_chunks = ceil(Int, length(Deltas) / delta_vals_per_row)
    for i in 1:n_chunks
        s = (i-1)*delta_vals_per_row + 1
        e = min(i*delta_vals_per_row, length(Deltas))
        push!(delta_chunks, Deltas[s:e])
    end
else
    push!(delta_chunks, Deltas)
end

# t1s chunks
t1_chunks = []
if t1_vals_per_row < length(t1s)
    n_chunks = ceil(Int, length(t1s) / t1_vals_per_row)
    for i in 1:n_chunks
        s = (i-1)*t1_vals_per_row + 1
        e = min(i*t1_vals_per_row, length(t1s))
        push!(t1_chunks, t1s[s:e])
    end
else
    push!(t1_chunks, t1s)
end

# t2s chunks
t2_chunks = []
if t2_vals_per_row < length(t2s)
    n_chunks = ceil(Int, length(t2s) / t2_vals_per_row)
    for i in 1:n_chunks
        s = (i-1)*t2_vals_per_row + 1
        e = min(i*t2_vals_per_row, length(t2s))
        push!(t2_chunks, t2s[s:e])
    end
else
    push!(t2_chunks, t2s)
end

# t3s chunks
t3_chunks = []
if t3_vals_per_row < length(t3s)
    n_chunks = ceil(Int, length(t3s) / t3_vals_per_row)
    for i in 1:n_chunks
        s = (i-1)*t3_vals_per_row + 1
        e = min(i*t3_vals_per_row, length(t3s))
        push!(t3_chunks, t3s[s:e])
    end
else
    push!(t3_chunks, t3s)
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

    # Pre-calculate global N for fixed N mode (if needed)
    global_N = 0
    if mu_range_mode == :auto_mu
        max_phi_global = maximum(phis)
        max_t2_global = isempty(t2s) ? 1.0 : maximum(t2s)
        mu_max_global = max_t2_global * (1.5 + max_phi_global)
        mu_min_global = auto_mu_full_range ? -mu_max_global : 0.0
        range_width_global = mu_max_global - mu_min_global
        
        # Calculate N corresponding to the resolution param on the largest range
        global_N = ceil(Int, range_width_global / auto_mu_resolution_param) + 1
    end

    # Iterate over sequence groups first to allow per-group mu ranges
    for (g_idx, group) in enumerate(seq_groups)
        if g_idx == 1
             println("DEBUG: Processing group 1")
        end
        group_seqs = group.seqs
        group_phis = group.phis
        group_phasons = group.phasons

        # Determine mus for this group
        if mu_range_mode == :auto_mu
            # Calculate max slope in this group to determine required bandwidth
            # Bandwidth B(r) = t2 * (1.5 + r) where r = slope
            # We use the maximum t2 and maximum slope in the group to ensure coverage
            max_phi = maximum(group_phis)
            max_t2 = isempty(t2s) ? 1.0 : maximum(t2s) # Default to 1.0 if t2s empty (unlikely)
            
            mu_max = max_t2 * (1.5 + max_phi)
            mu_min = auto_mu_full_range ? -mu_max : 0.0
            
            if auto_mu_preserve_resolution
                # Use resolution param as step size
                n_points = ceil(Int, (mu_max - mu_min) / auto_mu_resolution_param) + 1
                group_mus = collect(range(mu_min, mu_max, length=n_points))
            else
                # Use fixed N derived from global max range
                n_points = global_N
                group_mus = collect(range(mu_min, mu_max, length=n_points))
            end
        else
            group_mus = mus_all
        end

        # Chunk group_mus
        local_mu_chunks = []
        if mu_vals_per_row < length(group_mus)
            n_chunks = ceil(Int, length(group_mus) / mu_vals_per_row)
            for i in 1:n_chunks
                s = (i-1)*mu_vals_per_row + 1
                e = min(i*mu_vals_per_row, length(group_mus))
                push!(local_mu_chunks, group_mus[s:e])
            end
        else
            push!(local_mu_chunks, group_mus)
        end
        
        if g_idx == 1
            println("DEBUG: group_mus length: $(length(group_mus))")
            println("DEBUG: local_mu_chunks length: $(length(local_mu_chunks))")
        end

        for (chunk_idx, mu_chunk) in enumerate(local_mu_chunks)
            
            # Sub-chunk Ns
            local_N_chunks = []
            if N_vals_per_row < length(Ns)
                n_n_chunks = ceil(Int, length(Ns) / N_vals_per_row)
                for i in 1:n_n_chunks
                    s = (i-1)*N_vals_per_row + 1
                    e = min(i*N_vals_per_row, length(Ns))
                    push!(local_N_chunks, Ns[s:e])
                end
            else
                push!(local_N_chunks, Ns)
            end

            for (n_chunk_idx, N_chunk) in enumerate(local_N_chunks)

                if g_idx == 1 && chunk_idx == 1 && n_chunk_idx == 1
                    println("DEBUG: Processing first chunk. N_chunk: $N_chunk")
                end
                
                # Truncate sequences to max(Ns) for this row
                max_N = maximum(N_chunk)
                truncated_sequences = Vector{Vector{Int}}()
                for seq in group_seqs
                    if length(seq) < max_N
                        error("Sequence length $(length(seq)) < max_N=$max_N for row $nrows")
                    end
                    push!(truncated_sequences, seq[1:max_N])
                end

                if g_idx == 1 && chunk_idx == 1 && n_chunk_idx == 1
                    println("DEBUG: Sequences truncated.")
                end

                # Collect sequence data for this group
                sequence_names = fill(sequence_name, length(group_seqs))  # Vector{String}
                sequence_ids = [(group_phis[i], group_phasons[i]) for i in 1:length(group_seqs)]  # Vector{Tuple{Float64,Float64}}
                sequences = truncated_sequences  # Vector{Vector{Int}} (truncated to max_N)

                for delta_chunk in delta_chunks
                    for t1_chunk in t1_chunks
                        for t2_chunk in t2_chunks
                            for t3_chunk in t3_chunks

                                # Use chunks for this row
                                row_vals = (
                                    literal(N_chunk),
                                    literal(mu_chunk),
                                    literal(delta_chunk),
                                    literal(t1_chunk),
                                    literal(t2_chunk),
                                    literal(t3_chunk),
                                    literal(sequence_names),
                                    literal(sequence_ids),
                                    literal(sequences)
                                )
                                
                                if g_idx == 1 && chunk_idx == 1 && n_chunk_idx == 1
                                    println("DEBUG: Row vals prepared. Writing...")
                                end
                                
                                write(io, join(row_vals, '\t') * "\n")
                                nrows += 1
                            end
                        end
                    end
                end
            end
        end
    end
    println("DEBUG: Finished loops. Internal nrows: $nrows")
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
println("Auto-mu info:")
if mu_range_mode == :auto_mu
    println("  Mode: :auto_mu")
    println("  Full range: $(auto_mu_full_range)")
    println("  Preserve resolution: $(auto_mu_preserve_resolution)")
    println("  Resolution param: $(auto_mu_resolution_param)")
    
    # Calculate min/max mu points across all groups for reporting
    min_pts = typemax(Int)
    max_pts = 0
    min_mu_range = Inf
    max_mu_range = -Inf
    
    # Stats for smart grouping
    min_combs = typemax(Int)
    max_combs = 0
    total_combs = 0
    
    # Re-calculate global N/step for reporting context
    global_step_size_report = 0.0
    global_N_report = 0
    
    max_phi_global = maximum(phis)
    min_phi_global = minimum(phis)
    max_t2_global = isempty(t2s) ? 1.0 : maximum(t2s)
    
    # Global max range calculation
    mu_max_global = max_t2_global * (1.5 + max_phi_global)
    mu_min_global = auto_mu_full_range ? -mu_max_global : 0.0
    range_width_global = mu_max_global - mu_min_global
    
    if auto_mu_preserve_resolution
        if auto_mu_resolution_param > 1
             # This branch seems unused based on current logic (resolution param is step size)
             # but keeping for consistency if logic changes back
             global_step_size_report = range_width_global / (auto_mu_resolution_param - 1)
        else
             global_step_size_report = auto_mu_resolution_param
        end
    else
        global_N_report = ceil(Int, range_width_global / auto_mu_resolution_param) + 1
    end

    # Iterate groups to find min/max points
    for group in seq_groups
        max_phi = maximum(group.phis)
        mu_max = max_t2_global * (1.5 + max_phi)
        mu_min = auto_mu_full_range ? -mu_max : 0.0
        
        if auto_mu_preserve_resolution
            # Variable N, fixed step
            n_points = ceil(Int, (mu_max - mu_min) / auto_mu_resolution_param) + 1
        else
            # Fixed N (global_N_report), variable step
            n_points = global_N_report
        end
        
        global min_pts = min(min_pts, n_points)
        global max_pts = max(max_pts, n_points)
        global min_mu_range = min(min_mu_range, mu_max - mu_min)
        global max_mu_range = max(max_mu_range, mu_max - mu_min)
        
        # Combinations stats
        n_seqs = length(group.seqs)
        combs = n_seqs * n_points
        global min_combs = min(min_combs, combs)
        global max_combs = max(max_combs, combs)
        global total_combs += combs
    end
    
    println("  Global max range width: $(range_width_global)")
    println("  Min points per group: $min_pts")
    println("  Max points per group: $max_pts")
    println("  Min mu range width: $min_mu_range")
    println("  Max mu range width: $max_mu_range")
    
    if smart_sequence_grouping
        avg_combs = total_combs / length(seq_groups)
        println("  Smart Grouping Stats:")
        println("    Target rows: $target_number_of_rows")
        println("    Generated Groups (Logical): $(length(seq_groups))")
        println("    Target combinations/row: $target_combinations_per_row")
        println("    Min combinations: $min_combs")
        
        # if min_combs < target_combinations_per_row * 0.8
println("  Ns: $N_vals_per_row per row, mode: $N_split_mode")
        #     println("      (Note: Min combinations is low due to the final 'remainder' group - this is normal/harmless)")
        # end
            
        println("    Max combinations: $max_combs")
        println("    Avg combinations: $(round(avg_combs, digits=2))")
        
        if max_pts > mu_vals_per_row
            println("\n  WARNING: Max points ($max_pts) > mu_vals_per_row ($mu_vals_per_row).")
            println("  Some groups will be split into multiple rows, breaking the load balancing.")
            println("  Increase mu_vals_per_row to fix this.")
        end
    end
else
    println("  Mode: :user_defined (using mus_all)")
end
println("______________________________________________")
println("Chunking info:")
if smart_sequence_grouping && mu_range_mode == :auto_mu
    println("  mus: (ignored by Smart Grouping - dynamic)")
else
    println("  mus: $mu_vals_per_row per row, mode: $mu_split_mode")
end
println("  Deltas: $delta_vals_per_row per row, mode: $delta_split_mode")
println("  t1s: $t1_vals_per_row per row, mode: $t1_split_mode")
println("  t2s: $t2_vals_per_row per row, mode: $t2_split_mode")
println("  t3s: $t3_vals_per_row per row, mode: $t3_split_mode")

if smart_sequence_grouping && mu_range_mode == :auto_mu
    println("  slopes: (ignored by Smart Grouping)")
    println("  phasons: (ignored by Smart Grouping)")
    println("  sequences: (ignored by Smart Grouping - dynamic)")
else
    println("  slopes: $slope_vals_per_row per row")
    println("  phasons: $phason_vals_per_row per row, mode: $phason_split_mode")
    println("  sequences: $sequences_per_row per row")
end
println("______________________________________________")
println("Number of rows: $(countlines(out_file) - 1)")
if smart_sequence_grouping && mu_range_mode == :auto_mu
    println("Combinations per row: Variable (see Smart Grouping Stats above)")
else
    perms = N_vals_per_row * mu_vals_per_row * delta_vals_per_row * t1_vals_per_row * t2_vals_per_row * t3_vals_per_row * sequences_per_row
    println("Combinations per row: $perms")
end
