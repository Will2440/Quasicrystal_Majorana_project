module SeqGen

using ProgressMeter
using Random: AbstractRNG, MersenneTwister, shuffle!
using DataFrames
using BSON


function normal_SeqGen(N::Int)
    return fill(1, N)
end

function golden_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "AB", "B" => "A")
    end
    # println(sequence)
    
    number_sequence = [ch == 'A' ? 1 : 2 for ch in sequence]

    golden_phi = (1 + sqrt(5)) / 2
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(number_sequence, golden_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $golden_phi (tol=$tol)."
    end
    
    return number_sequence
end

function silver_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "BAA", "B" => "A")
    end
    # println(sequence)
    
    number_sequence = [ch == 'A' ? 1 : 2 for ch in sequence]

    silver_phi = 1 + sqrt(2)
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(number_sequence, silver_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $silver_phi (tol=$tol)."
    end
    
    return number_sequence
end

function thue_morse_SeqGen(n::Int)
    if n <= 0 || (n & (n - 1)) != 0
        throw(ArgumentError("n must be a power of 2."))
    end
    
    sequence = Int[]
    for i in 0:n-1
        num_ones = count(c -> c == '1', string(i, base=2))
        push!(sequence, num_ones % 2 + 1)
    end

    thue_morese_phi = 1.0
    tol = 1e-2
    ratio, is_valid = check_PV_ratio(sequence, thue_morese_phi, tol)
    if !is_valid
        @warn "Ratio $ratio is outside the tolerance of φ ≈ $thue_morese_phi (tol=$tol)."
    end
    
    return sequence
end

function plastic_SeqGen(N::Int)
    sequence = "A"
    while length(sequence) < N
        sequence = replace(sequence, "A" => "B", "B" => "AC", "C" => "A") # Van der Laan word generation
        # sequence = replace(sequence, "A" => "AB", "B" => "AC", "C" => "A") # ChatGPT
    end
    
    number_sequence = [ch == 'A' ? 1 : ch == 'B' ? 2 : 3 for ch in sequence]

    plastic_rho = 1.3247179572447460259609088
    tol = 1e-2
    ratio_AB, ratio_BC, AB_is_valid, BC_is_valid = check_plastic_ratio(number_sequence, plastic_rho, tol)
    if !AB_is_valid
        error("Ratio $ratio_AB is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    end
    if !BC_is_valid
        error("Ratio $ratio_BC is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    end

    # if !AB_is_valid || !BC_is_valid
    #     error("Ratio $(ifelse(!AB_is_valid, ratio_AB, ratio_BC)) is outside the tolerance of ρ ≈ $plastic_rho (tol=$tol).")
    # end
    
    return number_sequence
end

function check_PV_ratio(number_sequence::Vector{Int}, phi::Float64, tol::Float64)
    count_A = count(x -> x == 1, number_sequence)
    count_B = count(x -> x == 2, number_sequence)
    
    ratio = count_A / count_B
    bool = abs(ratio - phi) ≤ tol
    
    return ratio, bool
end

function check_plastic_ratio(number_sequence::Vector{Int}, rho::Float64, tol::Float64)
    count_A = count(x -> x == 1, number_sequence)
    count_B = count(x -> x == 2, number_sequence)
    count_C = count(x -> x == 3, number_sequence)

    ratio_AB = count_A / count_B
    ratio_BC = count_B / count_C

    bool_AB = abs(ratio_AB - rho) ≤ tol
    bool_BC = abs(ratio_AB - rho) ≤ tol

    return ratio_AB, ratio_BC, bool_AB, bool_BC
end

function cut_to_length(sequence::Vector{Int}, sequence_length::Int)
    cut_sequence = Vector{Int}(undef, sequence_length)

    for i in 1:sequence_length
        cut_sequence[i] = sequence[i]
    end

    return cut_sequence
end

function golden_LengthCalc(exp_range::Vector{Int})
    phi = (1 + sqrt(5) ) / 2
    phi_conj = (1 - sqrt(5) ) / 2
    sequence_lengths = Int[]

    for g in exp_range
        length = real((phi^(g+2) - phi_conj^(g+2)) / sqrt(5))
        length = floor.(Int, length)
        push!(sequence_lengths, length)
    end

    return sequence_lengths
end

function silver_LengthCalc(n_range::Vector{Int})
    lengths = Int[]
    for n in n_range
        if n == 0
            push!(lengths, 1)
        elseif n == 1
            push!(lengths, 3)
        else
            L_prev_2, L_prev_1 = 1, 3
            for _ in 2:n
                L_current = 2 * L_prev_1 + L_prev_2
                L_prev_2, L_prev_1 = L_prev_1, L_current
            end
            push!(lengths, L_prev_1)
        end
    end
    return lengths
end

function thue_morse_LengthCalc(n_range::Vector{Int})
    lengths = Int[]
    for n in n_range
        length = 2^n
        push!(lengths, length)
    end
    return lengths
end

function plastic_LengthCalc(x_range::Vector{Int})
    lengths = Int[]
    for n in x_range
        sequence = "A"
        for _ in 1:n
            sequence = replace(sequence, "A" => "B", "B" => "AC", "C" => "A") # Van der Laan word generation
        end
        # println(typeof(sequence))
        number_sequence = [ch == 'A' ? 1 : ch == 'B' ? 2 : 3 for ch in sequence]
        length = size(number_sequence,1)
        push!(lengths, length)
    end
    return lengths
end



# # # Example usage

# # # Create pairing parameter sequence for N sites
# sequence_length = 1000

# # Normal Crystal
# N_normal = 1000
# normal_sequence = normal_crystal_gen(N_normal)
# normal_sequence = cut_to_length(normal_sequence, sequence_length)

# # Golden Ratio
# N_gold = 1000
# golden_sequence = golden_ratio_sequence_gen(N_gold)
# golden_sequence = cut_to_length(golden_sequence, sequence_length)

# # Silver Ratio
# N_silver = 1000
# silver_sequence = silver_ratio_sequence_gen(N_silver)
# silver_sequence = cut_to_length(silver_sequence, sequence_length)

# # Thue-Morse 
# N_thue_morse = 1024
# thue_morse_sequence = thue_morse_sequence_gen(N_thue_morse)
# thue_morse_sequence = cut_to_length(thue_morse_sequence, sequence_length)

# # Plastic Ratio
# N_plastic = 1000
# plastic_sequence = plastic_ratio_sequence_gen(N_plastic)
# plastic_sequence = cut_to_length(plastic_sequence, sequence_length)

# println("finished")





#############################################
##### Sturmian Slope Sequence Generation #####
#############################################

# Continued-fraction slope evaluator
function cf_prefix_to_slope(prefix::Vector{Int}; tail_repeats::Int=30)
    @assert tail_repeats >= 1
    x = 1.0
    # approximate infinite tail of ones
    for _ in 2:tail_repeats
        x = 1.0 + 1.0 / x
    end
    # fold the prefix elements from back to front
    for a in reverse(prefix)
        x = float(a) + 1.0 / x
    end
    return 1.0 / x
end

# Streaming enumeration of prefixes up to length L, digits 1:K
function enumerate_prefixes_stream(K::Int, L::Int; measure::Bool=false)
    t0 = time_ns()
    total = sum(K^k for k in 1:L)
    println("Stage 1/2: Enumerating prefixes (total ≈ $total)...")
    prefixes = Channel{Vector{Int}}(Inf) do ch
        function build_prefix(prefix::Vector{Int}, depth::Int)
            if depth > L
                return
            end
            Threads.@threads for a in 1:K
                newprefix = [prefix; a]
                put!(ch, newprefix)
                build_prefix(newprefix, depth + 1)
            end
        end
        build_prefix(Int[], 1)
    end
    if measure
        println("enumerate_prefixes_stream elapsed = $( (time_ns() - t0)/1e9 ) s")
    end
    return prefixes
end

## bypasses use of channel and makes prefixes straight into vector (more efficicent)
function enumerate_prefixes_stream_vector(K::Int, L::Int; measure::Bool=false)
    t0 = time_ns()
    total = sum(K^k for k in 1:L)
    println("Stage 1/2: Enumerating prefixes (total ≈ $total)...")
    
    prefixes = Vector{Vector{Int}}()
    sizehint!(prefixes, total)  # Pre-allocate for efficiency
    
    # Start with length-1 prefixes
    current = [[a] for a in 1:K]
    append!(prefixes, current)
    
    # Iteratively build longer prefixes
    for depth in 2:L
        new_current = Vector{Vector{Int}}()
        for p in current
            for a in 1:K
                new_p = [p; a]
                push!(new_current, new_p)
                push!(prefixes, new_p)
            end
        end
        current = new_current
    end
    
    if measure
        println("enumerate_prefixes_stream elapsed = $((time_ns() - t0)/1e9) s")
    end
    return prefixes
end

# Main function: sample dense set of Sturmian slopes
function sample_sturmian_slopes(K::Int, L::Int; tail_repeats::Int=30, measure::Bool=false)
    t0 = time_ns()
    prefixes_ch = enumerate_prefixes_stream(K, L; measure=true)

    # collect channel into a Vector so it is indexable/has length for threading
    prefixes = collect(prefixes_ch)

    n = length(prefixes)
    println("Stage 2/2: Computing slopes for $n prefixes...")

    slopes = Vector{Float64}(undef, n)
    pref_out = Vector{Vector{Int}}(undef, n)

    Threads.@threads for i in 1:n
        p = prefixes[i]
        s = cf_prefix_to_slope(p; tail_repeats=tail_repeats)
        slopes[i] = s
        pref_out[i] = p
    end

    # sort by slope and return DataFrame in same format as original
    order = sortperm(slopes)
    out = DataFrame(prefix = pref_out[order],
                    slope  = slopes[order])
    if measure
        println("sample_sturmian_slopes elapsed = $( (time_ns() - t0)/1e9 ) s")
    end
    return out
end

# Downsample slopes to rebalance distribution across [0,1]
function rebalance_slopes(slopes_df::DataFrame; bins::Int=200, max_per_bin::Int=5, rng::AbstractRNG=MersenneTwister(42))
    """
    Downsample a (possibly uneven) collection of Sturmian slopes to produce a more
    evenly distributed set across [0,1].

    Input:
      - slopes_df: DataFrame with columns :prefix (Vector{Int}) and :slope (Float64),
                   exactly the same format produced by `sample_sturmian_slopes`.
      - bins: number of equal-width bins across [0,1] used for stratification.
      - max_per_bin: maximum number of samples to keep from any single bin.
      - rng: optional RNG for reproducible random downsampling.

    Output:
      - DataFrame(prefix, slope) with a subset of rows of `slopes_df`, sorted by :slope,
        and with the same column names / types as `sample_sturmian_slopes`.
    """
    # @assert (:prefix in names(slopes_df)) && (:slope in names(slopes_df)) "Expect columns :prefix and :slope"

    # clamp slopes into [0,1] for robust binning
    slopes = clamp.(Float64.(slopes_df.slope), 0.0, 1.0)

    # bin edges and assignments
    edges = collect(range(0.0, stop=1.0, length=bins+1))
    bin_idx = [min(searchsortedlast(edges, s), bins) for s in slopes]  # 1..bins

    # collect row indices per bin
    groups = Dict{Int, Vector{Int}}()
    for (i, b) in enumerate(bin_idx)
        push!(get!(groups, b, Int[]), i)
    end

    # stratified downsample: keep at most max_per_bin items per bin
    selected = Int[]
    for idxs in values(groups)
        if length(idxs) <= max_per_bin
            append!(selected, idxs)
        else
            # random shuffle then take first max_per_bin (no extra dependencies)
            perm = sort(idxs, by = _ -> rand(rng))
            append!(selected, perm[1:max_per_bin])
        end
    end

    if isempty(selected)
        return DataFrame(prefix = Vector{Vector{Int}}(), slope = Float64[])
    end

    # build output DataFrame in the same format as sample_sturmian_slopes (prefix, slope)
    out = slopes_df[selected, [:prefix, :slope]]
    sort!(out, :slope)
    return out
end

function load_sturmian_seq_bson(path::AbstractString)
    @assert isfile(path) "BSON file not found: $path"
    raw = BSON.load(path)

    # try symbol or string keys
    phis = haskey(raw, :phis) ? raw[:phis] : (haskey(raw, "phis") ? raw["phis"] : nothing)
    seqs = haskey(raw, :seqs) ? raw[:seqs] : (haskey(raw, "seqs") ? raw["seqs"] : nothing)
    phasons = haskey(raw, :phasons) ? raw[:phasons] : (haskey(raw, "phasons") ? raw["phasons"] : nothing)

    phis === nothing && error("Key 'phis' not found in BSON file: $path. Keys: $(collect(keys(raw)))")
    seqs === nothing && error("Key 'seqs' not found in BSON file: $path. Keys: $(collect(keys(raw)))")
    phasons === nothing && error("Key 'phasons' not found in BSON file: $path. Keys: $(collect(keys(raw)))")

    # @assert length(phis) == length(seqs) "Length mismatch: phis has $(length(phis)) but seqs has $(length(seqs))"

    # build DataFrame with canonical column names
    return DataFrame(phi = deepcopy(phis), seq = deepcopy(seqs), phason = deepcopy(phasons))
end


# # # Example Usage
# K = 3
# L = 8
# sturmian_slopes_df = sample_sturmian_slopes(K, L; tail_repeats=30)
# balanced_sturm_df = rebalance_slopes(sturmian_slopes_df; bins=200, max_per_bin=5, rng=MersenneTwister(42))




#############################################
##### Generalised Sturmian Slope Sequence Generation #####
#############################################

# Modified continued-fraction slope evaluator with customizable tail
function cf_prefix_to_slope_custom(prefix::Vector{Int}; tail_option::Union{Nothing, String, Int}=nothing, tail_length::Int=30, rng::AbstractRNG=MersenneTwister(42))
    x = 1.0
    
    if tail_option === nothing
        # No tail: approximate with infinite 1's (x starts at 1.0)
        for a in reverse(prefix)
            x = float(a) + 1.0 / x
        end
        return 1.0 / x
    elseif tail_option == "all_ones"
        tail = fill(1, tail_length)
    elseif tail_option == "all_zeros"
        error("Tail option 'all_zeros' not supported: continued fraction with zeros is invalid.")
    elseif tail_option == "random"
        tail = rand(rng, 1:maximum(prefix), tail_length)  # Random 1s and 2s for binary Sturmian
    elseif isa(tail_option, Int)
        if tail_option < 0
            error("Tail option Int must be positive.")
        end
        tail = fill(tail_option, tail_length)
    else
        error("Invalid tail_option: $tail_option. Use nothing, 'all_ones', 'random', or a positive Int.")
    end
    
    # Fold tail from innermost to outermost (end to start)
    for i in tail_length:-1:1
        a = tail[i]
        x = float(a) + 1.0 / x
    end
    
    # Fold prefix from innermost to outermost (end to start)
    for a in reverse(prefix)
        x = float(a) + 1.0 / x
    end
    
    return 1.0 / x
end

# New function: sample Sturmian slopes with customizable tail options, outputting DataFrame like sample_sturmian_slopes
function sample_sturmian_slopes_custom(K::Int, L::Int; tail_option::Union{Nothing, String, Int}=nothing, tail_length::Int=30, measure::Bool=false, rng::AbstractRNG=MersenneTwister(42))
    t0 = time_ns()
    prefixes_ch = enumerate_prefixes_stream(K, L; measure=true)
    
    # Collect channel into a Vector
    prefixes = collect(prefixes_ch)
    
    n = length(prefixes)
    if measure
        println("Stage 2/2: Computing slopes for $n prefixes with custom tail $tail_option...")
    end
    
    slopes = Vector{Float64}(undef, n)
    pref_out = Vector{Vector{Int}}(undef, n)
    
    Threads.@threads for i in 1:n
        p = prefixes[i]
        s = cf_prefix_to_slope_custom(p; tail_option=tail_option, tail_length=tail_length, rng=rng)
        slopes[i] = s
        pref_out[i] = p
    end
    
    # Sort by slope and return DataFrame
    order = sortperm(slopes)
    out = DataFrame(prefix = pref_out[order], slope = slopes[order])
    if measure
        println("sample_sturmian_slopes_custom elapsed = $((time_ns() - t0)/1e9) s")
    end
    return out
end



#################################
####### Version which can balance slopes before genrated -- more efficicent ###########
#################################

# New function: sample Sturmian slopes with customizable tail options, outputting DataFrame like sample_sturmian_slopes
function sample_sturmian_slopes_custom_optBalance(K::Int, L::Int; tail_option::Union{Nothing, String, Int}=nothing, tail_length::Int=30, measure::Bool=false, rng::AbstractRNG=MersenneTwister(42), balanced::Bool=false, bins::Int=200, max_per_bin::Int=5)
    t0 = time_ns()
    prefixes = enumerate_prefixes_stream_vector(K, L; measure=true)
    
    n = length(prefixes)
    if measure
        if balanced
            println("Stage 2/2: Balancing slopes nd computing for up to $(bins*max_per_bin) prefixes with custom tail...")
        else
            println("Stage 2/2: Computing slopes for $n prefixes with custom tail...")
        end
    end
    
    if balanced
        # Prepare bins
        edges = collect(range(0.0, stop=1.0, length=bins+1))
        bin_counts = Dict{Int, Int}()  # Track how many per bin
        selected_prefixes = Vector{Vector{Int}}()
        selected_slopes = Vector{Float64}()
        bin_lock = ReentrantLock()  # Declare the lock for thread safety
        
        Threads.@threads for p in prefixes
            s = cf_prefix_to_slope_custom(p; tail_option=tail_option, tail_length=tail_length, rng=rng)
            s_clamped = clamp(s, 0.0, 1.0)
            bin_idx = min(searchsortedlast(edges, s_clamped), bins)
            
            lock(bin_lock)
            try
                if get(bin_counts, bin_idx, 0) < max_per_bin
                    push!(selected_prefixes, p)
                    push!(selected_slopes, s)
                    bin_counts[bin_idx] = get(bin_counts, bin_idx, 0) + 1
                end
            finally
                unlock(bin_lock)
            end
        end
        
        # Sort by slope
        order = sortperm(selected_slopes)
        out = DataFrame(prefix = selected_prefixes[order], slope = selected_slopes[order])
    else
        # Original method: compute all slopes in parallel
        slopes = Vector{Float64}(undef, n)
        pref_out = Vector{Vector{Int}}(undef, n)
        
        Threads.@threads for i in 1:n
            p = prefixes[i]
            s = cf_prefix_to_slope_custom(p; tail_option=tail_option, tail_length=tail_length, rng=rng)
            slopes[i] = s
            pref_out[i] = p
        end
        
        # Sort by slope and return DataFrame
        order = sortperm(slopes)
        out = DataFrame(prefix = pref_out[order], slope = slopes[order])
    end
    
    if measure
        println("sample_sturmian_slopes_custom elapsed = $((time_ns() - t0)/1e9) s")
    end
    return out
end

end # end module
