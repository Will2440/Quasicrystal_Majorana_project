module ParamCombGen
export log_range, t_ranges_combs, split_sequence

using IterTools

function log_range(base::Float64, start_power::Float64, end_power::Float64, length::Int)
    powers = range(start_power, stop=end_power, length=length)
    values = base .^ powers
    return values
end

function t_ranges_combs(t_ranges::Vector{Vector{Float64}})
    t_combinations_matrix = collect(IterTools.product(t_ranges...))
    t_combinations = [collect(t) for t in vec(t_combinations_matrix)]
    return t_combinations
end

# function split_sequence(sequence::Vector{Int}, m::Int, N::Int)
#     total_length = m * N
    
#     if total_length > length(sequence)
#         error("Requested $m chunks of length $N, but sequence length is only $(length(sequence)).")
#     end
    
#     return [sequence[(i-1)*N+1 : i*N] for i in 1:m]
# end

function split_sequence(sequence::Vector{Int}, N::Int, start_indices::Vector{Int})
    m = length(start_indices)
    chunks = Vector{Vector{Int}}(undef, m)

    for (i, start) in enumerate(start_indices)
        if start + N - 1 > length(sequence)
            error("Starting index $start with chunk length $N exceeds sequence length $(length(sequence)).")
        end
        chunks[i] = sequence[start:start+N-1]
    end
    
    return chunks
end


end