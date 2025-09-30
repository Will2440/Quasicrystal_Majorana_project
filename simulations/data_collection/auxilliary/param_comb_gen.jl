module ParamCombGen
export log_range, t_ranges_combs, split_sequence, restrict_range, angle_to_gradient

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

function restrict_range(xs::Vector{Float64}, ys::Vector{Float64}, cuts::Vector{Dict{Symbol, Any}})
    """
    restrict_range(xs, ys, cuts)

    Filter points in the 2D parameter space defined by `xs` and `ys` according to 
    a list of cutting rules.

    # Arguments
    - `xs::Vector{Float64}`: x values
    - `ys::Vector{Float64}`: y values
    - `cuts::Vector{Dict}`: list of cutting rules. Each rule is a Dict with keys:
        - `:gradient` (slope of line)
        - `:y_intercept` (intercept of line)
        - `:x_range` (Tuple `(xmin, xmax)` where the rule applies)
        - `:y_range` (Tuple `(ymin, ymax)` where the rule applies)
        - `:cut_which_side` (`"above"` or `"below"`)

    # Returns
    - A vector of tuples `(x, y)` that survive all cuts.
    """

    points = [(x, y) for x in xs, y in ys]
    survivors = Vector{Tuple{Float64, Float64}}()

    for (x, y) in points
        keep = true
        for cut in cuts
            m  = cut[:gradient]
            b  = cut[:y_intercept]
            xr = cut[:x_range]
            yr = cut[:y_range]
            side = cut[:cut_which_side]

            # only apply cut if point is inside the region where the cut is active
            if xr[1] ≤ x ≤ xr[2] && yr[1] ≤ y ≤ yr[2]
                y_line = m * x + b
                if side == "below"
                    # cut points BELOW the line
                    if y < y_line
                        keep = false
                        break
                    end
                elseif side == "above"
                    # cut points ABOVE the line
                    if y > y_line
                        keep = false
                        break
                    end
                else
                    error("cut_which_side must be \"left\" or \"right\"")
                end
            end
        end

        if !keep
            push!(survivors, (x, y))
        end
    end

    return survivors::Vector{Tuple{Float64, Float64}}
end

function angle_to_gradient(angle_deg::Float64)
    angle_rad = deg2rad(angle_deg)
    return tan(angle_rad)
end

end