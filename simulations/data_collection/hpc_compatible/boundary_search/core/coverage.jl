module BoundaryCoverage

export plan_grid

function plan_grid(bounds, count::Integer, exclusions, audit_count::Integer)
    lo, hi = Float64.(bounds)
    isfinite(lo) && isfinite(hi) && lo < hi || error("u_range must have finite increasing bounds")
    count >= 3 || error("u_baseline_points must be at least 3")
    audit_count >= 0 || error("audit_points_per_interval cannot be negative")
    intervals = sort([(max(lo, Float64(interval[1])), min(hi, Float64(interval[2])))
                      for interval in exclusions if interval[1] < interval[2]]; by=first)
    merged = Tuple{Float64,Float64}[]
    for (left, right) in intervals
        left < right || continue
        if !isempty(merged) && left <= merged[end][2]
            merged[end] = (merged[end][1], max(merged[end][2], right))
        else
            push!(merged, (left, right))
        end
    end
    lengths = Float64[]
    regions = Tuple{Float64,Float64}[]
    cursor = lo
    for (left, right) in merged
        left > cursor && push!(regions, (cursor, left))
        cursor = right
    end
    cursor < hi && push!(regions, (cursor, hi))
    append!(lengths, [right - left for (left, right) in regions])
    total = sum(lengths)
    points = Float64[lo, hi]
    if total > 0
        for ((left, right), length_region) in zip(regions, lengths)
            intervals_for_region = max(1, ceil(Int, (count - 1) * length_region / total))
            append!(points, range(left, right; length=intervals_for_region + 1))
        end
    end
    for (left, right) in merged
        audit_count > 0 && append!(points, range(left, right; length=audit_count + 2)[2:end-1])
    end
    sort!(unique!(points))
    return (points=points, exclusions=merged, regions=regions)
end

end