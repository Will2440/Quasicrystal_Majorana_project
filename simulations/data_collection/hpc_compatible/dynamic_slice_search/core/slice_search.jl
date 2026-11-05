module SliceSearch

export build_rho_grid,
       baseline_mu_grid,
       detect_transitions,
       refine_transition_bracket,
       refine_low_gap_intervals!

function build_rho_grid(rho_start::Real, rho_max::Real, n_slices::Integer, spacing::Symbol=:linear)
    rho_start > 0 || error("rho_start must be positive")
    rho_max >= rho_start || error("rho_max must be greater than or equal to rho_start")
    n_slices >= 1 || error("n_rho_slices must be at least one")

    if spacing == :linear
        return collect(range(Float64(rho_start), Float64(rho_max); length=n_slices))
    elseif spacing == :log
        return exp.(range(log(Float64(rho_start)), log(Float64(rho_max)); length=n_slices))
    end
    error("Unknown rho spacing mode: $spacing (expected :linear or :log)")
end

function baseline_mu_grid(mu_lo::Real, mu_hi::Real, n_points::Integer)
    mu_hi > mu_lo || error("mu_hi must be greater than mu_lo")
    n_points >= 2 || error("mu_baseline_points must be at least two")
    return collect(range(Float64(mu_lo), Float64(mu_hi); length=n_points))
end

function detect_transitions(samples)
    ordered = sort(collect(samples); by=s -> s.mu)
    transitions = NamedTuple[]
    length(ordered) < 2 && return transitions
    for i in 1:(length(ordered) - 1)
        left, right = ordered[i], ordered[i + 1]
        left.phase_label == right.phase_label && continue
        push!(transitions, (
            mu_lo=Float64(left.mu),
            mu_hi=Float64(right.mu),
            mu=0.5 * (Float64(left.mu) + Float64(right.mu)),
            phase_lo=Int(left.phase_label),
            phase_hi=Int(right.phase_label),
            gap_min=min(Float64(left.gap_min), Float64(right.gap_min)),
        ))
    end
    return transitions
end

function refine_transition_bracket(eval_at_mu, left, right, mu_tol::Real, max_iters::Integer)
    left.phase_label != right.phase_label || error("Transition refinement requires opposite endpoint labels")
    mu_tol > 0 || error("transition_mu_tol must be positive")

    lo, hi = left, right
    iterations = 0
    while (Float64(hi.mu) - Float64(lo.mu) > mu_tol) && iterations < max_iters
        mid_mu = 0.5 * (Float64(lo.mu) + Float64(hi.mu))
        mid = eval_at_mu(mid_mu)
        if mid.phase_label == lo.phase_label
            lo = mid
        else
            hi = mid
        end
        iterations += 1
    end

    return (
        mu_lo=Float64(lo.mu),
        mu_hi=Float64(hi.mu),
        mu=0.5 * (Float64(lo.mu) + Float64(hi.mu)),
        phase_lo=Int(lo.phase_label),
        phase_hi=Int(hi.phase_label),
        gap_min=min(Float64(lo.gap_min), Float64(hi.gap_min)),
        refinement_iters=iterations,
        converged=(Float64(hi.mu) - Float64(lo.mu) <= mu_tol),
    )
end

function refine_low_gap_intervals!(eval_at_mu, samples, gap_tol::Real, mu_tol::Real, max_depth::Integer)
    ordered = sort(collect(samples); by=s -> s.mu)
    added = NamedTuple[]
    stack = [(ordered[i], ordered[i + 1], 0) for i in 1:(length(ordered) - 1)]

    while !isempty(stack)
        left, right, depth = pop!(stack)
        width = Float64(right.mu) - Float64(left.mu)
        width <= mu_tol && continue
        depth >= max_depth && continue
        min(Float64(left.gap_min), Float64(right.gap_min)) <= gap_tol || continue

        mid = eval_at_mu(0.5 * (Float64(left.mu) + Float64(right.mu)))
        push!(added, mid)
        push!(stack, (left, mid, depth + 1), (mid, right, depth + 1))
    end

    return added
end

end # module