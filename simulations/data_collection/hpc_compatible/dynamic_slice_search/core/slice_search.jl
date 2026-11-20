module SliceSearch

export build_rho_grid,
    build_delta_grid,
    rho_progress,
       baseline_mu_grid,
    scaled_baseline_point_count,
    boundary_points_per_prediction,
       detect_transitions,
       refine_transition_bracket,
       refine_low_gap_intervals!

function build_rho_grid(rho_start::Real, rho_max::Real, n_slices::Integer, spacing::Symbol=:linear;
                        concentration::Real=4.0)
    rho_start > 0 || error("rho_start must be positive")
    rho_max >= rho_start || error("rho_max must be greater than or equal to rho_start")
    n_slices >= 1 || error("n_rho_slices must be at least one")

    if spacing == :linear
        n_slices == 1 && return [Float64(rho_start)]
        return collect(range(Float64(rho_start), Float64(rho_max); length=n_slices))
    elseif spacing == :log
        isfinite(concentration) && 0 < concentration <= 700 ||
            error("rho log concentration must be finite and in (0, 700]")
        n_slices == 1 && return [Float64(rho_start)]
        rho_lo, rho_hi = Float64(rho_start), Float64(rho_max)
        positions = range(0.0, 1.0; length=n_slices)
        bias_to_low_rho = rho_lo * rho_hi <= 1.0
        log_fractions = if bias_to_low_rho
            expm1.(Float64(concentration) .* positions) ./ expm1(Float64(concentration))
        else
            1 .- expm1.(Float64(concentration) .* (1 .- positions)) ./ expm1(Float64(concentration))
        end
        return exp.(log(rho_lo) .+ log(rho_hi / rho_lo) .* log_fractions)
    end
    error("Unknown rho spacing mode: $spacing (expected :linear or :log)")
end

function build_delta_grid(delta_start::Real, delta_max::Real, n_slices::Integer, spacing::Symbol=:linear; concentration::Real=4.0)
    isfinite(delta_start) && isfinite(delta_max) && delta_start >= 0 && delta_max >= 0 ||
        error("Delta endpoints must be finite and nonnegative")
    n_slices >= 1 || error("n_Delta_slices must be at least one")

    if spacing == :linear
        n_slices == 1 && return [Float64(delta_start)]
        return collect(range(Float64(delta_start), Float64(delta_max); length=n_slices))
    elseif spacing == :log
        delta_start > 0 && delta_max > 0 || error("Delta endpoints must be positive for log spacing")
        n_slices == 1 && return [Float64(delta_start)]
        return exp.(range(log(Float64(delta_start)), log(Float64(delta_max)); length=n_slices))
    elseif spacing == :exponential
        isfinite(concentration) && 0 < concentration <= 700 || error("Delta exponential concentration must be finite and in (0, 700]")
        n_slices == 1 && return [Float64(delta_start)]
        lo, hi = minmax(Float64(delta_start), Float64(delta_max))
        positions = range(delta_start <= delta_max ? 0.0 : 1.0, delta_start <= delta_max ? 1.0 : 0.0; length=n_slices)
        return [lo + (hi - lo) * expm1(Float64(concentration) * position) / expm1(Float64(concentration)) for position in positions]
    end
    error("Unknown Delta spacing mode: $spacing (expected :linear, :log, or :exponential)")
end

function rho_progress(rho::Real, rho_start::Real, rho_max::Real)
    rho_max >= rho_start || error("rho_max must be greater than or equal to rho_start")
    rho_max == rho_start && return 0.0
    return clamp((Float64(rho) - Float64(rho_start)) / (Float64(rho_max) - Float64(rho_start)), 0.0, 1.0)
end

function baseline_mu_grid(mu_lo::Real, mu_hi::Real, n_points::Integer)
    mu_hi > mu_lo || error("mu_hi must be greater than mu_lo")
    n_points >= 2 || error("mu_baseline_points must be at least two")
    return collect(range(Float64(mu_lo), Float64(mu_hi); length=n_points))
end

function scaled_baseline_point_count(
    current_mu_span::Real,
    reference_mu_span::Real,
    reference_points::Integer;
    rho_fraction::Real=0.0,
    rho_growth_factor::Real=1.0,
    rho_growth_power::Real=1.0,
)
    current_mu_span > 0 || error("current mu span must be positive")
    reference_mu_span > 0 || error("reference mu span must be positive")
    reference_points >= 2 || error("reference baseline point count must be at least two")
    rho_growth_factor >= 1 || error("rho discovery growth factor must be at least one")
    rho_growth_power > 0 || error("rho discovery growth power must be positive")

    reference_intervals = reference_points - 1
    growth_progress = clamp(Float64(rho_fraction), 0.0, 1.0)^Float64(rho_growth_power)
    rho_growth = 1.0 + (Float64(rho_growth_factor) - 1.0) * growth_progress
    scaled_intervals = ceil(Int,
        reference_intervals * Float64(current_mu_span) / Float64(reference_mu_span) * rho_growth)
    return max(2, scaled_intervals + 1)
end

function boundary_points_per_prediction(start_points::Integer, max_points::Integer, rho_fraction::Real, growth_power::Real=1.0)
    start_points >= 3 || error("boundary points at rho_start must be at least three")
    max_points >= start_points || error("maximum boundary points must be at least the starting count")
    growth_power > 0 || error("boundary growth power must be positive")

    progress = clamp(Float64(rho_fraction), 0.0, 1.0)^Float64(growth_power)
    return round(Int, start_points + (max_points - start_points) * progress)
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