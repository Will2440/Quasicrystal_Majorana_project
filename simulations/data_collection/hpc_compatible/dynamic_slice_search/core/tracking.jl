module TransitionTracking

export EdgeTrack, predict_edges, update_tracks!

mutable struct EdgeTrack
    id::Int
    orientation::Tuple{Int,Int}
    history::Vector{NamedTuple}
    last_seen::Int
end

function predict_edges(tracks::Vector{EdgeTrack}, next_value::Real, current_slice::Int, mu_bounds, cfg::Dict{Symbol,Any})
    predictions = NamedTuple[]
    min_history = Int(get(cfg, :momentum_min_history, 3))
    momentum_enabled = Bool(get(cfg, :momentum_enabled, true))
    damping = Float64(get(cfg, :momentum_damping, 0.2))
    max_shift_fraction = Float64(get(cfg, :momentum_max_shift_fraction, 0.02))
    min_band_fraction = Float64(get(cfg, :prediction_band_min_fraction,
        get(cfg, :prediction_band_half_width_fraction, 0.03)))
    max_band_fraction = Float64(get(cfg, :prediction_band_max_fraction, 0.08))
    motion_multiplier = Float64(get(cfg, :prediction_band_motion_multiplier, 2.5))
    curvature_multiplier = Float64(get(cfg, :prediction_band_curvature_multiplier, 1.5))
    active_gap = Int(get(cfg, :track_max_missed_slices, 1))
    mu_span = Float64(mu_bounds[2] - mu_bounds[1])

    for track in tracks
        isempty(track.history) && continue
        current_slice - track.last_seen > active_gap && continue
        last = track.history[end]
        predicted_mu = Float64(last.mu)
        predicted_motion = 0.0
        curvature_uncertainty = 0.0

        if momentum_enabled && length(track.history) >= min_history
            prior = track.history[end - 2]
            delta_value = Float64(last.sweep_value) - Float64(prior.sweep_value)
            next_delta = Float64(next_value) - Float64(last.sweep_value)
            if delta_value != 0 && next_delta * delta_value > 0
                slope = (Float64(last.mu) - Float64(prior.mu)) / delta_value
                predicted_motion = abs(slope * next_delta)
                raw_shift = damping * slope * next_delta
                shift_cap = max_shift_fraction * mu_span
                predicted_mu += clamp(raw_shift, -shift_cap, shift_cap)

                if length(track.history) >= 4
                    earlier = track.history[end - 2]
                    older = track.history[end - 3]
                    older_delta = Float64(earlier.sweep_value) - Float64(older.sweep_value)
                    recent_delta = Float64(last.sweep_value) - Float64(earlier.sweep_value)
                    if older_delta * recent_delta > 0
                        older_slope = (Float64(earlier.mu) - Float64(older.mu)) / older_delta
                        recent_slope = (Float64(last.mu) - Float64(earlier.mu)) / recent_delta
                        curvature_uncertainty = abs(recent_slope - older_slope) * abs(next_delta)
                    end
                end
            end
        end

        predicted_mu = clamp(predicted_mu, Float64(mu_bounds[1]), Float64(mu_bounds[2]))
        min_band = min_band_fraction * mu_span
        max_band = max(max_band_fraction * mu_span, min_band)
        band_half_width = clamp(
            max(min_band, motion_multiplier * predicted_motion + curvature_multiplier * curvature_uncertainty),
            min_band,
            max_band,
        )
        push!(predictions, (
            track_id=track.id,
            mu=predicted_mu,
            orientation=track.orientation,
            band_half_width=band_half_width,
            predicted_motion=predicted_motion,
            curvature_uncertainty=curvature_uncertainty,
        ))
    end
    return predictions
end

function ordered_match_cost(active_tracks, transitions, predictions, mu_span::Float64, cfg)
    n_tracks = length(active_tracks)
    n_transitions = length(transitions)
    skip_penalty = Float64(get(cfg, :track_unmatched_penalty, 0.75))
    base_tolerance = Float64(get(cfg, :track_match_base_tolerance_fraction, 0.01)) * mu_span
    max_tolerance = Float64(get(cfg, :track_match_max_tolerance_fraction, 0.05)) * mu_span
    step_scale = Float64(get(cfg, :track_match_step_scale, 1.5))

    costs = fill(Inf, n_tracks + 1, n_transitions + 1)
    actions = fill(UInt8(0), n_tracks + 1, n_transitions + 1)
    costs[1, 1] = 0.0

    for i in 0:n_tracks, j in 0:n_transitions
        current = costs[i + 1, j + 1]
        isfinite(current) || continue

        if i < n_tracks
            candidate = current + skip_penalty
            if candidate < costs[i + 2, j + 1]
                costs[i + 2, j + 1] = candidate
                actions[i + 2, j + 1] = 0x01
            end
        end
        if j < n_transitions
            candidate = current + skip_penalty
            if candidate < costs[i + 1, j + 2]
                costs[i + 1, j + 2] = candidate
                actions[i + 1, j + 2] = 0x02
            end
        end
        if i < n_tracks && j < n_transitions
            track = active_tracks[i + 1]
            transition = transitions[j + 1]
            predicted_mu = get(predictions, track.id, Float64(track.history[end].mu))
            last_mu = Float64(track.history[end].mu)
            tolerance = min(max_tolerance, base_tolerance + step_scale * abs(predicted_mu - last_mu))
            orientation = (Int(transition.phase_lo), Int(transition.phase_hi))
            distance = abs(predicted_mu - Float64(transition.mu))

            if track.orientation == orientation && distance <= tolerance
                match_cost = distance / max(tolerance, eps(Float64))
                candidate = current + match_cost
                if candidate < costs[i + 2, j + 2]
                    costs[i + 2, j + 2] = candidate
                    actions[i + 2, j + 2] = 0x03
                end
            end
        end
    end

    matches = Tuple{Int,Int}[]
    i, j = n_tracks + 1, n_transitions + 1
    while i > 1 || j > 1
        action = actions[i, j]
        if action == 0x03
            push!(matches, (i - 1, j - 1))
            i -= 1
            j -= 1
        elseif action == 0x01
            i -= 1
        elseif action == 0x02
            j -= 1
        else
            error("Could not reconstruct ordered transition matching")
        end
    end
    reverse!(matches)
    return matches
end

function update_tracks!(tracks::Vector{EdgeTrack}, transitions, value::Real, slice_index::Int, mu_bounds, cfg::Dict{Symbol,Any}; predictions=NamedTuple[])
    mu_span = Float64(mu_bounds[2] - mu_bounds[1])
    active_gap = Int(get(cfg, :track_max_missed_slices, 1))
    active_tracks = filter(track -> slice_index - track.last_seen <= active_gap, tracks)
    sort!(active_tracks; by=track -> track.history[end].mu)
    original_transitions = collect(transitions)
    transition_order = sortperm(eachindex(original_transitions); by=i -> original_transitions[i].mu)
    ordered_transitions = original_transitions[transition_order]
    next_id = isempty(tracks) ? 1 : maximum(t.id for t in tracks) + 1
    prediction_map = Dict{Int,Float64}(p.track_id => Float64(p.mu) for p in predictions)
    matches = ordered_match_cost(active_tracks, ordered_transitions, prediction_map, mu_span, cfg)
    transition_ids = zeros(Int, length(ordered_transitions))

    for (track_index, transition_index) in matches
        track = active_tracks[track_index]
        transition = ordered_transitions[transition_index]
        push!(track.history, (sweep_value=Float64(value), mu=Float64(transition.mu)))
        track.last_seen = slice_index
        transition_ids[transition_index] = track.id
    end

    matched_transitions = Set(last(pair) for pair in matches)
    for transition_index in eachindex(ordered_transitions)
        transition_index in matched_transitions && continue
        transition = ordered_transitions[transition_index]
        orientation = (Int(transition.phase_lo), Int(transition.phase_hi))
        point = (sweep_value=Float64(value), mu=Float64(transition.mu))
        push!(tracks, EdgeTrack(next_id, orientation, [point], slice_index))
        transition_ids[transition_index] = next_id
        next_id += 1
    end

    input_order_ids = similar(transition_ids)
    for (sorted_index, original_index) in enumerate(transition_order)
        input_order_ids[original_index] = transition_ids[sorted_index]
    end
    return input_order_ids
end

end # module