module TransitionTracking

using Statistics: median

export EdgeTrack, predict_edges, update_tracks!, predicted_exclusions

mutable struct EdgeTrack
    id::Int
    orientation::Tuple{Int,Int}
    history::Vector{NamedTuple}
    last_seen::Int
end

function predict_edges(tracks::Vector{EdgeTrack}, next_value::Real, current_slice::Int, u_bounds, cfg::Dict{Symbol,Any})
    predictions = NamedTuple[]
    lower, upper = Float64.(u_bounds)
    span = upper - lower
    span > 0 || error("u bounds must be increasing")
    min_band = Float64(get(cfg, :boundary_padding_min, 0.03))
    max_band = Float64(get(cfg, :boundary_padding_max, 0.08))
    0 <= min_band <= max_band || error("boundary padding must satisfy 0 <= min <= max")
    active_gap = Int(get(cfg, :track_max_missed_slices, 1))
    min_history = max(1, Int(get(cfg, :momentum_min_history, 3)))
    confidence_history = max(1, Int(get(cfg, :confidence_min_history, 3)))
    damping = Float64(get(cfg, :momentum_damping, 0.2))
    shift_cap = max(0.0, Float64(get(cfg, :momentum_max_shift_fraction, 0.02))) * span

    for track in tracks
        isempty(track.history) && continue
        gap = current_slice - track.last_seen
        0 < gap <= active_gap || continue
        last = track.history[end]
        next_delta = Float64(next_value) - Float64(last.sweep_value)
        slopes = Float64[]
        if Bool(get(cfg, :momentum_enabled, true)) && length(track.history) >= min_history
            first_index = max(2, length(track.history) - 4)
            for index in length(track.history):-1:first_index
                previous, point = track.history[index - 1], track.history[index]
                get(previous, :converged, true) && get(point, :converged, true) || break
                delta = Float64(point.sweep_value) - Float64(previous.sweep_value)
                delta * next_delta > 0 || break
                push!(slopes, (Float64(point.u) - Float64(previous.u)) / delta)
            end
        end
        slope = isempty(slopes) ? 0.0 : median(slopes)
        raw_shift = damping * slope * next_delta
        shift = clamp(raw_shift, -shift_cap, shift_cap)
        predicted_u = clamp(Float64(last.u) + shift, lower, upper)
        predicted_motion = abs(predicted_u - Float64(last.u))
        slope_uncertainty = isempty(slopes) ? 0.0 : median(abs.(slopes .- slope)) * abs(damping * next_delta)
        motion_uncertainty = max(slope_uncertainty, abs(raw_shift - (predicted_u - Float64(last.u))))
        confidence = gap == 1 && get(last, :validated_matches, 0) >= confidence_history &&
            get(last, :converged, true)
        band = confidence ? clamp(max(min_band, motion_uncertainty), min_band, max_band) : max_band
        push!(predictions, (
            track_id=track.id, u=predicted_u, orientation=track.orientation,
            band_half_width=band, confidence=confidence, predicted_motion=predicted_motion,
            motion_uncertainty=motion_uncertainty, sweep_value=Float64(next_value),
            slice_index=current_slice, source_last_seen=track.last_seen,
            source_history_length=length(track.history),
        ))
    end
    return predictions
end

function match_tolerance(track, predicted_u, span, cfg)
    base = Float64(get(cfg, :track_match_base_tolerance_fraction, 0.01)) * span
    maximum_tolerance = Float64(get(cfg, :track_match_max_tolerance_fraction, 0.05)) * span
    step_scale = Float64(get(cfg, :track_match_step_scale, 1.5))
    return min(maximum_tolerance, base + step_scale * abs(predicted_u - Float64(track.history[end].u)))
end

function ordered_match_cost(active_tracks, transitions, predictions, span::Float64, cfg)
    n_tracks, n_transitions = length(active_tracks), length(transitions)
    skip_penalty = Float64(get(cfg, :track_unmatched_penalty, 0.75))
    eligible = falses(n_tracks, n_transitions)
    match_costs = fill(Inf, n_tracks, n_transitions)
    for track_index in 1:n_tracks, transition_index in 1:n_transitions
        track, transition = active_tracks[track_index], transitions[transition_index]
        predicted_u = get(predictions, track.id, Float64(track.history[end].u))
        tolerance = match_tolerance(track, predicted_u, span, cfg)
        distance = abs(predicted_u - Float64(transition.u))
        if track.orientation == (Int(transition.phase_lo), Int(transition.phase_hi)) && distance <= tolerance
            eligible[track_index, transition_index] = true
            match_costs[track_index, transition_index] = distance / max(tolerance, eps(Float64))
        end
    end
    costs = fill(Inf, n_tracks + 1, n_transitions + 1)
    actions = fill(UInt8(0), n_tracks + 1, n_transitions + 1)
    costs[1, 1] = 0.0
    for track_index in 0:n_tracks, transition_index in 0:n_transitions
        current = costs[track_index + 1, transition_index + 1]
        isfinite(current) || continue
        if track_index < n_tracks && current + skip_penalty < costs[track_index + 2, transition_index + 1]
            costs[track_index + 2, transition_index + 1] = current + skip_penalty
            actions[track_index + 2, transition_index + 1] = 0x01
        end
        if transition_index < n_transitions && current + skip_penalty < costs[track_index + 1, transition_index + 2]
            costs[track_index + 1, transition_index + 2] = current + skip_penalty
            actions[track_index + 1, transition_index + 2] = 0x02
        end
        if track_index < n_tracks && transition_index < n_transitions
            candidate = current + match_costs[track_index + 1, transition_index + 1]
            if candidate < costs[track_index + 2, transition_index + 2]
                costs[track_index + 2, transition_index + 2] = candidate
                actions[track_index + 2, transition_index + 2] = 0x03
            end
        end
    end
    matches = Tuple{Int,Int}[]
    track_index, transition_index = n_tracks + 1, n_transitions + 1
    while track_index > 1 || transition_index > 1
        action = actions[track_index, transition_index]
        if action == 0x03
            push!(matches, (track_index - 1, transition_index - 1))
            track_index -= 1
            transition_index -= 1
        elseif action == 0x01
            track_index -= 1
        elseif action == 0x02
            transition_index -= 1
        else
            error("Could not reconstruct ordered transition matching")
        end
    end
    reverse!(matches)
    return matches, eligible
end

function update_tracks!(tracks::Vector{EdgeTrack}, transitions, value::Real, slice_index::Int, u_bounds, cfg::Dict{Symbol,Any}; predictions=NamedTuple[])
    span = Float64(u_bounds[2] - u_bounds[1])
    span > 0 || error("u bounds must be increasing")
    active_gap = Int(get(cfg, :track_max_missed_slices, 1))
    active_tracks = filter(track -> !isempty(track.history) && 0 < slice_index - track.last_seen <= active_gap, tracks)
    sort!(active_tracks; by=track -> track.history[end].u)
    original_transitions = collect(transitions)
    transition_order = sortperm(eachindex(original_transitions); by=index -> original_transitions[index].u)
    ordered_transitions = original_transitions[transition_order]
    valid_predictions = Dict(p.track_id => p for p in predictions if
        hasproperty(p, :sweep_value) && p.sweep_value == Float64(value) &&
        hasproperty(p, :slice_index) && p.slice_index == slice_index &&
        hasproperty(p, :source_last_seen) && hasproperty(p, :source_history_length))
    filter!(pair -> any(track -> track.id == pair.first && track.last_seen == pair.second.source_last_seen &&
        length(track.history) == pair.second.source_history_length, active_tracks), valid_predictions)
    prediction_map = Dict(id => Float64(prediction.u) for (id, prediction) in valid_predictions)
    matches, eligible = ordered_match_cost(active_tracks, ordered_transitions, prediction_map, span, cfg)
    transition_ids = zeros(Int, length(ordered_transitions))
    matched_ids = Set{Int}()
    max_residual = Float64(get(cfg, :confidence_max_residual, 0.01))

    for (track_index, transition_index) in matches
        track, transition = active_tracks[track_index], ordered_transitions[transition_index]
        converged = Bool(get(transition, :converged, true))
        unambiguous = count(eligible[track_index, :]) == 1 && count(eligible[:, transition_index]) == 1
        prediction = get(valid_predictions, track.id, nothing)
        residual = prediction === nothing ? Inf : abs(Float64(prediction.u) - Float64(transition.u))
        validated = converged && unambiguous && slice_index - track.last_seen == 1 && residual <= max_residual
        streak = validated ? get(track.history[end], :validated_matches, 0) + 1 : 0
        push!(track.history, (sweep_value=Float64(value), u=Float64(transition.u),
            converged=converged, validated_matches=streak, residual=residual))
        track.last_seen = slice_index
        transition_ids[transition_index] = track.id
        push!(matched_ids, track.id)
    end
    for track in tracks
        isempty(track.history) && continue
        if !(track.id in matched_ids)
            track.history[end] = merge(track.history[end], (validated_matches=0,))
        end
    end
    next_id = isempty(tracks) ? 1 : maximum(track.id for track in tracks) + 1
    for transition_index in eachindex(ordered_transitions)
        transition_ids[transition_index] != 0 && continue
        transition = ordered_transitions[transition_index]
        point = (sweep_value=Float64(value), u=Float64(transition.u),
            converged=Bool(get(transition, :converged, true)), validated_matches=0, residual=Inf)
        push!(tracks, EdgeTrack(next_id, (Int(transition.phase_lo), Int(transition.phase_hi)), [point], slice_index))
        transition_ids[transition_index] = next_id
        next_id += 1
    end
    input_order_ids = similar(transition_ids)
    for (sorted_index, original_index) in enumerate(transition_order)
        input_order_ids[original_index] = transition_ids[sorted_index]
    end
    return input_order_ids
end

function predicted_exclusions(predictions, latest_transitions, ids, cfg, bounds)
    transitions = collect(latest_transitions)
    length(transitions) == length(ids) || throw(ArgumentError("transition IDs must align with transitions"))
    lower, upper = Float64.(bounds)
    upper > lower || error("u bounds must be increasing")
    order = sortperm(eachindex(transitions); by=index -> transitions[index].u)
    prediction_map = Dict(prediction.track_id => prediction for prediction in predictions)
    exclusions = Tuple{Float64,Float64}[]
    for index in 1:(length(order) - 1)
        left_index, right_index = order[index], order[index + 1]
        left, right = transitions[left_index], transitions[right_index]
        (Int(left.phase_lo), Int(left.phase_hi)) == (1, 0) || continue
        (Int(right.phase_lo), Int(right.phase_hi)) == (0, 1) || continue
        Bool(get(left, :converged, true)) && Bool(get(right, :converged, true)) || continue
        left_id, right_id = ids[left_index], ids[right_index]
        left_id > 0 && right_id > 0 && left_id != right_id || continue
        left_prediction = get(prediction_map, left_id, nothing)
        right_prediction = get(prediction_map, right_id, nothing)
        left_prediction === nothing && continue
        right_prediction === nothing && continue
        get(left_prediction, :confidence, false) && get(right_prediction, :confidence, false) || continue
        left_prediction.orientation == (1, 0) && right_prediction.orientation == (0, 1) || continue
        min_padding = Float64(get(cfg, :boundary_padding_min, 0.03))
        left_pad = max(min_padding, Float64(left_prediction.band_half_width))
        right_pad = max(min_padding, Float64(right_prediction.band_half_width))
        start = max(lower, Float64(left.u), Float64(left_prediction.u) + left_pad)
        stop = min(upper, Float64(right.u), Float64(right_prediction.u) - right_pad)
        start < stop && push!(exclusions, (start, stop))
    end
    return exclusions
end

end