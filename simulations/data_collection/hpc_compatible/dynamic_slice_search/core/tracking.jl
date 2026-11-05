module TransitionTracking

export EdgeTrack, predict_edges, update_tracks!

mutable struct EdgeTrack
    id::Int
    orientation::Tuple{Int,Int}
    history::Vector{NamedTuple}
    last_seen::Int
end

function predict_edges(tracks::Vector{EdgeTrack}, next_rho::Real, current_slice::Int, mu_bounds, cfg::Dict{Symbol,Any})
    predictions = NamedTuple[]
    min_history = Int(get(cfg, :momentum_min_history, 3))
    momentum_enabled = Bool(get(cfg, :momentum_enabled, true))
    damping = Float64(get(cfg, :momentum_damping, 0.2))
    max_shift_fraction = Float64(get(cfg, :momentum_max_shift_fraction, 0.02))
    active_gap = Int(get(cfg, :track_max_missed_slices, 1))
    mu_span = Float64(mu_bounds[2] - mu_bounds[1])

    for track in tracks
        isempty(track.history) && continue
        current_slice - track.last_seen > active_gap && continue
        last = track.history[end]
        predicted_mu = Float64(last.mu)

        if momentum_enabled && length(track.history) >= min_history
            prior = track.history[end - 2]
            delta_rho = Float64(last.rho) - Float64(prior.rho)
            if delta_rho > 0
                slope = (Float64(last.mu) - Float64(prior.mu)) / delta_rho
                raw_shift = damping * slope * (Float64(next_rho) - Float64(last.rho))
                shift_cap = max_shift_fraction * mu_span
                predicted_mu += clamp(raw_shift, -shift_cap, shift_cap)
            end
        end

        predicted_mu = clamp(predicted_mu, Float64(mu_bounds[1]), Float64(mu_bounds[2]))
        push!(predictions, (track_id=track.id, mu=predicted_mu, orientation=track.orientation))
    end
    return predictions
end

function update_tracks!(tracks::Vector{EdgeTrack}, transitions, rho::Real, slice_index::Int, mu_bounds, cfg::Dict{Symbol,Any})
    mu_span = Float64(mu_bounds[2] - mu_bounds[1])
    max_match = Float64(get(cfg, :track_match_max_fraction, 0.15)) * mu_span
    unmatched_tracks = Set(eachindex(tracks))
    next_id = isempty(tracks) ? 1 : maximum(t.id for t in tracks) + 1
    transition_ids = Int[]

    for transition in sort(collect(transitions); by=t -> t.mu)
        orientation = (Int(transition.phase_lo), Int(transition.phase_hi))
        best_idx = nothing
        best_distance = Inf

        for idx in unmatched_tracks
            track = tracks[idx]
            track.orientation == orientation || continue
            distance = abs(Float64(track.history[end].mu) - Float64(transition.mu))
            if distance < best_distance && distance <= max_match
                best_idx = idx
                best_distance = distance
            end
        end

        point = (rho=Float64(rho), mu=Float64(transition.mu))
        if best_idx === nothing
            push!(tracks, EdgeTrack(next_id, orientation, [point], slice_index))
            push!(transition_ids, next_id)
            next_id += 1
        else
            track = tracks[best_idx]
            push!(track.history, point)
            track.last_seen = slice_index
            push!(transition_ids, track.id)
            delete!(unmatched_tracks, best_idx)
        end
    end

    return transition_ids
end

end # module