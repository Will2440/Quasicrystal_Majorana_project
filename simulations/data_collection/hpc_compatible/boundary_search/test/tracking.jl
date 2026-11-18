using Test

include(joinpath(@__DIR__, "..", "core", "tracking.jl"))
using .TransitionTracking

transition(u, lo, hi; converged=true) = (u=u, phase_lo=lo, phase_hi=hi, converged=converged)
cfg = Dict{Symbol,Any}(:confidence_min_history => 2, :confidence_max_residual => 0.011,
    :boundary_padding_min => 0.02, :boundary_padding_max => 0.08,
    :momentum_min_history => 2, :momentum_damping => 1.0,
    :momentum_max_shift_fraction => 0.1, :track_match_base_tolerance_fraction => 0.05,
    :track_match_max_tolerance_fraction => 0.2, :track_max_missed_slices => 2)
bounds = (0.0, 1.0)

@testset "normalized predictions and confidence" begin
    tracks = EdgeTrack[]
    @test update_tracks!(tracks, [(u=0.2, phase_lo=1, phase_hi=0)], 0.0, 1, bounds, cfg) == [1]
    first_prediction = only(predict_edges(tracks, 1.0, 2, bounds, cfg))
    @test first_prediction.u == 0.2
    @test !first_prediction.confidence && first_prediction.band_half_width == 0.08
    update_tracks!(tracks, [transition(0.21, 1, 0)], 1.0, 2, bounds, cfg; predictions=[first_prediction])
    second_prediction = only(predict_edges(tracks, 2.0, 3, bounds, cfg))
    @test second_prediction.u ≈ 0.22
    @test !second_prediction.confidence
    update_tracks!(tracks, [transition(0.22, 1, 0)], 2.0, 3, bounds, cfg; predictions=[second_prediction])
    third_prediction = only(predict_edges(tracks, 4.0, 4, bounds, cfg))
    @test third_prediction.confidence && third_prediction.band_half_width ≈ 0.02
    @test third_prediction.u ≈ 0.24
    @test third_prediction.predicted_motion ≈ 0.02
    capped = only(predict_edges(tracks, 100.0, 4, bounds, cfg))
    @test capped.u ≈ 0.32
    @test capped.band_half_width == 0.08
    @test capped.motion_uncertainty > 0
    update_tracks!(tracks, [transition(0.24, 1, 0; converged=false)], 4.0, 4, bounds, cfg; predictions=[third_prediction])
    @test !only(predict_edges(tracks, 5.0, 5, bounds, cfg)).confidence
end

@testset "misses, residuals and stale predictions" begin
    tracks = EdgeTrack[]
    update_tracks!(tracks, [transition(0.3, 1, 0)], 0.0, 1, bounds, cfg)
    for slice in 2:3
        prediction = predict_edges(tracks, slice - 1, slice, bounds, cfg)
        update_tracks!(tracks, [transition(0.3, 1, 0)], slice - 1, slice, bounds, cfg; predictions=prediction)
    end
    @test only(predict_edges(tracks, 3.0, 4, bounds, cfg)).confidence
    update_tracks!(tracks, [], 3.0, 4, bounds, cfg)
    @test !only(predict_edges(tracks, 4.0, 5, bounds, cfg)).confidence
    @test isempty(predict_edges(tracks, 5.0, 6, bounds, cfg))
    stale = predict_edges(tracks, 4.0, 5, bounds, cfg)
    update_tracks!(tracks, [transition(0.3, 1, 0)], 4.0, 5, bounds, cfg; predictions=stale)
    update_tracks!(tracks, [transition(0.3, 1, 0)], 5.0, 6, bounds, cfg; predictions=stale)
    @test !only(predict_edges(tracks, 6.0, 7, bounds, cfg)).confidence
    prediction = predict_edges(tracks, 6.0, 7, bounds, cfg)
    update_tracks!(tracks, [transition(0.32, 1, 0)], 6.0, 7, bounds, cfg; predictions=prediction)
    @test !only(predict_edges(tracks, 7.0, 8, bounds, cfg)).confidence
end

@testset "confidence requires validated observations" begin
    tracks = EdgeTrack[]
    update_tracks!(tracks, [transition(0.3, 1, 0)], 0.0, 1, bounds, cfg)
    for slice in 2:5
        update_tracks!(tracks, [transition(0.3, 1, 0)], slice - 1, slice, bounds, cfg)
    end
    @test !only(predict_edges(tracks, 5.0, 6, bounds, cfg)).confidence
    for slice in 6:7
        predictions = predict_edges(tracks, slice - 1, slice, bounds, cfg)
        update_tracks!(tracks, [transition(0.3, 1, 0)], slice - 1, slice, bounds, cfg; predictions=predictions)
    end
    @test only(predict_edges(tracks, 7.0, 8, bounds, cfg)).confidence
    predictions = predict_edges(tracks, 7.0, 8, bounds, cfg)
    update_tracks!(tracks, [transition(0.32, 1, 0)], 7.0, 8, bounds, cfg; predictions=predictions)
    @test !only(predict_edges(tracks, 8.0, 9, bounds, cfg)).confidence
end

@testset "ordered association and ambiguity" begin
    tracks = EdgeTrack[]
    original = [transition(0.2, 1, 0), transition(0.8, 0, 1)]
    @test update_tracks!(tracks, original, 0.0, 1, bounds, cfg) == [1, 2]
    prediction = predict_edges(tracks, 1.0, 2, bounds, cfg)
    inserted = [transition(0.8, 0, 1), transition(0.5, 1, 0), transition(0.2, 1, 0)]
    @test update_tracks!(tracks, inserted, 1.0, 2, bounds, cfg; predictions=prediction) == [2, 3, 1]
    prediction = predict_edges(tracks, 2.0, 3, bounds, cfg)
    ambiguous = [transition(0.2, 1, 0), transition(0.21, 1, 0), transition(0.8, 0, 1)]
    update_tracks!(tracks, ambiguous, 2.0, 3, bounds, cfg; predictions=prediction)
    @test !only(filter(pred -> pred.track_id == 1, predict_edges(tracks, 3.0, 4, bounds, cfg))).confidence
end

@testset "ambiguity resets confidence" begin
    tracks = EdgeTrack[]
    update_tracks!(tracks, [transition(0.2, 1, 0)], 0.0, 1, bounds, cfg)
    for slice in 2:3
        predictions = predict_edges(tracks, slice - 1, slice, bounds, cfg)
        update_tracks!(tracks, [transition(0.2, 1, 0)], slice - 1, slice, bounds, cfg; predictions=predictions)
    end
    predictions = predict_edges(tracks, 3.0, 4, bounds, cfg)
    @test only(predictions).confidence
    ids = update_tracks!(tracks, [transition(0.2, 1, 0), transition(0.21, 1, 0)], 3.0, 4, bounds, cfg; predictions=predictions)
    @test 1 in ids
    @test !only(filter(prediction -> prediction.track_id == 1, predict_edges(tracks, 4.0, 5, bounds, cfg))).confidence
end

@testset "descending variable steps" begin
    tracks = EdgeTrack[]
    update_tracks!(tracks, [transition(0.2, 1, 0)], 3.0, 1, bounds, cfg)
    prediction = predict_edges(tracks, 2.0, 2, bounds, cfg)
    update_tracks!(tracks, [transition(0.21, 1, 0)], 2.0, 2, bounds, cfg; predictions=prediction)
    prediction = only(predict_edges(tracks, 0.0, 3, bounds, cfg))
    @test prediction.u ≈ 0.23
    @test prediction.predicted_motion ≈ 0.02
    @test only(predict_edges(tracks, 3.0, 3, bounds, cfg)).u ≈ 0.21
end

@testset "robust recent slope" begin
    tracks = [EdgeTrack(1, (1, 0), NamedTuple[(sweep_value=0.0, u=0.2),
        (sweep_value=1.0, u=0.21), (sweep_value=2.0, u=0.31),
        (sweep_value=3.0, u=0.32)], 4)]
    @test only(predict_edges(tracks, 4.0, 5, bounds, cfg)).u ≈ 0.33
    push!(tracks[1].history, (sweep_value=2.0, u=0.30))
    tracks[1].last_seen = 5
    @test only(predict_edges(tracks, 2.5, 6, bounds, cfg)).u ≈ 0.30
    @test only(predict_edges(tracks, 1.0, 6, bounds, cfg)).u ≈ 0.28
end

@testset "predicted exclusions" begin
    transitions = [transition(0.8, 0, 1), transition(0.2, 1, 0)]
    predictions = [(track_id=1, u=0.25, orientation=(1, 0), band_half_width=0.03, confidence=true),
        (track_id=2, u=0.75, orientation=(0, 1), band_half_width=0.04, confidence=true)]
    interval = only(predicted_exclusions(predictions, transitions, [2, 1], cfg, bounds))
    @test interval[1] ≈ 0.28 && interval[2] ≈ 0.71
    outward = [merge(predictions[1], (u=0.1,)), merge(predictions[2], (u=0.9,))]
    @test only(predicted_exclusions(outward, transitions, [2, 1], cfg, bounds)) == (0.2, 0.8)
    overlap = [merge(predictions[1], (u=0.6, band_half_width=0.1)), merge(predictions[2], (u=0.7, band_half_width=0.1))]
    @test isempty(predicted_exclusions(overlap, transitions, [2, 1], cfg, bounds))
    @test isempty(predicted_exclusions([merge(predictions[1], (confidence=false,)), predictions[2]], transitions, [2, 1], cfg, bounds))
    @test isempty(predicted_exclusions(predictions[1:1], transitions, [2, 1], cfg, bounds))
    @test isempty(predicted_exclusions(predictions, [transitions[1], transition(0.5, 0, 1), transitions[2]], [2, 3, 1], cfg, bounds))
    @test isempty(predicted_exclusions(predictions, [transitions[1], transition(0.2, 1, 0; converged=false)], [2, 1], cfg, bounds))
    @test isempty(predicted_exclusions(predictions, transitions, [2, 0], cfg, bounds))
    @test isempty(predicted_exclusions(predictions, [], Int[], cfg, bounds))
    @test only(predicted_exclusions(predictions, transitions, [2, 1], cfg, (0.4, 0.6))) == (0.4, 0.6)
    @test_throws ArgumentError predicted_exclusions(predictions, transitions, [1], cfg, bounds)
end