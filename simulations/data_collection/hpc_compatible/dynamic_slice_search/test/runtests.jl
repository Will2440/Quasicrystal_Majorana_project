using Test
using Statistics

include(joinpath(@__DIR__, "..", "core", "slice_search.jl"))
include(joinpath(@__DIR__, "..", "core", "tracking.jl"))

using .SliceSearch
using .TransitionTracking

sample(mu, phase; gap=0.1) = (mu=Float64(mu), phase_label=Int(phase), gap_min=Float64(gap))

@testset "Slice search" begin
    linear = SliceSearch.build_rho_grid(1.0, 3.0, 5, :linear)
    @test linear ≈ [1.0, 1.5, 2.0, 2.5, 3.0]
    logarithmic = SliceSearch.build_rho_grid(1.0, 9.0, 3, :log)
    @test logarithmic ≈ [1.0, 3.0, 9.0]
    @test_throws ErrorException SliceSearch.build_rho_grid(0.0, 3.0, 5, :log)

    grid = SliceSearch.baseline_mu_grid(-1.0, 1.0, 5)
    @test grid ≈ [-1.0, -0.5, 0.0, 0.5, 1.0]

    samples = [sample(-1, 0), sample(0, 1), sample(1, 0), sample(2, 0), sample(3, 1)]
    transitions = SliceSearch.detect_transitions(samples)
    @test length(transitions) == 3
    @test [(t.phase_lo, t.phase_hi) for t in transitions] == [(0, 1), (1, 0), (0, 1)]
    @test isempty(SliceSearch.detect_transitions([sample(0, 0)]))

    eval_at_mu = mu -> sample(mu, mu < 0.25 ? 0 : 1; gap=abs(mu - 0.25))
    bracket = SliceSearch.refine_transition_bracket(
        eval_at_mu, sample(0.0, 0), sample(1.0, 1), 1e-3, 30)
    @test bracket.mu_lo <= 0.25 <= bracket.mu_hi
    @test bracket.mu_hi - bracket.mu_lo <= 1e-3
    @test bracket.converged

    low_gap_eval = mu -> sample(mu, 1; gap=abs(mu - 0.5))
    added = SliceSearch.refine_low_gap_intervals!(
        low_gap_eval, [sample(0.0, 1; gap=0.5), sample(1.0, 1; gap=0.5)],
        0.51, 0.1, 2)
    @test !isempty(added)
    @test all(s -> s.mu == 0.5, added[1:1])
end

@testset "Conservative transition tracking" begin
    cfg = Dict{Symbol,Any}(
        :mu_min => 0.0,
        :mu_max => 10.0,
        :momentum_min_history => 3,
        :momentum_damping => 0.2,
        :momentum_max_shift_fraction => 0.02,
        :prediction_band_half_width_fraction => 0.03,
        :track_match_max_fraction => 0.5,
        :track_max_missed_slices => 1,
        :momentum_enabled => true,
    )
    tracks = TransitionTracking.EdgeTrack[]
    transition(mu) = (mu=mu, phase_lo=0, phase_hi=1)

    ids1 = TransitionTracking.update_tracks!(tracks, [transition(2.0)], 1.0, 1, (0.0, 10.0), cfg)
    ids2 = TransitionTracking.update_tracks!(tracks, [transition(2.5)], 1.1, 2, (0.0, 10.0), cfg)
    @test ids1 == ids2

    # Before three observed positions, prediction stays at the latest position.
    pred = TransitionTracking.predict_edges(tracks, 1.2, 3, (0.0, 10.0), cfg)
    @test pred[1].mu == 2.5

    TransitionTracking.update_tracks!(tracks, [transition(3.0)], 1.2, 3, (0.0, 10.0), cfg)
    pred = TransitionTracking.predict_edges(tracks, 1.3, 4, (0.0, 10.0), cfg)
    @test length(pred) == 1
    @test 3.0 < pred[1].mu <= 3.2

    cfg[:momentum_enabled] = false
    no_momentum = TransitionTracking.predict_edges(tracks, 1.4, 4, (0.0, 10.0), cfg)
    @test no_momentum[1].mu == 3.0
end