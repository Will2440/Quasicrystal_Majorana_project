using Test
using Statistics

include(joinpath(@__DIR__, "..", "core", "slice_search.jl"))
include(joinpath(@__DIR__, "..", "core", "tracking.jl"))
include(joinpath(@__DIR__, "..", "core", "bandwidth.jl"))
include(joinpath(@__DIR__, "..", "core", "sweep.jl"))

using .SliceSearch
using .TransitionTracking

sample(mu, phase; gap=0.1) = (mu=Float64(mu), phase_label=Int(phase), gap_min=Float64(gap))

@testset "Sweep parameters" begin
    cfg = Dict{Symbol,Any}(
        :sweep_mode => :vary_t2, :t1_fixed => 1.0, :t2_fixed => 1.5,
        :Delta => 0.05, :rho_start => 1.0, :rho_max => 2.0,
        :n_rho_slices => 3, :rho_spacing_mode => :linear,
        :Delta_start => 0.0, :Delta_max => 0.2, :n_Delta_slices => 3,
        :Delta_spacing_mode => :linear)
    @test SliceSweep.sweep_grid(cfg) ≈ [1.0, 1.5, 2.0]
    rho_params = SliceSweep.physical_parameters(2.0, cfg, 0.5)
    @test (rho_params.rho, rho_params.Delta, rho_params.t1, rho_params.t2) == (2.0, 0.05, 1.0, 2.0)
    cfg[:sweep_parameter] = :Delta
    @test SliceSweep.sweep_grid(cfg) ≈ [0.0, 0.1, 0.2]
    cfg[:n_Delta_slices] = 1
    @test SliceSweep.sweep_grid(cfg) == [0.0]
    cfg[:n_Delta_slices] = 3
    critical = SliceSweep.fixed_mu_critical(cfg, 0.5)
    first = SliceSweep.physical_parameters(0.0, cfg, 0.5, critical)
    last = SliceSweep.physical_parameters(0.2, cfg, 0.5, critical)
    @test (first.rho, first.t1, first.t2, first.Delta) == (1.5, 1.0, 1.5, 0.0)
    @test last.Delta == 0.2
    @test first.mu_critical == last.mu_critical == 2.5
    cfg[:Delta_spacing_mode] = :log
    @test_throws ErrorException SliceSweep.sweep_grid(cfg)
    cfg[:Delta_start] = 0.01
    @test SliceSweep.sweep_grid(cfg) ≈ [0.01, sqrt(0.002), 0.2]
    cfg[:t1_fixed] = 0.0
    @test_throws ErrorException SliceSweep.fixed_mu_critical(cfg, 0.5)
end

@testset "Slice search" begin
    linear = SliceSearch.build_rho_grid(1.0, 3.0, 5, :linear)
    @test linear ≈ [1.0, 1.5, 2.0, 2.5, 3.0]
    @test SliceSearch.build_rho_grid(1.0, 3.0, 1, :linear) == [1.0]
    @test SliceSearch.build_rho_grid(1.0, 3.0, 1, :log) == [1.0]
    logarithmic = SliceSearch.build_rho_grid(1.0, 9.0, 3, :log)
    @test logarithmic ≈ [1.0, 3.0, 9.0]
    @test_throws ErrorException SliceSearch.build_rho_grid(0.0, 3.0, 5, :log)
    @test SliceSearch.build_delta_grid(0.0, 0.2, 3, :linear) ≈ [0.0, 0.1, 0.2]
    @test SliceSearch.build_delta_grid(0.01, 1.0, 3, :log) ≈ [0.01, 0.1, 1.0]
    @test_throws ErrorException SliceSearch.build_delta_grid(0.0, 1.0, 3, :log)
    @test_throws ErrorException SliceSearch.build_delta_grid(-0.1, 1.0, 3, :linear)
    @test SliceSearch.rho_progress(1.0, 1.0, 3.0) == 0.0
    @test SliceSearch.rho_progress(2.0, 1.0, 3.0) == 0.5
    @test SliceSearch.rho_progress(4.0, 1.0, 3.0) == 1.0

    grid = SliceSearch.baseline_mu_grid(-1.0, 1.0, 5)
    @test grid ≈ [-1.0, -0.5, 0.0, 0.5, 1.0]
    @test SliceSearch.scaled_baseline_point_count(2.0, 2.0, 201) == 201
    @test SliceSearch.scaled_baseline_point_count(4.0, 2.0, 201) == 401
    @test SliceSearch.scaled_baseline_point_count(2.5, 2.0, 201) == 251
    @test SliceSearch.scaled_baseline_point_count(
        4.0, 2.0, 201; rho_fraction=0.0, rho_growth_factor=1.25, rho_growth_power=2.0) == 401
    @test SliceSearch.scaled_baseline_point_count(
        4.0, 2.0, 201; rho_fraction=1.0, rho_growth_factor=1.25, rho_growth_power=2.0) == 501
    @test SliceSearch.boundary_points_per_prediction(9, 25, 0.0, 1.5) == 9
    @test SliceSearch.boundary_points_per_prediction(9, 25, 1.0, 1.5) == 25
    @test 9 < SliceSearch.boundary_points_per_prediction(9, 25, 0.5, 1.5) < 25
    @test_throws ErrorException SliceSearch.scaled_baseline_point_count(2.0, 0.0, 201)

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
        :momentum_min_history => 3,
        :momentum_damping => 0.2,
        :momentum_max_shift_fraction => 0.02,
        :prediction_band_min_fraction => 0.03,
        :prediction_band_max_fraction => 0.08,
        :prediction_band_motion_multiplier => 2.5,
        :prediction_band_curvature_multiplier => 1.5,
        :track_match_base_tolerance_fraction => 0.01,
        :track_match_max_tolerance_fraction => 0.05,
        :track_match_step_scale => 1.5,
        :track_unmatched_penalty => 0.75,
        :track_max_missed_slices => 1,
        :momentum_enabled => true,
    )
    tracks = TransitionTracking.EdgeTrack[]
    transition(mu, phase_lo=0, phase_hi=1) = (mu=mu, phase_lo=phase_lo, phase_hi=phase_hi)

    ids1 = TransitionTracking.update_tracks!(tracks, [transition(2.0)], 1.0, 1, (0.0, 10.0), cfg)
    predictions = TransitionTracking.predict_edges(tracks, 1.1, 2, (0.0, 10.0), cfg)
    ids2 = TransitionTracking.update_tracks!(tracks, [transition(2.05)], 1.1, 2, (0.0, 10.0), cfg; predictions)
    @test ids1 == ids2

    # Before three observed positions, prediction stays at the latest position.
    pred = TransitionTracking.predict_edges(tracks, 1.2, 3, (0.0, 10.0), cfg)
    @test pred[1].mu == 2.05
    @test pred[1].band_half_width == 0.3

    TransitionTracking.update_tracks!(tracks, [transition(2.1)], 1.2, 3, (0.0, 10.0), cfg; predictions=pred)
    pred = TransitionTracking.predict_edges(tracks, 1.3, 4, (0.0, 10.0), cfg)
    @test length(pred) == 1
    @test 2.1 < pred[1].mu <= 2.12

    cfg[:momentum_enabled] = false
    no_momentum = TransitionTracking.predict_edges(tracks, 1.4, 4, (0.0, 10.0), cfg)
    @test no_momentum[1].mu == 2.1
end

@testset "Motion-informed prediction bands" begin
    cfg = Dict{Symbol,Any}(
        :momentum_min_history => 3,
        :momentum_damping => 0.2,
        :momentum_max_shift_fraction => 0.02,
        :prediction_band_min_fraction => 0.03,
        :prediction_band_max_fraction => 0.08,
        :prediction_band_motion_multiplier => 2.5,
        :prediction_band_curvature_multiplier => 1.5,
        :track_max_missed_slices => 1,
        :momentum_enabled => true,
    )
    track = TransitionTracking.EdgeTrack(
        1,
        (0, 1),
        [(rho=1.0, mu=1.0), (rho=1.1, mu=1.1), (rho=1.2, mu=1.2), (rho=1.3, mu=1.3)],
        4,
    )
    prediction = only(TransitionTracking.predict_edges([track], 1.4, 5, (0.0, 10.0), cfg))
    @test prediction.mu ≈ 1.32
    @test prediction.predicted_motion ≈ 0.1
    @test prediction.band_half_width ≈ 0.3
    @test prediction.curvature_uncertainty ≈ 0.0
end

@testset "Ordered matching protects established boundaries" begin
    cfg = Dict{Symbol,Any}(
        :momentum_min_history => 3,
        :momentum_damping => 0.2,
        :momentum_max_shift_fraction => 0.02,
        :prediction_band_min_fraction => 0.03,
        :prediction_band_max_fraction => 0.08,
        :prediction_band_motion_multiplier => 2.5,
        :prediction_band_curvature_multiplier => 1.5,
        :track_match_base_tolerance_fraction => 0.01,
        :track_match_max_tolerance_fraction => 0.05,
        :track_match_step_scale => 1.5,
        :track_unmatched_penalty => 0.75,
        :track_max_missed_slices => 1,
        :momentum_enabled => true,
    )
    tracks = TransitionTracking.EdgeTrack[]
    rising(mu) = (mu=mu, phase_lo=0, phase_hi=1)
    falling(mu) = (mu=mu, phase_lo=1, phase_hi=0)

    TransitionTracking.update_tracks!(tracks, [rising(2.0), falling(8.0)], 1.0, 1, (0.0, 10.0), cfg)
    TransitionTracking.update_tracks!(tracks, [rising(2.05), falling(8.05)], 1.1, 2, (0.0, 10.0), cfg)
    predictions = TransitionTracking.predict_edges(tracks, 1.2, 3, (0.0, 10.0), cfg)
    ids = TransitionTracking.update_tracks!(
        tracks,
        [rising(6.5), falling(8.1)],
        1.2,
        3,
        (0.0, 10.0),
        cfg;
        predictions,
    )

    @test ids[2] == 2
    @test ids[1] != 2
    @test tracks[2].history[end].mu == 8.1
    @test tracks[2].history[end].rho == 1.2
    @test tracks[1].last_seen == 2
end