using Test
include(joinpath(@__DIR__, "..", "core", "coverage.jl"))
using .BoundaryCoverage

@testset "Physical spacing target" begin
    @test baseline_point_count((0.0, 1.0), 11, 2.0, nothing) == 11
    @test baseline_point_count((0.0, 1.0), 11, 2.0, 0.1) == 21
    @test baseline_point_count((0.0, 1.0), 11, 4.0, 0.1) == 41
    @test baseline_point_count((0.0, 1.0), 51, 2.0, 0.1) == 51
    @test_throws ErrorException baseline_point_count((0.0, 1.0), 11, 2.0, 0.0)
end

@testset "Normalized sampling coverage" begin
    baseline = plan_grid((0.0, 1.0), 11, [], 0)
    @test baseline.points ≈ collect(0.0:0.1:1.0)
    focused = plan_grid((0.0, 1.0), 11, [(0.2, 0.6)], 0)
    @test 0.2 in focused.points && 0.6 in focused.points
    @test !any(point -> 0.2 < point < 0.6, focused.points)
    @test length(focused.points) >= length(baseline.points)
    @test maximum(diff(filter(point -> point <= 0.2, focused.points))) < 0.1
    @test maximum(diff(filter(point -> point >= 0.6, focused.points))) < 0.1
    audited = plan_grid((0.0, 1.0), 11, [(0.2, 0.6)], 1)
    @test any(point -> 0.2 < point < 0.6, audited.points)
    @test plan_grid((0.0, 1.0), 11, [(0.2, 0.5), (0.4, 0.7)], 0).exclusions == [(0.2, 0.7)]
    shifted = plan_grid((0.0, 1.0), 11, [(0.23, 0.57)], 0)
    @test 0.23 in shifted.points && 0.57 in shifted.points
    @test !any(point -> 0.23 < point < 0.57, shifted.points)
end