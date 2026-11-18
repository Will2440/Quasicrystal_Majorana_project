using Test
include(joinpath(@__DIR__, "..", "core", "coverage.jl"))
using .BoundaryCoverage

@testset "Normalized sampling coverage" begin
    baseline = plan_grid((0.0, 1.0), 11, [], 0)
    @test baseline.points ≈ collect(0.0:0.1:1.0)
    focused = plan_grid((0.0, 1.0), 11, [(0.2, 0.6)], 0)
    @test 0.2 in focused.points && 0.6 in focused.points
    @test !any(point -> 0.2 < point < 0.6, focused.points)
    @test length(focused.points) >= 11
    audited = plan_grid((0.0, 1.0), 11, [(0.2, 0.6)], 1)
    @test any(point -> 0.2 < point < 0.6, audited.points)
    @test plan_grid((0.0, 1.0), 11, [(0.2, 0.5), (0.4, 0.7)], 0).exclusions == [(0.2, 0.7)]
end