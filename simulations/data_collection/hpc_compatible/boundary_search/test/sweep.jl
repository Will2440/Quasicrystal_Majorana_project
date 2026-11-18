using Test

include(joinpath(@__DIR__, "..", "core", "slice_search.jl"))
module DynamicBandwidth
end
include(joinpath(@__DIR__, "..", "core", "sweep.jl"))
using .SliceSearch
using .SliceSweep

@testset "rho grids" begin
    for (start, stop) in ((1.1, 4.0), (0.9, 0.1), (1.0, 4.0), (1.0, 0.1)), spacing in (:linear, :log)
        grid = build_rho_grid(start, stop, 11, spacing)
        @test length(grid) == 11
        @test first(grid) ≈ start
        @test last(grid) ≈ stop
        @test all(diff(grid) .* sign(stop - start) .> 0)
        @test build_rho_grid(start, start, 1, spacing) ≈ [start]
        @test build_rho_grid(start, stop, 1, spacing) ≈ [start]
    end
    @test build_rho_grid(0.9, 0.1, 5, :linear) ≈ [0.9, 0.7, 0.5, 0.3, 0.1]
    positions = range(0.0, 1.0; length=11)
    fractions = expm1.(4.0 .* positions) ./ expm1(4.0)
    expected_low = exp.(log(0.1) .+ (log(0.9) - log(0.1)) .* fractions)
    expected_high = exp.(log(1.1) .+ (log(4.0) - log(1.1)) .* (1 .- reverse(fractions)))
    @test build_rho_grid(0.9, 0.1, 11, :log) ≈ reverse(expected_low)
    @test build_rho_grid(1.1, 4.0, 11, :log) ≈ expected_high
    for (start, stop, count) in ((0.1, 0.9, 2), (4.0, 1.1, 2), (0.9, 1.1, 2),
                                  (1.1, 0.9, 2), (0.9, 1.0, 2),
                                  (0.1, 0.0, 2), (1.1, Inf, 2),
                                  (NaN, 2.0, 2), (1.1, 1.1, 2), (1.1, 2.0, 0))
        @test_throws ErrorException build_rho_grid(start, stop, count)
    end
    @test build_rho_grid(1.0, 1.0, 1) == [1.0]
    @test_throws ErrorException build_rho_grid(1.1, 2.0, 3, :log; concentration=0)
    @test_throws ErrorException build_rho_grid(1.1, 2.0, 3, :unknown)
    @test rho_progress(0.5, 0.9, 0.1) ≈ 0.5
    @test rho_progress(1.1, 1.1, 1.1) == 0.0
end

@testset "Delta grids" begin
    for spacing in (:linear, :log, :exponential)
        grid = build_delta_grid(1.0, 0.01, 11, spacing)
        @test length(grid) == 11
        @test first(grid) ≈ 1.0
        @test last(grid) ≈ 0.01
        @test all(diff(grid) .< 0)
        @test build_delta_grid(0.1, 0.1, 1, spacing) ≈ [0.1]
        @test build_delta_grid(1.0, 0.1, 1, spacing) ≈ [1.0]
    end
        @test build_delta_grid(1.0, 0.1, 5, :linear) ≈ collect(range(1.0, 0.1; length=5))
        @test build_delta_grid(1.0, 0.1, 5, :log) ≈ exp.(range(log(1.0), log(0.1); length=5))
        positions = range(1.0, 0.0; length=5)
        @test build_delta_grid(1.0, 0.1, 5, :exponential) ≈
                    0.1 .+ 0.9 .* expm1.(4.0 .* positions) ./ expm1(4.0)
    for (start, stop, count) in ((0.1, 1.0, 2), (1.0, 0.0, 2), (1.0, -0.1, 2),
                                  (Inf, 0.1, 2), (1.0, NaN, 2), (0.1, 0.1, 2), (1.0, 0.1, 0))
        @test_throws ErrorException build_delta_grid(start, stop, count)
    end
    @test_throws ErrorException build_delta_grid(1.0, 0.1, 3, :exponential; concentration=0)
    @test_throws ErrorException build_delta_grid(1.0, 0.1, 3, :unknown)
end

@testset "configuration and sweep dispatch" begin
    cfg = include(joinpath(@__DIR__, "..", "config", "default_config.jl"))
    @test isabspath(cfg[:sequence_bson_path])
    @test cfg[:skip_trivial_interiors] === false
    @test sweep_grid(cfg) == build_rho_grid(cfg[:rho_start], cfg[:rho_max], cfg[:n_rho_slices], cfg[:rho_spacing_mode])
    @test sweep_progress(cfg[:rho_start], cfg) == 0.0
    @test sweep_progress(cfg[:rho_max], cfg) == 1.0
    cfg[:rho_start], cfg[:rho_max] = 0.9, 0.1
    @test sweep_progress(0.5, cfg) ≈ 0.5
    @test issorted(sweep_grid(cfg); rev=true)
    cfg[:sweep_parameter] = :Delta
    @test issorted(sweep_grid(cfg); rev=true)
    @test sweep_progress(cfg[:Delta_start], cfg) == 0.0
    @test sweep_progress(cfg[:Delta_max], cfg) == 1.0
    cfg[:sweep_parameter] = :unknown
    @test_throws ErrorException sweep_grid(cfg)
end