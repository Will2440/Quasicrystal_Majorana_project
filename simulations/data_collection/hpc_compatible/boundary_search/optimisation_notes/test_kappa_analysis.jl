using Test
include("kappa_analysis.jl")
using .KappaAnalysis

@testset "Kappa windows and label agreement" begin
    windows = positive_windows([0, 1, 1, 0, 1, 1])
    @test length(windows) == 2
    @test windows[1] == (first=2, last=3, lower_censored=false, upper_censored=false)
    @test windows[2].upper_censored
    @test isempty(positive_windows([0, 0]))
    @test only(positive_windows([1, 1])).lower_censored
    @test agreement_counts([1, 0, 1], [0, 1, 1]).topological_lost == 1
    @test agreement_counts([1, 0, 1], [0, 1, 1]).trivial_changed == 1
    @test_throws ArgumentError agreement_counts([1], [1, 0])
    @test aligned_changes([0.0, 1.0], [1, 0], [0.0, 1.0], [0, 1]) == [1, 2]
    @test_throws ArgumentError aligned_changes([0.0, 1.0], [1, 0], [0.0, 0.9], [1, 0])
    @test validate_comparable_configs(Dict(:N => 500), Dict(:N => 500)) === nothing
    @test_throws ArgumentError validate_comparable_configs(Dict(:N => 500), Dict(:N => 300))
end