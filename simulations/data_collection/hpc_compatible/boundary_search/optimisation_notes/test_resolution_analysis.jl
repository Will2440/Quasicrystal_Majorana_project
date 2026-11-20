using Test
include("resolution_analysis.jl")
using .ResolutionAnalysis

@testset "Nested grid and offset diagnostics" begin
    labels = [1, 1, 0, 1, 1, 1, 0, 1, 1]
    @test grid_indices(9, 4) == [1, 5, 9]
    @test grid_indices(9, 4, 2) == [1, 3, 7, 9]
    @test length(phase_runs(labels)) == 5
    @test grid_metrics(labels, 4, 0, 0.1).phase0_runs == 0
    @test grid_metrics(labels, 4, 0, 0.1).missed_fine_runs == 2
    @test grid_metrics(labels, 4, 2, 0.1).phase0_runs == 1
    @test grid_metrics(labels, 4, 2, 0.1).missed_fine_runs == 1
    @test grid_metrics(labels, 1, 0, 0.1).missed_fine_runs == 0
    @test_throws ErrorException grid_indices(8, 4)
    rows = [merge((stride=stride, offset=offset), grid_metrics(labels, stride, offset, 0.1))
        for (stride, offset) in ((4, 0), (4, 2), (2, 0), (2, 1), (1, 0))]
    @test coarsest_passing_stride(rows) === nothing
    @test coarsest_passing_stride([merge(row, (missed_fine_runs=0, phase0_runs=2,
        phase1_runs=3, transitions=4)) for row in rows]) == 4
end