using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "processing.jl"))
using .IDOPProcessing

@testset "Outer normalization reference mu" begin
    mktempdir() do directory
        cache_path = joinpath(directory, "bounds.csv")
        CSV.write(cache_path, DataFrame(
            phi=[0.25, 0.25, 0.75, 0.75],
            mu=[4.0, 0.0, 0.0, 3.0],
            global_min=[-8.0, -4.0, -2.0, -5.0],
            global_max=[8.0, 4.0, 2.0, 5.0],
        ))
        data = DataFrame(phi=[0.25, 0.75], mu_range=[[-4.0, 0.0, 4.0], [-2.0, 0.0, 2.0]])

        IDOPProcessing.apply_outer_normalisation!(
            data, cache_path; target_col=:mu_range, out_col=:normalized, reference_mu=0.0,
        )
        @test data.normalized[1] ≈ [-1.0, 0.0, 1.0]
        @test data.normalized[2] ≈ [-1.0, 0.0, 1.0]

        IDOPProcessing.apply_outer_normalisation!(
            data, cache_path; target_col=:mu_range, out_col=:legacy, verbose_warnings=false,
        )
        @test data.legacy[1] ≈ [-0.5, 0.0, 0.5]
        @test_throws ErrorException IDOPProcessing.apply_outer_normalisation!(
            data, cache_path; target_col=:mu_range, reference_mu=1.0,
        )
    end
end