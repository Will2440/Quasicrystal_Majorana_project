using Test
using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "slice_data.jl"))
using .SliceBoundaryData

@testset "Slice boundary data readers" begin
    mktempdir() do run_dir
        slices_dir = joinpath(run_dir, "slices")
        slice_dir = joinpath(slices_dir, "slice_0001")
        mkpath(slice_dir)

        CSV.write(joinpath(run_dir, "transition_boundaries.csv"), DataFrame(
            slice_index=[1], rho=[1.0], track_id=[1], transition_index=[1],
            mu=[0.5], mu_lo=[0.49], mu_hi=[0.51], phase_lo=[0], phase_hi=[1],
            gap_min=[1e-5], refinement_iters=[8], converged=[true],
        ))
        CSV.write(joinpath(run_dir, "slice_summary.csv"), DataFrame(
            slice_index=[1], rho=[1.0], gap_count=[1], n_transitions=[1],
            n_points=[20], n_baseline=[15], n_prediction_points=[5],
            n_gap_refinement_points=[0], eval_seconds=[1.0],
        ))
        CSV.write(joinpath(slice_dir, "point_results.csv"), DataFrame(
            point_key=["1|0.4", "1|0.6"], rho=[1.0, 1.0], mu=[0.4, 0.6],
            u=[0.4, 0.6], t1=[1.0, 1.0], t2=[1.0, 1.0],
            phase_label=[0, 1], topo_fraction=[0.0, 1.0],
            invariant_mean=[0.0, 1.0], invariant_min=[0.0, 1.0],
            invariant_max=[0.0, 1.0], gap_min=[0.1, 0.1],
            gap_mean=[0.1, 0.1], gap_max=[0.1, 0.1], n_samples=[5, 5],
        ))

        boundaries = SliceBoundaryData.load_transitions(run_dir)
        summary = SliceBoundaryData.load_slice_summary(run_dir)
        samples = SliceBoundaryData.load_samples(run_dir)

        @test nrow(boundaries) == 1
        @test boundaries.mu[1] == 0.5
        @test summary.gap_count[1] == 1
        @test nrow(samples) == 2
        @test all(samples.slice_index .== 1)
        @test first(SliceBoundaryData.available_slices(run_dir))[1] == 1
    end
end