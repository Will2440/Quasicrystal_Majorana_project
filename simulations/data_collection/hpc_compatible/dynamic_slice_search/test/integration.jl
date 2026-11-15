using Test
using BSON
using CSV
using DataFrames
using JLD2

include(joinpath(@__DIR__, "..", "master_slice.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "..", "data_processing", "boundary_processing", "slice_boundary_processing", "slice_data.jl"))
module SliceMergeTest
include(joinpath(@__DIR__, "..", "merge_slice_runs.jl"))
end

@testset "Real localizer smoke test" begin
    cfg = Dict{Symbol,Any}(:sweep_parameter => :Delta, :t1_fixed => 1.0, :t2_fixed => 1.5)
    ctx = (N=6, gamma=0.5, fixed_mu_critical=SliceSweep.fixed_mu_critical(cfg, 0.5),
           sequence=[1, 2, 1, 2, 1, 2], x0s=[3.0], kappas=[0.01],
           n_lowest_evals=2, kk_tol=1e-8, kk_maxiter=100,
           topological_threshold=0.5, topological_fraction=0.5)
    first = evaluate_point(0.0, 0.4, cfg, ctx)
    second = evaluate_point(0.3, 0.4, cfg, ctx)
    @test first.Delta == 0.0 && second.Delta == 0.3
    @test first.rho == second.rho == 1.5
    @test isfinite(first.gap_min) && isfinite(second.gap_min)
    @test !isapprox(first.gap_min, second.gap_min; atol=1e-5)
end

@eval DynamicPhysics function evaluate_specloc_point(
    N::Int, t_n::Vector{Float64}, mu::Float64, Delta::Float64,
    sequence::Vector{Int}, x0s::Vector{Float64}, kappas::Vector{Float64}; kwargs...)
    phase = mu > 2.0 || (Delta >= 0.1 && 0.5 <= mu <= 1.5) ? 0 : 1
    return (phase_label=phase, topo_fraction=Float64(phase), invariant_mean=Float64(phase),
            invariant_min=Float64(phase), invariant_max=Float64(phase), gap_min=0.1,
            gap_mean=0.1, gap_max=0.1, n_samples=1)
end

@testset "Full slice searches" begin
    mktempdir() do directory
        sequence_path = joinpath(directory, "sequence.bson")
        BSON.bson(sequence_path, Dict(:seqs => [[1, 2, 1, 2]], :phis => [0.5], :phasons => [0.0]))
        cfg = include(joinpath(@__DIR__, "..", "config", "default_slice_config.jl"))
        merge!(cfg, Dict{Symbol,Any}(
            :sequence_bson_path => sequence_path, :target_slope => 0.5, :N => 4,
            :output_root => directory, :run_root_pointer_file => joinpath(directory, "latest.txt"),
            :global_stop_file => joinpath(directory, "stop"), :t1_fixed => 1.0, :t2_fixed => 1.5,
            :Delta_start => 0.2, :Delta_max => 0.0, :n_Delta_slices => 2,
            :Delta_spacing_mode => :linear,
            :Delta_onset_tol => 0.026, :Delta_onset_max_extra_slices => 10,
            :rho_start => 1.0, :rho_max => 2.0, :n_rho_slices => 2,
            :mu_baseline_points => 9, :gap_refine_max_depth => 0, :transition_mu_tol => 0.02,
            :n_eval_workers => 1, :progress_log_every_points => 1000))
        original_args = copy(ARGS)
        run_roots = Dict{Symbol,String}()
        try
            for parameter in (:Delta, :rho)
                cfg[:sweep_parameter] = parameter
                cfg[:run_name] = "test_$(parameter)"
                config_path = joinpath(directory, "config_$(parameter).jl")
                open(config_path, "w") do io
                    print(io, repr(cfg))
                end
                empty!(ARGS)
                push!(ARGS, config_path)
                main()
                run_root = strip(read(cfg[:run_root_pointer_file], String))
                run_roots[parameter] = run_root
                summary = CSV.read(joinpath(run_root, "slice_summary.csv"), DataFrame)
                @test SliceBoundaryData.sweep_axis(summary) == parameter
                transitions = CSV.read(joinpath(run_root, "transition_boundaries.csv"), DataFrame)
                points = CSV.read(joinpath(run_root, "slices", "slice_0001", "point_results.csv"), DataFrame)
                onsets = CSV.read(joinpath(run_root, "onset_brackets.csv"), DataFrame)
                metadata = JLD2.load(joinpath(run_root, "run_summary.jld2"))
                @test all(==(String(parameter)), summary.sweep_parameter)
                @test nrow(transitions) > 0
                @test all(==(String(parameter)), points.sweep_parameter)
                @test metadata["sweep_parameter"] == parameter
                if parameter == :Delta
                    @test all(==(1.5), summary.rho)
                    @test summary.Delta[1:2] ≈ [0.2, 0.0]
                    @test summary.n_baseline[1:2] == [9, 17]
                    @test all(summary.Delta .== summary.sweep_value)
                    @test length(unique(summary.mu_critical)) == 1
                    @test unique(summary.mu_critical) ≈ [2.5]
                    @test all(==(1.0), points.t1) && all(==(1.5), points.t2)
                    @test nrow(summary) > 2
                    @test nrow(onsets) == 1
                    @test onsets.Delta_lo[1] < 0.1 <= onsets.Delta_hi[1]
                    @test onsets.Delta_hi[1] - onsets.Delta_lo[1] <= cfg[:Delta_onset_tol]
                    @test metadata["stop_reason"] == "completed"
                else
                    @test summary.rho == summary.sweep_value
                    @test all(==(0.05), summary.Delta)
                    @test :rho_lo in propertynames(onsets)
                end
            end
            @test_throws ErrorException SliceMergeTest.load_run(run_roots[:Delta])
            legacy = joinpath(directory, "legacy")
            priority = joinpath(directory, "priority")
            cp(run_roots[:rho], legacy)
            cp(run_roots[:rho], priority)
            legacy_summary = CSV.read(joinpath(legacy, "slice_summary.csv"), DataFrame)
            legacy_transitions = CSV.read(joinpath(legacy, "transition_boundaries.csv"), DataFrame)
            filter!(row -> row.slice_index == 1, legacy_summary)
            filter!(row -> row.slice_index == 1, legacy_transitions)
            select!(legacy_summary, Not([:Delta, :sweep_value, :sweep_parameter, :mu_min, :mu_max, :mu_critical]))
            @test SliceBoundaryData.sweep_axis(legacy_summary) == :rho
            select!(legacy_transitions, Not([:Delta, :sweep_value, :sweep_parameter]))
            CSV.write(joinpath(legacy, "slice_summary.csv"), legacy_summary)
            CSV.write(joinpath(legacy, "transition_boundaries.csv"), legacy_transitions)
            priority_summary = CSV.read(joinpath(priority, "slice_summary.csv"), DataFrame)
            priority_transitions = CSV.read(joinpath(priority, "transition_boundaries.csv"), DataFrame)
            filter!(row -> row.slice_index == 2, priority_summary)
            filter!(row -> row.slice_index == 2, priority_transitions)
            CSV.write(joinpath(priority, "slice_summary.csv"), priority_summary)
            CSV.write(joinpath(priority, "transition_boundaries.csv"), priority_transitions)
            merged = joinpath(directory, "merged")
            SliceMergeTest.merge_runs(legacy, priority, merged)
            merged_run = SliceMergeTest.load_run(merged)
            @test merged_run.summary.rho == [1.0, 2.0]
            @test ismissing(merged_run.summary.sweep_parameter[1])
            @test merged_run.summary.sweep_parameter[2] == "rho"
        finally
            empty!(ARGS)
            append!(ARGS, original_args)
        end
    end
end