using Test, BSON, CSV, DataFrames, JLD2

include(joinpath(@__DIR__, "..", "master.jl"))
using .BoundarySearch

function fixture(directory)
    path = joinpath(directory, "sequence.bson")
    BSON.bson(path, Dict(:seqs => [[1, 2, 1, 2, 1, 2]], :phis => [0.5], :phasons => [0.0]))
    cfg = include(joinpath(@__DIR__, "..", "config", "default_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(
        :sequence_bson_path => path, :target_slope => 0.5, :target_phason => 0.0,
        :N => 6, :output_root => directory, :run_root_pointer_file => joinpath(directory, "latest.txt"),
        :global_stop_file => joinpath(directory, "stop"), :run_name => "integration",
        :u_range => (0.0, 1.0), :u_baseline_points => 41, :transition_u_tol => 0.002,
        :sweep_parameter => :rho, :skip_trivial_interiors => false,
        :rho_start => 1.1, :rho_max => 1.5, :n_rho_slices => 3, :rho_spacing_mode => :linear,
        :Delta_start => 0.5, :Delta_max => 0.1, :n_Delta_slices => 3,
        :Delta_spacing_mode => :linear, :rho_onset_tol => 0.025, :Delta_onset_tol => 0.025,
        :rho_onset_max_extra_slices => 40, :Delta_onset_max_extra_slices => 40,
        :n_eval_workers => 2, :blas_threads => 1, :n_lowest_evals => 2,
        :confidence_min_history => 1, :confidence_max_residual => 0.1,
        :momentum_min_history => 2, :momentum_damping => 1.0, :momentum_max_shift_fraction => 0.2,
        :track_match_base_tolerance_fraction => 0.1, :track_match_max_tolerance_fraction => 0.3,
        :u_refine_points_per_prediction => 5))
    return cfg
end

@testset "real localizer and CLI" begin
    mktempdir() do directory
        cfg = fixture(directory)
        cfg[:n_rho_slices] = 1
        cfg[:u_baseline_points] = 3
        cfg[:transition_max_refine_iters] = 2
        cfg[:kappa_mode] = :delta_over_L
        config_path = joinpath(directory, "config.jl")
        open(config_path, "w") do io
            print(io, repr(cfg))
        end
        run = BoundarySearch.main([config_path])
        @test run.stop_reason == "completed"
        @test occursin(r"^integration_\d{2}-\d{2}-\d{2}-\d{2}-\d{2}-\d{2}_(\d+|local)(_\d+)?$", basename(run.run_root))
        points = CSV.read(joinpath(run.run_root, "slices", "slice_0001", "point_results.csv"), DataFrame)
        @test nrow(points) >= 3
        @test all(isfinite, points.gap_min)
        @test all(points.mu .≈ points.u .* points.mu_critical)
        @test all(points.kappa_values .≈ cfg[:Delta] / cfg[:N])
        @test JLD2.load(joinpath(run.run_root, "run_summary.jld2"), "stop_reason") == "completed"
        project = normpath(joinpath(@__DIR__, ".."))
        Base.run(`$(Base.julia_cmd()) --startup-file=no --project=$project $(joinpath(project, "master.jl")) $config_path`)
        persisted = strip(read(cfg[:run_root_pointer_file], String))
        @test isfile(joinpath(persisted, "run_summary.jld2"))
    end
end

const classifier = Ref{Function}((rho, delta, u) -> 1)
const received_kappas = Ref{Vector{Float64}}(Float64[])

@eval BoundarySearch.DynamicPhysics function evaluate_specloc_point(
    N::Int, hoppings::Vector{Float64}, mu::Float64, delta::Float64,
    sequence::Vector{Int}, x0s::Vector{Float64}, kappas::Vector{Float64}; kwargs...)
    Main.received_kappas[] = kappas
    critical = hoppings[1] + hoppings[2]
    phase = Main.classifier[](hoppings[2] / hoppings[1], delta, mu / critical)
    return (phase_label=phase, topo_fraction=Float64(phase), invariant_mean=Float64(phase),
        invariant_min=Float64(phase), invariant_max=Float64(phase), gap_min=0.1,
        gap_mean=0.1, gap_max=0.1, n_samples=length(x0s) * length(kappas))
end

@testset "Physical baseline spacing in slices" begin
    mktempdir() do directory
        cfg = fixture(directory)
        cfg[:n_rho_slices] = 2
        cfg[:rho_start], cfg[:rho_max] = 1.1, 1.5
        cfg[:u_baseline_points] = 3
        cfg[:mu_baseline_spacing_max] = 0.1
        cfg[:u_refine_points_per_prediction] = 0
        classifier[] = (rho, delta, u) -> 1
        run = BoundarySearch.run_search(cfg)
        @test run.stop_reason == "completed"
        @test length(run.results) == 2
        for row in run.results
            @test row.n_baseline == max(3, ceil(Int, row.mu_critical / 0.1) + 1)
            @test row.mu_critical / (row.n_baseline - 1) <= 0.1
        end
        @test last(run.results).n_baseline > first(run.results).n_baseline
        cfg[:mu_baseline_spacing_max] = 0.0
        @test_throws ErrorException BoundarySearch.validate_config(cfg)
    end
end

@testset "ordered refinement for both sweep directions" begin
    for mode in (:rho_up, :rho_down, :delta_down)
        mktempdir() do directory
            cfg = fixture(directory)
            cfg[:n_rho_slices] = 4
            cfg[:n_Delta_slices] = 4
            if mode == :rho_down
                cfg[:rho_start], cfg[:rho_max] = 0.9, 0.5
            elseif mode == :delta_down
                cfg[:sweep_parameter] = :Delta
                cfg[:kappa_mode] = :delta_over_L
            end
            progress(rho, delta) = mode == :rho_up ? (rho - 1.1) / 0.4 :
                mode == :rho_down ? (0.9 - rho) / 0.4 : (0.5 - delta) / 0.4
            classifier[] = (rho, delta, u) -> begin
                fraction = progress(rho, delta)
                left = 0.18 + 0.05 * fraction
                right = 0.35 + 0.05 * fraction
                left <= u <= right || (fraction >= 0.56 && 0.62 <= u <= 0.77) ? 0 : 1
            end
            run = BoundarySearch.run_search(cfg)
            @test run.stop_reason == "completed"
            @test !isempty(run.onsets)
            @test all(row.status == "resolved" && row.value_hi-row.value_lo <= 0.025 for row in run.onsets)
            @test isempty(run.pending)
            @test run.extra_slices > 0
            progress_values = [BoundarySearch.SliceSweep.sweep_progress(row.sweep_value, cfg) for row in run.results]
            @test issorted(progress_values)
            @test any(row.kind == "onset_probe" for row in run.results)
            @test any(row.execution_id != index for (index, row) in enumerate(run.results))
            @test all(row.track_id > 0 for row in run.boundaries)
            @test all(issorted([BoundarySearch.SliceSweep.sweep_progress(point.sweep_value, cfg) for point in track.history]) for track in run.tracks)
            @test any(length(track.history) > cfg[:n_rho_slices] for track in run.tracks)
            final_directory = joinpath(run.run_root, "slices", "slice_$(lpad(last(run.results).execution_id, 4, '0'))")
            final_predictions = CSV.read(joinpath(final_directory, "predictions.csv"), DataFrame)
            @test any(final_predictions.source_history_length .> 3)
            first_ids = Set(row.track_id for row in run.boundaries if row.execution_id == first(run.results).execution_id)
            final_ids = Set(row.track_id for row in run.boundaries if row.execution_id == last(run.results).execution_id)
            @test first_ids ⊆ final_ids
            @test all(row.mu ≈ row.u * only(filter(result -> result.execution_id == row.execution_id, run.results)).mu_critical for row in run.boundaries)
            table = CSV.read(joinpath(run.run_root, "transition_boundaries.csv"), DataFrame)
            @test all(table.track_id .> 0)
            @test nrow(table) == length(run.boundaries)
            @test sum(row.n_points for row in run.results) == run.total_points
            @test first(run.results).sweep_value ≈ (mode == :delta_down ? 0.5 : cfg[:rho_start])
            @test last(run.results).sweep_value ≈ (mode == :delta_down ? 0.1 : cfg[:rho_max])
            if mode == :delta_down
                BoundarySearch.evaluate_point(last(run.results).sweep_value, 0.5, cfg, BoundarySearch.make_context(cfg))
                @test received_kappas[] ≈ [last(run.results).Delta / cfg[:N]]
            end
        end
    end
end

@testset "exclusions and fallback" begin
    for scenario in (:stable, :audit, :confirmation)
        mktempdir() do directory
            cfg = fixture(directory)
            cfg[:n_rho_slices] = 5
            cfg[:skip_trivial_interiors] = true
            cfg[:audit_points_per_interval] = 3
            cfg[:rho_onset_tol] = 1.0
            classifier[] = (rho, delta, u) -> begin
                changed = rho > 1.49
                left = scenario == :confirmation && changed ? 0.08 : 0.2
                right = scenario == :confirmation && changed ? 0.92 : 0.8
                phase = left <= u <= right ? 0 : 1
                scenario == :audit && changed && 0.46 <= u <= 0.54 ? 1 : phase
            end
            run = BoundarySearch.run_search(cfg)
            @test run.stop_reason == "completed"
            @test any(row.n_exclusions > 0 for row in run.results)
            last_result = last(run.results)
            @test last_result.fallback == (scenario == :stable ? "" : scenario == :audit ? "audit_phase_1" : "boundary_confirmation_failed")
            points = CSV.read(joinpath(run.run_root, "slices", "slice_0005", "point_results.csv"), DataFrame)
            exclusions = CSV.read(joinpath(run.run_root, "slices", "slice_0005", "excluded_intervals.csv"), DataFrame)
            @test nrow(exclusions) > 0
            @test all(points.mu .≈ points.u .* points.mu_critical)
            @test nrow(points) == last_result.n_points
            if scenario == :audit
                failed_audit = findfirst(index -> points.stage[index] == "audit" && points.phase_label[index] == 1, 1:nrow(points))
                @test failed_audit !== nothing
                @test points.stage[failed_audit + 1] == "fallback_baseline"
                @test count(==("audit"), points.stage) < cfg[:audit_points_per_interval]
            else
                @test count(==("audit"), points.stage) == cfg[:audit_points_per_interval]
            end
            for boundary in filter(row -> row.execution_id == last_result.execution_id, run.boundaries)
                @test any(row.u == boundary.u_lo && row.phase_label == boundary.phase_lo for row in eachrow(points))
                @test any(row.u == boundary.u_hi && row.phase_label == boundary.phase_hi for row in eachrow(points))
            end
            if scenario != :stable
                @test last_result.n_exclusions == 0
                @test count(==("fallback_baseline"), points.stage) == cfg[:u_baseline_points]
                @test all(any(isapprox(point, observed; atol=1e-12) for observed in points.u) for point in range(0.0, 1.0; length=41))
            else
                @test last_result.gap_count == 1
                @test length(last_result.transitions) == 2
            end
        end
    end
end

@testset "onset probe exclusions" begin
    for mode in (:rho_up, :delta_down), skipping in (false, true)
        mktempdir() do directory
            cfg = fixture(directory)
            cfg[:sweep_parameter] = mode == :rho_up ? :rho : :Delta
            cfg[:n_rho_slices] = 5
            cfg[:n_Delta_slices] = 5
            cfg[:skip_trivial_interiors] = skipping
            cfg[:audit_points_per_interval] = 3
            progress(rho, delta) = mode == :rho_up ? (rho - 1.1) / 0.4 : (0.5 - delta) / 0.4
            classifier[] = (rho, delta, u) ->
                (0.2 <= u <= 0.5 || (progress(rho, delta) >= 0.875 && 0.7 <= u <= 0.8)) ? 0 : 1
            run = BoundarySearch.run_search(cfg)
            probes = filter(row -> row.kind == "onset_probe", run.results)
            @test run.stop_reason == "completed"
            @test !isempty(probes)
            @test all(row -> skipping ? row.n_exclusions > 0 : row.n_exclusions == 0, probes)
            @test all(row -> isempty(row.fallback), probes)
            @test length(run.onsets) == 1
            for probe in probes
                directory_path = joinpath(run.run_root, "slices", "slice_$(lpad(probe.execution_id, 4, '0'))")
                points = CSV.read(joinpath(directory_path, "point_results.csv"), DataFrame)
                exclusions = CSV.read(joinpath(directory_path, "excluded_intervals.csv"), DataFrame)
                @test nrow(exclusions) == (skipping ? 1 : 0)
                @test count(==("audit"), points.stage) == (skipping ? cfg[:audit_points_per_interval] : 0)
            end
        end
    end
    mktempdir() do directory
        cfg = fixture(directory)
        cfg[:n_rho_slices] = 5
        cfg[:skip_trivial_interiors] = true
        cfg[:audit_points_per_interval] = 3
        classifier[] = (rho, delta, u) -> begin
            phase = (0.2 <= u <= 0.5 || (rho >= 1.45 && 0.7 <= u <= 0.8)) ? 0 : 1
            1.44 <= rho <= 1.46 && 0.34 <= u <= 0.36 ? 1 : phase
        end
        run = BoundarySearch.run_search(cfg)
        probes = filter(row -> row.kind == "onset_probe", run.results)
        failed = filter(row -> row.fallback == "audit_phase_1", probes)
        @test !isempty(failed)
        for probe in failed
            points = CSV.read(joinpath(run.run_root, "slices", "slice_$(lpad(probe.execution_id, 4, '0'))", "point_results.csv"), DataFrame)
            @test probe.n_exclusions == 0
            @test any(==("fallback_baseline"), points.stage)
        end
    end
end

@testset "budgets, stop and validation" begin
    mktempdir() do directory
        cfg = fixture(directory)
        classifier[] = (rho, delta, u) -> rho >= 1.3 && 0.3 <= u <= 0.7 ? 0 : 1
        cfg[:rho_onset_max_extra_slices] = 1
        run = BoundarySearch.run_search(cfg)
        @test run.stop_reason == "onset_refinement_limit"
        @test run.extra_slices == 1
        @test length(run.pending) == 1
        @test maximum(row.sweep_value for row in run.results) ≈ 1.3
        @test JLD2.load(joinpath(run.run_root, "run_checkpoint.jld2"), "stop_reason") == "onset_refinement_limit"
        @test nrow(CSV.read(joinpath(run.run_root, "onset_brackets.csv"), DataFrame)) >= 1
        cfg[:rho_onset_max_extra_slices] = 40
        cfg[:max_total_points] = 45
        run = BoundarySearch.run_search(cfg)
        @test run.stop_reason == "point_budget"
        @test run.total_points == 45
        @test length(run.results) == 1
        @test isfile(joinpath(run.run_root, "slices", "slice_0002", "point_results.csv"))
        cfg[:max_total_points] = 100000
        cfg[:n_eval_workers] = 1
        classifier[] = (rho, delta, u) -> begin
            rho > 1.1 && write(cfg[:global_stop_file], "stop")
            1
        end
        run = BoundarySearch.run_search(cfg)
        @test run.stop_reason == "stop_requested"
        @test length(run.results) == 1
        @test run.total_points == cfg[:u_baseline_points] + 1
        @test nrow(CSV.read(joinpath(run.run_root, "slices", "slice_0002", "point_results.csv"), DataFrame)) == 1
        rm(cfg[:global_stop_file])
        write(cfg[:global_stop_file], "stop")
        run = BoundarySearch.run_search(cfg)
        @test run.stop_reason == "stop_requested"
        @test isempty(run.results)
        rm(cfg[:global_stop_file])
        cfg[:kappas] = [-1.0]
        @test_throws ErrorException BoundarySearch.validate_config(cfg)
        cfg[:kappas] = [0.01]
        cfg[:transition_u_tol] = NaN
        @test_throws ErrorException BoundarySearch.validate_config(cfg)
        cfg[:transition_u_tol] = 0.001
        cfg[:skip_trivial_interiors] = true
        cfg[:audit_points_per_interval] = 0
        @test_throws ErrorException BoundarySearch.validate_config(cfg)
    end
end