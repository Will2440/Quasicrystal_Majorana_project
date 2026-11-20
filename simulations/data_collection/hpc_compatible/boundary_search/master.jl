#!/usr/bin/env julia

module BoundarySearch

using CSV, DataFrames, Dates, JLD2, LinearAlgebra, Printf
using Base.Threads

const ROOT = @__DIR__
include(joinpath(ROOT, "core", "slice_search.jl"))
include(joinpath(ROOT, "core", "tracking.jl"))
include(joinpath(ROOT, "core", "coverage.jl"))
include(joinpath(ROOT, "core", "physics.jl"))
include(joinpath(ROOT, "core", "sequences.jl"))
include(joinpath(ROOT, "core", "bandwidth.jl"))
include(joinpath(ROOT, "core", "sweep.jl"))
include(joinpath(ROOT, "core", "io_utils.jl"))

using .TransitionTracking: EdgeTrack, predict_edges, update_tracks!, predicted_exclusions
using .BoundaryCoverage: plan_grid, baseline_point_count

struct SearchStop <: Exception
    reason::String
end

function validate_config(cfg::Dict{Symbol,Any})
    values = SliceSweep.sweep_grid(cfg)
    bounds = Float64.(cfg[:u_range])
    length(bounds) == 2 && all(isfinite, bounds) && bounds[1] < bounds[2] || error("u_range must have finite increasing bounds")
    for key in (:N, :x0_count, :n_lowest_evals, :kk_maxiter, :max_points_per_slice,
                :max_total_points, :n_eval_workers, :progress_log_every_points, :blas_threads,
                :momentum_min_history, :confidence_min_history, :track_max_missed_slices)
        cfg[key] isa Integer && cfg[key] > 0 || error("$key must be a positive integer")
    end
    for key in (:transition_max_refine_iters, :audit_points_per_interval,
                :rho_onset_max_extra_slices, :Delta_onset_max_extra_slices)
        cfg[key] isa Integer && cfg[key] >= 0 || error("$key must be a nonnegative integer")
    end
    cfg[:u_baseline_points] isa Integer && cfg[:u_baseline_points] >= 3 || error("u_baseline_points must be at least 3")
    spacing = get(cfg, :mu_baseline_spacing_max, nothing)
    spacing === nothing || (spacing isa Real && isfinite(spacing) && spacing > 0) || error("mu_baseline_spacing_max must be nothing or finite and positive")
    prediction_count = get(cfg, :u_refine_points_per_prediction, 0)
    prediction_count isa Integer && (prediction_count == 0 || prediction_count >= 3) || error("u_refine_points_per_prediction must be 0 or at least 3")
    for key in (:t1_fixed, :t2_fixed, :Delta, :transition_u_tol, :rho_onset_tol, :Delta_onset_tol, :kk_tol)
        isfinite(cfg[key]) && cfg[key] > 0 || error("$key must be finite and positive")
    end
    for key in (:x0_alpha, :slope_tolerance, :phason_tolerance, :confidence_max_residual,
                :boundary_padding_min, :boundary_padding_max, :momentum_damping,
                :momentum_max_shift_fraction, :track_match_base_tolerance_fraction,
                :track_match_max_tolerance_fraction, :track_match_step_scale, :track_unmatched_penalty)
        isfinite(cfg[key]) && cfg[key] >= 0 || error("$key must be finite and nonnegative")
    end
    cfg[:boundary_padding_min] <= cfg[:boundary_padding_max] || error("boundary padding min exceeds max")
    cfg[:track_match_base_tolerance_fraction] <= cfg[:track_match_max_tolerance_fraction] || error("track match base exceeds maximum")
    isfinite(cfg[:topological_threshold]) || error("topological_threshold must be finite")
    0 < cfg[:topological_fraction] <= 1 || error("topological_fraction must be in (0, 1]")
    for key in (:target_slope, :target_phason)
        cfg[key] === nothing || isfinite(cfg[key]) || error("$key must be finite or nothing")
    end
    cfg[:sequence_selection_mode] in (:first, :closest) || error("invalid sequence_selection_mode")
    cfg[:sweep_mode] in (:vary_t1, :vary_t2) || error("invalid sweep_mode")
    cfg[:skip_trivial_interiors] isa Bool && cfg[:momentum_enabled] isa Bool || error("skip and momentum settings must be boolean")
    !cfg[:skip_trivial_interiors] || cfg[:audit_points_per_interval] > 0 || error("skipping requires at least one audit point per interval")
    delta = SliceSweep.sweep_parameter(cfg) == :Delta ? first(values) : cfg[:Delta]
    SliceSweep.kappas_for_delta(cfg, delta)
    return cfg
end

function make_context(cfg)
    sequences, slopes, phasons = DynamicSequences.load_sequences_from_bson(String(cfg[:sequence_bson_path]))
    _, sequence, gamma, phason = DynamicSequences.select_sequence(sequences, slopes, phasons;
        target_slope=isnothing(cfg[:target_slope]) ? nothing : Float64(cfg[:target_slope]), slope_tolerance=Float64(cfg[:slope_tolerance]),
        target_phason=isnothing(cfg[:target_phason]) ? nothing : Float64(cfg[:target_phason]), phason_tolerance=Float64(cfg[:phason_tolerance]),
        selection_mode=cfg[:sequence_selection_mode])
    isfinite(gamma) && 0 <= gamma <= 1 || error("selected sequence slope must be in [0, 1]")
    isfinite(phason) || error("selected phason must be finite")
    length(sequence) >= cfg[:N] || error("selected sequence is shorter than N")
    all(index -> index in (1, 2), sequence[1:cfg[:N]]) || error("sequence entries must be 1 or 2")
    x0s = cfg[:x0_count] == 1 ? [Float64(cfg[:N] / 2)] :
        collect(range(cfg[:N] / 2 - cfg[:x0_alpha] * cfg[:N],
                      cfg[:N] / 2 + cfg[:x0_alpha] * cfg[:N]; length=cfg[:x0_count]))
    return (sequence=sequence, gamma=gamma, phason=phason, x0s=x0s,
            fixed_mu_critical=SliceSweep.fixed_mu_critical(cfg, gamma))
end

function evaluate_point(value::Float64, u::Float64, cfg, ctx)
    params = SliceSweep.physical_parameters(value, cfg, ctx.gamma, ctx.fixed_mu_critical)
    isfinite(params.mu_critical) && params.mu_critical > 0 || error("mu_critical must be finite and positive")
    mu = u * params.mu_critical
    kappas = SliceSweep.kappas_for_delta(cfg, params.Delta)
    result = DynamicPhysics.evaluate_specloc_point(Int(cfg[:N]), [params.t1, params.t2], mu,
        params.Delta, ctx.sequence, ctx.x0s, kappas; n_lowest_evals=Int(cfg[:n_lowest_evals]),
        kk_tol=cfg[:kk_tol], kk_maxiter=Int(cfg[:kk_maxiter]),
        topological_threshold=Float64(cfg[:topological_threshold]), topological_fraction=Float64(cfg[:topological_fraction]))
    result.phase_label in (0, 1) || error("classifier must return phase 0 or 1")
    return merge((sweep_parameter=String(SliceSweep.sweep_parameter(cfg)), sweep_value=value,
        rho=params.rho, Delta=params.Delta, t1=params.t1, t2=params.t2, mu_critical=params.mu_critical,
        u=u, mu=mu, kappa_mode=String(cfg[:kappa_mode]), kappa_values=join(kappas, ";")), result)
end

function write_rows(path, rows, columns)
    frame = isempty(rows) ? DataFrame([name => Any[] for name in columns]) : DataFrame(rows)
    CSV.write(path, frame)
end

function check_stop(stop_files)
    any(isfile, stop_files) && throw(SearchStop("stop_requested"))
end

function create_run_root(output_root, run_name)
    timestamp = Dates.format(now(), "yy-mm-dd-HH-MM-SS")
    raw_job_id = get(ENV, "SLURM_JOB_ID", "local")
    job_id = occursin(r"^\d+$", raw_job_id) ? raw_job_id : "local"
    base = joinpath(output_root, "$(run_name)_$(timestamp)_$(job_id)")
    candidate = base
    suffix = 1
    while true
        try
            mkdir(candidate)
            return candidate
        catch error
            ispath(candidate) || rethrow(error)
            candidate = "$(base)_$(suffix)"
            suffix += 1
        end
    end
end

function scan_slice(value, execution_id, cfg, ctx, predictions, exclusions, total, stop_files, run_root; kind="frontier")
    directory = joinpath(run_root, "slices", @sprintf("slice_%04d", execution_id))
    mkpath(directory)
    params = SliceSweep.physical_parameters(value, cfg, ctx.gamma, ctx.fixed_mu_critical)
    isfinite(params.mu_critical) && params.mu_critical > 0 || error("mu_critical must be finite and positive")
    baseline_count = baseline_point_count(cfg[:u_range], cfg[:u_baseline_points], params.mu_critical,
        get(cfg, :mu_baseline_spacing_max, nothing))
    plan = plan_grid(cfg[:u_range], baseline_count, exclusions, cfg[:audit_points_per_interval])
    samples = Dict{Float64,NamedTuple}()
    observations = NamedTuple[]
    fallback = ""
    active_exclusions = plan.exclusions
    started = time()
    prediction_count = get(cfg, :u_refine_points_per_prediction, 0)
    predicted_points = Float64[]
    for prediction in predictions
        prediction_count == 0 && continue
        lower = max(cfg[:u_range][1], prediction.u - prediction.band_half_width)
        upper = min(cfg[:u_range][2], prediction.u + prediction.band_half_width)
        append!(predicted_points, range(lower, upper; length=prediction_count))
    end
    filter!(point -> !any(interval -> interval[1] < point < interval[2], active_exclusions), predicted_points)
    prediction_rows = [(execution_id=execution_id, sweep_value=value, track_id=p.track_id, u=p.u,
        mu=p.u * params.mu_critical, phase_lo=p.orientation[1], phase_hi=p.orientation[2],
        band_half_width=p.band_half_width, confidence=p.confidence,
        source_last_seen=p.source_last_seen, source_history_length=p.source_history_length) for p in predictions]
    exclusion_rows = [(execution_id=execution_id, sweep_value=value, u_lo=interval[1], u_hi=interval[2],
        mu_lo=interval[1] * params.mu_critical, mu_hi=interval[2] * params.mu_critical) for interval in plan.exclusions]

    function evaluate_many(points, stage; force=false)
        remaining = sort!(unique!(Float64[point for point in points if force || !haskey(samples, point)]))
        workers = min(Int(cfg[:n_eval_workers]), Threads.nthreads())
        for offset in 1:workers:length(remaining)
            check_stop(stop_files)
            batch = remaining[offset:min(offset + workers - 1, end)]
            capacity = min(Int(cfg[:max_points_per_slice]) - length(observations), Int(cfg[:max_total_points]) - total[])
            capacity > 0 || throw(SearchStop("point_budget"))
            batch = batch[1:min(length(batch), capacity)]
            jobs = [Threads.@spawn try
                (result=evaluate_point(value, point, cfg, ctx), error=nothing)
            catch err
                (result=nothing, error=err)
            end for point in batch]
            failure = nothing
            for job in jobs
                outcome = fetch(job)
                if outcome.error === nothing
                    point_stage = stage == "planned" && any(interval -> interval[1] < outcome.result.u < interval[2], active_exclusions) ? "audit" : stage
                    row = merge(outcome.result, (execution_id=execution_id, stage=point_stage))
                    push!(observations, row)
                    total[] += 1
                    if haskey(samples, row.u) && samples[row.u].phase_label != row.phase_label
                        failure = SearchStop("inconsistent_measurement")
                    else
                        samples[row.u] = row
                    end
                    if length(observations) % Int(cfg[:progress_log_every_points]) == 0
                        DynamicIO.logmsg("$(SliceSweep.sweep_parameter(cfg))=$value evaluated $(length(observations)) points")
                    end
                else
                    failure = outcome.error
                end
            end
            failure === nothing || throw(failure)
            length(batch) < min(workers, length(remaining) - offset + 1) && throw(SearchStop("point_budget"))
        end
    end

    function coarse_transitions()
        ordered = sort!(collect(values(samples)); by=point -> point.u)
        return [(left, right) for (left, right) in zip(ordered[1:end-1], ordered[2:end]) if
            left.phase_label != right.phase_label && !any(interval -> left.u < interval[2] && right.u > interval[1], active_exclusions)]
    end

    try
        audit_points = [point for point in plan.points if any(interval -> interval[1] <= point <= interval[2], active_exclusions)]
        for point in audit_points
            evaluate_many([point], "planned")
            if samples[point].phase_label == 1
                fallback = "audit_phase_1"
                break
            end
        end
        if !isempty(fallback)
            active_exclusions = Tuple{Float64,Float64}[]
            evaluate_many(plan_grid(cfg[:u_range], baseline_count, [], 0).points, "fallback_baseline"; force=true)
            evaluate_many(predicted_points, "planned")
        else
            evaluate_many(vcat(plan.points, predicted_points), "planned")
        end
        if !isempty(active_exclusions)
            coarse = coarse_transitions()
            for (lower, upper) in active_exclusions
                left_predictions = filter(p -> p.confidence && p.orientation == (1, 0) && p.u <= lower, predictions)
                right_predictions = filter(p -> p.confidence && p.orientation == (0, 1) && p.u >= upper, predictions)
                left_prediction = isempty(left_predictions) ? nothing : left_predictions[argmax(p.u for p in left_predictions)]
                right_prediction = isempty(right_predictions) ? nothing : right_predictions[argmin(p.u for p in right_predictions)]
                left_ok = any(pair -> pair[1].phase_label == 1 && pair[2].phase_label == 0 &&
                    left_prediction !== nothing && pair[1].u >= left_prediction.u - left_prediction.band_half_width &&
                    pair[2].u <= left_prediction.u + left_prediction.band_half_width, coarse)
                right_ok = any(pair -> pair[1].phase_label == 0 && pair[2].phase_label == 1 &&
                    right_prediction !== nothing && pair[1].u >= right_prediction.u - right_prediction.band_half_width &&
                    pair[2].u <= right_prediction.u + right_prediction.band_half_width, coarse)
                if !(left_ok && right_ok)
                    fallback = "boundary_confirmation_failed"
                    break
                end
            end
            if !isempty(fallback)
                active_exclusions = Tuple{Float64,Float64}[]
                evaluate_many(plan_grid(cfg[:u_range], baseline_count, [], 0).points, "fallback_baseline"; force=true)
            end
        end
        transitions = NamedTuple[]
        for (left, right) in coarse_transitions()
            iterations = 0
            while right.u - left.u > cfg[:transition_u_tol] && iterations < cfg[:transition_max_refine_iters]
                midpoint = left.u + (right.u - left.u) / 2
                midpoint == left.u || midpoint == right.u ? break : nothing
                evaluate_many([midpoint], "transition_refinement")
                middle = samples[midpoint]
                if middle.phase_label == left.phase_label
                    left = middle
                else
                    right = middle
                end
                iterations += 1
            end
            center = (left.u + right.u) / 2
            push!(transitions, (u=center, u_lo=left.u, u_hi=right.u, mu=center * params.mu_critical,
                mu_lo=left.mu, mu_hi=right.mu, phase_lo=left.phase_label, phase_hi=right.phase_label,
                gap_min=min(left.gap_min, right.gap_min), refinement_iters=iterations,
                converged=right.u-left.u <= cfg[:transition_u_tol]))
        end
        ordered = sort!(collect(values(samples)); by=point -> point.u)
        gap_count = count(index -> ordered[index].phase_label == 0 &&
            (index == 1 || ordered[index-1].phase_label != 0), eachindex(ordered))
        return (execution_id=execution_id, kind=kind, sweep_value=value, rho=params.rho, Delta=params.Delta,
            sweep_parameter=String(SliceSweep.sweep_parameter(cfg)), mu_critical=params.mu_critical,
            mu_min=cfg[:u_range][1] * params.mu_critical, mu_max=cfg[:u_range][2] * params.mu_critical,
            kappa_mode=String(cfg[:kappa_mode]), kappa_values=join(SliceSweep.kappas_for_delta(cfg, params.Delta), ";"),
            gap_count=gap_count, transitions=transitions, n_points=length(observations),
            n_baseline=length(plan.points), n_predictions=length(predictions),
            n_exclusions=length(active_exclusions), fallback=fallback, elapsed_seconds=time()-started)
    finally
        write_rows(joinpath(directory, "point_results.csv"), observations, [:execution_id, :u, :mu, :phase_label, :stage])
        write_rows(joinpath(directory, "predictions.csv"), prediction_rows, [:execution_id, :sweep_value, :track_id, :u, :mu])
        write_rows(joinpath(directory, "excluded_intervals.csv"), [merge(row, (fallback=fallback,)) for row in exclusion_rows],
            [:execution_id, :sweep_value, :u_lo, :u_hi, :mu_lo, :mu_hi, :fallback])
    end
end

function save_run(run_root, cfg, ctx, results, boundary_rows, tracks, onsets, pending, total, extra, status, config_path)
    summaries = [merge((sequence_index=index,), NamedTuple{Tuple(filter(!=(:transitions), collect(keys(row))))}(
        Tuple(getproperty(row, key) for key in keys(row) if key != :transitions)), (n_transitions=length(row.transitions),)) for (index, row) in enumerate(results)]
    write_rows(joinpath(run_root, "slice_summary.csv"), summaries, [:sequence_index, :execution_id, :sweep_value, :gap_count])
    write_rows(joinpath(run_root, "transition_boundaries.csv"), boundary_rows, [:sequence_index, :execution_id, :track_id, :u, :mu])
    write_rows(joinpath(run_root, "onset_brackets.csv"), vcat(onsets, pending),
        [:sweep_parameter, :value_lo, :value_hi, :gap_count_lo, :gap_count_hi, :status])
    metadata = (config=cfg, config_path=config_path, run_root=run_root, stop_reason=status,
        timestamp=string(now()), total_points=total, extra_slices=extra, sweep_parameter=SliceSweep.sweep_parameter(cfg),
        results=results, boundaries=boundary_rows, tracks=tracks, onset_brackets=onsets, pending_brackets=pending,
        gamma=ctx.gamma, phason=ctx.phason, julia_threads=Threads.nthreads(), blas_threads=BLAS.get_num_threads(),
        provenance="Point rows are actual evaluations; boundaries are brackets of measured samples. Exclusions are provisional, and equal gap counts do not exclude onsets.")
    JLD2.jldsave(joinpath(run_root, "run_checkpoint.jld2"); metadata...)
    JLD2.jldsave(joinpath(run_root, "run_summary.jld2"); metadata...)
end

function bracket_row(left, right, parameter, status)
    lower, upper = left.sweep_value <= right.sweep_value ? (left, right) : (right, left)
    return (sweep_parameter=String(parameter), value_lo=lower.sweep_value, value_hi=upper.sweep_value,
        gap_count_lo=lower.gap_count, gap_count_hi=upper.gap_count, status=status,
        execution_id_lo=lower.execution_id, execution_id_hi=upper.execution_id)
end

function run_search(cfg::Dict{Symbol,Any}; config_path="")
    validate_config(cfg)
    BLAS.set_num_threads(Int(cfg[:blas_threads]))
    ctx = make_context(cfg)
    mkpath(String(cfg[:output_root]))
    run_root = create_run_root(String(cfg[:output_root]), String(cfg[:run_name]))
    pointer = String(cfg[:run_root_pointer_file])
    mkpath(dirname(pointer))
    open(pointer, "w") do io
        println(io, run_root)
    end
    stop_files = [joinpath(run_root, "STOP_REQUESTED"), String(cfg[:global_stop_file])]
    parameter = SliceSweep.sweep_parameter(cfg)
    tolerance = Float64(cfg[parameter == :rho ? :rho_onset_tol : :Delta_onset_tol])
    budget = Int(cfg[parameter == :rho ? :rho_onset_max_extra_slices : :Delta_onset_max_extra_slices])
    tracks = EdgeTrack[]
    results, boundaries, onsets, pending = NamedTuple[], NamedTuple[], NamedTuple[], NamedTuple[]
    total, extra, execution_id = Ref(0), 0, 0
    status = "completed"
    previous = nothing
    latest_ids = Int[]
    buffered = NamedTuple[]
    queue = Tuple{NamedTuple,NamedTuple}[]

    function replay!(row)
        sequence_index = length(results) + 1
        predictions = predict_edges(tracks, row.sweep_value, sequence_index, cfg[:u_range], cfg)
        ids = update_tracks!(tracks, row.transitions, row.sweep_value, sequence_index, cfg[:u_range], cfg; predictions=predictions)
        push!(results, row)
        for (index, transition) in enumerate(row.transitions)
            push!(boundaries, merge((sequence_index=sequence_index, execution_id=row.execution_id,
                sweep_parameter=row.sweep_parameter, sweep_value=row.sweep_value, rho=row.rho,
                Delta=row.Delta, track_id=ids[index], transition_index=index), transition))
        end
        return ids
    end

    try
        for value in SliceSweep.sweep_grid(cfg)
            check_stop(stop_files)
            execution_id += 1
            predictions = predict_edges(tracks, value, length(results) + 1, cfg[:u_range], cfg)
            exclusions = previous !== nothing && cfg[:skip_trivial_interiors] ?
                predicted_exclusions(predictions, previous.transitions, latest_ids, cfg, cfg[:u_range]) : []
            endpoint = scan_slice(value, execution_id, cfg, ctx, predictions, exclusions, total, stop_files, run_root)
            empty!(buffered)
            push!(buffered, endpoint)
            empty!(queue)
            if previous !== nothing && previous.gap_count != endpoint.gap_count
                push!(queue, (previous, endpoint))
            end
            while !isempty(queue)
                check_stop(stop_files)
                left, right = first(queue)
                if abs(right.sweep_value - left.sweep_value) <= tolerance
                    push!(onsets, bracket_row(left, right, parameter, "resolved"))
                    popfirst!(queue)
                    continue
                end
                extra < budget || throw(SearchStop("onset_refinement_limit"))
                midpoint = left.sweep_value + (right.sweep_value - left.sweep_value) / 2
                midpoint != left.sweep_value && midpoint != right.sweep_value || throw(SearchStop("onset_refinement_limit"))
                execution_id += 1
                probe_predictions = predict_edges(tracks, midpoint, length(results) + 1, cfg[:u_range], cfg)
                probe_exclusions = cfg[:skip_trivial_interiors] ?
                    predicted_exclusions(probe_predictions, previous.transitions, latest_ids, cfg, cfg[:u_range]) : []
                probe = scan_slice(midpoint, execution_id, cfg, ctx, probe_predictions, probe_exclusions, total, stop_files, run_root; kind="onset_probe")
                extra += 1
                push!(buffered, probe)
                popfirst!(queue)
                left.gap_count != probe.gap_count && push!(queue, (left, probe))
                probe.gap_count != right.gap_count && push!(queue, (probe, right))
            end
            sort!(buffered; by=row -> SliceSweep.sweep_progress(row.sweep_value, cfg))
            for row in buffered
                latest_ids = replay!(row)
            end
            previous = endpoint
            empty!(buffered)
            save_run(run_root, cfg, ctx, results, boundaries, tracks, onsets, pending, total[], extra, "in_progress", config_path)
            DynamicIO.logmsg("$parameter=$value gaps=$(endpoint.gap_count), points=$(total[]), probes=$extra")
        end
    catch err
        if err isa SearchStop
            status = err.reason
        elseif err isa InterruptException
            status = "interrupt"
        else
            status = "error: " * sprint(showerror, err)
        end
        append!(pending, [bracket_row(left, right, parameter, "pending") for (left, right) in queue])
        sort!(buffered; by=row -> SliceSweep.sweep_progress(row.sweep_value, cfg))
        for row in buffered
            replay!(row)
        end
        save_run(run_root, cfg, ctx, results, boundaries, tracks, onsets, pending, total[], extra, status, config_path)
        err isa SearchStop || err isa InterruptException || rethrow()
    end
    save_run(run_root, cfg, ctx, results, boundaries, tracks, onsets, pending, total[], extra, status, config_path)
    return (run_root=run_root, stop_reason=status, results=results, boundaries=boundaries,
        tracks=tracks, onsets=onsets, pending=pending, total_points=total[], extra_slices=extra)
end

function main(args=ARGS)
    Base.exit_on_sigint(false)
    length(args) <= 1 || error("Usage: julia --project=boundary_search boundary_search/master.jl [config.jl]")
    config_path = isempty(args) ? joinpath(ROOT, "config", "default_config.jl") : abspath(only(args))
    cfg = DynamicIO.read_config(config_path)
    mkpath(String(cfg[:output_root]))
    result = run_search(cfg; config_path=config_path)
    DynamicIO.logmsg("Finished: $(result.stop_reason), output=$(result.run_root)")
    return result
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    BoundarySearch.main()
end