#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Statistics
using Base.Threads

project_root = @__DIR__
include(joinpath(project_root, "core", "slice_search.jl"))
include(joinpath(project_root, "core", "tracking.jl"))
include(joinpath(project_root, "core", "physics.jl"))
include(joinpath(project_root, "core", "sequences.jl"))
include(joinpath(project_root, "core", "bandwidth.jl"))
include(joinpath(project_root, "core", "sweep.jl"))
include(joinpath(project_root, "core", "io_utils.jl"))

using .SliceSearch
using .TransitionTracking
using .DynamicPhysics
using .DynamicSequences
using .DynamicBandwidth
using .SliceSweep
using .DynamicIO

struct PointBudgetExceeded <: Exception
    message::String
end
Base.showerror(io::IO, err::PointBudgetExceeded) = print(io, err.message)

function stop_requested(stop_files::Vector{String})
    return any(isfile, stop_files)
end

function build_x0s(N::Int, alpha::Float64, npoints::Int)
    return collect(range(N / 2 - alpha * N, N / 2 + alpha * N; length=npoints))
end

function make_eval_context(cfg::Dict{Symbol,Any})
    seqs, phis, phasons = DynamicSequences.load_sequences_from_bson(cfg[:sequence_bson_path])
    _, sequence, gamma, phason = DynamicSequences.select_sequence(
        seqs, phis, phasons;
        target_slope=cfg[:target_slope],
        slope_tolerance=cfg[:slope_tolerance],
        target_phason=cfg[:target_phason],
        phason_tolerance=cfg[:phason_tolerance],
        selection_mode=cfg[:sequence_selection_mode],
    )

    N = Int(cfg[:N])
    return (
        N=N,
        sequence=sequence,
        gamma=gamma,
        phason=phason,
        fixed_mu_critical=SliceSweep.fixed_mu_critical(cfg, gamma),
        x0s=build_x0s(N, Float64(cfg[:x0_alpha]), Int(cfg[:x0_count])),
        kappas=Vector{Float64}(cfg[:kappas]),
        n_lowest_evals=Int(cfg[:n_lowest_evals]),
        kk_tol=Float64(cfg[:kk_tol]),
        kk_maxiter=Int(cfg[:kk_maxiter]),
        topological_threshold=Float64(cfg[:topological_threshold]),
        topological_fraction=Float64(cfg[:topological_fraction]),
    )
end

function point_key(rho::Float64, mu::Float64)
    return @sprintf("%.12g|%.12g", rho, mu)
end

function physical_mu_bounds(value::Float64, cfg::Dict{Symbol,Any}, ctx)
    params = SliceSweep.physical_parameters(value, cfg, ctx.gamma, ctx.fixed_mu_critical)
    mu_critical = params.mu_critical
    u_lo, u_hi = cfg[:u_range]
    return (Float64(u_lo * mu_critical), Float64(u_hi * mu_critical), Float64(mu_critical))
end

function evaluate_point(value::Float64, mu::Float64, cfg::Dict{Symbol,Any}, ctx)
    params = SliceSweep.physical_parameters(value, cfg, ctx.gamma, ctx.fixed_mu_critical)
    u = mu / params.mu_critical
    result = DynamicPhysics.evaluate_specloc_point(
        ctx.N, [params.t1, params.t2], mu, params.Delta, ctx.sequence, ctx.x0s, ctx.kappas;
        n_lowest_evals=ctx.n_lowest_evals,
        kk_tol=ctx.kk_tol,
        kk_maxiter=ctx.kk_maxiter,
        topological_threshold=ctx.topological_threshold,
        topological_fraction=ctx.topological_fraction,
    )
    return (
        point_key=point_key(value, mu), sweep_parameter=String(SliceSweep.sweep_parameter(cfg)),
        sweep_value=value, rho=params.rho, Delta=params.Delta,
        mu=mu, u=u, t1=params.t1, t2=params.t2,
        phase_label=result.phase_label, topo_fraction=result.topo_fraction,
        invariant_mean=result.invariant_mean, invariant_min=result.invariant_min,
        invariant_max=result.invariant_max, gap_min=result.gap_min,
        gap_mean=result.gap_mean, gap_max=result.gap_max, n_samples=result.n_samples,
    )
end

function evaluate_mus_parallel(rho::Float64, mus, cfg::Dict{Symbol,Any}, ctx)
    n = length(mus)
    n == 0 && return NamedTuple[]
    n_workers = max(1, min(Int(get(cfg, :n_eval_workers, Threads.nthreads())), Threads.nthreads(), n))
    log_every = max(1, Int(get(cfg, :progress_log_every_points, 25)))
    results = Vector{Any}(undef, n)
    jobs = Channel{Int}(n)
    for i in eachindex(mus)
        put!(jobs, i)
    end
    close(jobs)
    done = Threads.Atomic{Int}(0)
    t0 = time()

    @sync for _ in 1:n_workers
        Threads.@spawn begin
            for idx in jobs
                results[idx] = evaluate_point(rho, Float64(mus[idx]), cfg, ctx)
                completed = Threads.atomic_add!(done, 1) + 1
                if completed % log_every == 0 || completed == n
                    elapsed = time() - t0
                    DynamicIO.logmsg("$(SliceSweep.sweep_parameter(cfg))=$rho evaluation $completed/$n, rate=$(round(completed / max(elapsed, 1e-9), digits=2)) points/s")
                end
            end
        end
    end
    return NamedTuple[results...]
end

function count_phase_runs(samples, phase::Int)
    ordered = sort(collect(samples); by=s -> s.mu)
    count_runs = 0
    in_run = false
    for sample in ordered
        if sample.phase_label == phase
            if !in_run
                count_runs += 1
                in_run = true
            end
        else
            in_run = false
        end
    end
    return count_runs
end

function evaluate_mu_cached!(cache::Dict{String,NamedTuple}, rho::Float64, mu::Float64,
                              cfg::Dict{Symbol,Any}, ctx, total_points::Base.RefValue{Int})
    key = point_key(rho, mu)
    haskey(cache, key) && return cache[key]
    if length(cache) >= Int(cfg[:max_points_per_slice])
        throw(PointBudgetExceeded("Reached max_points_per_slice=$(cfg[:max_points_per_slice]) at $(SliceSweep.sweep_parameter(cfg))=$rho"))
    end
    if total_points[] >= Int(cfg[:max_total_points])
        throw(PointBudgetExceeded("Reached max_total_points=$(cfg[:max_total_points])"))
    end
    result = evaluate_point(rho, mu, cfg, ctx)
    cache[key] = result
    total_points[] += 1
    return result
end

function scan_slice(rho::Float64, slice_index::Int, cfg::Dict{Symbol,Any}, ctx,
                    predictions, total_points::Base.RefValue{Int}, slice_dir::String)
    mu_lo, mu_hi, mu_critical = physical_mu_bounds(rho, cfg, ctx)
    cache = Dict{String,NamedTuple}()
    baseline = SliceSearch.baseline_mu_grid(mu_lo, mu_hi, Int(cfg[:mu_baseline_points]))
    candidate_mus = copy(baseline)
    prediction_rows = NamedTuple[]
    band_half_width = Float64(get(cfg, :prediction_band_half_width_fraction, 0.03)) * (mu_hi - mu_lo)
    n_prediction_points = max(3, Int(get(cfg, :mu_refine_points_per_prediction, 9)))

    for prediction in predictions
        band_lo = max(mu_lo, prediction.mu - band_half_width)
        band_hi = min(mu_hi, prediction.mu + band_half_width)
        append!(candidate_mus, range(band_lo, band_hi; length=n_prediction_points))
        push!(prediction_rows, (
            track_id=prediction.track_id,
            predicted_mu=prediction.mu,
            band_lo=band_lo,
            band_hi=band_hi,
            orientation_lo=prediction.orientation[1],
            orientation_hi=prediction.orientation[2],
        ))
    end

    sort!(candidate_mus)
    unique!(candidate_mus)
    length(candidate_mus) <= Int(cfg[:max_points_per_slice]) ||
        throw(PointBudgetExceeded("Initial candidates exceed max_points_per_slice at $(SliceSweep.sweep_parameter(cfg))=$rho"))
    length(candidate_mus) + total_points[] <= Int(cfg[:max_total_points]) ||
        throw(PointBudgetExceeded("Initial candidates exceed max_total_points at $(SliceSweep.sweep_parameter(cfg))=$rho"))

    started = time()
    rows = evaluate_mus_parallel(rho, candidate_mus, cfg, ctx)
    for row in rows
        cache[row.point_key] = row
    end
    total_points[] += length(rows)

    eval_at_mu = mu -> evaluate_mu_cached!(cache, rho, Float64(mu), cfg, ctx, total_points)
    samples = collect(values(cache))
    added_gap_samples = SliceSearch.refine_low_gap_intervals!(
        eval_at_mu, samples, Float64(cfg[:gap_refine_tol]),
        Float64(cfg[:gap_refine_mu_tol]), Int(cfg[:gap_refine_max_depth]))
    samples = collect(values(cache))

    coarse_transitions = SliceSearch.detect_transitions(samples)
    refined_transitions = NamedTuple[]
    for transition in coarse_transitions
        left = cache[point_key(rho, transition.mu_lo)]
        right = cache[point_key(rho, transition.mu_hi)]
        push!(refined_transitions, SliceSearch.refine_transition_bracket(
            eval_at_mu, left, right, Float64(cfg[:transition_mu_tol]),
            Int(cfg[:transition_max_refine_iters])))
    end

    samples = sort(collect(values(cache)); by=s -> s.mu)
    gap_count = count_phase_runs(samples, 0)
    mkpath(slice_dir)
    CSV.write(joinpath(slice_dir, "point_results.csv"), DataFrame(samples))
    CSV.write(joinpath(slice_dir, "predictions.csv"), DataFrame(prediction_rows))

    params = SliceSweep.physical_parameters(rho, cfg, ctx.gamma, ctx.fixed_mu_critical)
    return (
        slice_index=slice_index,
        sweep_parameter=String(SliceSweep.sweep_parameter(cfg)),
        sweep_value=rho,
        rho=params.rho,
        Delta=params.Delta,
        gap_count=gap_count,
        transitions=refined_transitions,
        samples=samples,
        n_points=length(cache),
        n_baseline=length(baseline),
        n_prediction_points=length(candidate_mus) - length(baseline),
        n_gap_refinement_points=length(added_gap_samples),
        eval_seconds=time() - started,
        mu_min=mu_lo,
        mu_max=mu_hi,
        mu_critical=mu_critical,
        slice_dir=slice_dir,
    )
end

function write_summary(summary_rows, transition_rows, run_root::String)
    summary_df = DataFrame(
        slice_index=Int[], rho=Float64[], Delta=Float64[], sweep_value=Float64[], sweep_parameter=String[],
        mu_min=Float64[], mu_max=Float64[], mu_critical=Float64[], gap_count=Int[], n_transitions=Int[],
        n_points=Int[], n_baseline=Int[], n_prediction_points=Int[],
        n_gap_refinement_points=Int[], eval_seconds=Float64[],
    )
    for row in summary_rows
        push!(summary_df, (
            row.slice_index, row.rho, row.Delta, row.sweep_value, row.sweep_parameter,
            row.mu_min, row.mu_max, row.mu_critical, row.gap_count, length(row.transitions),
            row.n_points, row.n_baseline, row.n_prediction_points,
            row.n_gap_refinement_points, row.eval_seconds,
        ))
    end
    CSV.write(joinpath(run_root, "slice_summary.csv"), summary_df)

    transition_df = DataFrame(
        slice_index=Int[], rho=Float64[], Delta=Float64[], sweep_value=Float64[], sweep_parameter=String[], track_id=Int[], transition_index=Int[],
        mu=Float64[], mu_lo=Float64[], mu_hi=Float64[], phase_lo=Int[],
        phase_hi=Int[], gap_min=Float64[], refinement_iters=Int[], converged=Bool[],
    )
    for row in transition_rows
        push!(transition_df, row)
    end
    CSV.write(joinpath(run_root, "transition_boundaries.csv"), transition_df)
    return summary_df, transition_df
end

function save_checkpoint(run_root, summary_rows, transition_rows, tracks, total_points, stop_reason, parameter)
    summary_df, transition_df = write_summary(summary_rows, transition_rows, run_root)
    DynamicIO.write_jld2(
        joinpath(run_root, "run_checkpoint.jld2");
        timestamp=string(Dates.now()),
        stop_reason=stop_reason,
        total_points=total_points,
        sweep_parameter=parameter,
        completed_slices=[r.sweep_value for r in summary_rows],
        gap_counts=[r.gap_count for r in summary_rows],
        summary_columns=Dict(String(n) => summary_df[!, n] for n in names(summary_df)),
        transition_columns=Dict(String(n) => transition_df[!, n] for n in names(transition_df)),
        tracks=tracks,
    )
end

function onset_brackets(summary_rows, tolerance::Float64, parameter::Symbol)
    ordered = sort(collect(summary_rows); by=s -> s.sweep_value)
    brackets = NamedTuple[]
    for i in 1:(length(ordered) - 1)
        left, right = ordered[i], ordered[i + 1]
        if left.gap_count != right.gap_count && right.sweep_value - left.sweep_value <= tolerance
            names = (Symbol("$(parameter)_lo"), Symbol("$(parameter)_hi"), :gap_count_lo, :gap_count_hi)
            push!(brackets, NamedTuple{names}((left.sweep_value, right.sweep_value,
                                               left.gap_count, right.gap_count)))
        end
    end
    return brackets
end

function validate_config(cfg::Dict{Symbol,Any})
    cfg[:mu_baseline_points] >= 3 || error("mu_baseline_points must be at least 3")
    cfg[:transition_mu_tol] > 0 || error("transition_mu_tol must be positive")
    cfg[:gap_refine_max_depth] >= 0 || error("gap_refine_max_depth cannot be negative")
    parameter = SliceSweep.sweep_parameter(cfg)
    SliceSweep.sweep_grid(cfg)
    tolerance_key = parameter == :rho ? :rho_onset_tol : :Delta_onset_tol
    budget_key = parameter == :rho ? :rho_onset_max_extra_slices : :Delta_onset_max_extra_slices
    cfg[tolerance_key] > 0 || error("$tolerance_key must be positive")
    cfg[budget_key] >= 0 || error("$budget_key cannot be negative")
    parameter == :Delta && SliceSweep.fixed_mu_critical(cfg, Float64(cfg[:target_slope]))
    return cfg
end

function main()
    Base.exit_on_sigint(false)
    config_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(project_root, "config", "default_slice_config.jl")
    cfg = validate_config(DynamicIO.read_config(config_path))
    BLAS.set_num_threads(max(1, Int(cfg[:blas_threads])))

    run_root = DynamicIO.ensure_dir(joinpath(cfg[:output_root], "$(cfg[:run_name])_$(DynamicIO.now_tag())"))
    DynamicIO.logmsg("Slice-search output: $run_root")
    pointer = String(cfg[:run_root_pointer_file])
    DynamicIO.ensure_dir(dirname(pointer))
    open(pointer, "w") do io
        println(io, run_root)
    end
    stop_files = [joinpath(run_root, "STOP_REQUESTED"), String(cfg[:global_stop_file])]

    ctx = make_eval_context(cfg)
    DynamicIO.logmsg("Selected sequence slope/phason: $(ctx.gamma), $(ctx.phason)")
    parameter = SliceSweep.sweep_parameter(cfg)
    sweep_values = SliceSweep.sweep_grid(cfg)
    onset_tol = Float64(cfg[parameter == :rho ? :rho_onset_tol : :Delta_onset_tol])
    onset_budget = Int(cfg[parameter == :rho ? :rho_onset_max_extra_slices : :Delta_onset_max_extra_slices])
    parameter == :Delta && DynamicIO.logmsg("Fixed mu_c estimate: $(ctx.fixed_mu_critical)")
    tracks = EdgeTrack[]
    summary_rows = NamedTuple[]
    transition_rows = NamedTuple[]
    total_points = Ref(0)
    slice_index = 0
    stop_reason = "completed"
    onset_queue = Tuple{NamedTuple,NamedTuple}[]

    try
        for (grid_index, rho) in enumerate(sweep_values)
            if stop_requested(stop_files)
                stop_reason = "stop_requested"
                break
            end
            slice_index += 1
            mu_lo, mu_hi, _ = physical_mu_bounds(Float64(rho), cfg, ctx)
            predictions = TransitionTracking.predict_edges(tracks, rho, slice_index, (mu_lo, mu_hi), cfg)
            slice_dir = DynamicIO.ensure_dir(joinpath(run_root, "slices", "slice_$(lpad(slice_index, 4, '0'))"))
            result = scan_slice(Float64(rho), slice_index, cfg, ctx, predictions, total_points, slice_dir)
            push!(summary_rows, result)

            transition_ids = TransitionTracking.update_tracks!(
                tracks, result.transitions, result.sweep_value, slice_index, (result.mu_min, result.mu_max), cfg)
            for (j, transition) in enumerate(result.transitions)
                push!(transition_rows, (
                    slice_index=slice_index,
                    rho=result.rho,
                    Delta=result.Delta,
                    sweep_value=result.sweep_value,
                    sweep_parameter=result.sweep_parameter,
                    track_id=transition_ids[j],
                    transition_index=j,
                    mu=transition.mu,
                    mu_lo=transition.mu_lo,
                    mu_hi=transition.mu_hi,
                    phase_lo=transition.phase_lo,
                    phase_hi=transition.phase_hi,
                    gap_min=transition.gap_min,
                    refinement_iters=transition.refinement_iters,
                    converged=transition.converged,
                ))
            end

            if length(summary_rows) > 1
                previous = summary_rows[end - 1]
                result.gap_count != previous.gap_count && push!(onset_queue, (previous, result))
            end
            save_checkpoint(run_root, summary_rows, transition_rows, tracks, total_points[], "in_progress", parameter)
            DynamicIO.logmsg("Slice $(grid_index)/$(length(sweep_values)): $parameter=$rho gaps=$(result.gap_count), transitions=$(length(result.transitions)), points=$(result.n_points), total=$(total_points[])")
        end

        extra_slices = 0
        while stop_reason == "completed" && !isempty(onset_queue)
            stop_requested(stop_files) && (stop_reason = "stop_requested_during_onset_refinement"; break)
            left, right = popfirst!(onset_queue)
            right.sweep_value - left.sweep_value <= onset_tol && continue
            if extra_slices >= onset_budget
                pushfirst!(onset_queue, (left, right))
                break
            end

            midpoint = 0.5 * (left.sweep_value + right.sweep_value)
            slice_index += 1
            slice_dir = DynamicIO.ensure_dir(joinpath(run_root, "slices", "slice_$(lpad(slice_index, 4, '0'))"))
            result = scan_slice(midpoint, slice_index, cfg, ctx, NamedTuple[], total_points, slice_dir)
            push!(summary_rows, result)
            for (j, transition) in enumerate(result.transitions)
                push!(transition_rows, (
                    slice_index=slice_index,
                    rho=result.rho,
                    Delta=result.Delta,
                    sweep_value=result.sweep_value,
                    sweep_parameter=result.sweep_parameter,
                    track_id=0,
                    transition_index=j,
                    mu=transition.mu,
                    mu_lo=transition.mu_lo,
                    mu_hi=transition.mu_hi,
                    phase_lo=transition.phase_lo,
                    phase_hi=transition.phase_hi,
                    gap_min=transition.gap_min,
                    refinement_iters=transition.refinement_iters,
                    converged=transition.converged,
                ))
            end
            extra_slices += 1

            result.gap_count != left.gap_count && push!(onset_queue, (left, result))
            right.gap_count != result.gap_count && push!(onset_queue, (result, right))
            save_checkpoint(run_root, summary_rows, transition_rows, tracks, total_points[], "onset_refinement", parameter)
            DynamicIO.logmsg("Onset probe $parameter=$midpoint: gaps=$(result.gap_count), bracket width=$(right.sweep_value-left.sweep_value)")
        end
        if !isempty(onset_queue) && stop_reason == "completed"
            stop_reason = "onset_refinement_limit"
        end
    catch err
        if err isa PointBudgetExceeded
            stop_reason = "point_budget"
            DynamicIO.logmsg("Stopping at point budget: $(err.message)")
        elseif err isa InterruptException
            stop_reason = "interrupt"
            DynamicIO.logmsg("Interrupt received; preserving completed slices")
        else
            rethrow(err)
        end
    end

    save_checkpoint(run_root, summary_rows, transition_rows, tracks, total_points[], stop_reason, parameter)
    brackets = onset_brackets(summary_rows, onset_tol, parameter)
    onset_df = DataFrame(
        Symbol("$(parameter)_lo") => Float64[], Symbol("$(parameter)_hi") => Float64[],
        :gap_count_lo => Int[], :gap_count_hi => Int[])
    for bracket in brackets
        push!(onset_df, bracket)
    end
    CSV.write(joinpath(run_root, "onset_brackets.csv"), onset_df)
    DynamicIO.write_jld2(
        joinpath(run_root, "run_summary.jld2");
        run_root=run_root,
        config_path=config_path,
        stop_reason=stop_reason,
        timestamp=string(Dates.now()),
        total_points=total_points[],
        slices_completed=length(summary_rows),
        onset_brackets=brackets,
        sweep_parameter=parameter,
        sweep_values=sort([r.sweep_value for r in summary_rows]),
        rho_values=[r.rho for r in sort(summary_rows; by=s -> s.sweep_value)],
        Delta_values=[r.Delta for r in sort(summary_rows; by=s -> s.sweep_value)],
        fixed_mu_critical=ctx.fixed_mu_critical,
        gap_counts=[r.gap_count for r in sort(summary_rows; by=s -> s.sweep_value)],
        julia_threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),
    )
    DynamicIO.logmsg("Slice search complete: stop_reason=$stop_reason, slices=$(length(summary_rows)), points=$(total_points[]), output=$run_root")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end