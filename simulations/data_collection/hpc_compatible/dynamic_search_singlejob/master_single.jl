#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using Statistics
using LinearAlgebra
using Base.Threads

project_root = @__DIR__
include(joinpath(project_root, "core", "adaptive.jl"))
include(joinpath(project_root, "core", "io_utils.jl"))
include(joinpath(project_root, "core", "physics.jl"))
include(joinpath(project_root, "core", "sequences.jl"))
include(joinpath(project_root, "core", "bandwidth.jl"))

using .DynamicAdaptive
using .DynamicIO
using .DynamicPhysics
using .DynamicSequences
using .DynamicBandwidth

function stop_requested(stop_files::Vector{String})
    return any(isfile(f) for f in stop_files)
end

function row_to_result_tuple(r)
    return (
        phase_label = Int(r.phase_label),
        topo_fraction = Float64(r.topo_fraction),
        invariant_mean = Float64(r.invariant_mean),
        invariant_min = Float64(r.invariant_min),
        invariant_max = Float64(r.invariant_max),
        gap_min = Float64(r.gap_min),
        gap_mean = Float64(r.gap_mean),
        gap_max = Float64(r.gap_max),
        n_samples = Int(r.n_samples),
    )
end

function write_point_results_csv(out_path::String, known_results::Dict{String,NamedTuple})
    point_df = DataFrame(
        point_key = String[],
        phase_label = Int[],
        topo_fraction = Float64[],
        invariant_mean = Float64[],
        invariant_min = Float64[],
        invariant_max = Float64[],
        gap_min = Float64[],
        gap_mean = Float64[],
        gap_max = Float64[],
        n_samples = Int[],
    )

    for (k, v) in known_results
        push!(point_df, (k, v.phase_label, v.topo_fraction, v.invariant_mean, v.invariant_min, v.invariant_max, v.gap_min, v.gap_mean, v.gap_max, v.n_samples))
    end

    CSV.write(out_path, point_df)
    return point_df
end

function save_cells_csv(cells, out_path::String)
    df = DataFrame(
        id = Int[],
        depth = Int[],
        rho_lo = Float64[],
        rho_hi = Float64[],
        u_lo = Float64[],
        u_hi = Float64[],
    )

    for c in cells
        push!(df, (c.id, c.depth, c.rho_lo, c.rho_hi, c.u_lo, c.u_hi))
    end

    CSV.write(out_path, df)
end

function save_resolved_cells_csv(cells, out_path::String)
    df = DataFrame(
        id = Int[],
        depth = Int[],
        rho_lo = Float64[],
        rho_hi = Float64[],
        u_lo = Float64[],
        u_hi = Float64[],
        phase_label = Int[],
        confidence_gap = Float64[],
        mixed_phase = Bool[],
    )

    for c in cells
        push!(df, (c.id, c.depth, c.rho_lo, c.rho_hi, c.u_lo, c.u_hi, c.phase_label, c.confidence_gap, c.mixed_phase))
    end

    CSV.write(out_path, df)
end

function save_checkpoint(
    run_root::String,
    known_results::Dict{String,NamedTuple},
    history_df::DataFrame,
    gap_history::Vector{Int},
    stop_reason::String,
    iter::Int,
)
    point_df = write_point_results_csv(joinpath(run_root, "point_results_compact.csv"), known_results)
    CSV.write(joinpath(run_root, "search_history.csv"), history_df)

    history_cols = Dict(String(n) => history_df[!, n] for n in names(history_df))
    point_cols = Dict(String(n) => point_df[!, n] for n in names(point_df))

    DynamicIO.write_jld2(
        joinpath(run_root, "run_checkpoint.jld2");
        iter=iter,
        stop_reason=stop_reason,
        timestamp=string(Dates.now()),
        gap_history=gap_history,
        history_columns=history_cols,
        point_columns=point_cols,
    )
end

function build_x0s(N::Int, alpha::Float64, npoints::Int)
    return collect(range(N / 2 - alpha * N, N / 2 + alpha * N, length=npoints))
end

function make_eval_context(cfg::Dict{Symbol,Any})
    seqs, phis, phasons = DynamicSequences.load_sequences_from_bson(cfg[:sequence_bson_path])
    _, sequence, gamma, phason = DynamicSequences.select_sequence(
        seqs,
        phis,
        phasons;
        target_slope=cfg[:target_slope],
        slope_tolerance=cfg[:slope_tolerance],
        target_phason=cfg[:target_phason],
        phason_tolerance=cfg[:phason_tolerance],
        selection_mode=cfg[:sequence_selection_mode],
    )

    N = Int(cfg[:N])
    return (
        N=N,
        Delta=Float64(cfg[:Delta]),
        sequence=sequence,
        gamma=gamma,
        phason=phason,
        x0s=build_x0s(N, Float64(cfg[:x0_alpha]), Int(cfg[:x0_count])),
        kappas=Vector{Float64}(cfg[:kappas]),
        n_lowest_evals=Int(cfg[:n_lowest_evals]),
        kk_tol=Float64(cfg[:kk_tol]),
        kk_maxiter=Int(cfg[:kk_maxiter]),
        topological_threshold=Float64(cfg[:topological_threshold]),
        topological_fraction=Float64(cfg[:topological_fraction]),
    )
end

function evaluate_candidate_point(p, cfg::Dict{Symbol,Any}, ctx)
    rho = Float64(p.rho)
    u = Float64(p.u)

    t1, t2 = DynamicBandwidth.tpair_from_rho(
        rho;
        sweep_mode=cfg[:sweep_mode],
        t1_fixed=cfg[:t1_fixed],
        t2_fixed=cfg[:t2_fixed],
    )

    mu_c_est = DynamicBandwidth.mu_critical_estimate(ctx.gamma, t1, t2)
    mu = u * mu_c_est

    eval_out = DynamicPhysics.evaluate_specloc_point(
        ctx.N,
        [t1, t2],
        mu,
        ctx.Delta,
        ctx.sequence,
        ctx.x0s,
        ctx.kappas;
        n_lowest_evals=ctx.n_lowest_evals,
        kk_tol=ctx.kk_tol,
        kk_maxiter=ctx.kk_maxiter,
        topological_threshold=ctx.topological_threshold,
        topological_fraction=ctx.topological_fraction,
    )

    return (
        point_key=String(p.point_key),
        rho=rho,
        u=u,
        mu=mu,
        t1=t1,
        t2=t2,
        phase_label=eval_out.phase_label,
        topo_fraction=eval_out.topo_fraction,
        invariant_mean=eval_out.invariant_mean,
        invariant_min=eval_out.invariant_min,
        invariant_max=eval_out.invariant_max,
        gap_min=eval_out.gap_min,
        gap_mean=eval_out.gap_mean,
        gap_max=eval_out.gap_max,
        n_samples=eval_out.n_samples,
    )
end

function evaluate_candidates_parallel(candidates, cfg::Dict{Symbol,Any}, ctx, stop_files::Vector{String})
    n = length(candidates)
    n == 0 && return NamedTuple[]

    requested = Int(get(cfg, :n_eval_workers, Threads.nthreads()))
    n_workers = max(1, min(requested, Threads.nthreads(), n))
    log_every = max(1, Int(get(cfg, :progress_log_every_points, 25)))

    DynamicIO.logmsg("Point evaluation: n_points=$n, workers=$n_workers, julia_threads=$(Threads.nthreads())")

    results = Vector{Any}(undef, n)
    jobs = Channel{Int}(n)
    for i in 1:n
        put!(jobs, i)
    end
    close(jobs)

    done_ctr = Threads.Atomic{Int}(0)
    t0 = time()

    @sync for _ in 1:n_workers
        Threads.@spawn begin
            for idx in jobs
                if stop_requested(stop_files)
                    continue
                end
                results[idx] = evaluate_candidate_point(candidates[idx], cfg, ctx)
                old = Threads.atomic_add!(done_ctr, 1)
                done = old + 1
                if (done % log_every == 0) || (done == n)
                    elapsed = time() - t0
                    rate = done / max(elapsed, 1e-9)
                    DynamicIO.logmsg("Evaluation progress: $done/$n points, elapsed=$(round(elapsed, digits=1))s, rate=$(round(rate, digits=2)) pts/s")
                end
            end
        end
    end

    rows = NamedTuple[]
    for r in results
        r === nothing && continue
        push!(rows, r)
    end
    return rows
end

function estimate_iteration_time(n_points::Int, cfg::Dict{Symbol,Any})
    est_rate = Float64(get(cfg, :estimated_points_per_second, 44.0))
    if est_rate <= 0
        return 0.0
    end
    return n_points / est_rate
end

function main()
    Base.exit_on_sigint(false)

    config_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(project_root, "config", "default_single_config.jl")
    cfg = DynamicIO.read_config(config_path)

    blas_threads = Int(get(cfg, :blas_threads, 1))
    BLAS.set_num_threads(max(1, blas_threads))

    run_root = DynamicIO.ensure_dir(joinpath(cfg[:output_root], "$(cfg[:run_name])_$(DynamicIO.now_tag())"))
    DynamicIO.logmsg("Run directory: $run_root")

    run_root_pointer = String(get(cfg, :run_root_pointer_file, joinpath(cfg[:output_root], "latest_run_root.txt")))
    DynamicIO.ensure_dir(dirname(run_root_pointer))
    open(run_root_pointer, "w") do io
        println(io, run_root)
    end

    local_stop_file = joinpath(run_root, "STOP_REQUESTED")
    global_stop_file = String(get(cfg, :global_stop_file, joinpath(cfg[:output_root], "STOP_REQUESTED")))
    stop_files = [local_stop_file, global_stop_file]
    DynamicIO.logmsg("Early-stop files: $(join(stop_files, ", "))")

    ctx = make_eval_context(cfg)
    DynamicIO.logmsg("Selected sequence slope/phason: $(ctx.gamma), $(ctx.phason)")

    cells = DynamicAdaptive.initial_cells(cfg)
    known_results = Dict{String,NamedTuple}()

    history_df = DataFrame(
        iter = Int[],
        timestamp = String[],
        n_cells = Int[],
        n_points_total = Int[],
        n_candidates = Int[],
        n_refined_children = Int[],
        n_resolved = Int[],
        rho_max_trivial_gaps = Int[],
        eval_seconds = Float64[],
        eval_points_per_second = Float64[],
    )

    gap_history = Int[]
    stop_reason = "completed"
    max_iters = Int(cfg[:max_iters])
    max_points = Int(cfg[:max_points])

    total_eval_points = 0
    total_eval_seconds = 0.0

    try
        for iter in 1:max_iters
            if stop_requested(stop_files)
                stop_reason = "stop_requested_before_iteration"
                DynamicIO.logmsg("Stop requested before iteration $iter. Ending loop safely.")
                break
            end

            iter_dir = DynamicIO.ensure_dir(joinpath(run_root, "iter_$(lpad(iter, 3, '0'))"))
            points_dir = DynamicIO.ensure_dir(joinpath(iter_dir, "points"))
            results_dir = DynamicIO.ensure_dir(joinpath(iter_dir, "results"))

            candidates = DynamicAdaptive.collect_candidate_points(cells, known_results, cfg)
            n_candidates = length(candidates)

            if n_candidates == 0
                stop_reason = "no_new_candidates"
                DynamicIO.logmsg("No new candidate points. Terminating.")
                break
            end

            points_csv = joinpath(points_dir, "points_all.csv")
            DynamicIO.write_points_csv(candidates, points_csv)
            save_cells_csv(cells, joinpath(iter_dir, "cells_input.csv"))

            est_iter_seconds = estimate_iteration_time(n_candidates, cfg)
            if est_iter_seconds > 0
                DynamicIO.logmsg("Iter $iter: cells=$(length(cells)), new_points=$n_candidates, est_eval_time=$(round(est_iter_seconds, digits=1))s")
            else
                DynamicIO.logmsg("Iter $iter: cells=$(length(cells)), new_points=$n_candidates")
            end

            eval_t0 = time()
            rows = evaluate_candidates_parallel(candidates, cfg, ctx, stop_files)
            eval_seconds = time() - eval_t0

            total_eval_points += length(rows)
            total_eval_seconds += eval_seconds

            DynamicIO.write_worker_results(rows, joinpath(results_dir, "worker_results_local.csv"))
            worker_results_df = DynamicIO.read_worker_results(results_dir)

            n_result_rows = nrow(worker_results_df)
            n_unique_points = n_result_rows == 0 ? 0 : length(unique(worker_results_df.point_key))
            iter_rate = n_result_rows / max(eval_seconds, 1e-9)

            DynamicIO.logmsg("Iter $iter worker outputs: rows=$n_result_rows, unique_points=$n_unique_points/$n_candidates")
            DynamicIO.logmsg("Iter $iter evaluation timing: $(round(eval_seconds, digits=2))s, rate=$(round(iter_rate, digits=2)) pts/s")

            if n_result_rows == 0
                if startswith(stop_reason, "stop_requested") || stop_requested(stop_files)
                    stop_reason = "stop_requested_during_evaluation"
                    DynamicIO.logmsg("No worker results collected before stop; preserving state and exiting.")
                    save_checkpoint(run_root, known_results, history_df, gap_history, stop_reason, iter)
                    break
                end
                error("No local evaluation results found in $results_dir")
            end

            require_all = Bool(get(cfg, :require_all_worker_results, true))
            if require_all && n_unique_points < n_candidates && !stop_requested(stop_files)
                error("Incomplete local results in iter $iter: expected $n_candidates unique points, got $n_unique_points")
            end

            for r in eachrow(worker_results_df)
                known_results[String(r.point_key)] = row_to_result_tuple(r)
            end

            refine_out = DynamicAdaptive.classify_and_refine(cells, known_results, cfg)
            cells = refine_out.next_cells
            resolved_cells = refine_out.resolved_cells

            rho_max_gap_count = DynamicAdaptive.estimate_trivial_gap_count_at_rho_max(resolved_cells, cfg)
            push!(gap_history, rho_max_gap_count)

            save_cells_csv(cells, joinpath(iter_dir, "cells_next.csv"))
            save_resolved_cells_csv(resolved_cells, joinpath(iter_dir, "cells_resolved.csv"))

            push!(history_df, (
                iter,
                string(Dates.now()),
                length(cells),
                length(known_results),
                n_candidates,
                refine_out.n_refined,
                length(resolved_cells),
                rho_max_gap_count,
                eval_seconds,
                iter_rate,
            ))

            save_checkpoint(run_root, known_results, history_df, gap_history, "in_progress", iter)
            DynamicIO.logmsg("Iter $iter summary: refined_children=$(refine_out.n_refined), rho_max_trivial_gaps=$rho_max_gap_count")

            if length(known_results) >= max_points
                stop_reason = "max_points"
                DynamicIO.logmsg("Reached max_points=$(max_points). Terminating.")
                break
            end

            if refine_out.n_refined == 0
                stop_reason = "no_refinement"
                DynamicIO.logmsg("No cells refined in this iteration. Terminating.")
                break
            end

            stable_iters = Int(cfg[:gap_count_stable_iters])
            gap_tol = Int(cfg[:gap_count_tol])
            if length(gap_history) >= stable_iters
                window = gap_history[end-stable_iters+1:end]
                if (maximum(window) - minimum(window)) <= gap_tol
                    stop_reason = "stable_gap_count"
                    DynamicIO.logmsg("Gap count stable across last $(stable_iters) iterations. Terminating.")
                    break
                end
            end

            max_expected_gaps = Int(cfg[:max_expected_gaps])
            max_expected_gaps_tol = Int(cfg[:max_expected_gaps_tol])
            if max_expected_gaps > 0 && rho_max_gap_count >= (max_expected_gaps - max_expected_gaps_tol)
                stop_reason = "expected_gap_cap"
                DynamicIO.logmsg("Reached expected gap-count target ($(rho_max_gap_count) >= $(max_expected_gaps - max_expected_gaps_tol)). Terminating.")
                break
            end
        end
    catch err
        if err isa InterruptException
            stop_reason = "interrupt"
            DynamicIO.logmsg("Interrupt received. Preserving current state and exiting cleanly.")
        else
            rethrow(err)
        end
    end

    save_checkpoint(run_root, known_results, history_df, gap_history, stop_reason, nrow(history_df))

    overall_rate = total_eval_points / max(total_eval_seconds, 1e-9)
    DynamicIO.write_jld2(
        joinpath(run_root, "run_summary.jld2");
        run_root=run_root,
        config_path=config_path,
        stop_reason=stop_reason,
        timestamp=string(Dates.now()),
        n_points_total=length(known_results),
        n_iterations_completed=nrow(history_df),
        gap_history=gap_history,
        total_eval_points=total_eval_points,
        total_eval_seconds=total_eval_seconds,
        overall_points_per_second=overall_rate,
        n_eval_workers=min(Int(get(cfg, :n_eval_workers, Threads.nthreads())), Threads.nthreads()),
        julia_threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),
    )

    DynamicIO.logmsg("Search complete. stop_reason=$stop_reason")
    DynamicIO.logmsg("Total evaluated points: $total_eval_points in $(round(total_eval_seconds, digits=2))s (rate=$(round(overall_rate, digits=2)) pts/s)")
    DynamicIO.logmsg("Outputs in: $run_root")
end

main()
