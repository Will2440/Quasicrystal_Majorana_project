#!/usr/bin/env julia

using CSV
using DataFrames
using Dates
using Statistics

project_root = @__DIR__
include(joinpath(project_root, "core", "adaptive.jl"))
include(joinpath(project_root, "core", "io_utils.jl"))

using .DynamicAdaptive
using .DynamicIO

function effective_max_workers(n_chunks::Int, cfg::Dict{Symbol,Any})
    hard_cap = Int(get(cfg, :hard_max_workers, 199))
    user_cap = Int(get(cfg, :max_workers_per_array, hard_cap))
    user_cap > hard_cap && DynamicIO.logmsg("Requested max_workers_per_array=$user_cap exceeds hard_max_workers=$hard_cap; clamping.")
    return max(1, min(n_chunks, user_cap, hard_cap))
end

function write_active_jobs(path::String, jobids::Vector{String})
    open(path, "w") do io
        for j in jobids
            println(io, j)
        end
    end
end

function clear_active_jobs(path::String)
    open(path, "w") do io
        write(io, "")
    end
end

function worker_runtime_profile(cfg::Dict{Symbol,Any})
    fallback_time = String(get(cfg, :fallback_worker_time, "02:00:00"))
    fallback_mem_mb = Int(get(cfg, :fallback_worker_mem_mb, 4096))
    return Dict{Symbol,Any}(
        :time_str => fallback_time,
        :mem_mb => fallback_mem_mb,
        :source => "fallback",
        :calibrated => false,
        :elapsed_seconds => 0.0,
        :maxrss_mb => 0.0,
    )
end

function build_common_sbatch_args(cfg::Dict{Symbol,Any}; for_calibration::Bool=false)
    args = String[]

    partition_key = for_calibration ? :calibration_partition : :slurm_partition
    account_key = for_calibration ? :calibration_account : :slurm_account

    partition = strip(String(get(cfg, partition_key, "")))
    account = strip(String(get(cfg, account_key, "")))
    qos = strip(String(get(cfg, :slurm_qos, "")))
    constraint = strip(String(get(cfg, :slurm_constraint, "")))
    reservation = strip(String(get(cfg, :slurm_reservation, "")))

    !isempty(partition) && append!(args, ["--partition", partition])
    !isempty(account) && append!(args, ["--account", account])
    !isempty(qos) && append!(args, ["--qos", qos])
    !isempty(constraint) && append!(args, ["--constraint", constraint])
    !isempty(reservation) && append!(args, ["--reservation", reservation])

    return args
end

function build_worker_resource_args(cfg::Dict{Symbol,Any}, runtime_profile::Dict{Symbol,Any})
    args = build_common_sbatch_args(cfg)
    push!(args, "--time", String(runtime_profile[:time_str]))
    push!(args, "--mem", "$(Int(runtime_profile[:mem_mb]))M")
    return args
end

function parse_slurm_memory_mb(mem::AbstractString)
    s = strip(mem)
    isempty(s) && return 0.0
    s == "0" && return 0.0

    m = match(r"^([0-9]*\.?[0-9]+)([KMGTP])(?:i?B?)?(?:n|c)?$", s)
    m === nothing && return 0.0

    v = parse(Float64, m.captures[1])
    unit = m.captures[2]
    factor = unit == "K" ? (1.0 / 1024.0) :
             unit == "M" ? 1.0 :
             unit == "G" ? 1024.0 :
             unit == "T" ? 1024.0 * 1024.0 :
             unit == "P" ? 1024.0 * 1024.0 * 1024.0 : 0.0
    return v * factor
end

function estimate_job_elapsed_seconds(jobid::String)
    try
        out = readchomp(`sacct -n -X -j $jobid --format=ElapsedRaw`)
        vals = parse_elapsed_seconds(split(out, '\n'))
        return isempty(vals) ? 0.0 : maximum(vals)
    catch err
        DynamicIO.logmsg("Could not query elapsed runtime for $jobid ($err).")
        return 0.0
    end
end

function estimate_job_maxrss_mb(jobid::String)
    try
        out = readchomp(`sacct -n -X -j $jobid --format=MaxRSS`)
        vals = [parse_slurm_memory_mb(ln) for ln in split(out, '\n')]
        vals = filter(x -> x > 0, vals)
        return isempty(vals) ? 0.0 : maximum(vals)
    catch err
        DynamicIO.logmsg("Could not query MaxRSS for $jobid ($err).")
        return 0.0
    end
end

function wait_for_job_simple(jobid::String, poll_seconds::Int, stop_files::Vector{String})
    while true
        if stop_requested(stop_files)
            DynamicIO.logmsg("Stop requested during calibration; cancelling calibration job $jobid")
            try
                run(`scancel $jobid`)
            catch err
                DynamicIO.logmsg("scancel failed for calibration job $jobid: $err")
            end
            return :stopped
        end

        q = readchomp(`squeue -h -j $jobid`)
        isempty(strip(q)) && return :done
        sleep(poll_seconds)
    end
end

function build_calibration_points(cfg::Dict{Symbol,Any}, n_points::Int)
    n = max(1, n_points)
    rho_lo, rho_hi = cfg[:rho_range]
    u_lo, u_hi = cfg[:u_range]

    rhos = collect(range(rho_lo, rho_hi, length=n))
    us = collect(range(u_lo, u_hi, length=n))

    points = NamedTuple[]
    for i in 1:n
        rho = rhos[i]
        u = us[min(i, length(us))]
        key = string(round(rho, digits=12), "|", round(u, digits=12))
        push!(points, (point_key=key, rho=rho, u=u))
    end
    return points
end

function calibrate_worker_resources(
    cfg::Dict{Symbol,Any},
    config_path::String,
    run_root::String,
    stop_files::Vector{String},
)
    profile = worker_runtime_profile(cfg)

    if cfg[:execution_mode] != :submit_and_wait
        DynamicIO.logmsg("Skipping calibration in execution_mode=$(cfg[:execution_mode]).")
        return profile
    end

    if !Bool(get(cfg, :auto_calibrate_worker_resources, true))
        DynamicIO.logmsg("Auto calibration disabled. Using fallback worker resources.")
        return profile
    end

    calib_dir = DynamicIO.ensure_dir(joinpath(run_root, "calibration"))
    points_dir = DynamicIO.ensure_dir(joinpath(calib_dir, "points"))
    results_dir = DynamicIO.ensure_dir(joinpath(calib_dir, "results"))

    n_points = Int(get(cfg, :calibration_points, 3))
    points = build_calibration_points(cfg, n_points)
    points_csv = joinpath(points_dir, "points_all.csv")
    DynamicIO.write_points_csv(points, points_csv)
    chunk_size = max(1, n_points)
    DynamicIO.split_points_into_chunks(points_csv, points_dir, chunk_size)

    script = cfg[:sbatch_script]
    seed_time = String(get(cfg, :calibration_time_seed, "00:30:00"))
    seed_mem_mb = Int(get(cfg, :calibration_mem_seed_mb, Int(profile[:mem_mb])))
    poll_seconds = Int(get(cfg, :calibration_poll_seconds, 15))

    args = String["sbatch", "--array=1-1%1"]
    append!(args, build_common_sbatch_args(cfg; for_calibration=true))
    append!(args, ["--time", seed_time, "--mem", "$(seed_mem_mb)M"])

    extra_args = strip(String(get(cfg, :sbatch_extra_args, "")))
    if !isempty(extra_args)
        append!(args, split(extra_args))
    end

    append!(args, [script, config_path, calib_dir])

    DynamicIO.logmsg("Submitting calibration worker with n_points=$n_points, time=$seed_time, mem=$(seed_mem_mb)M")
    out = readchomp(Cmd(args))
    m = match(r"Submitted batch job (\d+)", out)
    m === nothing && error("Could not parse calibration sbatch output: $out")
    jobid = String(m.captures[1])
    DynamicIO.logmsg(out)

    st = wait_for_job_simple(jobid, poll_seconds, stop_files)
    st == :stopped && return profile

    elapsed_seconds = estimate_job_elapsed_seconds(jobid)
    maxrss_mb = estimate_job_maxrss_mb(jobid)

    if elapsed_seconds <= 0 || maxrss_mb <= 0
        msg = "Calibration metrics unavailable from sacct for job $jobid; using fallback worker resources."
        if Bool(get(cfg, :calibration_require_success, false))
            error(msg)
        end
        DynamicIO.logmsg(msg)
        return profile
    end

    points_per_job = Int(cfg[:points_per_job])
    sec_per_point = elapsed_seconds / max(n_points, 1)

    time_sf = Float64(get(cfg, :calibration_time_safety_factor, 1.8))
    time_buf = Int(get(cfg, :calibration_time_buffer_seconds, 120))
    min_time = Int(get(cfg, :min_worker_time_seconds, 300))

    est_time = ceil(Int, sec_per_point * points_per_job * time_sf + time_buf)
    est_time = max(est_time, min_time)

    mem_sf = Float64(get(cfg, :calibration_mem_safety_factor, 1.6))
    mem_buf = Int(get(cfg, :calibration_mem_buffer_mb, 512))
    min_mem = Int(get(cfg, :min_worker_mem_mb, 2048))

    est_mem = ceil(Int, maxrss_mb * mem_sf + mem_buf)
    est_mem = max(est_mem, min_mem)

    profile[:time_str] = format_seconds(est_time)
    profile[:mem_mb] = est_mem
    profile[:source] = "calibrated"
    profile[:calibrated] = true
    profile[:elapsed_seconds] = elapsed_seconds
    profile[:maxrss_mb] = maxrss_mb

    DynamicIO.write_jld2(
        joinpath(calib_dir, "calibration_summary.jld2");
        calibration_jobid=jobid,
        calibration_points=n_points,
        measured_elapsed_seconds=elapsed_seconds,
        measured_maxrss_mb=maxrss_mb,
        sec_per_point=sec_per_point,
        estimated_worker_time_seconds=est_time,
        estimated_worker_time_str=String(profile[:time_str]),
        estimated_worker_mem_mb=est_mem,
        points_per_job=points_per_job,
        time_safety_factor=time_sf,
        mem_safety_factor=mem_sf,
    )

    DynamicIO.logmsg(
        "Calibration complete: elapsed=$(round(elapsed_seconds, digits=2))s, maxrss=$(round(maxrss_mb, digits=1))MB, " *
        "recommended worker time=$(profile[:time_str]), mem=$(est_mem)M"
    )

    return profile
end

function submit_array_job(
    config_path::String,
    iter_dir::String,
    n_chunks::Int,
    cfg::Dict{Symbol,Any},
    runtime_profile::Dict{Symbol,Any},
)
    script = cfg[:sbatch_script]
    max_workers = effective_max_workers(n_chunks, cfg)
    array_spec = "1-$(n_chunks)%$(max_workers)"

    cmd = String["sbatch", "--array=$array_spec"]
    append!(cmd, build_worker_resource_args(cfg, runtime_profile))
    extra_args = strip(String(get(cfg, :sbatch_extra_args, "")))
    if !isempty(extra_args)
        append!(cmd, split(extra_args))
    end
    append!(cmd, [script, config_path, iter_dir])

    out = readchomp(Cmd(cmd))
    m = match(r"Submitted batch job (\d+)", out)
    m === nothing && error("Could not parse sbatch output: $out")

    return String(m.captures[1]), out, max_workers
end

function stop_requested(stop_files::Vector{String})
    return any(isfile(f) for f in stop_files)
end

function format_seconds(sec::Real)
    total = max(0, round(Int, sec))
    h = total ÷ 3600
    m = (total % 3600) ÷ 60
    s = total % 60
    return lpad(string(h), 2, '0') * ":" * lpad(string(m), 2, '0') * ":" * lpad(string(s), 2, '0')
end

function count_worker_result_files(results_dir::String)
    isdir(results_dir) || return 0
    return count(f -> endswith(f, ".csv") && occursin("worker_results_", f), readdir(results_dir))
end

function wait_for_job(
    jobid::String,
    poll_seconds::Int,
    results_dir::String,
    expected_chunks::Int,
    eta_per_chunk_seconds::Float64,
    max_workers::Int,
    stop_files::Vector{String},
)
    t0 = time()

    while true
        if stop_requested(stop_files)
            DynamicIO.logmsg("Stop requested while waiting; cancelling worker array job $jobid")
            try
                run(`scancel $jobid`)
            catch err
                DynamicIO.logmsg("scancel failed for $jobid: $err")
            end
            return :stopped
        end

        q = readchomp(`squeue -h -j $jobid -o %T`)
        states = isempty(strip(q)) ? String[] : split(strip(q), '\n')
        n_running = count(==("RUNNING"), states)
        n_pending = count(==("PENDING"), states)
        n_done = count_worker_result_files(results_dir)
        elapsed = time() - t0

        msg = "Worker job $jobid progress: done=$n_done/$expected_chunks running=$n_running pending=$n_pending elapsed=$(format_seconds(elapsed))"
        if eta_per_chunk_seconds > 0
            remaining = max(expected_chunks - n_done, 0)
            eta_seconds = remaining * eta_per_chunk_seconds / max(max_workers, 1)
            msg *= " eta≈$(format_seconds(eta_seconds))"
        end
        DynamicIO.logmsg(msg)

        isempty(states) && return :done
        sleep(poll_seconds)
    end
end

function parse_elapsed_seconds(lines)
    vals = Int[]
    for ln in lines
        s = strip(ln)
        isempty(s) && continue
        try
            v = parse(Int, s)
            v > 0 && push!(vals, v)
        catch
        end
    end
    return vals
end

function estimate_chunk_runtime_seconds(jobid::String)
    try
        out = readchomp(`sacct -n -X -j $jobid --format=ElapsedRaw`)
        lines = split(out, '\n')
        vals = parse_elapsed_seconds(lines)
        isempty(vals) && return 0.0
        return quantile(vals, 0.90)
    catch err
        DynamicIO.logmsg("Could not query sacct runtime stats for $jobid ($err).")
        return 0.0
    end
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

function main()
    Base.exit_on_sigint(false)

    config_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(project_root, "config", "default_config.jl")
    cfg = DynamicIO.read_config(config_path)

    # Accept partition and account from CLI (passed by master.sbatch from Slurm allocation)
    cli_partition = length(ARGS) >= 2 ? strip(String(ARGS[2])) : ""
    cli_account = length(ARGS) >= 3 ? strip(String(ARGS[3])) : ""

    if !isempty(cli_partition)
        cfg[:slurm_partition] = cli_partition
        DynamicIO.logmsg("CLI override: slurm_partition=$cli_partition")
    end
    if !isempty(cli_account)
        cfg[:slurm_account] = cli_account
        DynamicIO.logmsg("CLI override: slurm_account=$cli_account")
    end

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
    active_jobs_file = String(get(cfg, :active_jobs_file, joinpath(run_root, "active_worker_jobs.txt")))
    clear_active_jobs(active_jobs_file)

    DynamicIO.logmsg("Early-stop files: $(join(stop_files, ", "))")
    DynamicIO.logmsg("Active worker job registry: $active_jobs_file")

    runtime_profile = calibrate_worker_resources(cfg, config_path, run_root, stop_files)
    DynamicIO.logmsg(
        "Worker resource profile: source=$(runtime_profile[:source]), time=$(runtime_profile[:time_str]), mem=$(runtime_profile[:mem_mb])M"
    )

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
    )

    gap_history = Int[]
    stop_reason = "completed"
    eta_per_chunk_seconds = Float64(get(cfg, :initial_chunk_eta_seconds, 0.0))

    max_iters = cfg[:max_iters]
    max_points = cfg[:max_points]

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
            chunk_paths = DynamicIO.split_points_into_chunks(points_csv, points_dir, cfg[:points_per_job])
            n_chunks = length(chunk_paths)
            max_workers = effective_max_workers(n_chunks, cfg)

            DynamicIO.logmsg("Iter $iter: cells=$(length(cells)), new_points=$n_candidates, chunks=$n_chunks, max_parallel_workers=$max_workers")
            DynamicIO.logmsg("Iter $iter worker resources: time=$(runtime_profile[:time_str]), mem=$(runtime_profile[:mem_mb])M")

            if cfg[:execution_mode] == :prepare_only
                DynamicIO.logmsg("Prepared iteration only. Submit manually with:")
                part = strip(String(get(cfg, :slurm_partition, "")))
                acct = strip(String(get(cfg, :slurm_account, "")))
                part_arg = isempty(part) ? "" : " --partition=$part"
                acct_arg = isempty(acct) ? "" : " --account=$acct"
                DynamicIO.logmsg(
                    "  sbatch --array=1-$(n_chunks)%$(max_workers)$part_arg$acct_arg --time=$(runtime_profile[:time_str]) --mem=$(runtime_profile[:mem_mb])M $(cfg[:sbatch_script]) $(config_path) $(iter_dir)"
                )
                save_cells_csv(cells, joinpath(iter_dir, "cells_input.csv"))
                stop_reason = "prepare_only"
                break
            elseif cfg[:execution_mode] == :submit_and_wait
                jobid, sbatch_out, _ = submit_array_job(config_path, iter_dir, n_chunks, cfg, runtime_profile)
                write_active_jobs(active_jobs_file, [jobid])

                DynamicIO.logmsg(sbatch_out)
                DynamicIO.logmsg("Waiting for worker job $jobid ...")
                wait_status = wait_for_job(
                    jobid,
                    cfg[:poll_seconds],
                    results_dir,
                    n_chunks,
                    eta_per_chunk_seconds,
                    max_workers,
                    stop_files,
                )

                chunk_runtime_est = estimate_chunk_runtime_seconds(jobid)
                if chunk_runtime_est > 0
                    eta_per_chunk_seconds = chunk_runtime_est
                    DynamicIO.logmsg("Updated chunk runtime estimate (p90): $(format_seconds(eta_per_chunk_seconds))")
                end

                clear_active_jobs(active_jobs_file)
                if wait_status == :stopped
                    stop_reason = "stop_requested_during_wait"
                end
            else
                error("Unknown execution_mode: $(cfg[:execution_mode])")
            end

            worker_results_df = DynamicIO.read_worker_results(results_dir)
            n_result_rows = nrow(worker_results_df)
            n_unique_points = n_result_rows == 0 ? 0 : length(unique(worker_results_df.point_key))
            DynamicIO.logmsg("Iter $iter worker outputs: rows=$n_result_rows, unique_points=$n_unique_points/$n_candidates")

            if n_result_rows == 0
                if startswith(stop_reason, "stop_requested")
                    DynamicIO.logmsg("No worker results collected before stop; preserving state and exiting.")
                    save_checkpoint(run_root, known_results, history_df, gap_history, stop_reason, iter)
                    break
                end
                error("No worker results found in $results_dir")
            end

            require_all = Bool(get(cfg, :require_all_worker_results, true))
            if require_all && n_unique_points < n_candidates && !startswith(stop_reason, "stop_requested")
                error("Incomplete worker results in iter $iter: expected $n_candidates unique points, got $n_unique_points")
            end

            for r in eachrow(worker_results_df)
                known_results[String(r.point_key)] = row_to_result_tuple(r)
            end

            if startswith(stop_reason, "stop_requested")
                DynamicIO.logmsg("Stop requested after partial/complete worker collection; writing checkpoint and ending cleanly.")
                save_checkpoint(run_root, known_results, history_df, gap_history, stop_reason, iter)
                break
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
            ))
            save_checkpoint(run_root, known_results, history_df, gap_history, "in_progress", iter)

            DynamicIO.logmsg("Iter $iter summary: refined_children=$(refine_out.n_refined), rho_max_trivial_gaps=$rho_max_gap_count")

            # Stop criteria: point budget
            if length(known_results) >= max_points
                stop_reason = "max_points"
                DynamicIO.logmsg("Reached max_points=$(max_points). Terminating.")
                break
            end

            # Stop criteria: no further refinement
            if refine_out.n_refined == 0
                stop_reason = "no_refinement"
                DynamicIO.logmsg("No cells refined in this iteration. Terminating.")
                break
            end

            # Stop criteria: stable gap count at rho_max
            stable_iters = cfg[:gap_count_stable_iters]
            gap_tol = cfg[:gap_count_tol]
            if length(gap_history) >= stable_iters
                window = gap_history[end-stable_iters+1:end]
                if (maximum(window) - minimum(window)) <= gap_tol
                    stop_reason = "stable_gap_count"
                    DynamicIO.logmsg("Gap count stable across last $(stable_iters) iterations. Terminating.")
                    break
                end
            end

            # Optional stop criteria: expected maximum number of trivial gaps
            max_expected_gaps = cfg[:max_expected_gaps]
            max_expected_gaps_tol = cfg[:max_expected_gaps_tol]
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
            clear_active_jobs(active_jobs_file)
            rethrow(err)
        end
    end

    clear_active_jobs(active_jobs_file)
    save_checkpoint(run_root, known_results, history_df, gap_history, stop_reason, length(history_df.iter))

    DynamicIO.write_jld2(
        joinpath(run_root, "run_summary.jld2");
        run_root=run_root,
        config_path=config_path,
        stop_reason=stop_reason,
        timestamp=string(Dates.now()),
        n_points_total=length(known_results),
        n_iterations_completed=length(history_df.iter),
        gap_history=gap_history,
        worker_time_str=String(runtime_profile[:time_str]),
        worker_mem_mb=Int(runtime_profile[:mem_mb]),
        worker_profile_source=String(runtime_profile[:source]),
    )

    DynamicIO.logmsg("Search complete. stop_reason=$stop_reason")
    DynamicIO.logmsg("Outputs in: $run_root")
end

main()
