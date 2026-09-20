#!/usr/bin/env julia

using CSV
using DataFrames
using Dates

project_root = @__DIR__
include(joinpath(project_root, "core", "physics.jl"))
include(joinpath(project_root, "core", "sequences.jl"))
include(joinpath(project_root, "core", "bandwidth.jl"))
include(joinpath(project_root, "core", "io_utils.jl"))

using .DynamicPhysics
using .DynamicSequences
using .DynamicBandwidth
using .DynamicIO

function logmsg(msg::AbstractString)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    println("[$ts] [worker] $msg")
end

function build_x0s(N::Int, alpha::Float64, npoints::Int)
    return collect(range(N / 2 - alpha * N, N / 2 + alpha * N, length=npoints))
end

function main()
    if length(ARGS) < 3
        error("Usage: julia worker.jl <config_path> <points_chunk.csv> <worker_results.csv>")
    end

    config_path = ARGS[1]
    points_path = ARGS[2]
    out_path = ARGS[3]

    cfg = DynamicIO.read_config(config_path)
    worker_log_every = Int(get(cfg, :worker_log_every, 10))

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

    N = cfg[:N]
    Delta = cfg[:Delta]
    x0s = build_x0s(N, cfg[:x0_alpha], cfg[:x0_count])
    kappas = cfg[:kappas]

    points_df = CSV.read(points_path, DataFrame)
    n_total = nrow(points_df)
    rows = NamedTuple[]

    logmsg("Starting chunk: $points_path")
    logmsg("Total points in chunk: $n_total")
    t0 = time()

    for (i, r) in enumerate(eachrow(points_df))
        rho = Float64(r.rho)
        u = Float64(r.u)

        t1, t2 = DynamicBandwidth.tpair_from_rho(
            rho;
            sweep_mode=cfg[:sweep_mode],
            t1_fixed=cfg[:t1_fixed],
            t2_fixed=cfg[:t2_fixed],
        )

        mu_c_est = DynamicBandwidth.mu_critical_estimate(gamma, t1, t2)
        mu = u * mu_c_est

        eval_out = DynamicPhysics.evaluate_specloc_point(
            N,
            [t1, t2],
            mu,
            Delta,
            sequence,
            x0s,
            kappas;
            n_lowest_evals=cfg[:n_lowest_evals],
            kk_tol=cfg[:kk_tol],
            kk_maxiter=cfg[:kk_maxiter],
            topological_threshold=cfg[:topological_threshold],
            topological_fraction=cfg[:topological_fraction],
        )

        push!(rows, (
            point_key=String(r.point_key),
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
        ))

        if (i % worker_log_every == 0) || (i == n_total)
            elapsed = time() - t0
            avg = elapsed / max(i, 1)
            remaining = n_total - i
            eta = remaining * avg
            logmsg("Progress: $i/$n_total points, elapsed=$(round(elapsed, digits=1))s, eta=$(round(eta, digits=1))s")
        end
    end

    DynamicIO.write_worker_results(rows, out_path)

    elapsed_total = time() - t0
    logmsg("Wrote worker results: $out_path")
    logmsg("Sequence slope/phason used: $(gamma), $(phason)")
    logmsg("Chunk complete in $(round(elapsed_total, digits=2))s")
end

main()
