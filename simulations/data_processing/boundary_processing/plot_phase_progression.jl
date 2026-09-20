#!/usr/bin/env julia
"""
    plot_phase_progression.jl

Animate phase boundary refinement across iterations.

Usage:
    julia plot_phase_progression.jl <run_directory> [--save <output_dir>] [--interval <N>]

Creates individual plots for each iteration (or every N iterations) showing boundary evolution.
Useful for understanding how the algorithm refines the boundary over time.
"""

using CSV
using DataFrames
using Plots
using Statistics
using Glob

const GAMMA = 2 / (1 + sqrt(5))
const T1_FIXED = 1.0
const T2_FIXED = 1.0
const SWEEP_MODE = :vary_t2

function parse_point_key(key::String)
    parts = split(key, "|")
    return parse(Float64, parts[1]), parse(Float64, parts[2])
end

function calculate_mu_critical(rho::Float64, gamma::Float64)
    if SWEEP_MODE == :vary_t2
        t1 = T1_FIXED
        t2 = rho * t1
    elseif SWEEP_MODE == :vary_t1
        t2 = T2_FIXED
        t1 = t2 / rho
    else
        error("Unknown sweep_mode: $SWEEP_MODE")
    end
    tbar = (1.0 - gamma) * t1 + gamma * t2
    bandwidth = 4.0 * tbar
    return 0.5 * bandwidth
end

function load_iteration_data(run_dir::String, iter::Int)
    """Load data accumulated up to and including specified iteration."""
    
    results_file = joinpath(run_dir, "point_results_compact.csv")
    if !isfile(results_file)
        return nothing
    end
    
    df = CSV.read(results_file, DataFrame)
    
    rhos = Float64[]
    mus = Float64[]
    phases = Int[]
    
    for row in eachrow(df)
        rho, u = parse_point_key(row.point_key)
        mu_c = calculate_mu_critical(rho, GAMMA)
        mu = u * mu_c
        push!(rhos, rho)
        push!(mus, mu)
        push!(phases, row.phase_label)
    end
    
    return (rhos=rhos, mus=mus, phases=phases)
end

function create_plot(data::NamedTuple, iter::Int, n_points::Int)
    """Create phase boundary plot for specific iteration."""
    
    rhos = data.rhos
    mus = data.mus
    phases = data.phases
    
    topo_mask = phases .== 1
    trivial_mask = phases .== 0
    
    p = scatter(
        rhos[topo_mask], mus[topo_mask],
        label="Topological",
        color=:red,
        alpha=0.6,
        markersize=4,
        legend=:topright,
        xlabel="ρ = t₂/t₁",
        ylabel="μ (chemical potential)",
        title="Iteration $iter (n=$(n_points) points)",
        size=(900, 600),
        dpi=100
    )
    
    scatter!(p,
        rhos[trivial_mask], mus[trivial_mask],
        label="Trivial",
        color=:blue,
        alpha=0.6,
        markersize=4
    )
    
    return p
end

function get_available_iterations(run_dir::String)
    """Find all iteration directories in run."""
    iter_dirs = glob("iter_*", run_dir)
    iters = Int[]
    for dir in iter_dirs
        m = match(r"iter_(\d+)", basename(dir))
        if m !== nothing
            push!(iters, parse(Int, m.captures[1]))
        end
    end
    return sort(iters)
end

function main()
    if length(ARGS) < 1
        error("Usage: julia plot_phase_progression.jl <run_directory> [--save <output_dir>] [--interval <N>]")
    end
    
    run_dir = ARGS[1]
    if !isdir(run_dir)
        error("Run directory not found: $run_dir")
    end
    
    save_dir = nothing
    interval = 1
    
    i = 2
    while i <= length(ARGS)
        if ARGS[i] == "--save"
            i += 1
            save_dir = ARGS[i]
            mkpath(save_dir)
        elseif ARGS[i] == "--interval"
            i += 1
            interval = parse(Int, ARGS[i])
        end
        i += 1
    end
    
    # Get available iterations
    iters = get_available_iterations(run_dir)
    println("[plot_phase_progression] Found iterations: $iters")
    
    # Generate plots
    results_file = joinpath(run_dir, "point_results_compact.csv")
    initial_df = CSV.read(results_file, DataFrame)
    n_total = nrow(initial_df)
    
    for iter in iters
        if mod(iter - 1, interval) != 0 && iter != iters[end]
            continue
        end
        
        println("[plot_phase_progression] Processing iteration $iter...")
        data = load_iteration_data(run_dir, iter)
        
        if data !== nothing
            p = create_plot(data, iter, length(data.rhos))
            
            if save_dir !== nothing
                out_file = joinpath(save_dir, "phase_iter_$(lpad(iter, 3, '0')).png")
                println("  → Saving to $out_file")
                savefig(p, out_file)
            else
                display(p)
            end
        end
    end
    
    println("[plot_phase_progression] Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
