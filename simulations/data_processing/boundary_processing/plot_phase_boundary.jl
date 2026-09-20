#!/usr/bin/env julia
"""
    plot_phase_boundary.jl

Plot the topological phase boundary from dynamic search results on a mu vs rho diagram.

Usage:
    julia plot_phase_boundary.jl <run_directory> [--iter <iteration>] [--save <output_path>] [--show]

Arguments:
    run_directory     Path to the dynamic search run output directory
    --iter N          Use data up to iteration N (default: use all available)
    --save PATH       Save plot to file (default: show only)
    --show            Display plot interactively (default if no --save)

Example:
    julia plot_phase_boundary.jl ../data_collection/hpc_compatible/dynamic_search_singlejob/runs/dynamic_mb_search_singlejob_20260920_171150 --iter 5 --save phase_diagram.png
"""

using CSV
using DataFrames
using Plots
using Statistics

# Configuration parameters (match the run configuration)
const GAMMA = 2 / (1 + sqrt(5))  # Golden ratio conjugate ≈ 0.618
const T1_FIXED = 1.0
const T2_FIXED = 1.0
const SWEEP_MODE = :vary_t2

function parse_point_key(key::String)
    """Parse point_key 'rho|u' into (rho, u)"""
    parts = split(key, "|")
    return parse(Float64, parts[1]), parse(Float64, parts[2])
end

function calculate_mu_critical(rho::Float64, gamma::Float64)
    """Calculate mu_c from rho, given sweep mode and fixed t values."""
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

function load_data(run_dir::String, iter_limit::Union{Int, Nothing}=nothing)
    """Load point results from run directory, optionally limiting to specified iteration."""
    
    results_file = joinpath(run_dir, "point_results_compact.csv")
    if !isfile(results_file)
        error("Results file not found: $results_file")
    end
    
    df = CSV.read(results_file, DataFrame)
    
    # Parse point keys and calculate mu values
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

function create_plot(data::NamedTuple, title::String="")
    """Create phase boundary plot."""
    
    rhos = data.rhos
    mus = data.mus
    phases = data.phases
    
    # Separate topological (phase=1) and trivial (phase=0) points
    topo_mask = phases .== 1
    trivial_mask = phases .== 0
    
    p = scatter(
        rhos[topo_mask], mus[topo_mask],
        label="Topological (ν=1)",
        color=:red,
        alpha=0.6,
        markersize=3,
        legend=:topright,
        xlabel="ρ = t₂/t₁",
        ylabel="μ (chemical potential)",
        title=title,
        size=(900, 600),
        dpi=150
    )
    
    scatter!(p,
        rhos[trivial_mask], mus[trivial_mask],
        label="Trivial (ν=0)",
        color=:blue,
        alpha=0.6,
        markersize=3
    )
    
    return p
end

function main()
    if length(ARGS) < 1
        println(__doc__)
        error("Run directory argument required")
    end
    
    run_dir = ARGS[1]
    if !isdir(run_dir)
        error("Run directory not found: $run_dir")
    end
    
    # Parse optional arguments
    iter_limit = nothing
    save_path = nothing
    show_plot = true
    
    i = 2
    while i <= length(ARGS)
        if ARGS[i] == "--iter"
            i += 1
            iter_limit = parse(Int, ARGS[i])
        elseif ARGS[i] == "--save"
            i += 1
            save_path = ARGS[i]
            show_plot = false
        elseif ARGS[i] == "--show"
            show_plot = true
        end
        i += 1
    end
    
    # Load data
    println("[plot_phase_boundary] Loading results from: $run_dir")
    data = load_data(run_dir, iter_limit)
    println("[plot_phase_boundary] Loaded $(length(data.rhos)) points")
    
    # Create plot title
    title_str = "Topological Phase Diagram (ρ vs μ)"
    if iter_limit !== nothing
        title_str *= " - Up to iteration $iter_limit"
    end
    
    # Create plot
    println("[plot_phase_boundary] Creating plot...")
    p = create_plot(data, title_str)
    
    # Save or show
    if save_path !== nothing
        println("[plot_phase_boundary] Saving plot to: $save_path")
        savefig(p, save_path)
    end
    
    if show_plot
        println("[plot_phase_boundary] Displaying plot...")
        display(p)
    end
    
    println("[plot_phase_boundary] Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
