#!/usr/bin/env julia
"""
    analyze_boundary.jl

Analyze boundary statistics and quality metrics from a dynamic search run.

Usage:
    julia analyze_boundary.jl <run_directory> [--report FILE]

Outputs:
    - Total points and phase distribution
    - Boundary sharpness metrics
    - Point density by phase
    - Suggested confidence regions
"""

using CSV
using DataFrames
using Statistics

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

function analyze_run(run_dir::String)
    """Compute boundary statistics from run."""
    
    results_file = joinpath(run_dir, "point_results_compact.csv")
    if !isfile(results_file)
        error("Results file not found: $results_file")
    end
    
    df = CSV.read(results_file, DataFrame)
    
    # Parse and compute
    rhos = Float64[]
    mus = Float64[]
    phases = Int[]
    gaps = Float64[]
    topos = Float64[]
    
    for row in eachrow(df)
        rho, u = parse_point_key(row.point_key)
        mu_c = calculate_mu_critical(rho, GAMMA)
        mu = u * mu_c
        
        push!(rhos, rho)
        push!(mus, mu)
        push!(phases, row.phase_label)
        push!(gaps, row.gap_min)
        push!(topos, row.topo_fraction)
    end
    
    n_total = length(phases)
    n_topo = count(==(1), phases)
    n_trivial = count(==(0), phases)
    
    # Boundary analysis: find points with mixed classification
    boundary_mask = (topos .> 0.0) .& (topos .< 1.0)
    n_boundary = count(boundary_mask)
    
    # Gap statistics
    mean_gap = mean(gaps)
    min_gap = minimum(gaps)
    median_gap = median(gaps)
    
    # Spatial stats
    rho_range = (minimum(rhos), maximum(rhos))
    mu_range = (minimum(mus), maximum(mus))
    
    # Point density (points per unit area)
    rho_span = rho_range[2] - rho_range[1]
    mu_span = mu_range[2] - mu_range[1]
    density = n_total / (rho_span * mu_span)
    
    return (
        n_total=n_total,
        n_topological=n_topo,
        n_trivial=n_trivial,
        n_boundary_ambiguous=n_boundary,
        mean_gap=mean_gap,
        min_gap=min_gap,
        median_gap=median_gap,
        rho_range=rho_range,
        mu_range=mu_range,
        point_density=density,
        topo_fraction=n_topo / n_total,
        boundary_fraction=n_boundary / n_total
    )
end

function load_run_history(run_dir::String)
    """Load search_history.csv for iteration statistics."""
    
    history_file = joinpath(run_dir, "search_history.csv")
    if !isfile(history_file)
        return nothing
    end
    
    df = CSV.read(history_file, DataFrame)
    return df
end

function print_report(stats::NamedTuple, history::Union{DataFrame, Nothing})
    """Pretty-print analysis report."""
    
    println("\n" * "="^70)
    println("TOPOLOGICAL PHASE BOUNDARY ANALYSIS")
    println("="^70)
    
    println("\n📊 OVERALL STATISTICS:")
    println("  Total points evaluated:        $(stats.n_total)")
    println("  Topological phase (ν=1):       $(stats.n_topological) ($(round(stats.topo_fraction*100; digits=1))%)")
    println("  Trivial phase (ν=0):           $(stats.n_trivial) ($(round((1-stats.topo_fraction)*100; digits=1))%)")
    println("  Ambiguous/mixed points:        $(stats.n_boundary_ambiguous) ($(round(stats.boundary_fraction*100; digits=1))%)")
    
    println("\n🎯 BOUNDARY QUALITY:")
    println("  Mean gap (boundary detection):  $(round(stats.mean_gap; sigdigits=4))")
    println("  Median gap:                     $(round(stats.median_gap; sigdigits=4))")
    println("  Minimum gap:                    $(round(stats.min_gap; sigdigits=4))")
    if stats.min_gap < 1e-4
        println("    ⚠️  Very small gaps detected - high refinement resolution")
    end
    
    println("\n📐 PARAMETER SPACE COVERAGE:")
    println("  ρ range [t₂/t₁]:                [$(round(stats.rho_range[1]; digits=3)), $(round(stats.rho_range[2]; digits=3))]")
    println("  μ range:                        [$(round(stats.mu_range[1]; digits=3)), $(round(stats.mu_range[2]; digits=3))]")
    println("  Point density:                  $(round(stats.point_density; sigdigits=3)) points/unit²")
    
    println("\n📈 ITERATION HISTORY:")
    if history !== nothing
        n_iters = nrow(history)
        println("  Total iterations:               $n_iters")
        println("  Final cell count:               $(history[end, :n_cells])")
        println("  Points per iteration (avg):     $(round(mean(history.n_candidates); digits=0))")
        
        # Find when boundary was found
        if n_iters > 1
            first_refined = findfirst(>(0), history.n_refined_children)
            if first_refined !== nothing
                println("  Boundary found at iteration:    $first_refined")
            end
        end
    else
        println("  (search_history.csv not found)")
    end
    
    println("\n" * "="^70 * "\n")
    
    return nothing
end

function main()
    if length(ARGS) < 1
        error("Usage: julia analyze_boundary.jl <run_directory> [--report FILE]")
    end
    
    run_dir = ARGS[1]
    if !isdir(run_dir)
        error("Run directory not found: $run_dir")
    end
    
    report_file = nothing
    if length(ARGS) >= 3 && ARGS[2] == "--report"
        report_file = ARGS[3]
    end
    
    # Run analysis
    println("Analyzing boundary from: $run_dir")
    stats = analyze_run(run_dir)
    history = load_run_history(run_dir)
    
    # Generate report
    if report_file !== nothing
        open(report_file, "w") do io
            redirect_stdout(io) do
                print_report(stats, history)
            end
        end
        println("Report saved to: $report_file")
    else
        print_report(stats, history)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
