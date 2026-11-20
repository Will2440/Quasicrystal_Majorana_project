using CSV, DataFrames, LinearAlgebra
root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))
include(joinpath(@__DIR__, "resolution_analysis.jl"))
using .ResolutionAnalysis

all_cases = (
    ("rho35_delta005_N500", "rho35_grid1601_k1e4_latest.txt", 3.5, 0.05, 500),
    ("rho5_delta005_N500", "rho5_grid3201_k1e4_latest.txt", 5.0, 0.05, 500),
    ("rho35_delta002_N500", "rho35_delta002_latest.txt", 3.5, 0.02, 500),
    ("rho5_delta005_N300", "rho5_grid3201_k1e4_latest.txt", 5.0, 0.05, 300),
)
ARGS in (String[], ["fine"]) || error("Expected no arguments or fine")
cases = isempty(ARGS) ? all_cases : (all_cases[2], all_cases[3])
output = joinpath(@__DIR__, "runs", "kappa_map_" * get(ENV, "SLURM_JOB_ID", "local"))
mkpath(output)
for (label, pointer, rho, delta, size) in cases
    kappas = if isempty(ARGS)
        10.0 .^ range(-7.0, -3.0; length=17)
    else
        lower = delta == 0.02 ? 5e-5 : 2e-5
        values = exp.(range(log(lower), log(2e-4); length=33))
        sort!(vcat(filter(value -> !isapprox(value, 1e-4; rtol=1e-10), values), 1e-4))
    end
    reference = strip(read(joinpath(@__DIR__, "runs", pointer), String))
    points = CSV.read(joinpath(reference, "slices", "slice_0001", "point_results.csv"), DataFrame)
    points = points[points.stage .== "planned", :]
    sort!(points, :u)
    labels = Int.(points.phase_label)
    runs = phase_runs(labels)
    selected = sort!(unique!(vcat([round(Int, (run.first + run.last) / 2) for run in runs],
        [round(Int, index) for index in range(1, nrow(points); length=33)])))
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(:N => size, :Delta => delta, :rho_start => rho,
        :rho_max => rho, :n_eval_workers => 4, :blas_threads => 1))
    BoundarySearch.validate_config(cfg)
    BLAS.set_num_threads(1)
    ctx = BoundarySearch.make_context(cfg)
    rows = NamedTuple[]
    for kappa in kappas
        cfg[:kappas] = [kappa]
        for offset in 1:4:length(selected)
            batch = selected[offset:min(offset + 3, end)]
            jobs = [Threads.@spawn BoundarySearch.evaluate_point(rho, Float64(points.u[index]), cfg, ctx) for index in batch]
            for (index, job) in zip(batch, jobs)
                result = fetch(job)
                push!(rows, (case=label, N=size, rho=rho, Delta=delta, u=result.u, mu=result.mu,
                    reference_kappa=1e-4, reference_phase=labels[index], kappa=kappa,
                    phase_label=result.phase_label, invariant=result.invariant_mean,
                    localizer_gap=result.gap_min, reference_path=reference))
            end
        end
        CSV.write(joinpath(output, label * ".csv"), DataFrame(rows))
        println("KAPPA ", label, " kappa=", kappa, " points=", length(selected))
        flush(stdout)
    end
end