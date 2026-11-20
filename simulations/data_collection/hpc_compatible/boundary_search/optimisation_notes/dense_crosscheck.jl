using CSV, DataFrames, LinearAlgebra, SparseArrays
root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(root, "master.jl"))
physics = BoundarySearch.DynamicPhysics
BLAS.set_num_threads(4)
rows = NamedTuple[]
for (rho, delta, u, kappa) in ((5.0, 0.05, 0.064640625, 1e-4),
    (5.0, 0.05, 0.065625, 1e-4), (5.0, 0.05, 0.065625, 1e-5),
    (3.5, 0.02, 0.6601875, 1e-5), (3.5, 0.02, 0.6601875, 4e-5),
    (3.5, 0.02, 0.6601875, 1e-4))
    cfg = include(joinpath(@__DIR__, "calibration_config.jl"))
    merge!(cfg, Dict{Symbol,Any}(:rho_start => rho, :rho_max => rho, :Delta => delta,
        :N => 500, :kappas => [kappa]))
    ctx = BoundarySearch.make_context(cfg)
    params = BoundarySearch.SliceSweep.physical_parameters(rho, cfg, ctx.gamma, ctx.fixed_mu_critical)
    result = BoundarySearch.evaluate_point(rho, u, cfg, ctx)
    hamiltonian = physics.np_create_bdg_hamiltonian(500, [params.t1, params.t2],
        u * params.mu_critical, delta, ctx.sequence)
    position, gamma = physics.build_operators_block_basis(500)
    localizer = sparse(real.(hamiltonian)) + kappa * (gamma * sparse(Diagonal(position.diag .- 250.0)))
    eigenvalues = eigvals(Symmetric(Matrix(localizer)))
    dense_gap = minimum(abs, eigenvalues)
    dense_invariant = (count(>(0), eigenvalues) - count(<(0), eigenvalues)) / 2
    push!(rows, (rho=rho, Delta=delta, u=u, mu=result.mu, kappa=kappa,
        sparse_invariant=result.invariant_mean, dense_invariant=dense_invariant,
        sparse_gap=result.gap_min, dense_gap=dense_gap,
        relative_gap_error=abs(result.gap_min - dense_gap) / dense_gap))
    println("CROSSCHECK ", last(rows))
    flush(stdout)
end
path = joinpath(@__DIR__, "runs", "dense_crosscheck_" * get(ENV, "SLURM_JOB_ID", "local") * ".csv")
CSV.write(path, DataFrame(rows))