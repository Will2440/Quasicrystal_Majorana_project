module DynamicPhysics

using LinearAlgebra
using SparseArrays
using KrylovKit
using Statistics

export evaluate_specloc_point

function np_create_bdg_hamiltonian(
    N::Int,
    t_n::Vector,
    mu::Float64,
    Delta::Float64,
    sequence::Vector
)
    H0 = zeros(Complex{Float64}, N, N)
    for i in 1:N-1
        H0[i, i] = -mu
        hopping_index = sequence[i]
        H0[i, i+1] = -t_n[hopping_index]
        H0[i+1, i] = -t_n[hopping_index]
    end
    H0[N, N] = -mu

    Delta_matrix = zeros(Complex{Float64}, N, N)
    for i in 1:N-1
        Delta_matrix[i, i+1] = Delta
        Delta_matrix[i+1, i] = -Delta
    end

    BdG = zeros(Complex{Float64}, 2N, 2N)
    BdG[1:N, 1:N] = H0
    BdG[N+1:end, N+1:end] = -conj(H0)
    BdG[1:N, N+1:end] = Delta_matrix
    BdG[N+1:end, 1:N] = Delta_matrix'

    return BdG
end

function build_operators_block_basis(N::Int)
    dim = 2 * N

    xvec = zeros(Float64, dim)
    for i in 1:N
        xvec[i] = i
        xvec[N + i] = i
    end
    X = Diagonal(xvec)

    Gamma = spzeros(Float64, dim, dim)
    for i in 1:N
        Gamma[i, N + i] = 1.0
        Gamma[N + i, i] = 1.0
    end

    return X, Gamma
end

function np_calc_specloc_signature_sparse_ldlt(
    L::SparseMatrixCSC{Float64, Int};
    max_shift_tries::Int=6,
    base_shift::Float64=1e-12
)::Int
    M = Symmetric(L)
    shift = 0.0

    for _ in 0:max_shift_tries
        try
            F = ldlt(M; shift=shift)
            LDmat = sparse(F.LD)
            Dvals = diag(LDmat)
            n_pos = count(x -> x > 0.0, Dvals)
            n_neg = count(x -> x < 0.0, Dvals)
            return n_pos - n_neg
        catch err
            estr = string(typeof(err))
            if err isa LinearAlgebra.ZeroPivotException || occursin("CHOLMOD", estr) || occursin("Pivot", estr)
                shift = shift == 0.0 ? base_shift : shift * 10.0
                continue
            else
                rethrow(err)
            end
        end
    end

    vals = eigvals(Symmetric(Matrix(L)))
    return count(vals .> 0.0) - count(vals .< 0.0)
end

function np_calc_specloc_low_lying_shift_invert(
    L::SparseMatrixCSC{Float64, Int};
    n_eigenpairs::Int=4,
    kk_tol::Real=1e-8,
    kk_maxiter::Int=300,
    regularization::Float64=1e-10
)::Vector{Float64}
    n = size(L, 1)

    try
        F = lu(L)
        L_inv_action(x) = F \ x
        evals_inv, _, _ = eigsolve(
            L_inv_action,
            rand(Float64, n),
            n_eigenpairs,
            :LM;
            ishermitian=true,
            tol=kk_tol,
            maxiter=kk_maxiter
        )

        low_vals = real.(1.0 ./ evals_inv)
        sort!(low_vals, by=abs)
        return low_vals[1:min(n_eigenpairs, length(low_vals))]
    catch
        L_reg = L + regularization * I(n)
        try
            F = lu(L_reg)
            L_inv_action(x) = F \ x
            evals_inv, _, _ = eigsolve(
                L_inv_action,
                rand(Float64, n),
                n_eigenpairs,
                :LM;
                ishermitian=true,
                tol=kk_tol,
                maxiter=kk_maxiter
            )

            low_vals = real.(1.0 ./ evals_inv) .- regularization
            sort!(low_vals, by=abs)
            return low_vals[1:min(n_eigenpairs, length(low_vals))]
        catch
            vals = eigvals(Symmetric(Matrix(L)))
            sort!(vals, by=abs)
            return vals[1:min(n_eigenpairs, length(vals))]
        end
    end
end

function evaluate_specloc_point(
    N::Int,
    t_n::Vector{Float64},
    mu::Float64,
    Delta::Float64,
    sequence::Vector{Int},
    x0s::Vector{Float64},
    kappas::Vector{Float64};
    n_lowest_evals::Int=4,
    kk_tol::Real=1e-8,
    kk_maxiter::Int=300,
    topological_threshold::Float64=0.5,
    topological_fraction::Float64=0.5
)
    truncated_sequence = Vector(sequence[1:N])
    H_BdG = np_create_bdg_hamiltonian(N, t_n, mu, Delta, truncated_sequence)

    X, Gamma = build_operators_block_basis(N)

    invariants = Float64[]
    localiser_gaps = Float64[]

    for x0 in x0s
        X_diff = sparse(Diagonal(X.diag .- x0))
        for kappa in kappas
            L = sparse(real.(H_BdG)) + kappa .* (Gamma * X_diff)
            low_lying = np_calc_specloc_low_lying_shift_invert(
                L;
                n_eigenpairs=n_lowest_evals,
                kk_tol=kk_tol,
                kk_maxiter=kk_maxiter
            )
            specloc_gap = minimum(abs.(low_lying))

            sig = np_calc_specloc_signature_sparse_ldlt(L)
            specloc_invariant = sig / 2.0

            push!(invariants, specloc_invariant)
            push!(localiser_gaps, specloc_gap)
        end
    end

    topo_votes = count(inv -> inv > topological_threshold, invariants)
    topo_frac = topo_votes / max(length(invariants), 1)
    phase_label = topo_frac >= topological_fraction ? 1 : 0

    return (
        phase_label = phase_label,
        topo_fraction = topo_frac,
        invariant_mean = mean(invariants),
        invariant_min = minimum(invariants),
        invariant_max = maximum(invariants),
        gap_min = minimum(localiser_gaps),
        gap_mean = mean(localiser_gaps),
        gap_max = maximum(localiser_gaps),
        n_samples = length(invariants),
    )
end

end # module
