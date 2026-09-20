module SpecFunc

using LinearAlgebra
using Base.Threads
using DataFrames
using ProgressMeter

function cut_and_project_sequence(N::Int, alpha::Float64; gamma::Float64=0.0, map_to::Tuple{Int,Int}=(1,2))
    @assert N >= 1 "N must be >= 1"
    # @assert 0.0 < alpha < 1.0 "alpha expected in (0,1) for standard Sturmian words"
    gamma = gamma - floor(gamma) # wrap into [0,1)
    seq = Vector{Int}(undef, N)
    frac = gamma
    for n in 1:N
        frac += alpha
        if frac >= 1.0
            frac -= 1.0
            # crossing an integer -> difference = 1
            seq[n] = map_to[1]   # "A"
        else
            seq[n] = map_to[2]   # "B"
        end
    end
    return seq
end

function regen_sequence(
    N::Int,
    seq_id::Tuple{Float64, Float64}
)::Vector{Int}
    phi, phason = seq_id
    map_to = (2, 1)
    seq = cut_and_project_sequence(N, phi; gamma=phason, map_to=map_to)

    return seq
end

function create_bdg_hamiltonian(
    N::Int, 
    t_n::Vector{Float64}, 
    mu::Float64, 
    Delta::Float64, 
    sequence::Vector{Int}
)::Matrix{ComplexF64}
    H0 = zeros(ComplexF64, N, N)
    for i in 1:N
        H0[i, i] = -mu
    end
    for i in 1:N-1
        hopping_index = sequence[i]
        hopping_term = t_n[hopping_index]
        H0[i, i+1] = -hopping_term
        H0[i+1, i] = -hopping_term
    end

    Delta_matrix = zeros(ComplexF64, N, N)
    for i in 1:N-1
        Delta_matrix[i, i+1] = Delta
        Delta_matrix[i+1, i] = -Delta
    end

    BdG = zeros(ComplexF64, 2N, 2N)
    BdG[1:N, 1:N] = H0
    BdG[N+1:end, N+1:end] = -conj(H0)
    BdG[1:N, N+1:end] = Delta_matrix
    BdG[N+1:end, 1:N] = Delta_matrix'
    return BdG
end

function calc_spectral_func_for_all!(
    df::DataFrame,
    eta::Float64,
    omega_range::Vector{Float64}
)::DataFrame

    n = nrow(df)
    spectral_funcs_col = Vector{Matrix{Float64}}(undef, n)
    
    # prog = Progress(n; desc="Spectral function", dt=1.0)
    # prog_lock = SpinLock()
    @showprogress for i in 1:n
        row = df[i, :]
        sequence = regen_sequence(row.N, row.sequence_id)
        bdg_matrix = create_bdg_hamiltonian(row.N, row.t_n, row.mu, row.Delta, sequence)
        spectral_values = zeros(row.N, length(omega_range))
        
        for j in 1:length(omega_range)
            omega = omega_range[j]
            G = inv((omega+ im * eta) * I - bdg_matrix)
            
            for site_idx in 1:row.N
                particle_contribution = -imag(G[site_idx, site_idx]) / π
                hole_contribution = -imag(G[row.N + site_idx, row.N + site_idx]) / π
                spectral_values[site_idx, j] = (particle_contribution + hole_contribution) / 2
            end
        end
        spectral_funcs_col[i] = spectral_values
        
        # lock(prog_lock)
        # next!(prog)
        # unlock(prog_lock)
    end

    df[!, :spectral_function] = spectral_funcs_col
    return df
end

function calc_spectral_func_for_all_optimised!(
    df::DataFrame,
    eta::Float64,
    omega_range::Vector{Float64}
)::DataFrame
    n = nrow(df)
    spectral_funcs_col = Vector{Matrix{Float64}}(undef, n)

    BLAS.set_num_threads(1)  # Important: avoid oversubscription

    @showprogress Threads.@threads for i in 1:n
        row = df[i, :]
        sequence = regen_sequence(row.N, row.sequence_id)
        bdg_matrix = create_bdg_hamiltonian(row.N, row.t_n, row.mu, row.Delta, sequence)
        spectral_values = zeros(row.N, length(omega_range))

        for j in 1:length(omega_range)
            omega = omega_range[j]
            G = inv((omega + im * eta) * I - bdg_matrix)

            for site_idx in 1:row.N
                particle_contribution = -imag(G[site_idx, site_idx]) / π
                hole_contribution = -imag(G[row.N + site_idx, row.N + site_idx]) / π
                spectral_values[site_idx, j] = (particle_contribution + hole_contribution) / 2
            end
        end
        spectral_funcs_col[i] = spectral_values
    end

    df[!, :spectral_function] = spectral_funcs_col
    return df
end

end # module SpecFunc