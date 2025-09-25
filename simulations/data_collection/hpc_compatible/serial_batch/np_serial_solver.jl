module npSerialSolver

export np_serial_solver_mu_loop, np_serial_solver_N_loop, np_serial_solver_N_loop_threaded

using LinearAlgebra
using Base.Threads
using DataFrames
using BSON: @save, @load


function np_create_bdg_hamiltonian_generalised(
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

function np_calc_maj_mp(
    eigenvectors::Matrix{Float64} # N.B. the eigen function gives only real eigenvectors since Hermitian is explicit
)
    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2
    middle_state_index = N

    u = eigenvectors[1:N, middle_state_index]
    v = eigenvectors[(N+1):2N, middle_state_index]

    numerator_L = BigFloat(sum(u[i] * conj(v[i]) + v[i] * conj(u[i]) for i in 1:(N÷2)))
    denominator_L = BigFloat(sum(abs2(u[i]) + abs2(v[i]) for i in 1:(N÷2)))
    MP_L = numerator_L / denominator_L

    numerator_R = BigFloat(sum(u[i] * conj(v[i]) + v[i] * conj(u[i]) for i in ((N÷2)+1):N))
    denominator_R = BigFloat(sum(abs2(u[i]) + abs2(v[i]) for i in ((N÷2)+1):N))
    MP_R = numerator_R / denominator_R

    overall_MP = MP_L * MP_R
    # overall_MP = Float64(overall_MP)

    return overall_MP
end



function np_serial_solver_N_loop_threaded(
    N_range::Vector{Int},
    t_n::Vector{Float64},
    mu_range::Vector{Float64},
    Delta::Float64,
    sequence::Vector{Int},
    sequence_name::String,
    filepath::String
)
    truncated_sequences = [Vector(sequence[1:N]) for N in N_range]

    # Collect results in a thread-safe manner by accumulating them separately
    results_chunks = Vector{DataFrame}(undef, Threads.nthreads())

    @Threads.threads for thread_id in 1:Threads.nthreads()
        local_results = DataFrame(
            N = Int[],
            t_n = Vector{Float64}[],
            mu = Float64[],
            Delta = Float64[],
            mp = Float64[]
            # ipr = Float64[],
            # loc_len = Float64[]
        )

        for (i, mu) in enumerate(mu_range)
            if i % Threads.nthreads() != thread_id - 1
                continue  # Each thread only handles a portion of mu_range
            end

            for (j, N) in enumerate(N_range)
                BdG = real(np_create_bdg_hamiltonian_generalised(N, t_n, mu, Delta, truncated_sequences[j]))
                eigenvalues, eigenvectors = LinearAlgebra.eigen(Hermitian(BdG))

                maj_mp = np_calc_maj_mp(eigenvectors)
                # maj_ipr = np_calc_maj_ipr(eigenvectors)
                # loc_length = calc_maj_loc_length(eigenvectors)

                push!(local_results, (
                    N = N,
                    t_n = t_n,
                    mu = mu,
                    Delta = Delta,
                    mp = maj_mp
                    # ipr = maj_ipr,
                    # loc_len = loc_length
                ))
            end
        end

        results_chunks[thread_id] = local_results
    end

    # Combine results from all threads
    results_df = vcat(results_chunks...)

    m1, m2 = minimum(mu_range), maximum(mu_range)
    filename = "$(sequence_name)_Nmin$(minimum(N_range))_Nmax$(maximum(N_range))_tn$(t_n[1])_$(t_n[2])_Delta$(Delta)_mu$(m1)-$(m2).bson"
    full_filepath = "$(filepath)$(filename)"

    @save full_filepath results_df

    return nothing
end



end

# N = 20
# t_n = [1.0, 2.0]
# Delta = 1.0
# mu_range = collect(range(0.0, 4.0, 41))
# sequence = golden_sequence
# seq_name = "GQC"
# precision = 512
# filepath = "/Users/Will/Documents/FINAL_PROJECT/simulations/np_results/serial_test"


# np_serial_solver(N, t_n, mu_range, Delta, sequence, seq_name, precision, filepath)
# println("completed run")