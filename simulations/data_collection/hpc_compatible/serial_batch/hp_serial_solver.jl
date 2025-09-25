module hpSerialSolver

export hp_serial_solver_mu_loop, hp_serial_solver_N_loop, hp_serial_solver_N_loop_threaded

using GenericLinearAlgebra
using LinearAlgebra
using Base.Threads
using DataFrames
using BSON: @save, @load
# using ProgressMeter


# function golden_ratio_sequence_gen(N::Int)
#     sequence = "A"
#     while length(sequence) < N
#         sequence = replace(sequence, "A" => "AB", "B" => "A")
#     end
#     # println(sequence)
    
#     number_sequence = [ch == 'A' ? 1 : 2 for ch in sequence]
    
#     return number_sequence
# end

# N_gold = 1000
# golden_sequence = golden_ratio_sequence_gen(N_gold)
# # println("Golden ratio sequence: ", golden_sequence)
# println("Sequence length: ", length(golden_sequence))


function hp_create_bdg_hamiltonian_generalised(
    N::Int, 
    t_n::Vector,
    mu::Float64, 
    Delta::Float64,
    sequence::Vector,
    precision::Int
)

    setprecision(BigFloat, precision)
    
    t_n = BigFloat.(t_n)
    mu = BigFloat(mu)
    Delta = BigFloat(Delta)
    
    H0 = zeros(Complex{BigFloat}, N, N)
    for i in 1:N-1
        H0[i, i] = -mu
        hopping_index = sequence[i]
        H0[i, i+1] = -t_n[hopping_index]
        H0[i+1, i] = -t_n[hopping_index]
    end
    H0[N, N] = -mu

    Delta_matrix = zeros(Complex{BigFloat}, N, N)
    for i in 1:N-1
        Delta_matrix[i, i+1] = Delta
        Delta_matrix[i+1, i] = -Delta
    end

    BdG = zeros(Complex{BigFloat}, 2N, 2N)
    BdG[1:N, 1:N] = H0
    BdG[N+1:end, N+1:end] = -conj(H0)
    BdG[1:N, N+1:end] = Delta_matrix
    BdG[N+1:end, 1:N] = Delta_matrix'

    return BdG
end

function hp_calc_maj_mp(
    eigenvectors::Matrix{BigFloat} # N.B. the eigen function gives only real eigenvectors since Hermitian is explicit
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
    overall_MP = Float64(overall_MP)

    return overall_MP
end

function hp_calc_mp_sitewise(
    eigenvectors::Matrix{BigFloat}
)
    num_cols = size(eigenvectors, 2)
    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2

    mp_sitewise = zeros(BigFloat, N, num_cols)
    nums = zeros(ComplexF64, N, num_cols)
    dens = zeros(ComplexF64, N, num_cols)
    
    for n in 1:num_cols
        u = eigenvectors[1:N, n]
        v = eigenvectors[(N+1):2N, n]

        for i in 1:N
            numerator = BigFloat(u[i] * conj(v[i]) + v[i] * conj(u[i]))
            denominator = BigFloat(abs2(u[i]) + abs2(v[i]))
            mp_sitewise[i, n] = BigFloat(numerator / denominator)
            nums[i] = numerator
            dens[i] = denominator
        end
    end
end

function hp_calc_maj_ipr(
    eigenvectors::Matrix{BigFloat}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]
    numerator = sum(site_probs .^ 2)
    denominator = sum(site_probs)^2

    ipr_value = numerator / denominator
    ipr_value = Float64(ipr_value)

    return ipr_value
end

function hp_calc_mbs_energy_gap(
    eigenvalues::Vector{BigFloat}
)
   
    middle_index = ceil(Int, length(eigenvalues) / 2)
    mbs_energy_gap = eigenvalues[middle_index] - eigenvalues[middle_index - 1]
    
    return mbs_energy_gap
end

function calc_maj_loc_length(
    eigenvectors::Matrix{BigFloat}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]

    i0 = argmax(site_probs)

    numerator = sum(site_probs .* abs.(collect(1:num_sites) .- i0))
    denominator = sum(site_probs)

    localisation_length = numerator / denominator
    return Float64(localisation_length)
end


function hp_serial_solver_mu_loop(
    N::Int,
    t_n::Vector{Float64},
    mu_range::Vector{Float64},
    Delta::Float64,
    sequence::Vector{Int},
    sequence_name::String,
    precision::Int,
    filepath::String
)

    # Initialize DataFrame for results
    results_df = DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        mp = Float64[],
        ipr = Float64[],
        mbs_gap = Float64[]
    )

    for mu in mu_range
        truncated_sequence = Vector(sequence[1:N])
        BdG = hp_create_bdg_hamiltonian_generalised(N, t_n, mu, Delta, truncated_sequence, precision)
        eigenvalues, eigenvectors = GenericLinearAlgebra.eigen(Hermitian(BdG))

        # eigenvalues = BigFloat.(eigenvalues)
        # eigenvectors = BigFloat.(eigenvectors)

        eigenvalues = real.(eigenvalues)
        eigenvectors = real.(eigenvectors)

        mp = hp_calc_maj_mp(eigenvectors)
        maj_ipr = hp_calc_maj_ipr(eigenvectors)
        gap = hp_calc_mbs_energy_gap(eigenvalues)
        # println("pushing mu value: $mu")

        push!(results_df, (
            N = N,
            t_n = t_n,
            mu = mu,
            Delta = Delta,
            sequence_name = sequence_name,
            mp = mp,
            ipr = maj_ipr,
            mbs_gap = gap
        ))
    end

    m1 = minimum(mu_range)
    m2 = maximum(mu_range)
    l = length(mu_range)

    t1 = t_n[1]
    t2 = t_n[2]

    filename = "$(sequence_name)_N$(N)_tn$(t1)_$(t2)_Delta$(Delta)_mu$(m1)-$(m2)_$(l).bson"
    filepath = "$(filepath)$(filename)"

    @save filepath results_df

    return nothing
end

function hp_serial_solver_N_loop(
    N_range::Vector{Int},
    t_n::Vector{Float64},
    mu_range::Vector{Float64},
    Delta::Float64,
    sequence::Vector{Int},
    sequence_name::String,
    precision::Int,
    filepath::String
)
    for mu in mu_range
        results_df = DataFrame(
            N = Int[],
            t_n = Vector{Float64}[],
            mu = Float64[],
            Delta = Float64[],
            sequence_name = String[],
            # mp = Float64[],
            ipr = Float64[],
            # mbs_gap = Float64[],
            loc_len = Float64[]
        )

        for N in N_range
            truncated_sequence = Vector(sequence[1:N])
            BdG = real(hp_create_bdg_hamiltonian_generalised(N, t_n, mu, Delta, truncated_sequence, precision))
            eigenvalues, eigenvectors = GenericLinearAlgebra.eigen(Hermitian(BdG))

            # eigenvalues = BigFloat.(eigenvalues)
            # eigenvectors = BigFloat.(eigenvectors)

            # mp = hp_calc_maj_mp(eigenvectors)
            maj_ipr = hp_calc_maj_ipr(eigenvectors)
            # gap = hp_calc_mbs_energy_gap(eigenvalues)
            loc_length = calc_maj_loc_length(eigenvectors)

            push!(results_df, (
                N = N,
                t_n = t_n,
                mu = mu,
                Delta = Delta,
                sequence_name = sequence_name,
                # mp = mp,
                ipr = maj_ipr,
                # mbs_gap = gap,
                loc_len = loc_length
            ))
        end

        m1 = minimum(mu_range)
        m2 = maximum(mu_range)
        l = length(mu_range)

        t1 = t_n[1]
        t2 = t_n[2]

        # Filename for this specific mu iteration
        filename = "$(sequence_name)_Nmin$(minimum(N_range))_Nmax$(maximum(N_range))_tn$(t1)_$(t2)_Delta$(Delta)_range$(m1)-$(m2)_$(l)_mu$(mu).bson"
        full_filepath = "$(filepath)$(filename)"

        @save full_filepath results_df
    end

    return nothing
end

function hp_serial_solver_N_loop_threaded(
    N_range::Vector{Int},
    t_n::Vector{Float64},
    mu_range::Vector{Float64},
    Delta::Float64,
    sequence::Vector{Int},
    sequence_name::String,
    precision::Int,
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
            mp = Float64[],
            ipr = Float64[],
            loc_len = Float64[]
        )

        for (i, mu) in enumerate(mu_range)
            if i % Threads.nthreads() != thread_id - 1
                continue  # Each thread only handles a portion of mu_range
            end

            for (j, N) in enumerate(N_range)
                BdG = real(hp_create_bdg_hamiltonian_generalised(N, t_n, mu, Delta, truncated_sequences[j], precision))
                eigenvalues, eigenvectors = GenericLinearAlgebra.eigen(Hermitian(BdG))

                maj_mp = hp_calc_maj_mp(eigenvectors)
                maj_ipr = hp_calc_maj_ipr(eigenvectors)
                loc_length = calc_maj_loc_length(eigenvectors)

                push!(local_results, (
                    N = N,
                    t_n = t_n,
                    mu = mu,
                    Delta = Delta,
                    mp = maj_mp,
                    ipr = maj_ipr,
                    loc_len = loc_length
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
# filepath = "/Users/Will/Documents/FINAL_PROJECT/simulations/hp_results/serial_test"


# hp_serial_solver(N, t_n, mu_range, Delta, sequence, seq_name, precision, filepath)
# println("completed run")