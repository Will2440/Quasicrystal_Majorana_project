"""
    file name:   solver.jl
    created:     24/09/2025
    last edited: 30/09/2025

    overview:
        This file generates an importable module containing all of the calculations relevant to the analysing Majoranas in the Kitaev chain.
        Only the solving functions of Sec 3 are exported for use outside of this module.
    
    structure:
        Sec 1:  Constructing the Hamiltonian
                    Contains the functions needed to create the Kitaev chain Hamiltonian with an arbitrary hopping sequence in normal and high precision.
        Sec 2:  Calculations
                    Contains the calculation functions of the majorana polarisation (mp), majorana energy gap size (mbs gap), inverse participation ratio (ipr) and the localisation length all in normal and high precision.
        Sec 3:  Solving Functions
                    [IN PROGRESS] Contains solving functions which leverage the construction and calculations of Secs 1 & 2 with different applications and saving methods.
                    Understanding the different solver types:
                        - :generic
                            This solver type iterates over all parameter combinations in a single disordered for loop, saving results in chunks to avoid memory issues.
                            It is the most flexible but may not be the most efficient for specific parameter sweeps.
                        - :restricted
                            This solver type is designed for sweeping over a 2D parameter space (e.g., mu vs rho) with restrictions applied to limit the parameter combinations based on user-defined rules.
                            It is useful to reduce computational cost when there are large regions of parameter space that are known to be trivial.
                            This solver type is optimized for sweeping over chemical potential (mu) while keeping other parameters fixed.
                            It is useful for studying phase transitions as a function of mu. (N.B. this was initially designed for HPC batch jobs where a second parameter was varied between each job in the batch.)
                        - :N_loop
                            This solver type is designed for varying the system size (N) and mu while keeping other parameters constant.
                            It consists of a nexted for loop where the outer loop iterates over mu values and the inner loop iterates over N values.
                            This is useful for maintaining consistent compeltion times for each step of the outer for loop over a range of varying N.
                            (N.B. this is ready for paralleisation over mu_range on a local machine.)
                        - :N_loop_threaded
                            This solver type is similar to :N_loop but leverages Julia's multithreading capabilities to parallelise the inner loop.
                            (N.B. this was originally designed for HPC where multiple nodes of multiple cores could be accessed, so a balance between speed of each job in the batch and overall number of jobs could be found.)
"""

module LocalSolv
export hp_generic_solver, np_generic_solver, np_mu_rho_restricted_solver, UserOptions

using GenericLinearAlgebra
using LinearAlgebra
# using Arpack
# using SparseArrays
using Base.Threads
using DataFrames
using BSON: @save, @load
using ProgressMeter


struct UserOptions
    calc_mp::Bool
    calc_ipr::Bool
    calc_mbs_energy_gap::Bool
    calc_loc_len::Bool
    calc_precision::Symbol # :hp or :np
    save_evecs::Symbol # :all_np, :all_hp, :maj_np, :maj_hp or :none
    save_evals::Symbol # :all_np, :all_hp, :maj_np, :maj_hp or :none
    solver_type::Symbol # :generic, :mu_loop or :N_loop
end


###########################################################
########### Sec 1: Constructing the Hamiltonian ###########
###########################################################

function hp_create_bdg_hamiltonian(
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

###########################################################
################### Sec 2: Calculations ###################
###########################################################

function hp_calc_maj_mp(
    eigenvectors::Matrix{Complex{BigFloat}} # N.B. the eigen function gives only real eigenvectors since Hermitian is explicit
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

    return overall_MP::Float64
end

function np_calc_maj_mp(
    eigenvectors::Matrix{ComplexF64} # N.B. the eigen function gives only real eigenvectors since Hermitian is explicit
)
    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2
    middle_state_index = N

    u = eigenvectors[1:N, middle_state_index]
    v = eigenvectors[(N+1):2N, middle_state_index]

    numerator_L = Float64(sum(u[i] * conj(v[i]) + v[i] * conj(u[i]) for i in 1:(N÷2)))
    denominator_L = Float64(sum(abs2(u[i]) + abs2(v[i]) for i in 1:(N÷2)))
    MP_L = numerator_L / denominator_L

    numerator_R = Float64(sum(u[i] * conj(v[i]) + v[i] * conj(u[i]) for i in ((N÷2)+1):N))
    denominator_R = Float64(sum(abs2(u[i]) + abs2(v[i]) for i in ((N÷2)+1):N))
    MP_R = numerator_R / denominator_R

    overall_MP = MP_L * MP_R

    return overall_MP::Float64
end

function hp_calc_maj_ipr(
    eigenvectors::Matrix{Complex{BigFloat}}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]
    numerator = sum(site_probs .^ 2)
    denominator = sum(site_probs)^2

    ipr_value = numerator / denominator
    ipr_value = Float64(ipr_value)

    return ipr_value::Float64
end

function np_calc_maj_ipr(
    eigenvectors::Matrix{ComplexF64}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]
    numerator = sum(site_probs .^ 2)
    denominator = sum(site_probs)^2

    ipr_value = numerator / denominator
    ipr_value = Float64(ipr_value)

    return ipr_value::Float64
end

function hp_calc_maj_loc_len(
    eigenvectors::Matrix{Complex{BigFloat}}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]

    i0 = argmax(site_probs)

    numerator = sum(site_probs .* abs.(collect(1:num_sites) .- i0))
    denominator = sum(site_probs)
    localisation_length = Float64(numerator / denominator)
    
    return localisation_length::Float64
end

function np_calc_maj_loc_len(
    eigenvectors::Matrix{ComplexF64}
)
    num_sites = div(size(eigenvectors, 1), 2)
    middle_state_index = num_sites

    psi = eigenvectors[:, middle_state_index]
    site_probs = [abs(psi[j])^2 + abs(psi[num_sites + j])^2 for j in 1:num_sites]

    i0 = argmax(site_probs)

    numerator = sum(site_probs .* abs.(collect(1:num_sites) .- i0))
    denominator = sum(site_probs)
    localisation_length = Float64(numerator / denominator)
    
    return localisation_length::Float64
end

function hp_mbs_gap_size(
    eigenvalues::Vector{BigFloat}
)
    middle_index = ceil(Int, length(eigenvalues) / 2)
    gap_size = Float64(abs(eigenvalues[middle_index]) + abs(eigenvalues[middle_index + 1]))

    return gap_size::Float64
end

function np_mbs_gap_size(
    eigenvalues::Vector{Float64}
)
    middle_index = ceil(Int, length(eigenvalues) / 2)
    gap_size = abs(eigenvalues[middle_index]) + abs(eigenvalues[middle_index + 1])

    return gap_size::Float64
end

###########################################################
################ Sec 3: Solving Functions #################
###########################################################


function hp_generic_solver(
    N_range::Vector{Int},
    t_n_range::Vector{Vector{Float64}},
    mu_range::Vector{Float64},
    Delta_range::Vector{Float64},
    sequences::Vector{Vector{Int}},
    sequence_name::String,
    precision::Int,
    chunk_size::Int,
    filepath::String,
    opts::UserOptions
)
    """
    Notes:
        - hp_generic_solver iterates over all possible parameter ranges indiscriminately, hence 'generic'.
        - The use of thread IDs to distribute tasks and data storage should only be used on local machines where such ID-ing is known
        - This can be used when the optimal loop method is not known or not needed.
    CAUTION: 
        - This ProgressMeter @showprogress may not give accurate representation of time remaining if disordered parameter looping results in variation in loop time over runtime.
        - The chunk_size must be sufficiently small to not exceed the memory allocated to this task (again paying attention to indiscriminate loop orders).
    """

    # Thread-local data store
    thread_local_results = Dict(Threads.threadid() => DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Missing}[]
    ))

    # Thread-local chunk counters
    thread_local_chunks = Dict(Threads.threadid() => 1)

    # Iterate over parameter combinations in parallel
    @showprogress Threads.@threads for idx in CartesianIndices((length(N_range), length(t_n_range), length(mu_range), length(Delta_range), length(sequences)))
        thread_id = Threads.threadid()

        # Initialize thread-local storage for this thread (if not already initialized)
        if !haskey(thread_local_results, thread_id)
            thread_local_results[thread_id] = DataFrame(
                N = Int[],
                t_n = Vector{Float64}[],
                mu = Float64[],
                Delta = Float64[],
                sequence_name = String[],
                mp = Float64[],
                maj_gap = Float64[],
                ipr = Float64[],
                loc_len = Float64[],
                eigenvalues = Union{Vector{Float64}, Missing}[],
                eigenvectors = Union{Matrix{Float64}, Missing}[]
            )
            thread_local_chunks[thread_id] = 1
        end

        results_df = thread_local_results[thread_id]
        chunk_idx = thread_local_chunks[thread_id]

        # Extract parameters
        N = N_range[idx[1]]
        t_n = t_n_range[idx[2]]
        mu = mu_range[idx[3]]
        Delta = Delta_range[idx[4]]
        sequence = sequences[idx[5]]

        # Perform computations
        truncated_sequence = Vector(sequence[1:N])
        BdG = hp_create_bdg_hamiltonian(N, t_n, mu, Delta, truncated_sequence, precision)
        evals, evecs = GenericLinearAlgebra.eigen(Hermitian(BdG))

        # Calculation options (use opts)
        mp = opts.calc_mp ? hp_calc_maj_mp(evecs) : NaN
        gap = opts.calc_mbs_energy_gap ? hp_mbs_gap_size(evals) : NaN
        ipr = opts.calc_ipr ? hp_calc_maj_ipr(evecs) : NaN
        loc_len = opts.calc_loc_len ? hp_calc_maj_loc_len(evecs) : NaN

        # Eigenvalue saving options (symbol-based)
        eigenvalues_to_save = begin
            if opts.save_evals == :all_hp
                evals
            elseif opts.save_evals == :all_np
                Float64.(evals)
            elseif opts.save_evals == :maj_hp
                mid = length(evals) ÷ 2
                evals[mid:mid+1]
            elseif opts.save_evals == :maj_np
                mid = length(evals) ÷ 2
                (Float64.(evals))[mid:mid+1]
            else
                missing
            end
        end

        # Eigenvector saving options (symbol-based)
        eigenvectors_to_save = begin
            if opts.save_evecs == :all_hp
                evecs
            elseif opts.save_evecs == :all_np
                Float64.(evecs)
            elseif opts.save_evecs == :maj_hp
                mid = size(evecs, 2) ÷ 2
                evecs[:, mid:mid+1]
            elseif opts.save_evecs == :maj_np
                mid = size(evecs, 2) ÷ 2
                (Float64.(evecs))[:, mid:mid+1]
            else
                missing
            end
        end

        # Append results to the thread's local DataFrame
        push!(results_df, (
            N = N,
            t_n = t_n,
            mu = mu,
            Delta = Delta,
            sequence_name = sequence_name,
            mp = mp,
            maj_gap = gap,
            ipr = ipr,
            loc_len = loc_len,
            eigenvalues = eigenvalues_to_save,
            eigenvectors = eigenvectors_to_save
        ))

        # Save chunk if the DataFrame reaches the chunk size
        if nrow(results_df) >= chunk_size
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
            empty!(results_df)
            thread_local_chunks[thread_id] += 1
        end
    end

    # Save any remaining rows in each thread's DataFrame
    for thread_id in keys(thread_local_results)
        results_df = thread_local_results[thread_id]
        if nrow(results_df) > 0
            chunk_idx = thread_local_chunks[thread_id]
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
        end
    end

    return nothing
end

function np_generic_solver(
    N_range::Vector{Int},
    t_n_range::Vector{Vector{Float64}},
    mu_range::Vector{Float64},
    Delta_range::Vector{Float64},
    sequences::Vector{Vector{Int}},
    sequence_name::String,
    chunk_size::Int,
    filepath::String,
    opts::UserOptions
)
    """
        Notes:
            - np_generic_solver iterates over all possible parameter ranges idnescriminately, hence 'generic'.
            - The use of thread IDs to distribute tasks and data storage shoul donly be used on local machines where such ID-ing is known
            - This can be used when the optimal loop method is not known or not needed.
        CAUTION: 
            - This ProgressMeter @showprogress may not give accurate representation of time remaining if disordered parameter looping results in variation in loop time over runtime.
            - The chunk_size must be sufficiently small to not exceed the memory allocated to this task (again paying attention to idnescriminate loop orders).
    """

    # Thread-local data store
    thread_local_results = Dict(Threads.threadid() => DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Missing}[]
    ))

    # Thread-local chunk counters
    thread_local_chunks = Dict(Threads.threadid() => 1)

    # Iterate over parameter combinations in parallel
    @showprogress Threads.@threads for idx in CartesianIndices((length(N_range), length(t_n_range), length(mu_range), length(Delta_range), length(sequences)))
        thread_id = Threads.threadid()

        # Initialize thread-local storage for this thread (if not already initialized)
        if !haskey(thread_local_results, thread_id)
            thread_local_results[thread_id] = DataFrame(
                N = Int[],
                t_n = Vector{Float64}[],
                mu = Float64[],
                Delta = Float64[],
                sequence_name = String[],
                mp = Float64[],
                maj_gap = Float64[],
                ipr = Float64[],
                loc_len = Float64[],
                eigenvalues = Union{Vector{Float64}, Missing}[],
                eigenvectors = Union{Matrix{Float64}, Missing}[]
            )
            thread_local_chunks[thread_id] = 1
        end

        results_df = thread_local_results[thread_id]
        chunk_idx = thread_local_chunks[thread_id]

        # Extract parameters
        N = N_range[idx[1]]
        t_n = t_n_range[idx[2]]
        mu = mu_range[idx[3]]
        Delta = Delta_range[idx[4]]
        sequence = sequences[idx[5]]

        # Perform computations to solve Hamiltonian
        truncated_sequence = Vector(sequence[1:N])
        BdG = np_create_bdg_hamiltonian(N, t_n, mu, Delta, truncated_sequence)
        evals, evecs = LinearAlgebra.eigen(Hermitian(BdG))

        # Calculation options
        mp = opts.calc_mp ? np_calc_maj_mp(evecs) : NaN
        gap = opts.calc_mbs_energy_gap ? np_mbs_gap_size(evals) : NaN
        ipr = opts.calc_ipr ? np_calc_maj_ipr(evecs) : NaN
        loc_len = opts.calc_loc_len ? np_calc_maj_loc_len(evecs) : NaN


        # Eigenvalue saving options (symbol-based)
        eigenvalues_to_save = begin
            if opts.save_evals == :all_np
                evals
            elseif opts.save_evals == :maj_np
                mid = length(evals) ÷ 2
                evals[mid:mid+1]
            else
                missing
            end
        end

        # Eigenvector saving options (symbol-based)
        eigenvectors_to_save = begin
            if opts.save_evecs == :all_np
                evecs
            elseif opts.save_evecs == :maj_np
                mid = size(evecs, 2) ÷ 2
                evecs[:, mid:mid+1]
            else
                missing
            end
        end

        # Append results to the thread's local DataFrame
        push!(results_df, (
            N = N,
            t_n = t_n,
            mu = mu,
            Delta = Delta,
            sequence_name = sequence_name,
            mp = mp,
            maj_gap = gap,
            ipr = ipr,
            loc_len = loc_len,
            eigenvalues = eigenvalues_to_save,
            eigenvectors = eigenvectors_to_save
        ))

        # Save chunk if the DataFrame reaches the chunk size
        if nrow(results_df) >= chunk_size
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
            empty!(results_df)
            thread_local_chunks[thread_id] += 1
        end
    end

    # Save any remaining rows in each thread's DataFrame
    for thread_id in keys(thread_local_results)
        results_df = thread_local_results[thread_id]
        if nrow(results_df) > 0
            chunk_idx = thread_local_chunks[thread_id]
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
        end
    end

    return nothing
end

function np_mu_rho_restricted_solver(
    N_range::Vector{Int},
    t_n_range::Vector{Vector{Float64}},
    mu_range::Vector{Float64},
    unrestricted_points::Vector{Tuple{Float64, Float64}},
    Delta_range::Vector{Float64},
    sequences::Vector{Vector{Int}},
    sequence_name::String,
    chunk_size::Int,
    filepath::String,
    opts::UserOptions
)
    """
        Notes:
            - np_generic_solver iterates over all possible parameter ranges idnescriminately, hence 'generic'.
            - The use of thread IDs to distribute tasks and data storage shoul donly be used on local machines where such ID-ing is known
            - This can be used when the optimal loop method is not known or not needed.
        CAUTION: 
            - This ProgressMeter @showprogress may not give accurate representation of time remaining if disordered parameter looping results in variation in loop time over runtime.
            - The chunk_size must be sufficiently small to not exceed the memory allocated to this task (again paying attention to idnescriminate loop orders).
    """

    # Thread-local data store
    thread_local_results = Dict(Threads.threadid() => DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Missing}[]
    ))

    # Thread-local chunk counters
    thread_local_chunks = Dict(Threads.threadid() => 1)

    # Iterate over parameter combinations in parallel
    @showprogress Threads.@threads for idx in CartesianIndices((length(N_range), length(t_n_range), length(mu_range), length(Delta_range), length(sequences)))
        thread_id = Threads.threadid()

        # Initialize thread-local storage for this thread (if not already initialized)
        if !haskey(thread_local_results, thread_id)
            thread_local_results[thread_id] = DataFrame(
                N = Int[],
                t_n = Vector{Float64}[],
                mu = Float64[],
                Delta = Float64[],
                sequence_name = String[],
                mp = Float64[],
                maj_gap = Float64[],
                ipr = Float64[],
                loc_len = Float64[],
                eigenvalues = Union{Vector{Float64}, Missing}[],
                eigenvectors = Union{Matrix{Float64}, Missing}[]
            )
            thread_local_chunks[thread_id] = 1
        end

        results_df = thread_local_results[thread_id]
        chunk_idx = thread_local_chunks[thread_id]

        # Extract parameters
        N = N_range[idx[1]]
        t_n = t_n_range[idx[2]]
        mu = mu_range[idx[3]]
        Delta = Delta_range[idx[4]]
        sequence = sequences[idx[5]]

        # mu = mu_rho_point[1]
        # t_n = [1.0, mu_rho_point[2]]  # Based on t1=1.0 always!

        rho = t_n[2] / t_n[1] 

        is_unrestricted = !any(p -> p[1] == mu && p[2] == rho, unrestricted_points)

        if is_unrestricted
            # Perform computations to solve Hamiltonian
            truncated_sequence = Vector(sequence[1:N])
            BdG = np_create_bdg_hamiltonian(N, t_n, mu, Delta, truncated_sequence)
            evals, evecs = LinearAlgebra.eigen(Hermitian(BdG))

            # Calculation options
            mp = opts.calc_mp ? np_calc_maj_mp(evecs) : NaN
            gap = opts.calc_mbs_energy_gap ? np_mbs_gap_size(evals) : NaN
            ipr = opts.calc_ipr ? np_calc_maj_ipr(evecs) : NaN
            loc_len = opts.calc_loc_len ? np_calc_maj_loc_len(evecs) : NaN


            # Eigenvalue saving options (symbol-based)
            eigenvalues_to_save = begin
                if opts.save_evals == :all_np
                    evals
                elseif opts.save_evals == :maj_np
                    mid = length(evals) ÷ 2
                    evals[mid:mid+1]
                else
                    missing
                end
            end

            # Eigenvector saving options (symbol-based)
            eigenvectors_to_save = begin
                if opts.save_evecs == :all_np
                    evecs
                elseif opts.save_evecs == :maj_np
                    mid = size(evecs, 2) ÷ 2
                    evecs[:, mid:mid+1]
                else
                    missing
                end
            end

            # Append results to the thread's local DataFrame
            push!(results_df, (
                N = N,
                t_n = t_n,
                mu = mu,
                Delta = Delta,
                sequence_name = sequence_name,
                mp = mp,
                maj_gap = gap,
                ipr = ipr,
                loc_len = loc_len,
                eigenvalues = eigenvalues_to_save,
                eigenvectors = eigenvectors_to_save
            ))
        
        else
            push!(results_df, (
                N = N,
                t_n = t_n,
                mu = mu,
                Delta = Delta,
                sequence_name = sequence_name,
                mp = NaN,
                maj_gap = NaN,
                ipr = NaN,
                loc_len = NaN,
                eigenvalues = missing,
                eigenvectors = missing
            ))
        end

        # Save chunk if the DataFrame reaches the chunk size
        if nrow(results_df) >= chunk_size
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
            empty!(results_df)
            thread_local_chunks[thread_id] += 1
        end
    end

    # Save any remaining rows in each thread's DataFrame
    for thread_id in keys(thread_local_results)
        results_df = thread_local_results[thread_id]
        if nrow(results_df) > 0
            chunk_idx = thread_local_chunks[thread_id]
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
        end
    end

    return nothing
end

function hp_mu_rho_restricted_solver(
    N_range::Vector{Int},
    t_n_range::Vector{Vector{Float64}},
    mu_range::Vector{Float64},
    unrestricted_points::Vector{Tuple{Float64, Float64}},
    Delta_range::Vector{Float64},
    sequences::Vector{Vector{Int}},
    sequence_name::String,
    precision::Int,
    chunk_size::Int,
    filepath::String,
    opts::UserOptions
)
    """
        Notes:
            - np_generic_solver iterates over all possible parameter ranges idnescriminately, hence 'generic'.
            - The use of thread IDs to distribute tasks and data storage shoul donly be used on local machines where such ID-ing is known
            - This can be used when the optimal loop method is not known or not needed.
        CAUTION: 
            - This ProgressMeter @showprogress may not give accurate representation of time remaining if disordered parameter looping results in variation in loop time over runtime.
            - The chunk_size must be sufficiently small to not exceed the memory allocated to this task (again paying attention to idnescriminate loop orders).
    """

    # Thread-local data store
    thread_local_results = Dict(Threads.threadid() => DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Missing}[]
    ))

    # Thread-local chunk counters
    thread_local_chunks = Dict(Threads.threadid() => 1)

    # Iterate over parameter combinations in parallel
    @showprogress Threads.@threads for idx in CartesianIndices((length(N_range), length(t_n_range), length(mu_range), length(Delta_range), length(sequences)))
        thread_id = Threads.threadid()

        # Initialize thread-local storage for this thread (if not already initialized)
        if !haskey(thread_local_results, thread_id)
            thread_local_results[thread_id] = DataFrame(
                N = Int[],
                t_n = Vector{Float64}[],
                mu = Float64[],
                Delta = Float64[],
                sequence_name = String[],
                mp = Float64[],
                maj_gap = Float64[],
                ipr = Float64[],
                loc_len = Float64[],
                eigenvalues = Union{Vector{Float64}, Missing}[],
                eigenvectors = Union{Matrix{Float64}, Missing}[]
            )
            thread_local_chunks[thread_id] = 1
        end

        results_df = thread_local_results[thread_id]
        chunk_idx = thread_local_chunks[thread_id]

        # Extract parameters
        N = N_range[idx[1]]
        t_n = t_n_range[idx[2]]
        mu = mu_range[idx[3]]
        Delta = Delta_range[idx[4]]
        sequence = sequences[idx[5]]

        # mu = mu_rho_point[1]
        # t_n = [1.0, mu_rho_point[2]]  # Based on t1=1.0 always!

        rho = t_n[2] / t_n[1] 

        is_unrestricted = !any(p -> p[1] == mu && p[2] == rho, unrestricted_points)

        if is_unrestricted
            # Perform computations to solve Hamiltonian
            truncated_sequence = Vector(sequence[1:N])
            BdG = np_create_bdg_hamiltonian(N, t_n, mu, Delta, truncated_sequence)
            evals, evecs = GenericLinearAlgebra.eigen(Hermitian(BdG))

            # Calculation options
            mp = opts.calc_mp ? np_calc_maj_mp(evecs) : NaN
            gap = opts.calc_mbs_energy_gap ? np_mbs_gap_size(evals) : NaN
            ipr = opts.calc_ipr ? np_calc_maj_ipr(evecs) : NaN
            loc_len = opts.calc_loc_len ? np_calc_maj_loc_len(evecs) : NaN


            # Eigenvalue saving options (symbol-based)
            eigenvalues_to_save = begin
                if opts.save_evals == :all_hp
                    evals
                elseif opts.save_evals == :all_np
                    Float64.(evals)
                elseif opts.save_evals == :maj_hp
                    mid = length(evals) ÷ 2
                    evals[mid:mid+1]
                elseif opts.save_evals == :maj_np
                    mid = length(evals) ÷ 2
                    (Float64.(evals))[mid:mid+1]
                else
                    missing
                end
            end

            # Eigenvector saving options (symbol-based)
            eigenvectors_to_save = begin
                if opts.save_evecs == :all_hp
                    evecs
                elseif opts.save_evecs == :all_np
                    Float64.(evecs)
                elseif opts.save_evecs == :maj_hp
                    mid = size(evecs, 2) ÷ 2
                    evecs[:, mid:mid+1]
                elseif opts.save_evecs == :maj_np
                    mid = size(evecs, 2) ÷ 2
                    (Float64.(evecs))[:, mid:mid+1]
                else
                    missing
                end
            end

            # Append results to the thread's local DataFrame
            push!(results_df, (
                N = N,
                t_n = t_n,
                mu = mu,
                Delta = Delta,
                sequence_name = sequence_name,
                mp = mp,
                maj_gap = gap,
                ipr = ipr,
                loc_len = loc_len,
                eigenvalues = eigenvalues_to_save,
                eigenvectors = eigenvectors_to_save
            ))
        
        else
            push!(results_df, (
                N = N,
                t_n = t_n,
                mu = mu,
                Delta = Delta,
                sequence_name = sequence_name,
                mp = NaN,
                maj_gap = NaN,
                ipr = NaN,
                loc_len = NaN,
                eigenvalues = missing,
                eigenvectors = missing
            ))
        end

        # Save chunk if the DataFrame reaches the chunk size
        if nrow(results_df) >= chunk_size
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
            empty!(results_df)
            thread_local_chunks[thread_id] += 1
        end
    end

    # Save any remaining rows in each thread's DataFrame
    for thread_id in keys(thread_local_results)
        results_df = thread_local_results[thread_id]
        if nrow(results_df) > 0
            chunk_idx = thread_local_chunks[thread_id]
            file_name = "$(filepath)_thread_$(thread_id)_chunk_$(chunk_idx).bson"
            @save file_name results_df
        end
    end

    return nothing
end

# not integrated
function hp_mu_loop_solver(
    N::Int,
    t_n::Vector{Float64},
    mu_range::Vector{Float64},
    Delta::Float64,
    sequence::Vector{Int},
    sequence_name::String,
    precision::Int,
    filepath::String
)
    """
        Notes:
            - hp_mu_loop_solver is optimised for fixed N, t_n, Delta and sequence parameters, varying only a range of mu.
            - Originaly written for serial HPC batch jobs where each job would then vary a secondary parameter.
        
        CAUTION:
            - Must ensure parameters and density of mu_range are such that allocated memory is not exceeded before the end of single loop to allow data to be safely saved to .bson
    """

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

# not integrated
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

# not integrated
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


end # end module
