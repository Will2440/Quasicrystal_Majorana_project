module LocMPScales

using LinearAlgebra
using ProgressMeter
using Base.Threads


"""
    calc_local_scales_mp(eigenvectors::Matrix{Float64})

Given BdG eigenvectors (2Nx2N), computes the local Majorana polarisation (MP) for all possible local chains at all even scales from 2 up to N.
Returns a vector of results, each entry is a NamedTuple with all numerators, denominators, MPs, and indices.

Each result contains:
    :scale         => scale (even integer, 2, 4, ..., N)
    :start_site    => starting site of the local chain (1-based)
    :end_site      => ending site of the local chain
    :numerator_L, :denominator_L, :MP_L
    :numerator_R, :denominator_R, :MP_R
    :overall_MP    => MP_L * MP_R
"""
function calc_local_scales_mp(
    eigenvectors::Union{Matrix{Float64}, Matrix{BigFloat}, Matrix{Complex{Float64}}, Matrix{Complex{BigFloat}}}
)::Vector{NamedTuple}
    
    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2
    middle_state_index = N  # Use the middle eigenstate (candidate MBS)
    middle_plus_one_state_index = N + 1
    u_mid = eigenvectors[1:N, middle_state_index]
    v_mid = eigenvectors[N+1:end, middle_state_index]
    u_mpo = eigenvectors[1:N, middle_plus_one_state_index]
    v_mpo = eigenvectors[N+1:end, middle_plus_one_state_index]

    results = NamedTuple[]
    @showprogress for scale in 2:2:N # Only even scales to ensure equal left/right splits
        for start_site in 1:(N - scale + 1)
            end_site = start_site + scale - 1
            left_inds = start_site:(start_site + scale ÷ 2 - 1)
            right_inds = (start_site + scale ÷ 2):end_site

            numerator_L = sum(u_mid[i] * conj(v_mid[i]) + v_mid[i] * conj(u_mid[i]) for i in left_inds)
            denominator_L = sum(abs2(u_mid[i]) + abs2(v_mid[i]) for i in left_inds)
            # MP_L = denominator_L ≈ 0 ? 0.0 : numerator_L / denominator_L
            MP_L = real(numerator_L / denominator_L)

            numerator_R = sum(u_mpo[i] * conj(v_mpo[i]) + v_mpo[i] * conj(u_mpo[i]) for i in right_inds)
            denominator_R = sum(abs2(u_mpo[i]) + abs2(v_mpo[i]) for i in right_inds)
            # MP_R = denominator_R ≈ 0 ? 0.0 : numerator_R / denominator_R
            MP_R = real(numerator_R / denominator_R)

            overall_MP = MP_L * MP_R
            overall_MP = Float64(overall_MP)

            push!(results, (
                scale=scale,
                start_site=start_site,
                end_site=end_site,
                numerator_L=numerator_L,
                denominator_L=denominator_L,
                MP_L=MP_L,
                numerator_R=numerator_R,
                denominator_R=denominator_R,
                MP_R=MP_R,
                overall_MP=overall_MP
            ))
        end
    end
    return results
end

function calc_local_scales_mp_optimised(
    eigenvectors::Union{Matrix{Float64}, Matrix{BigFloat}, Matrix{Complex{Float64}}, Matrix{Complex{BigFloat}}}
)::Vector{NamedTuple}

    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2
    middle_state_index = N
    middle_plus_one_state_index = N + 1
    u_mid = eigenvectors[1:N, middle_state_index]
    v_mid = eigenvectors[N+1:end, middle_state_index]
    u_mpo = eigenvectors[1:N, middle_plus_one_state_index]
    v_mpo = eigenvectors[N+1:end, middle_plus_one_state_index]

    # Pre-compute all (scale, start_site) pairs to flatten the workload
    tasks = NamedTuple{(:scale, :start_site), Tuple{Int, Int}}[]
    for scale in 2:2:N
        for start_site in 1:(N - scale + 1)
            push!(tasks, (scale=scale, start_site=start_site))
        end
    end

    # Pre-allocate results vector
    n_tasks = length(tasks)
    results = Vector{NamedTuple}(undef, n_tasks)

    # Parallel computation over all tasks
    Threads.@threads for i in 1:n_tasks
        scale, start_site = tasks[i]
        end_site = start_site + scale - 1
        left_inds = start_site:(start_site + scale ÷ 2 - 1)
        right_inds = (start_site + scale ÷ 2):end_site

        # Vectorised sums for slight speedup (Julia broadcasts efficiently)
        numerator_L = sum(u_mid[left_inds] .* conj.(v_mid[left_inds]) .+ v_mid[left_inds] .* conj.(u_mid[left_inds]))
        denominator_L = sum(abs2.(u_mid[left_inds]) .+ abs2.(v_mid[left_inds]))
        MP_L = numerator_L / denominator_L

        numerator_R = sum(u_mpo[right_inds] .* conj.(v_mpo[right_inds]) .+ v_mpo[right_inds] .* conj.(u_mpo[right_inds]))
        denominator_R = sum(abs2.(u_mpo[right_inds]) .+ abs2.(v_mpo[right_inds]))
        MP_R = numerator_R / denominator_R

        overall_MP = conj(MP_L) * MP_R
        try
            overall_MP = Float64(overall_MP)
        catch err
            @warn "overall_MP is complex, taking real part. Details:" value=overall_MP error=err
            overall_MP = Float64(real(overall_MP))
        end

        results[i] = (
            scale=scale,
            start_site=start_site,
            end_site=end_site,
            numerator_L=numerator_L,
            denominator_L=denominator_L,
            MP_L=MP_L,
            numerator_R=numerator_R,
            denominator_R=denominator_R,
            MP_R=MP_R,
            overall_MP=overall_MP
        )
    end

    return results
end

function calc_local_scales_mp_allscales_optimised(
    eigenvectors::Union{Matrix{Float64}, Matrix{BigFloat}, Matrix{Complex{Float64}}, Matrix{Complex{BigFloat}}}
)::Vector{NamedTuple}

    num_rows = size(eigenvectors, 1)
    N = num_rows ÷ 2
    middle_state_index = N
    middle_plus_one_state_index = N + 1
    u_mid = eigenvectors[1:N, middle_state_index]
    v_mid = eigenvectors[N+1:end, middle_state_index]
    u_mpo = eigenvectors[1:N, middle_plus_one_state_index]
    v_mpo = eigenvectors[N+1:end, middle_plus_one_state_index]

    # Pre-compute all (scale, start_site) pairs for all scales
    tasks = NamedTuple{(:scale, :start_site), Tuple{Int, Int}}[]
    for scale in 1:N
        for start_site in 1:(N - scale + 1)
            push!(tasks, (scale=scale, start_site=start_site))
        end
    end

    n_tasks = length(tasks)
    num_threads = Threads.nthreads()
    results = Vector{NamedTuple}(undef, n_tasks)

    Threads.@threads for thread_idx in 1:num_threads
        for i in thread_idx:num_threads:n_tasks
            scale, start_site = tasks[i]
            end_site = start_site + scale - 1

            # For odd scales, include the middle site in both left and right
            if isodd(scale)
                mid = start_site + (scale - 1) ÷ 2
                left_inds = start_site:mid
                right_inds = mid:end_site
            else
                left_inds = start_site:(start_site + scale ÷ 2 - 1)
                right_inds = (start_site + scale ÷ 2):end_site
            end

            numerator_L = sum(u_mid[left_inds] .* conj.(v_mid[left_inds]) .+ v_mid[left_inds] .* conj.(u_mid[left_inds]))
            denominator_L = sum(abs2.(u_mid[left_inds]) .+ abs2.(v_mid[left_inds]))
            MP_L = numerator_L / denominator_L

            numerator_R = sum(u_mpo[right_inds] .* conj.(v_mpo[right_inds]) .+ v_mpo[right_inds] .* conj.(u_mpo[right_inds]))
            denominator_R = sum(abs2.(u_mpo[right_inds]) .+ abs2.(v_mpo[right_inds]))
            MP_R = numerator_R / denominator_R

            overall_MP = conj(MP_L) * MP_R
            try
                overall_MP = Float64(overall_MP)
            catch err
                @warn "overall_MP is complex, taking real part. Details:" value=overall_MP error=err
                overall_MP = Float64(real(overall_MP))
            end

            results[i] = (
                scale=scale,
                start_site=start_site,
                end_site=end_site,
                numerator_L=numerator_L,
                denominator_L=denominator_L,
                MP_L=MP_L,
                numerator_R=numerator_R,
                denominator_R=denominator_R,
                MP_R=MP_R,
                overall_MP=overall_MP
            )
        end
    end

    return results
end


end # module LocMPScales