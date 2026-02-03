module HofstadterButterfly

using LinearAlgebra
using SparseArrays
# using CairoMakie
using ProgressMeter

"""
    harper_hamiltonian(q, alpha, kx, ky; tx=1.0, ty=1.0)

Constructs the q×q Harper-Hofstadter Hamiltonian matrix for a specific magnetic flux and k-point.

# Arguments
- `q::Int`: Denominator of the flux fraction α = p/q. Defines the size of the magnetic unit cell.
- `alpha::Float64`: Magnetic flux per unit cell (in units of flux quantum).
- `kx, ky::Float64`: Crystal momenta in the magnetic Brillouin zone.
- `tx, ty::Float64`: Hopping amplitudes in x and y directions (default 1.0).

# Returns
- `H::Matrix{ComplexF64}`: The Hamiltonian matrix.
"""
function harper_hamiltonian(q::Int, alpha::Float64, kx::Float64, ky::Float64; tx::Float64=1.0, ty::Float64=1.0)
    H = zeros(ComplexF64, q, q)
    for n in 1:q
        # Diagonal terms: Potential modulation due to magnetic field (or hopping in y-direction after Fourier transform)
        # Corresponds to 2 * ty * cos(2π*α*n + kx)
        H[n, n] = 2 * ty * cos(2π * alpha * (n - 1) + kx)
        
        # Off-diagonal terms: Nearest-neighbor hopping within the magnetic unit cell
        # The phase exp(±iky) comes from the Bloch condition in the y-direction
        # We use mod1 to handle the periodic boundary condition of the unit cell (wrapping q -> 1)
        H[n, mod1(n + 1, q)] = -tx * exp(1im * ky)
        H[n, mod1(n - 1, q)] = -tx * exp(-1im * ky)
    end
    return H
end

"""
    butterfly_spectrum(qmax=200; kx=0.0, tx=1.0, ty=1.0)

Computes the energy spectrum of the Harper model for a range of rational fluxes p/q.

# Arguments
- `qmax::Int`: Maximum denominator q to scan (fluxes will be p/q for q in 2:qmax).
- `kx::Float64`: Fixed momentum kx (usually 0.0 for the butterfly plot).
- `tx, ty`: Hopping parameters.

# Returns
- `flux_vals`: Vector of alpha values.
- `energies`: Vector of corresponding energy eigenvalues.
"""
function butterfly_spectrum(qmax=200; kx=0.0, tx=1.0, ty=1.0)
    flux_vals = Float64[]
    energies  = Float64[]

    # Iterate over all denominators q
    @showprogress for q in 2:qmax
        # Iterate over all numerators p coprime to q (or all p for full spectrum)
        for p in 0:q
            alpha = p / q
            
            # Construct Hamiltonian at ky=0 (spectrum is usually plotted vs flux, independent of ky for the bulk bands)
            # Note: For accurate band edges, one might need to scan ky, but for the butterfly visual, ky=0 is standard.
            H = harper_hamiltonian(q, alpha, kx, 0.0; tx=tx, ty=ty)
            
            eigs = eigen(H).values
            
            # Store results
            append!(flux_vals, fill(alpha, length(eigs)))
            append!(energies, eigs)
        end
    end

    return flux_vals, energies
end

"""
    chern_number(Hk)

Calculates the Chern number of a band using the discretized Berry curvature (Fukui-Hatsugai-Suzuki method).

# Arguments
- `Hk`: A 2D array of eigenvectors on the discretized Brillouin zone mesh. 
        Hk[i,j] is the eigenvector at (kx_i, ky_j).
"""
function chern_number(Hk)
    nk = size(Hk, 1)
    Ux = zeros(ComplexF64, nk, nk)
    Uy = similar(Ux)

    # 1. Compute Link Variables U_μ(k) = <n(k) | n(k+μ)> / |<n(k) | n(k+μ)>|
    # These represent the parallel transport of the phase between neighboring k-points.
    for i in 1:nk, j in 1:nk
        ψ = Hk[i, j]
        ψx = Hk[mod1(i + 1, nk), j] # Neighbor in x
        ψy = Hk[i, mod1(j + 1, nk)] # Neighbor in y
        
        # Normalize the overlap to get the U(1) link variable
        overlap_x = ψ' * ψx
        Ux[i, j] = abs(overlap_x) < 1e-12 ? 1.0 + 0.0im : overlap_x / abs(overlap_x)
        
        overlap_y = ψ' * ψy
        Uy[i, j] = abs(overlap_y) < 1e-12 ? 1.0 + 0.0im : overlap_y / abs(overlap_y)
    end

    # 2. Compute Lattice Field Strength F_xy = ln( U_x(k) U_y(k+x) U_x(k+y)^† U_y(k)^† )
    # This is the Berry flux through the plaquette at k.
    F = 0.0
    for i in 1:nk, j in 1:nk
        # Plaquette product around the elementary square in k-space
        U1 = Ux[i, j] * Uy[mod1(i + 1, nk), j]
        U2 = Ux[i, mod1(j + 1, nk)] * Uy[i, j]
        
        # The flux is the phase of the product (branch cut handled by angle/log)
        F += angle(U1 / U2)
    end

    # 3. Sum over BZ and divide by 2π to get the integer invariant
    return F / (2π)
end

"""
    compute_chern_numbers(q, p; nk=15, tx=1.0, ty=1.0)

Computes the Chern numbers for all q bands of the Harper model with flux p/q.

# Arguments
- `nk`: Grid density for the Brillouin zone discretization (nk x nk mesh).
- `tx, ty`: Hopping parameters passed to the Hamiltonian.

Computes the Chern numbers for all q bands of the Harper model with flux p/q.
Optimized to perform eigen decomposition only once per k-point.
"""
function compute_chern_numbers(q, p; nk=15, tx=1.0, ty=1.0)
    alpha = p / q
    ch = zeros(Int, q)

    # 1. Pre-calculate eigenvectors on the k-grid
    # We store the full matrix of eigenvectors for each k-point
    eig_vecs_grid = Array{Matrix{ComplexF64}}(undef, nk, nk)
    
    for i in 1:nk, j in 1:nk
        kx = 2π * (i - 1) / nk
        ky = 2π * (j - 1) / nk
        
        H = harper_hamiltonian(q, alpha, kx, ky; tx=tx, ty=ty)
        # Use eigen! to save allocation since H is temporary
        _, vecs = eigen!(H)
        eig_vecs_grid[i, j] = vecs
    end
    
    # 2. Compute Chern number for each band reusing the pre-calculated eigenvectors
    # Re-use a buffer for the specific band's eigenvectors
    current_band_vecs = Array{Vector{ComplexF64}}(undef, nk, nk)

    for band in 1:q
        for i in 1:nk, j in 1:nk
            # Extract the specific band's eigenvector (column)
            current_band_vecs[i, j] = eig_vecs_grid[i, j][:, band]
        end
        
        ch[band] = round(Int, chern_number(current_band_vecs))
    end

    return ch
end


"""
    coloured_butterfly(qmax=40; tx=1.0, ty=1.0)

Generates the data for the "Coloured Hofstadter Butterfly", where gaps are labelled by the 
integrated Chern number (Hall conductance).

# Arguments
- `qmax`: Maximum denominator for flux generation.
- `tx, ty`: Hopping parameters.

# Returns
- `flux`: Vector of alpha values.
- `energy`: Vector of energy eigenvalues.
- `colour`: Vector of gap labels (integers) corresponding to each energy state.

Optimized with multithreading and load balancing.
"""
function coloured_butterfly(qmax; nk=15, tx=1.0, ty=1.0)
    flux_vals = Float64[]
    energies  = Float64[]
    colours   = Int[]

    # 1. Flatten loop for load balancing
    tasks = Tuple{Int, Int}[]
    for q in 1:qmax
        for p in 0:q
            # Optimization: Skip reducible fractions (e.g. 2/4 is same as 1/2).
            # This avoids redundant calculations and overplotting.
            if gcd(p, q) == 1
                push!(tasks, (q, p))
            end
        end
    end

    # 2. Sort tasks by q descending (Longest Processing Time first)
    # This prevents thread starvation at the end of the run.
    sort!(tasks, by = x -> x[1], rev=true)

    # Weighting for progress bar: Calculation scales roughly as q^3
    total_weight = sum(t[1]^3 for t in tasks)
    prog = Progress(total_weight, desc="Computing coloured butterfly: ", showspeed=true)
    
    io_lock = ReentrantLock()

    Threads.@threads :dynamic for (q, p) in tasks
        
        # Thread-local buffers
        local_flux   = Float64[]
        local_energy = Float64[]
        local_colour = Int[]

        alpha = p / q
        
        # 1. Solve Hamiltonian
        H = harper_hamiltonian(q, alpha, 0.0, 0.0; tx=tx, ty=ty)
        vals = eigen!(H).values 
        
        # 2. Compute Chern numbers
        # (Relies on the optimized compute_chern_numbers function)
        ch = compute_chern_numbers(q, p; nk=nk, tx=tx, ty=ty)
        
        # 3. Store results
        # For a given flux p/q, we have q energy levels and q Chern numbers
        for i in 1:q
            push!(local_flux, alpha)
            push!(local_energy, vals[i])
            push!(local_colour, ch[i])
        end
        
        # 4. Thread-safe write to main arrays
        lock(io_lock) do
            append!(flux_vals, local_flux)
            append!(energies, local_energy)
            append!(colours, local_colour)
            next!(prog; step=q^3)
        end
    end

    return flux_vals, energies, colours
end

"""
    coloured_butterfly_gaps(qmax=40; tx=1.0, ty=1.0, threshold=1e-4)

Generates data to colour the *gaps* of the Hofstadter butterfly.
Calculates eigenvalues and Chern numbers, then identifies open gaps and assigns 
them the integrated Chern number (Hall conductance).

# Arguments
- `qmax`: Maximum denominator for flux generation.
- `tx, ty`: Hopping parameters for the Harper Hamiltonian.
- `threshold`: Minimum energy difference to consider a gap "open".

# Returns
- `flux`: Vector of alpha values.
- `gap_min`: Vector of gap lower bounds (top of band n).
- `gap_max`: Vector of gap upper bounds (bottom of band n+1).
- `gap_label`: Vector of integers (Hall conductance) for each gap.
"""
function coloured_butterfly_gaps_unthreaded(qmax=40; nk=15, tx=1.0, ty=1.0, threshold=1e-4)
    flux_vals = Float64[]
    gap_min   = Float64[]
    gap_max   = Float64[]
    gap_label = Int[]

    q_range = 3:qmax
    # Weighting for progress bar: Calculation scales roughly as q^4
    weights = q_range .^ 4
    total_weight = sum(weights)
    
    prog = Progress(total_weight, desc="Computing butterfly gaps: ", showspeed=true)

    for (idx, q) in enumerate(q_range)
        for p in 0:q
            alpha = p / q
            
            # 1. Solve Hamiltonian at a reference point (k=0) to get band edges
            # Note: Strictly, one should scan k-space to find the global min/max of bands,
            # but for the Harper model, the gaps are usually open at k=0 or consistent enough for visualization.
            H = harper_hamiltonian(q, alpha, 0.0, 0.0; tx=tx, ty=ty)
            vals = eigen(H).values # Eigenvalues are sorted ascending
            
            # 2. Compute Chern numbers for all bands
            # This requires diagonalizing H(k) on a grid in the Brillouin zone
            ch = compute_chern_numbers(q, p; nk=nk, tx=tx, ty=ty)
            
            # 3. Calculate Gap Labels (Integrated Chern Number)
            # The label for the gap *after* band n is the sum of Chern numbers 1..n
            integ_ch = cumsum(ch)
            
            # 4. Extract Gaps
            # Iterate through the spaces between consecutive eigenvalues
            for n in 1:(q-1)
                E_n = vals[n]       # Top of band n (approx)
                E_next = vals[n+1]  # Bottom of band n+1 (approx)
                
                # Check if the gap is open
                if E_next - E_n > threshold
                    push!(flux_vals, alpha)
                    push!(gap_min, E_n)
                    push!(gap_max, E_next)
                    push!(gap_label, integ_ch[n])
                end
            end
        end
        next!(prog; step=weights[idx])
    end

    return flux_vals, gap_min, gap_max, gap_label
end

function coloured_butterfly_gaps_threaded(qmax=40; nk=15, tx=1.0, ty=1.0, threshold=1e-4)
    flux_vals = Float64[]
    gap_min   = Float64[]
    gap_max   = Float64[]
    gap_label = Int[]

    # Flatten the loop to improve load balancing.
    # Each task is a tuple (q, p).
    tasks = Tuple{Int, Int}[]
    for q in 3:qmax
        for p in 0:q
            push!(tasks, (q, p))
        end
    end

    # Sort tasks by q descending (Longest Processing Time first)
    # This ensures the most expensive calculations start immediately, preventing
    # a situation where one thread gets stuck with a huge q at the very end.
    sort!(tasks, by = x -> x[1], rev=true)

    # Weighting for progress bar: Calculation scales roughly as q^3 (dominated by eigen)
    total_weight = sum(t[1]^3 for t in tasks)
    
    prog = Progress(total_weight, desc="Computing butterfly gaps: ", showspeed=true)
    
    # Lock for thread-safe writing
    io_lock = ReentrantLock()

    Threads.@threads :dynamic for (q, p) in tasks
        
        # Use local buffers to avoid locking inside the inner loop
        local_flux   = Float64[]
        local_gmin   = Float64[]
        local_gmax   = Float64[]
        local_glabel = Int[]

        alpha = p / q
        
        # 1. Solve Hamiltonian
        H = harper_hamiltonian(q, alpha, 0.0, 0.0; tx=tx, ty=ty)
        vals = eigen!(H).values # Use eigen! to reduce memory allocation
        
        # 2. Compute Chern numbers
        ch = compute_chern_numbers(q, p; nk=nk, tx=tx, ty=ty)
        
        # 3. Calculate Gap Labels
        integ_ch = cumsum(ch)
        
        # 4. Extract Gaps
        for n in 1:(q-1)
            E_n = vals[n]       
            E_next = vals[n+1]  
            
            if E_next - E_n > threshold
                push!(local_flux, alpha)
                push!(local_gmin, E_n)
                push!(local_gmax, E_next)
                push!(local_glabel, integ_ch[n])
            end
        end
        
        # Safely merge local results into the main arrays
        lock(io_lock) do
            append!(flux_vals, local_flux)
            append!(gap_min, local_gmin)
            append!(gap_max, local_gmax)
            append!(gap_label, local_glabel)
            next!(prog; step=q^3)
        end
    end

    return flux_vals, gap_min, gap_max, gap_label
end

"""
    coloured_butterfly_diophantine_gaps(qmax=40; tx=1.0, ty=1.0, threshold=1e-4)

Computes the Hofstadter butterfly gaps.
Uses the Diophantine equation (TKNN formula) to label gaps, which is exact and much faster than integration.
"""
function coloured_butterfly_diophantine_gaps(qmax=40; tx=1.0, ty=1.0, threshold=1e-4)
    flux_vals = Float64[]
    gap_min   = Float64[]
    gap_max   = Float64[]
    gap_label = Int[]

    # 1. Generate irreducible flux fractions
    tasks = Tuple{Int, Int}[]
    for q in 3:qmax
        for p in 1:q-1
            # Only calculate for irreducible fractions.
            # Reducible ones (e.g. 2/4) are physically identical to 1/2 but computationally redundant.
            if gcd(p, q) == 1
                push!(tasks, (q, p))
            end
        end
    end

    # Sort by q descending for load balancing
    sort!(tasks, by = x -> x[1], rev=true)

    # Weighting: Eigenvalues scale as q^3
    total_weight = sum(t[1]^3 for t in tasks)
    prog = Progress(total_weight, desc="Computing butterfly gaps (Diophantine): ", showspeed=true)
    
    io_lock = ReentrantLock()

    Threads.@threads :dynamic for (q, p) in tasks
        
        local_flux   = Float64[]
        local_gmin   = Float64[]
        local_gmax   = Float64[]
        local_glabel = Int[]

        alpha = p / q
        
        # 1. Solve Hamiltonian for energies only
        H = harper_hamiltonian(q, alpha, 0.0, 0.0; tx=tx, ty=ty)
        vals = eigen!(H).values 
        
        # 2. Identify Gaps and Label them algebraically
        for n in 1:(q-1)
            E_n = vals[n]       
            E_next = vals[n+1]  
            
            if E_next - E_n > threshold
                # --- THE FIX: DIOPHANTINE LABELLING ---
                # We solve: n ≡ p * sigma (mod q)
                # sigma = n * inverse(p) (mod q)
                
                # invmod(p, q) finds x such that p*x ≡ 1 (mod q)
                sigma = (n * invmod(p, q)) % q
                
                # Wrap sigma to the range (-q/2, q/2] to get the physical Hall conductance
                if sigma > q/2
                    sigma -= q
                end

                push!(local_flux, alpha)
                push!(local_gmin, E_n)
                push!(local_gmax, E_next)
                push!(local_glabel, sigma)
            end
        end
        
        lock(io_lock) do
            append!(flux_vals, local_flux)
            append!(gap_min, local_gmin)
            append!(gap_max, local_gmax)
            append!(gap_label, local_glabel)
            next!(prog; step=q^3)
        end
    end

    return flux_vals, gap_min, gap_max, gap_label
end


end # module HofstadterButterfly