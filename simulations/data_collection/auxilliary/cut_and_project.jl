using BSON: @save, @load
using BSON
using DataFrames
using Base.Threads: @threads, nthreads, threadid, SpinLock
using ProgressMeter
using Plots
using StatsBase

const THREAD_SLOTS = try
    Threads.maxthreads()
catch
    Threads.nthreads()
end

function load_sturmian_grad_bson(path::AbstractString)
    @assert isfile(path) "BSON file not found: $path"
    raw = BSON.load(path)

    # normalize keys to Symbol -> value mapping
    norm = Dict{Symbol,Any}()
    for (k,v) in pairs(raw)
        # keys from BSON.load are usually strings; convert to Symbol safely
        kn = try
            Symbol(k)
        catch
            Symbol(string(k))
        end
        norm[kn] = v
    end

    # 1) If any value is already a DataFrame, return it
    for v in values(norm)
        if v isa DataFrame
            return deepcopy(v)   # return a copy to avoid accidental mutation of loaded object
        end
    end

    # 2) If top-level has prefix & slope arrays, build a DataFrame
    if (haskey(norm, :prefix) || haskey(norm, Symbol("prefix"))) &&
       (haskey(norm, :slope)  || haskey(norm, Symbol("slope")))
        pref = haskey(norm, :prefix) ? norm[:prefix] : norm[Symbol("prefix")]
        slp  = haskey(norm, :slope)  ? norm[:slope]  : norm[Symbol("slope")]
        return DataFrame(prefix = pref, slope = slp)
    end

    # 3) If the file contains a single value which is a Vector{NamedTuple}, convert it
    if length(norm) == 1
        v = first(values(norm))
        if isa(v, Vector) && !isempty(v) && v[1] isa NamedTuple
            return DataFrame(v)
        end
    end

    error("No DataFrame or recognizable (prefix,slope) data found in BSON file: $path. Keys: $(collect(keys(norm)))")
end

function add_rationals(n::Int; start_r::Int=1)
    @assert n >= 1 "n must be >= 1"
    @assert start_r >= 1 "start_r must be >= 1"

    rs = start_r:(start_r + n - 1)
    rats = [1.0 / Float64(r) for r in rs]
    comps = [1.0 - r for r in rats]
    return rats, comps
end

function add_extrema!(df::DataFrame)
    slopes = df.slope
    if !(0.0 in slopes)
        push!(df, (slope=0.0, origin="extrema", type="original", prefix=missing))
    end
    if !(1.0 in slopes)
        push!(df, (slope=1.0, origin="extrema", type="original", prefix=missing))
    end
    return df
end

function generate_full_slope_df(orig_df::DataFrame; n_rationals::Int=10, start_r::Int=2)
    # 1. Add rationals
    rats, _ = add_rationals(n_rationals; start_r=start_r)
    rat_df = DataFrame(slope=rats, origin="rational", type="original", prefix=fill(missing, length(rats)))
    # rat_df = add_extrema!(rat_df)
    # Mark extrema
    extrema_slopes = [0.0, 1.0]
    extrema_df = DataFrame(slope=extrema_slopes, origin="extrema", type="original", prefix=fill(missing, length(extrema_slopes)))

    # 2. Mark original irrationals
    irr_df = deepcopy(orig_df)
    irr_df.origin .= "irrational"
    irr_df.type .= "original"

    # 3. Combine
    base_df = vcat(irr_df, rat_df, extrema_df)

    # 4. Complementary (1 - slope)
    comp_df = deepcopy(base_df)
    comp_df.slope .= 1 .- base_df.slope
    comp_df.type .= "complementary"

    # 5. Reciprocal (1 / slope), avoid division by zero
    rec_df = deepcopy(base_df)
    rec_df.slope .= [s == 0.0 ? NaN : 1.0/s for s in base_df.slope]
    rec_df.type .= "reciprocal"

    # # 6. Negatives (for π/2 to π)
    # neg_base_df = deepcopy(base_df)
    # neg_base_df.slope .= -base_df.slope
    # neg_base_df.type .= "negative_original"

    # neg_comp_df = deepcopy(comp_df)
    # neg_comp_df.slope .= -comp_df.slope
    # neg_comp_df.type .= "negative_complementary"

    # neg_rec_df = deepcopy(rec_df)
    # neg_rec_df.slope .= -rec_df.slope
    # neg_rec_df.type .= "negative_reciprocal"

    # 7. Combine all, remove duplicates, sort
    all_df = vcat(base_df, comp_df, rec_df)#, neg_base_df, neg_comp_df, neg_rec_df)
    # Remove duplicates by slope and type
    all_df = unique(all_df, [:slope, :type])
    # Remove NaN slopes (from reciprocal of zero)
    all_df = all_df[.!isnan.(all_df.slope), :]
    # Sort by slope
    sort!(all_df, :slope)

    return all_df
end

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

function generalised_cut_and_project_sequence(N::Int, alpha::Float64; gamma::Float64=0.0, map_to::Tuple{Int,Int}=(1,2))
    alpha_floor = alpha - floor(alpha)
    gamma = gamma - floor(gamma) # wrap into [0,1)
    seq = Vector{Int}(undef, N)
    frac = gamma
    for n in 1:N
        frac += alpha_floor
        if frac >= 1.0
            frac -= 1.0
            seq[n] = map_to[1]
        elseif frac < 0.0
            frac += 1.0
            seq[n] = map_to[2]
        else
            seq[n] = map_to[2]
        end
    end
    return seq
end

function algorithmic_cap(N::Int, alpha::Float64; gamma::Float64=0.0, map_to::Tuple{Int,Int}=(1,2))
    @assert N >= 1 "N must be >= 1"
    gamma = gamma - floor(gamma) # wrap into [0,1)
    alpha = alpha - floor(alpha)
    seq = Vector{Int}(undef, N)
    for n in 1:N
        term_1 = floor((n + 1) * alpha + gamma)
        term_2 = floor(n * alpha + gamma)
        out_of_bounds_correction = floor(alpha)
        diff = term_1 - term_2 #- out_of_bounds_correction
        if diff == 1
            seq[n] = map_to[1]
        else
            seq[n] = map_to[2]
        end
    end
    return seq
end

# ##############################################################################
# ## Usage 1: slopes and 1-slopes in interval [0,1] with consistent map_to=(2,1)
# ##############################################################################
# slope_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"
# N=500
# map_to = (2,1) # (2,1) makes slope=0.0 1,1,1,1,... and slope=1.0 2,2,2,2,...
# phason_angle_range = collect(0.0:0.01:1.0)

# indir = joinpath(@__DIR__, "sturm_grad_sets")
# infile = joinpath(indir, slope_file)
# println("Loading Sturmian slopes from ", infile)
# sturmian_slopes_df = load_sturmian_grad_bson(infile)

# rats, _ = add_rationals(10; start_r=2)
# rats_df = DataFrame(slope=rats, origin="rational", type="original", prefix=fill(missing, length(rats)))
# rats_df = add_extrema!(rats_df)

# # Add all rational slopes to sturmian_slopes_df
# df_rats_and_irrats = vcat(
#     sturmian_slopes_df,
#     rats_df;
#     cols = :union
# )


# seqs = Vector{Vector{Int}}()
# phis = Vector{Float64}()
# phasons = Vector{Float64}()
# for phason in phason_angle_range
#     for (i, row) in enumerate(eachrow(df_rats_and_irrats))
#         phi = row.slope
#         comp_phi = 1.0 - phi
#         seq = cut_and_project_sequence(N, phi; gamma=phason, map_to=map_to)
#         comp_seq = cut_and_project_sequence(N, comp_phi; gamma=phason, map_to=map_to)
#         push!(seqs, seq)
#         push!(seqs, comp_seq)
#         push!(phis, phi)
#         push!(phis, comp_phi)
#         push!(phasons, phason)
#         push!(phasons, phason)
#     end
# end

# outdir = joinpath(@__DIR__, "sturm_seq_sets")
# outfile = joinpath(outdir, slope_file[1:end-5] * "_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end])_const_mapping.bson")
# @save outfile phis seqs phasons
# println("Saved $(length(phis)) phis, $(length(seqs)) sequences to ", outfile)
# ##############################################################################



# ##############################################################################
# ## Usage 2: For specific single slope value: slopes and 1-slopes in interval [0,1] with consistent map_to=(2,1)
# ##############################################################################
# slope_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"
# N=500
# map_to = (2,1) # (2,1) makes slope=0.0 1,1,1,1,... and slope=1.0 2,2,2,2,...
# phason_angle_range = collect(0.0:0.01:1.0)

# target_phi = 0.61803398875  # Example: Golden Mean (1/phi)
# tol = 1e-4                  # Tolerance

# indir = joinpath(@__DIR__, "sturm_grad_sets")
# infile = joinpath(indir, slope_file)
# println("Loading Sturmian slopes from ", infile)
# sturmian_slopes_df = load_sturmian_grad_bson(infile)

# # Filter for the specific slope
# mask = abs.(sturmian_slopes_df.slope .- target_phi) .< tol
# if count(mask) == 0
#     error("No slope found within tolerance $tol of $target_phi")
# end

# # Take the first match (or the closest one)
# subset_df = sturmian_slopes_df[mask, :]
# best_match_idx = argmin(abs.(subset_df.slope .- target_phi))
# selected_phi = subset_df.slope[best_match_idx]

# println("Selected slope: $selected_phi (Target: $target_phi)")

# seqs = Vector{Vector{Int}}()
# phis = Vector{Float64}()
# phasons = Vector{Float64}()
# for phason in phason_angle_range
#     # Generate sequence for the selected phi
#     seq = cut_and_project_sequence(N, selected_phi; gamma=phason, map_to=map_to)
#     push!(seqs, seq)
#     push!(phis, selected_phi)
#     push!(phasons, phason)

#     # # Optional: Generate complementary sequence (1 - phi) if needed
#     # comp_phi = 1.0 - selected_phi
#     # comp_seq = cut_and_project_sequence(N, comp_phi; gamma=phason, map_to=map_to)
#     # push!(seqs, comp_seq)
#     # push!(phis, comp_phi)
#     # push!(phasons, phason)
# end


# outdir = joinpath(@__DIR__, "sturm_seq_sets")
# outfile = joinpath(outdir, slope_file[1:end-5] * "_single_slope_$(round(selected_phi, digits=5))_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end])_const_mapping.bson")
# @save outfile phis seqs phasons
# println("Saved $(length(phis)) sequences for slope $selected_phi to ", outfile)
# ##############################################################################




# ##############################################################################
# ## Usage 2: slopes and 1-slopes in interval (0,1) with inverted mapping for complementary seqs
# ##############################################################################
# slope_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"
# N=1000

# indir = joinpath(@__DIR__, "sturm_grad_sets")
# infile = joinpath(indir, slope_file)
# println("Loading Sturmian slopes from ", infile)
# sturmian_slopes_df = load_sturmian_grad_bson(infile)

# rats, _ = add_rationals(10; start_r=2)
# rats = add_extrema!(DataFrame(slope=rats))

# ## Add rationals
# df_rats_and_irrats = vcat(
#     sturmian_slopes_df,
#     rats;
#     cols = :union
# )

# phason_angle_range = [0.0] #collect(0.0:0.1:1.0)

# seqs = Vector{Vector{Int}}()
# phis = Vector{Float64}()
# phasons = Vector{Float64}()
# for phason in phason_angle_range
#     for (i, row) in enumerate(eachrow(df_rats_and_irrats))
#         phi = row.slope
#         seq = cut_and_project_sequence(N, phi; gamma=phason, map_to=(1,2))
#         push!(seqs, seq)
#         push!(phis, phi)
#         push!(phasons, phason)
#     end
# end

# ## second cut-and-project with inverted mapping (t2 <-> t1)
# for phason in phason_angle_range
#     for (i, row) in enumerate(eachrow(df_rats_and_irrats))
#         phi = row.slope
#         comp_phi = 1.0 - phi
#         seq = cut_and_project_sequence(N, comp_phi; gamma=phason, map_to=(2,1))
#         push!(seqs, seq)
#         push!(phis, comp_phi)
#         push!(phasons, phason)
#     end
# end

# outdir = joinpath(@__DIR__, "sturm_seq_sets")
# outfile = joinpath(outdir, slope_file[1:end-5] * "_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end]).bson")
# @save outfile phis seqs phasons
# println("Saved $(length(phis)) phis, $(length(seqs)) sequences to ", outfile)
# ##############################################################################


# ##############################################################################
# ## Usage 3: slopes in (0,1) U (1, ∞) U (-∞, 0) with consistent map_to=(1,2)
# ##############################################################################
# slope_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"
# N = 2000

# indir = joinpath(@__DIR__, "sturm_grad_sets")
# infile = joinpath(indir, slope_file)
# println("Loading Sturmian slopes from ", infile)
# sturmian_slopes_df = load_sturmian_grad_bson(infile)

# # Generate full slope DataFrame with labels
# full_slope_df = generate_full_slope_df(sturmian_slopes_df; n_rationals=10, start_r=2)

# phason_angle_range = [0.0]  # collect(0.0:0.1:1.0)

# seqs = Vector{Vector{Int}}()
# phis = Vector{Float64}()
# phasons = Vector{Float64}()
# origins = Vector{String}()
# types = Vector{String}()
# for phason in phason_angle_range
#     for row in eachrow(full_slope_df)
#         phi = row.slope
#         # seq = generalised_cut_and_project_sequence(N, phi; gamma=phason, map_to=(1,2))
#         seq = algorithmic_cap(N, phi; gamma=phason, map_to=(1,2))
#         push!(seqs, seq)
#         push!(phis, phi)
#         push!(phasons, phason)
#         push!(origins, row.origin)
#         push!(types, row.type)
#     end
# end

# # # Print first few lines of data for verification
# # println("First 5 phis: ", phis[1:min(5, length(phis))])
# # println("First 5 origins: ", origins[1:min(5, length(origins))])
# # println("First 5 types: ", types[1:min(5, length(types))])
# # println("First sequence (first 10 elements): ", seqs[1][1:min(10, length(seqs[1]))])
# # println("Total phis: $(length(phis)), seqs: $(length(seqs)), phasons: $(length(phasons)), origins: $(length(origins)), types: $(length(types))")

# full_df = DataFrame(phi=phis, seq=seqs, phason=phasons, origin=origins, type=types)
# # println(full_df)
# # println("all sequences in full: ", seqs)

# outdir = joinpath(@__DIR__, "sturm_seq_sets")
# outfile = joinpath(outdir, slope_file[1:end-5] * "_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end])_algorithmic_generalised.bson")
# @save outfile phis seqs phasons origins types
# println("Saved $(length(phis)) phis with phason range $(phason_angle_range). Overall $(length(seqs)) sequences added to ", outfile)
# ##############################################################################




##############################################################################
## Usage 4: Generate slopes as in Hofstadter butterfly (all reduced p/q, plus irrationals)
##############################################################################
N = 100
map_to = (2,1)
qmin, qmax = 2, N

# 1. Threaded generation of all reduced rationals p/q in [0,1]
slopes = Float64[]
slope_lock = SpinLock()
slope_prog = Progress(qmax - qmin + 1; desc="Generating slopes")

@threads for q in qmin:qmax
    local_chunk = Float64[]
    sizehint!(local_chunk, q)
    for p in 1:q-1
        if gcd(p, q) == 1
            push!(local_chunk, p / q)
        end
    end
    lock(slope_lock)
    append!(slopes, local_chunk)
    next!(slope_prog)
    unlock(slope_lock)
end

push!(slopes, 0.0, 1.0)
sort!(slopes)
unique!(slopes)

phason_angle_range = collect(0.0:0.01:1.0)
n_phi = length(slopes)
n_phason = length(phason_angle_range)
slot_count = n_phi * n_phason

seqs = Vector{Vector{Int}}(undef, slot_count)
phis = Vector{Float64}(undef, slot_count)
phasons = Vector{Float64}(undef, slot_count)

seq_lock = SpinLock()
seq_prog = Progress(n_phason; desc="Cut-and-project")

@threads for ip in 1:n_phason
    phason = phason_angle_range[ip]
    base = (ip - 1) * n_phi
    for is in 1:n_phi
        idx = base + is
        phi = slopes[is]
        seqs[idx] = cut_and_project_sequence(N, phi; gamma=phason, map_to=map_to)
        phis[idx] = phi
        phasons[idx] = phason
    end
    lock(seq_lock); next!(seq_prog); unlock(seq_lock)
end
println("Generated $(length(phis)) slopes.")

## Rebalance distribution (keep ≤ nperbin unique φ per bin, retain all phasons per φ)
nbins = 1000
nperbin = 1
edges = range(0.0, stop=1.0, length=nbins + 1)

unique_phis = unique(phis)
phi_to_inds = Dict{Float64, Vector{Int}}()
for (i, phi) in enumerate(phis)
    push!(get!(phi_to_inds, phi, Int[]), i)
end

bin_idx = clamp.(searchsortedfirst.(Ref(edges), unique_phis) .- 1, 1, nbins)
selected_unique = falses(length(unique_phis))
counts = zeros(Int, nbins)
for i in eachindex(unique_phis)
    b = bin_idx[i]
    if counts[b] < nperbin
        selected_unique[i] = true
        counts[b] += 1
    end
end

selected_phis = unique_phis[selected_unique]
selected_inds = reduce(vcat, [phi_to_inds[phi] for phi in selected_phis])
phis = phis[selected_inds]
seqs = seqs[selected_inds]
phasons = phasons[selected_inds]
println("Rebalanced to $(length(selected_phis)) unique slopes -> $(length(phis)) sequences across $(length(phason_angle_range)) phasons (nbins=$nbins, nperbin=$nperbin).")

println("Plotting distribution of generated slopes...")
plt = histogram(
    phis;
    bins = nbins,
    xlabel = "Slope φ",
    ylabel = "Count",
    title = "Density distribution of generated slopes (N=$N)",
    legend = false,
    normalize = false
)
savefig_dir = joinpath(@__DIR__, "slope_density_dists")
isdir(savefig_dir) || mkpath(savefig_dir)
savefig(plt, joinpath(savefig_dir, "hof_style_slopes_N$(N)_slope_distribution.png"))

println("Saving $(length(phis)) phis, $(length(seqs)) sequences...")
outdir = joinpath(@__DIR__, "sturm_seq_sets")
outfile = joinpath(outdir, "hof_style_slopes_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end])_nbins$(nbins)_npb$(nperbin).bson")
@save outfile phis seqs phasons
println("Saved $(length(phis)) phis, $(length(seqs)) sequences to ", outfile)
##############################################################################