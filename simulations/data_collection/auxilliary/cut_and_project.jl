using BSON: @save, @load
using BSON
using DataFrames

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


## Usage
slope_file = "sturmian_slopes_K12_L6_balanced_bins500_mpb1.bson"

indir = joinpath(@__DIR__, "sturm_grad_sets")
infile = joinpath(indir, slope_file)
println("Loading Sturmian slopes from ", infile)
sturmian_slopes_df = load_sturmian_grad_bson(infile)
N = 500

seqs = Vector{Vector{Int}}()
phis = Vector{Float64}()
for (i, row) in enumerate(eachrow(sturmian_slopes_df))
    phi = row.slope
    seq = cut_and_project_sequence(N, phi)
    push!(seqs, seq)
    push!(phis, phi)
end

outdir = joinpath(@__DIR__, "sturm_seq_sets")
outfile = joinpath(outdir, slope_file[1:end-5] * "_N$(N).bson")
@save outfile phis seqs
println("Saved $(length(phis)) phis and $(length(seqs)) sequences to ", outfile)