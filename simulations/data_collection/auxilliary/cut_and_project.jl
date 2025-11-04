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

function add_rationals(n::Int; start_r::Int=1)
    @assert n >= 1 "n must be >= 1"
    @assert start_r >= 1 "start_r must be >= 1"

    rs = start_r:(start_r + n - 1)
    rats = [1.0 / Float64(r) for r in rs]
    comps = [1.0 - r for r in rats]
    return rats, comps
end

function add_extrema!(df::DataFrame)
    min_slope = 0.0
    max_slope = 1.0

    # add these slope values to the DataFrame
    slopes = collect(df.slope)
    new_slopes = Float64[]
    for s in (min_slope, max_slope)
        if !(s in slopes)
            push!(new_slopes, s)
        end
    end
    return df_ext
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
slope_file = "sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false.bson"

indir = joinpath(@__DIR__, "sturm_grad_sets")
infile = joinpath(indir, slope_file)
println("Loading Sturmian slopes from ", infile)
sturmian_slopes_df = load_sturmian_grad_bson(infile)

rats, comps = add_rationals(10; start_r=2)

## Add all rational slopes to sturmian_slopes_df straight away; without t1<->t2 inversion
# rational_slopes = vcat(rats, comps)
# df_rats_and_irrats = vcat(
#     sturmian_slopes_df,
#     DataFrame(slope = rational_slopes);
#     cols = :union
# )

## Add rationals to one df and comps to another
df_rats_and_irrats = vcat(
    sturmian_slopes_df,
    DataFrame(slope = rats);
    cols = :union
)

# df_complimentaries = vcat(
#     sturmian_slopes_df,
#     DataFrame(slope = comps);
#     cols = :union
# )


N = 1000

phason_angle_range = [0.0] #collect(0.0:0.1:1.0)

seqs = Vector{Vector{Int}}()
phis = Vector{Float64}()
phasons = Vector{Float64}()
for phason in phason_angle_range
    for (i, row) in enumerate(eachrow(df_rats_and_irrats))
        phi = row.slope
        seq = cut_and_project_sequence(N, phi; gamma=phason, map_to=(1,2))
        push!(seqs, seq)
        push!(phis, phi)
        push!(phasons, phason)
    end
end

## second cut-and-project with inverted mapping (t2 <-> t1)
for phason in phason_angle_range
    for (i, row) in enumerate(eachrow(df_rats_and_irrats))
        phi = row.slope
        comp_phi = 1.0 - phi
        seq = cut_and_project_sequence(N, comp_phi; gamma=phason, map_to=(2,1))
        push!(seqs, seq)
        push!(phis, comp_phi)
        push!(phasons, phason)
    end
end

outdir = joinpath(@__DIR__, "sturm_seq_sets")
outfile = joinpath(outdir, slope_file[1:end-5] * "_N$(N)_phason_$(phason_angle_range[1])-$(length(phason_angle_range))-$(phason_angle_range[end]).bson")
@save outfile phis seqs phasons
println("Saved $(length(phis)) phis, $(length(seqs)) sequences to ", outfile)