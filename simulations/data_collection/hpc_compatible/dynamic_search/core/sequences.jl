module DynamicSequences

using BSON

export load_sequences_from_bson, select_sequence

function load_sequences_from_bson(sequence_bson_path::String)
    @assert isfile(sequence_bson_path) "BSON file not found: $sequence_bson_path"

    raw = BSON.load(sequence_bson_path)
    @assert haskey(raw, :seqs) "BSON missing :seqs"
    @assert haskey(raw, :phis) "BSON missing :phis"

    seqs = Vector{Vector{Int}}(raw[:seqs])
    phis = Float64.(raw[:phis])
    phasons = haskey(raw, :phasons) ? Float64.(raw[:phasons]) : fill(0.0, length(phis))

    length(seqs) == length(phis) == length(phasons) || error("Mismatched lengths in BSON arrays")

    return seqs, phis, phasons
end

function select_sequence(
    seqs::Vector{Vector{Int}},
    phis::Vector{Float64},
    phasons::Vector{Float64};
    target_slope::Union{Nothing,Float64}=nothing,
    slope_tolerance::Float64=1e-6,
    target_phason::Union{Nothing,Float64}=nothing,
    phason_tolerance::Float64=1e-6,
    selection_mode::Symbol=:closest
)
    idxs = collect(1:length(seqs))

    if target_slope !== nothing
        slope_mask = abs.(phis .- target_slope) .<= slope_tolerance
        idxs = idxs[slope_mask]
    end

    if target_phason !== nothing
        phason_mask = abs.(phasons .- target_phason) .<= phason_tolerance
        idxs = idxs[phason_mask]
    end

    if isempty(idxs)
        if selection_mode == :closest && target_slope !== nothing
            # Fall back to closest slope if no tolerance match.
            i = argmin(abs.(phis .- target_slope))
            return i, seqs[i], phis[i], phasons[i]
        end
        error("No sequence matches the requested slope/phason filter")
    end

    if selection_mode == :first
        i = first(idxs)
        return i, seqs[i], phis[i], phasons[i]
    elseif selection_mode == :closest && target_slope !== nothing
        local_idx = argmin(abs.(phis[idxs] .- target_slope))
        i = idxs[local_idx]
        return i, seqs[i], phis[i], phasons[i]
    else
        i = first(idxs)
        return i, seqs[i], phis[i], phasons[i]
    end
end

end # module
