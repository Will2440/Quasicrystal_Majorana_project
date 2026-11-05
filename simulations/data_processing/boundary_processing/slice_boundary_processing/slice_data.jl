module SliceBoundaryData

using CSV
using DataFrames

export load_transitions, load_slice_summary, load_samples, available_slices

function require_run_dir(run_dir::AbstractString)
    isdir(run_dir) || error("Run directory not found: $run_dir")
    return String(run_dir)
end

function load_transitions(run_dir::AbstractString)
    path = joinpath(require_run_dir(run_dir), "transition_boundaries.csv")
    isfile(path) || error("Transition data not found: $path")
    df = CSV.read(path, DataFrame)
    required = (:slice_index, :rho, :track_id, :mu, :mu_lo, :mu_hi, :phase_lo, :phase_hi)
    all(name -> name in propertynames(df), required) || error("Unexpected transition table columns in $path")
    sort!(df, [:rho, :mu])
    return df
end

function load_slice_summary(run_dir::AbstractString)
    path = joinpath(require_run_dir(run_dir), "slice_summary.csv")
    isfile(path) || error("Slice summary not found: $path")
    df = CSV.read(path, DataFrame)
    required = (:slice_index, :rho, :gap_count, :n_points)
    all(name -> name in propertynames(df), required) || error("Unexpected slice summary columns in $path")
    sort!(df, [:rho])
    return df
end

function available_slices(run_dir::AbstractString)
    slices_dir = joinpath(require_run_dir(run_dir), "slices")
    isdir(slices_dir) || error("Slices directory not found: $slices_dir")
    found = Tuple{Int,String}[]
    for path in readdir(slices_dir; join=true)
        isdir(path) || continue
        matched = match(r"^slice_(\d+)$", basename(path))
        matched === nothing && continue
        push!(found, (parse(Int, matched.captures[1]), path))
    end
    sort!(found; by=first)
    return found
end

function load_samples(run_dir::AbstractString; max_slice_index::Union{Nothing,Int}=nothing)
    frames = DataFrame[]
    for (slice_index, path) in available_slices(run_dir)
        max_slice_index !== nothing && slice_index > max_slice_index && continue
        results_path = joinpath(path, "point_results.csv")
        isfile(results_path) || continue
        frame = CSV.read(results_path, DataFrame)
        isempty(frame) && continue
        if !(:slice_index in propertynames(frame))
            frame.slice_index = fill(slice_index, nrow(frame))
        end
        push!(frames, frame)
    end
    isempty(frames) && return DataFrame(
        slice_index=Int[], rho=Float64[], mu=Float64[], u=Float64[],
        phase_label=Int[], topo_fraction=Float64[], gap_min=Float64[])
    return vcat(frames...; cols=:union)
end

end # module