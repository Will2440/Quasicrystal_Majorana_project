module DynamicIO

using CSV
using DataFrames
using Dates
using JLD2

export now_tag,
       ensure_dir,
       logmsg,
       write_points_csv,
       read_worker_results,
       write_worker_results,
       write_jld2,
       read_config

function now_tag()
    return Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
end

function ensure_dir(path::String)
    isdir(path) || mkpath(path)
    return path
end

function logmsg(msg::AbstractString)
    ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    println("[$ts] $msg")
end

function write_points_csv(points, out_path::String)
    df = DataFrame(
        point_key = String[],
        rho = Float64[],
        u = Float64[],
    )

    for p in points
        push!(df, (p.point_key, p.rho, p.u))
    end

    CSV.write(out_path, df)
    return out_path
end

function write_worker_results(rows, out_path::String)
    df = DataFrame(
        point_key = String[],
        rho = Float64[],
        u = Float64[],
        mu = Float64[],
        t1 = Float64[],
        t2 = Float64[],
        phase_label = Int[],
        topo_fraction = Float64[],
        invariant_mean = Float64[],
        invariant_min = Float64[],
        invariant_max = Float64[],
        gap_min = Float64[],
        gap_mean = Float64[],
        gap_max = Float64[],
        n_samples = Int[],
    )

    for r in rows
        push!(df, (
            r.point_key,
            r.rho,
            r.u,
            r.mu,
            r.t1,
            r.t2,
            r.phase_label,
            r.topo_fraction,
            r.invariant_mean,
            r.invariant_min,
            r.invariant_max,
            r.gap_min,
            r.gap_mean,
            r.gap_max,
            r.n_samples,
        ))
    end

    CSV.write(out_path, df)
    return out_path
end

function write_jld2(out_path::String; kwargs...)
    jldsave(out_path; kwargs...)
    return out_path
end

function read_worker_results(results_dir::String)
    files = filter(f -> endswith(f, ".csv") && occursin("worker_results_", f), readdir(results_dir; join=true))
    if isempty(files)
        return DataFrame(
            point_key = String[], rho = Float64[], u = Float64[], mu = Float64[],
            t1 = Float64[], t2 = Float64[], phase_label = Int[], topo_fraction = Float64[],
            invariant_mean = Float64[], invariant_min = Float64[], invariant_max = Float64[],
            gap_min = Float64[], gap_mean = Float64[], gap_max = Float64[], n_samples = Int[]
        )
    end

    dfs = [CSV.read(f, DataFrame) for f in files]
    return vcat(dfs...)
end

function read_config(config_path::String)
    cfg = include(config_path)
    isa(cfg, Dict{Symbol,Any}) || error("Config file must return Dict{Symbol,Any}")
    return cfg
end

end # module
