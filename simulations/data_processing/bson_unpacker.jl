"""
    file name:   bson_unpacker.jl
    created:     24/09/2025
    last edited: 25/10/2025

    overview:


    structure:
        - Sec 1:  Auxilliary Functions
                    Contains helpful functions which unpack the .bson and do minor data processing
        - Sec 2:  Main Unpacker Function
                    Contains the main function, calling on auxilliary functions.
       
    
    structure:
        - Section 1: .bson unpacker into dataframe -- ready to use on the data folders provided
        - Section 2: Generalised plotting functions, can be used individually, or within the following function in Sec 3.
        - Section 3: Additional code to calculate the abundance integral and save this data
        - Section 4: Contains a single function which can be run on the dataframe created in Sec 1 to generate all the standard plots seen in the data folders

    latest edits:
        - Added sequence_id tracking in the .bson converter
        - Added specialised unpacking sequence for Hofstadter model
"""


using BSON
using DataFrames
using Glob
using Logging

###########################################################
############## Sec 1: Auxilliary Functions ################
###########################################################

function normalise_to_t(t_n::Vector{Float64})
    norm = t_n[1]
    return norm
end

function bson_to_dataframe(file_path::String)

    bson_data = BSON.load(file_path)
    dataframe = DataFrame(bson_data)
    
    return dataframe
end

function process_bson_files(folder_path::String)
    combined_dataframe = DataFrame(
        N = Int[],
        t_n = Vector{Float64}[],
        mu = Float64[],
        Delta = Float64[],
        sequence_name = String[],
        sequence_id = Any[],
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Vector{BigFloat}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Matrix{BigFloat}, Missing}[]
    )
    
    file_paths = glob("*.bson", folder_path)
    files_read = 0
    for file_path in file_paths
        bson_data = BSON.load(file_path)
        if haskey(bson_data, :results_df)
            # Extract the DataFrame stored in the `results_df` key
            results_df = DataFrame(bson_data[:results_df])
            
            append!(combined_dataframe, results_df)
            files_read += 1
        end
    end
    println("number of files read: $files_read")
    
    return combined_dataframe
end

function calc_norms_df!(df::DataFrame)
    df.rho = [row.t_n[2] / row.t_n[1] for row in eachrow(df)]
    df.mu_t = [row.mu / row.t_n[1] for row in eachrow(df)]
    df.Delta_t = [row.Delta / row.t_n[1] for row in eachrow(df)]

    if length(df.t_n[1]) == 3
        df.sigma = [row.t_n[3] / row.t_n[1] for row in eachrow(df)]
    else
        nothing
    end

    return df
end

function calc_mp_disc!(df::DataFrame, tol::Float64)
    df.mp_disc = [isapprox(row.mp, -1.0, atol=tol) ? -1.0 : 0.0 for row in eachrow(df)]
    return df
end

function mask_df!(df::DataFrame)
    df.maj_gap_masked = [row.mp_disc == -1.0 ? row.maj_gap : 0.0 for row in eachrow(df)]
    df.ipr_masked = [row.mp_disc == -1.0 ? row.ipr : 0.0 for row in eachrow(df)]
    df.loc_len_masked = [row.mp_disc == -1.0 ? row.loc_len : 0.0 for row in eachrow(df)]
    return df
end


function unpack_from_inspect(folder_path::String)
    files = glob("*.bson", folder_path)
    println("unpack_from_inspect: found $(length(files)) .bson files in ", folder_path)
    if isempty(files)
        return DataFrame()
    end

    dfs = DataFrame[]
    files_used = String[]

    for f in files
        println(" loading: ", f)
        raw = nothing
        try
            raw = BSON.load(f)
        catch e
            @warn "BSON.load failed" file=f exception=e
            continue
        end

        println("  top-level keys: ", collect(keys(raw)))

        # try common keys first (symbol and string)
        df = nothing
        for trykey in (:results_df, "results_df")
            if haskey(raw, trykey)
                val = raw[trykey]
                if val isa DataFrame
                    df = DataFrame(val)
                    println("   -> used key: ", trykey, " (DataFrame), rows=", nrow(df))
                    break
                elseif isa(val, AbstractVector)
                    # maybe a Vector{NamedTuple} or Vector{Dict} or vectorized rows
                    try
                        df = DataFrame(val)
                        println("   -> used key: ", trykey, " (Vector -> DataFrame), rows=", nrow(df))
                        break
                    catch e
                        @debug "failed to convert value for key" file=f key=trykey exception=e
                    end
                else
                    @debug "key found but value not table-like" file=f key=trykey typeof=typeof(val)
                end
            end
        end

        # fallback 1: first value that is a DataFrame
        if df === nothing
            for (k,v) in pairs(raw)
                if v isa DataFrame
                    df = DataFrame(v)
                    println("   -> using first DataFrame found under key: ", k, " rows=", nrow(df))
                    break
                end
            end
        end

        # fallback 2: single top-level vector of NamedTuples/Dicts
        if df === nothing && length(raw) == 1
            v = first(values(raw))
            if isa(v, AbstractVector) && !isempty(v) && (v[1] isa NamedTuple || v[1] isa Dict)
                try
                    df = DataFrame(v)
                    println("   -> converted single top-level vector to DataFrame, rows=", nrow(df))
                catch e
                    @debug "failed to build DataFrame from single top-level vector" file=f exception=e
                end
            end
        end

        # fallback 3: pick the largest group of vector-like entries with same length
        if df === nothing
            vec_keys = Dict{Any,Int}()
            for (k,v) in pairs(raw)
                if v isa AbstractVector
                    vec_keys[k] = length(v)
                end
            end
            if !isempty(vec_keys)
                # choose most-common non-zero length
                len_counts = Dict{Int,Int}()
                for l in values(vec_keys)
                    len_counts[l] = get(len_counts, l, 0) + 1
                end
                candidates = filter(x -> x[1] > 0, collect(len_counts))
                if !isempty(candidates)
                    chosen_len = sort(candidates, by = x -> -x[2])[1][1]
                    selected = [k for (k,l) in vec_keys if l == chosen_len]
                    tbl = Dict{Symbol,Any}()
                    for k in selected
                        tbl[Symbol(string(k))] = raw[k]
                    end
                    try
                        df = DataFrame(tbl)
                        println("   -> built DataFrame from selected vector-like keys, rows=", nrow(df))
                    catch e
                        @debug "failed to build DataFrame from selected keys" file=f keys=selected exception=e
                    end
                end
            end
        end

        if df === nothing
            @warn "Could not extract table from BSON file; skipping" file=f keys=collect(keys(raw))
            continue
        end

        push!(dfs, df)
        push!(files_used, string(f))
    end

    if isempty(dfs)
        @warn "No BSON file produced a DataFrame; returning empty DataFrame"
        return DataFrame()
    end

    combined = vcat(dfs...; cols = :union)
    println("unpack_from_inspect: used $(length(files_used)) files; combined rows=$(nrow(combined)) cols=$(ncol(combined))")
    return combined
end

function discretise_evals!(
    df::DataFrame, 
    target_value::Float64, 
    threshold::Float64
)
    """
    For each row take the two middle-indexed eigenvalues (middle, middle+1),
    discretise each w.r.t. target_value and store a Vector{Int} of length 2
    in column :disc_eigenvalues.
    """
    n = nrow(df)
    disc_col = Vector{Vector{Int}}(undef, n)

    for (i, row) in enumerate(eachrow(df))
        vals = get(row, :eigenvalues, nothing)
        if vals === nothing || length(vals) == 0
            disc_col[i] = [0, 0]
            continue
        end

        v = real.(vals)
        lenv = length(v)
        mid = lenv ÷ 2
        idx1 = clamp(mid, 1, lenv)
        idx2 = clamp(mid + 1, 1, lenv)

        d1 = abs(v[idx1] - target_value) < threshold ? 1 : 0
        d2 = abs(v[idx2] - target_value) < threshold ? 1 : 0

        disc_col[i] = [d1, d2]
    end

    df.disc_eigenvalues = disc_col
    return df
end

function discretise_mp!(
    df::DataFrame,
    target_value::Float64,
    threshold::Float64
)
    n = nrow(df)
    disc_col = Vector{Int}(undef, n)

    for (i, row) in enumerate(eachrow(df))
        mp_val = get(row, :mp, nothing)
        if mp_val === nothing
            disc_col[i] = 0
            continue
        end

        v = real(mp_val)
        disc_col[i] = abs(v - target_value) < threshold ? 1 : 0
    end

    df.disc_mp = disc_col
    return df
end

function inspect_bson_folder(folder_path::AbstractString)
    files = glob("*.bson", folder_path)
    if isempty(files)
        println("No .bson files found in: ", folder_path)
        return
    end

    for f in files
        println("=== ", f, " ===")
        raw = nothing
        try
            raw = BSON.load(f)
        catch e
            println("  ERROR loading BSON: ", e)
            continue
        end

        println("  top-level keys: ", collect(keys(raw)))
        for (k,v) in pairs(raw)
            # key name
            kn = string(k)
            vt = typeof(v)
            # summary for common types
            if v isa DataFrame
                df = v
                println("  key = $kn => DataFrame; nrow=$(nrow(df)), ncol=$(ncol(df))")
                println("    columns:")
                for col in names(df)
                    coldata = df[!, col]
                    println("      - $(col): eltype=$(eltype(coldata)), typeof column=$(typeof(coldata))")
                end
                println("    first rows:")
                show(first(df, min(5, nrow(df)))) ; println()
            elseif v isa AbstractVector
                len = length(v)
                eltp = if len>0 try eltype(v) catch _ typeof(v[1]) end else :Empty end
                println("  key = $kn => Vector (len=$len), eltype=$eltp, typeof=$(vt)")
                if len>0
                    # print sample of first up to 5 elements
                    println("    sample:", v[1:min(5,len)])
                end
            elseif v isa Dict
                println("  key = $kn => Dict (nkeys=$(length(keys(v))))")
            else
                println("  key = $kn => $(vt)")
            end
        end
        println()
    end
    return nothing
end


###########################################################
############# Sec 2: Main Unpacker Function ###############
###########################################################

function unpack_bason_standard(folder_path::String; mp_tol=1e-5)
    df = process_bson_files(folder_path)
    df = calc_norms_df!(df)
    df = calc_mp_disc!(df, mp_tol)
    df = mask_df!(df)

    println("DataFrame keynames: $(names(df))")
    return df
end

function unpack_bson_hofstadter(folder_path::String; mp_targ=-1.0, mp_tol=1e-1, eval_targ=0.0, eval_tol=1e-2)
    df = process_bson_files(folder_path)
    df = calc_norms_df!(df)
    df = discretise_evals!(df, eval_targ, eval_tol)
    df = discretise_mp!(df, mp_targ, mp_tol)

    println("DataFrame keynames: $(names(df))")
    return df
end

# # Example Usage
# folder_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/hp/abundance_data/SQC_N(50-50-1)_t1(1.0-1.0-1_t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(0.1-2.0-20)"
# mp_tol = 0.1
# df = unpack_bason_standard(folder_path; mp_tol=mp_tol)

# println("First 5 eigenvalues and their types:")
# for i in 1:min(5, nrow(df))
#     println("eigenvalues[$i]: ", df.eigenvalues[i], " (type: ", typeof(df.eigenvalues[i]), ")")
# end

# data_folder = "Sturmian_K12_L6_balanced_bins500_mpb1_N(100-100-1)_t1(1.0-1.0-1_t2(1.5-1.5-1)_mu(0.0-3.0-151)_Delta(0.05-0.05-1)"
# folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "simulations", "raw_data", "np", "sturmian_sweep", data_folder))
# println("Unpacking Hofstadter data from folder: ", folder_path_hof)

# df = unpack_bson_hofstadter(folder_path_hof)


