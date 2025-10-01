"""
    file name:   bson_unpacker.jl
    created:     24/09/2025
    last edited: 25/09/2025

    overview:


    structure:
        - Sec 1:  Auxilliary Functions
                    Contains helpful functions which unpack the .bson and do minor data processing
        - Sec 2:  Main Unpacker Function
                    Contains the main function, calling on auxilliary functions.
       
    
    structure:
        - Section 1: .bson upacker into dataframe -- ready to use on the data fodlers provided
        - Section 2: Generalised plotting functions, can be used individualy, or within the following function in Sec 3.
        - Section 3: Additional code to calculate the abundance integral and ssave this data
        - Section 4: Contains a single function which can be run on the dataframe created in Sec 1 to generate all the standard plots seen in the data folders
"""


using BSON
using DataFrames
using Glob

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
        mp = Float64[],
        maj_gap = Float64[],
        ipr = Float64[],
        loc_len = Float64[],
        eigenvalues = Union{Vector{Float64}, Missing}[],
        eigenvectors = Union{Matrix{Float64}, Missing}[]
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

#  # Example Usage
folder_path = "../../raw_data/np/all_crystal_grad_testruns/mu_vs_rho_mp_heatmaps/PQC_N(50-50-1)_t1(1.0-1.0-101__t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(0.0-2.0-21)"
mp_tol = 0.1
df = unpack_bason_standard(folder_path; mp_tol=mp_tol)

