"""
    file name:   main.jl
    created:     25/09/2025
    last edited: 25/10/2025

    # file name:   main.jl
    # created:     25/09/2025
    # last edited: 30/09/2025

    # overview:
    #     The main script to execute simulations on a local machine. Define parameter ranges and decide on calculation and datasaving requirements to call on the correct functions.
    
    # structure:
    #     - Sec 1:  Parameter Choice
    #                 Generate all sequence types and define parameter ranges for N, t_n, mu and Delta
    #     - Sec 2:  Data Save Path  
    #                 Ddefine data save filepath
    #     - Sec 3:  Tailoring
    #                 Answer list of options relating to exactly what values should be calculated and saved, in what precision and using which solver type.
    #     - Sec 4:  Cut Parameters
    #                 (Optional) Define cutting rules to restrict parameter space for :restricted solver type
    #     - Sec 5:  Run
    #                 Call on the dispatch function to run the selected solver with all of the above parameters and options
    
    usage instructions:
        1) 
            a) If running in VSCode or other persistent Julia REPL: Ensure the local_machine environment is initialised by running simulations/data_collection/local_machine/__init__.jl ONCE per Julia session. Comment OUT the inclusion preamble below.
            b) If running directly from terminal command line ensure the inclusion preamble below is commented IN. Do not use __init__.jl
        2) Complete Sec 1 to define your desired parameter ranges
        3) Check Sec 2 to ensure the data save path is correct (it will typically self-generate a folder name based on the parameters chosen)
        4) Complete Sec 3 to choose what calculations to perform, in what precision and using which solver type (for more details on solver types see the lead comments in solvers.jl)
        (5) (Optional) Complete Sec 4 to define cutting rules to restrict parameter space for :restricted solver type for less expensive 2D param space solving.)
        6) Run this script with the code in Sec 5 to execute the simulations and save the data
    
    Latest edits:
        - Added sequence_id tracking and saving throughout the simulation pipeline to allow identification of different sequence chunks in the output data. (This is needed for Sturmian sequence sweeps.)
"""

project_root = @__DIR__
######################################################################################################################
# # Leave this in if running main.jl directly from the terminal commandline 
# # Comment this out and use __init__.jl if running in VSCode
include(joinpath(project_root, "../../data_collection/hpc_compatible/serial_batch/solvers.jl"))
include(joinpath(project_root, "../../data_collection/hpc_compatible/serial_batch/compute.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/sequence_gen.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/param_comb_gen.jl"))
######################################################################################################################

using .SeqGen
using .ParamCombGen
using .HPCSolv: UserOptions
using BSON: @save, @load
using BSON
using DelimitedFiles

row_index = parse(Int64, ARGS[1]);


###########################################################
################# Sec 1: Parameter Choice #################
###########################################################
"""
    Read one parameter row from the latest param_prep .dat file
    .dat columns (tab-separated):
    Ns::Vector{Int}
    mus::Vector{Float64}
    Deltas::Vector{Float64}
    t1s::Vector{Float64}
    t2s::Vector{Float64}
    t3s::Vector{Float64}
    sequence_names::Vector{String}
    sequence_ids::Vector{Tuple{Float64,Float64}}
    sequences::Vector{Vector{Int}}
"""
function read_parameters(params_path::String, row_index::Int)
    raw = DelimitedFiles.readdlm(params_path, '\t', String; header=true)
    row = raw[1][row_index, :]

    # parse Julia literals safely (vectors/tuples)
    parse_lit(S::AbstractString) = Base.invokelatest(eval, Meta.parse(S))

    Ns              = Vector{Int}(parse_lit(row[1]))
    mus             = Vector{Float64}(parse_lit(row[2]))
    Deltas          = Vector{Float64}(parse_lit(row[3]))
    t1s             = Vector{Float64}(parse_lit(row[4]))
    t2s             = Vector{Float64}(parse_lit(row[5]))
    t3s             = Vector{Float64}(parse_lit(row[6]))
    sequence_nms    = Vector{String}(parse_lit(row[7]))
    seq_ids         = Vector{Tuple{Float64,Float64}}(parse_lit(row[8]))
    seqs            = Vector{Vector{Int}}(parse_lit(row[9]))

    return (
        Ns = Ns, mus = mus, Deltas = Deltas,
        t1s = t1s, t2s = t2s, t3s = t3s,
        sequence_names = sequence_nms, sequence_ids = seq_ids, sequences = seqs
    )
end

# Choose the .dat file and which row to run
params_filename = "params_sturmian_slopes_K8_L3_balanced_bins500_mpb1_r50_comp-false_N1000_phason_0.0-1-0.0_Ns_500-500-1_mus_0-3-601_Deltas_0.0-0.2-5_t1_1-1-1_t2_1p5-2p5-3_t3_none_nseq_1.dat"
params_dat_path = joinpath(project_root, "batch_params", "param_sets", params_filename)

p = read_parameters(params_dat_path, row_index)

# Map parsed params to variables used by run_selected_solver
N_range = p.Ns
mu_range = p.mus
Delta_range = p.Deltas

t1_range = p.t1s
t2_range = p.t2s
t3_range = p.t3s
t_ranges = isempty(t3_range) ? [t1_range, t2_range] : [t1_range, t2_range, t3_range]
t_combinations = ParamCombGen.t_ranges_combs(t_ranges)

sequences = [p.sequence]             # Vector{Vector{Int}}
sequence_ids = [p.sequence_id]          # Vector{Tuple{Float64,Float64}}
sequence_name = p.sequence_name

## Set Precision for hp calculations
BigFloat_precision = 512

## Set chunk size for data saving
chunk_size = 1000


###########################################################
################# Sec 2: Cut Parameters ###################
###########################################################

# This is specific to a 2-phopping system, adjust as needed for other cases
rho_min = minimum(t2_range) / maximum(t1_range)
rho_max = maximum(t2_range) / minimum(t1_range)
rho_step = (rho_max - rho_min) / length(t2_range)
rho_range = collect(range(rho_min, rho_max, length(t2_range)))

xs = mu_range
ys = rho_range
grad1 = ParamCombGen.angle_to_gradient(35.0)
grad2 = ParamCombGen.angle_to_gradient(45.0)

include(joinpath(project_root, "../../data_collection/auxilliary/param_restriction_cuts.jl"))
# cuts = TMQC_D20_cuts # see auxilliary/param_restriction_cuts.jl for predefined cutting rules
cuts = [
    Dict(
        :gradient => 0.0,
        :y_intercept => -0.1,
        :x_range => (2.75, maximum(xs)),
        :y_range => (minimum(ys), maximum(ys)),
        :cut_which_side => "above"
    )
]

mu_rho_restricted = ParamCombGen.restrict_range(xs, ys, cuts)



###########################################################
################ Sec 3: Sequence Scaling ##################
###########################################################

## Allow either t_1 or rho to be preserved when scaling hopping for each sequences
preserve_symbol = :rho  # options: :t1, :rho
rho_target = 1.5    # target rho value when preserve_symbol = :rho
tbar = 1.5         # target t_1 value when preserve_symbol = :t1




###########################################################
#################### Sec 4: Tailoring #####################
###########################################################


function get_user_options()
    return UserOptions(
        true,    # calc_mp
        false,   # calc_ipr
        false,   # calc_mbs_energy_gap
        false,   # calc_loc_len
        :np,     # calc_precision: :hp, :np
        :none, # save_evecs: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :all_np, # save_evals: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :generic # solver_type: :generic, seq_scaled, :mu_loop, :N_loop, :restricted
    )
end


###########################################################
################## Sec 5: Data Save Path ##################
###########################################################

# project_name = "sturmian_sweep_t1-t2_swap"

# root_path = joinpath(project_root, "../../../simulations/raw_data")
# if length(t_ranges) < 3
#     folder_name = "np/$project_name/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t1_range))_t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))/"
# else
#     folder_name = "np/$project_name/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t1_range))_t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_t3($(t3_range[1])-$(t3_range[end])-$(length(t3_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))/"
# end
# datasave_path = "$(root_path)/$(folder_name)/"

# # Create the folder if it doesn't exist
# datasave_path = normpath(datasave_path)
# isdir(datasave_path) || mkpath(datasave_path)
# println("Data will be saved to: $(datasave_path)")

float3(x) = replace(@sprintf("%.6g", x), "." => "p")

function range_info(vec)
    isempty(vec) && return "none"
    if eltype(vec) <: Integer
        return "$(vec[1])-$(vec[end])-$(length(vec))"
    else
        return string(float3(vec[1]), "-", float3(vec[end]), "-", length(vec))
    end
end

# Decide precision label from options
opts = get_user_options()
precision_label = opts.calc_precision == :hp ? "hp" :
                  opts.calc_precision == :np ? "np" : "arpack"

project_name = "sturmian_sweep_t1-t2_swap"

# Save under serial_batch/results/
root_path = joinpath(project_root, "results", precision_label, project_name)

n_info   = range_info(N_range)
t1_info  = range_info(t1_range)
t2_info  = range_info(t2_range)
t3_info  = range_info(t3_range)
mu_info  = range_info(mu_range)
d_info   = range_info(Delta_range)

if isempty(t3_range)
    folder_name = "$(sequence_name)_N($n_info)_t1($t1_info)_t2($t2_info)_mu($mu_info)_Delta($d_info)"
else
    folder_name = "$(sequence_name)_N($n_info)_t1($t1_info)_t2($t2_info)_t3($t3_info)_mu($mu_info)_Delta($d_info)"
end

datasave_path = normpath(joinpath(root_path, folder_name))
isdir(datasave_path) || mkpath(datasave_path)
println("Data will be saved to: $(datasave_path)")




###########################################################
####################### Sec 6: Run ########################
###########################################################

# ------------------------------
# --------- Single run ---------
opts = get_user_options()
run_selected_solver(
    opts, 
    N_range, 
    t_combinations, 
    mu_range, 
    Delta_range, 
    sequences, 
    sequence_name, 
    datasave_path; 
    precision=BigFloat_precision, 
    chunk_size=chunk_size, 
    param_restrictions=mu_rho_restricted, 
    sequence_ids=sequence_ids,
    rho_target=rho_target,
    tbar=tbar,
    preserve_symbol=preserve_symbol
)