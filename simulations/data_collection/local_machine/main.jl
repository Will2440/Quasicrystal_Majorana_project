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
include(joinpath(project_root, "../../data_collection/local_machine/solvers.jl"))
include(joinpath(project_root, "../../data_collection/local_machine/compute.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/sequence_gen.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/param_comb_gen.jl"))
######################################################################################################################

using .SeqGen
using .ParamCombGen
using .LocalSolv: UserOptions
# using BSON: @save, @load
using JLD2: @save, @load
using Plots
using Dates


###########################################################
################# Sec 1: Parameter Choice #################
###########################################################

## Sequence Generation (ensure seed is large enough for desired N_range)
N_seq_seed = 1024
normal_sequence = SeqGen.normal_SeqGen(N_seq_seed)
golden_sequence = SeqGen.golden_SeqGen(N_seq_seed)
silver_sequence = SeqGen.silver_SeqGen(N_seq_seed)
thue_morse_sequence = SeqGen.thue_morse_SeqGen(N_seq_seed)
plastic_sequence = SeqGen.plastic_SeqGen(N_seq_seed)

## Recover Sturmian sequences from .bson

## all slopes
# seq_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/local_mp_seq_sets/slopes_N1000_phason_0.61803-1-0.61803_nbins100_npb1_mapto(2, 1).bson"

## only fibonacci
seq_file = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/local_mp_seq_sets/slopes_N1000_phason_0.61803-1-0.61803_nbins100_npb1_mapto(2, 1).bson"
indir = joinpath(project_root, "..", "auxilliary", "local_mp_seq_sets", seq_file)
indir = normpath(indir)
sturm_df = SeqGen.load_sturmian_seq_bson(indir)

## N Range (even single value must be a Vector type)
N_range = [500] #collect(round.(Int, logrange(100, 3000, length=10))) # [200]

## t_n Range (combine any number of different hopping ranges)
t1_range = collect(range(1.0, 1.0, 1))
t2_range = collect(range(1.5, 1.5, 1))
t3_range = [2.0, 3.0, 4.0]
t_ranges = [t1_range, t2_range]#, t3_range]
t_combinations = ParamCombGen.t_ranges_combs(t_ranges)

## mu Range (even single value must be a Vector type)
mu_range = collect(range(0.0, 3.0, length=1201))

## Delta Range (even single value must be a vector type)
Delta_range = collect(range(0.01, 0.1, 10))


## Sequence Chunking (use this to generate sequence samples which can be compared)
# chunk_length = 100
# seq_start_indicies = [1]
# sequences = ParamCombGen.split_sequence(normal_sequence, chunk_length, seq_start_indicies)
# if maximum(N_range) > chunk_length
#     error("The largest number of sites (N = $(maximum(N_range))) is incompatible with the chosen sequence chunk sizes ($chunk_length)")
# end


# # Sequence Definition (choose which standard hopping sequence to use)
# sequences = [golden_sequence]
# sequence_name = "GQC"
# sequence_ids = [(2/(1 + sqrt(5)), 0.0)]

raw_seqs = collect(sturm_df.seq)
sequences = [Vector{Int}(s) for s in raw_seqs]
# sequence_name = String(seq_file[1:end-5])  # remove .bson
sequence_name = String(split(seq_file, '/')[end][1:end-5])
# raw_seq_ids = collect(sturm_df.phi)
# sequence_ids = [Float64(id) for id in raw_seq_ids]
phis = collect(sturm_df.phi)
phasons = collect(sturm_df.phason)
sequence_ids = [(phi, phason) for (phi, phason) in zip(phis, phasons)]
# println("Sequence ids (phi, phason): ", sequence_ids)


## Set Precision for hp calculations
BigFloat_precision = 128


## Set chunk size for data saving
chunk_size = 1000

## Set eta and omega_range for spectral function calculation
eta = 0.1
omega_range = collect(-3.0:0.1:3.0)


###########################################################
############ Sec 6: Spectral Localiser Terms ##############
###########################################################

# mid_x0s = [Float64(maximum(N_range)/2)]
# x0s = mid_x0s # [mid_x0]

alpha = 0.125
N = N_range[1]
# x0s = collect(range(N/2 - alpha * N, N/2 + alpha * N, length=10))
x0s=[N/2]


kappas =  [2.5e-4] #collect(logrange(1e-5, 1e-2, length=16)) # collect(range(1e-4, 5e-2, length=25))
n_lowest_evals = 4


###########################################################
################## Sec 2: Data Save Path ##################
###########################################################

project_name = "specloc_phason_winding_recalc" #"local_mp_data"

root_path = joinpath(project_root, "../../../simulations/raw_data")
if length(t_ranges) < 3
    folder_name = "$project_name/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t1_range)))_t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))_kappa($(kappas[1])-$(kappas[end])-$(length(kappas)))_x0s($(x0s[1])-$(x0s[end])-$(length(x0s)))/"
else
    folder_name = "$project_name/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t1_range)))_t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_t3($(t3_range[1])-$(t3_range[end])-$(length(t3_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))_kappa($(kappas[1])-$(kappas[end])-$(length(kappas)))_x0s($(x0s[1])-$(x0s[end])-$(length(x0s)))/"
end
datasave_path = "$(root_path)/$(folder_name)/"

# Create the folder if it doesn't exist
datasave_path = normpath(datasave_path)
isdir(datasave_path) || mkpath(datasave_path)
println("Data will be saved to: $(datasave_path)")


###########################################################
#################### Sec 3: Tailoring #####################
###########################################################


function get_user_options()
    return UserOptions(
        true,    # calc_mp
        false,   # calc_ipr
        false,   # calc_mbs_energy_gap
        false,   # calc_loc_len
        false,    # calc_proj_wind
        false,   # calc_spectral_func
        eta,     # eta for spectral function calculation
        omega_range, # omega range for spectral function calculation
        :np,     # calc_precision: :hp, :np
        :none, # save_evecs: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :all_np, # save_evals: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :generic # solver_type: :generic, seq_scaled, :mu_loop, :N_loop, :restricted, :phason_loop, :spec_loc, :spec_loc_fast
    )
end


###########################################################
################# Sec 4: Cut Parameters ###################
###########################################################

# This is specific to a 2-phopping system, adjust as needed for other cases
# rho_min = minimum(t2_range) / maximum(t1_range)
# rho_max = maximum(t2_range) / minimum(t1_range)
# rho_step = (rho_max - rho_min) / length(t2_range)
# rho_range = collect(range(rho_min, rho_max, length(t2_range)))
rho_range = collect(range(1.0, 10.0, 10))

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
    ),
    # Dict(
    #     :gradient => 0.0,
    #     :y_intercept => -0.01,
    #     :x_range => (minimum(xs), 1.0),
    #     :y_range => (minimum(ys), maximum(ys)),
    #     :cut_which_side => "above"
    # )
]

mu_rho_restricted = ParamCombGen.restrict_range(xs, ys, cuts)


# ## Function to check the output of cuts before proceeding with simulation
# function plot_mu_rho_restricted(mu_range, rho_range, mu_rho_restricted)

#     restricted_set = Set(mu_rho_restricted)

#     all_points = [(mu, rho) for mu in mu_range, rho in rho_range]

#     restricted_x = Float64[]
#     restricted_y = Float64[]
#     unrestricted_x = Float64[]
#     unrestricted_y = Float64[]
#     for (mu, rho) in all_points
#         if (mu, rho) in restricted_set
#             push!(restricted_x, mu)
#             push!(restricted_y, rho)
#         else
#             push!(unrestricted_x, mu)
#             push!(unrestricted_y, rho)
#         end
#     end

#     plt = scatter(unrestricted_x, unrestricted_y, color=:gray, label="Unrestricted", legend=:topright)
#     scatter!(restricted_x, restricted_y, color=:red, label="Restricted")
#     xlabel!("mu")
#     ylabel!("rho")
#     title!("mu-rho restricted points")

#     savefig(plt, "sample_cuts.png")
# end

# plot_mu_rho_restricted(mu_range, rho_range, mu_rho_restricted)


###########################################################
################ Sec 5: Sequence Scaling ##################
###########################################################

## Allow either t_1 or rho to be preserved when scaling hopping for each sequences
preserve_symbol = :rho  # options: :t1, :rho
rho_target = 1.5    # target rho value when preserve_symbol = :rho
tbar = 1.5         # target t_1 value when preserve_symbol = :t1



###########################################################
####################### Sec 7: Run ########################
###########################################################

# ------------------------------
# --------- Single run ---------
opts = get_user_options()
start_time = now()

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
    preserve_symbol=preserve_symbol,
    x0s=x0s,
    kappas=kappas,
    n_lowest_evals=n_lowest_evals
)

end_time = now()
elapsed = end_time - start_time
println("Simulation completed in $(Dates.value(elapsed)/1000) seconds ($(elapsed))")


# # ------------------------------------------------------------------------------
# # --------- Repeated runs for all sequences (comment out Sec 2 to use) ---------
# sequence_list = [normal_sequence, golden_sequence, silver_sequence, thue_morse_sequence, plastic_sequence]
# sequence_name_list = ["NC", "GQC", "SQC", "TMQC", "PQC"]

# for (i, sequence) in enumerate(sequence_list)
    
#     sequence_name = sequence_name_list[i]
#     sequences = [sequence]

#     # generate data save path
#     root_path = joinpath(project_root, "../../../simulations/raw_data")
#     folder_name = "np/abundance/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t1_range)))_t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))/"
#     datasave_path = "$(root_path)/$(folder_name)/"

#     # Create the folder if it doesn't exist
#     isdir(datasave_path) || mkpath(datasave_path)

#     # run simulation
#     opts = get_user_options()
#     run_selected_solver(opts, N_range, t_combinations, mu_range, Delta_range, sequences, sequence_name, datasave_path; precision=BigFloat_precision, chunk_size=chunk_size, param_restrictions=nothing)
# end