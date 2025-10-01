
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
    
    # usage instructions:
    #     1) Ensure the local_machine environment is initialised by running simulations/data_collection/local_machine/__init__.jl ONCE per Julia session.
    #     2) Complete Sec 1 to define your desired parameter ranges
    #     3) Check Sec 2 to ensure the data save path is correct (it will typically self-generate a folder name based on the parameters chosen)
    #     4) Complete Sec 3 to choose what calculations to perform, in what precision and using which solver type (for more details on solver types see the lead comments in solvers.jl)
    #     (5) (Optional) Complete Sec 4 to define cutting rules to restrict parameter space for :restricted solver type for less expensive 2D param space solving.)
    #     6) Run this script with the code in Sec 5 to execute the simulations and save the data
script_dir = @__DIR__              # e.g. folder1/folder1_1
include(joinpath(script_dir, "../../data_collection/local_machine/solvers.jl"))
include(joinpath(script_dir, "../../data_collection/local_machine/compute.jl"))
include(joinpath(script_dir, "../../data_collection/auxilliary/sequence_gen.jl"))
include(joinpath(script_dir, "../../data_collection/auxilliary/param_comb_gen.jl"))


using .SeqGen
using .ParamCombGen
using .LocalSolv: UserOptions
using BSON: @save, @load
using Plots


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


## N Range (even single value must be a Vector type)
N_range = [15]
#collect(Int, range(55,55,1))
#floor.(Int, log_range(2.0, 3.0, 7.0, 10))


## t_n Range (combine any number of different hopping ranges)
t1_range = collect(range(1.0, 1.0, 1))
t2_range = collect(range(0.0, 10.0, 101))
t3_range = [4.0] #[2.0, 3.0, 4.0]
t_ranges = [t1_range, t2_range, t3_range]
t_combinations_1 = ParamCombGen.t_ranges_combs(t_ranges)

# t1_range = collect(range(0.01, 0.99, 99))
# t2_range = collect(range(1.0, 1.0, 1))
# t_ranges = [t1_range, t2_range]
# t_combinations_2 = ParamCombGen.t_ranges_combs(t_ranges)
# # println(t_combinations_2)
# # println(length(t_combinations_2))

# t_combinations = vcat(t_combinations_1, t_combinations_2)
t_combinations = t_combinations_1


## mu Range (even single calue must be a Vector type)
mu_range = collect(range(0.0, 10.0, length=101))


## Delta Range (even single value must be a vector type)
# Delta_range = [0.1] 
Delta_range = collect(range(0.0, 2.0, 3)) 
#log_range(10.0, 0.0, 0.0, 1)


## Sequence Chunking (use this to generate sequence samples which can be compared)
# chunk_length = 100
# seq_start_indicies = [1]
# sequences = ParamCombGen.split_sequence(normal_sequence, chunk_length, seq_start_indicies)
# if maximum(N_range) > chunk_length
#     error("The largest number of sites (N = $(maximum(N_range))) is incompatible with the chosen sequence chunk sizes ($chunk_length)")
# end


## Sequence Definition (chosoe which standard hopping sequence to use)
sequences = [plastic_sequence]
sequence_name = "PQC"


## Set Precision for hp calculations
BigFloat_precision = 512


## Set chunk size for data saving
chunk_size = 1000



###########################################################
################## Sec 2: Data Save Path ##################
###########################################################


root_path = joinpath(script_dir, "../../../simulations/raw_data")
# root_path = "/Users/mas/Documents/Research/projects/Majoranas/Quasicrystal_Majorana_project/simulations/raw_data"
folder_name = "np/all_crystal_grad_testruns/restricted_mu_vs_rho_mp_heatmaps/$(sequence_name)_N$(N_range[1])-$(N_range[end])-$(length(N_range))_t1$(t1_range[1])-$(t1_range[end])-$(length(t1_range))__t2$(t2_range[1])-$(t2_range[end])-$(length(t2_range))_mu$(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta$(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range))"
path = "$(root_path)/$(folder_name)/"

# Create the folder if it doesn't exist
# println("Data will be saved to path: $path\n")
# println("The path exists: $(isdir(path))\n")
#isdir(path) || mkpath(path)
# Only create if it does not exist
run(`bash -c 'if [ ! -d "$1" ]; then mkdir "$1"; fi' _ $path`)


###########################################################
#################### Sec 3: Tailoring #####################
###########################################################


function get_user_options()
    return UserOptions(
        true,    # calc_mp
        false,   # calc_ipr
        false,   # calc_mbs_energy_gap
        false,   # calc_loc_len
        :np,     # calc_precision: :hp, :np
        :maj_np, # save_evecs: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :maj_np, # save_evals: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :restricted # solver_type: :generic, :mu_loop, :N_loop, :restricted
    )
end


###########################################################
################# Sec 4: Cut Parameters ###################
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
cuts = GQC_D01_cuts # see auxilliary/param_restriction_cuts.jl for predefined cutting rules

mu_rho_restricted = ParamCombGen.restrict_range(xs, ys, cuts)


# Function to check the output of cuts before proceeding with simulation
function plot_mu_rho_restricted(mu_range, rho_range, mu_rho_restricted)
    # Convert restricted tuples to a Set for fast lookup
    restricted_set = Set(mu_rho_restricted)
    # Prepare grid of all (mu, rho) pairs
    all_points = [(mu, rho) for mu in mu_range, rho in rho_range]
    # Separate restricted and unrestricted points
    restricted_x = Float64[]
    restricted_y = Float64[]
    unrestricted_x = Float64[]
    unrestricted_y = Float64[]
    for (mu, rho) in all_points
        if (mu, rho) in restricted_set
            push!(restricted_x, mu)
            push!(restricted_y, rho)
        else
            push!(unrestricted_x, mu)
            push!(unrestricted_y, rho)
        end
    end
    # Plot
    scatter(unrestricted_x, unrestricted_y, color=:gray, label="Unrestricted", legend=:topright)
    scatter!(restricted_x, restricted_y, color=:red, label="Restricted")
    xlabel!("mu")
    ylabel!("rho")
    title!("mu-rho restricted points")
end

plot_mu_rho_restricted(mu_range, rho_range, mu_rho_restricted)


###########################################################
####################### Sec 5: Run ########################
###########################################################

# ------------------------------
# --------- Single run ---------
opts = get_user_options()
run_selected_solver(opts, N_range, t_combinations, mu_range, Delta_range, sequences, sequence_name, path; precision=BigFloat_precision, chunk_size=chunk_size, param_restrictions=mu_rho_restricted)


# # ------------------------------------------------------------------------------
# # --------- Repeated runs for all sequences (comment out Sec 2 to use) ---------
# sequence_list = [normal_sequence, golden_sequence, silver_sequence, thue_morse_sequence, plastic_sequence]
# sequence_name_list = ["NC", "GQC", "SQC", "TMQC", "PQC"]

# for (i, sequence) in enumerate(sequence_list)
    
#     sequence_name = sequence_name_list[i]
#     sequences = [sequence]

#     # generate data save path
#     root_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data"
#     folder_name = "np/all_crystal_grad_testruns/mu_vs_rho_mp_heatmaps/$(sequence_name)_N($(N_range[1])-$(N_range[end])-$(length(N_range)))_t1($(t1_range[1])-$(t1_range[end])-$(length(t2_range))__t2($(t2_range[1])-$(t2_range[end])-$(length(t2_range)))_mu($(mu_range[1])-$(mu_range[end])-$(length(mu_range)))_Delta($(Delta_range[1])-$(Delta_range[end])-$(length(Delta_range)))/"
#     path = "$(root_path)/$(folder_name)/"

#     # Create the folder if it doesn't exist
#     isdir(path) || mkpath(path)

#     # run simulation
#     opts = get_user_options()
#     run_selected_solver(opts, N_range, t_combinations, mu_range, Delta_range, sequences, sequence_name, path; precision=BigFloat_precision, chunk_size=chunk_size)
# end