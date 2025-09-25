"""
    file name:   main.jl
    created:     25/09/2025
    last edited: 25/09/2025

    overview:
        The main script to execute simulations on a local machine. Define parameter ranges and decide on calculation and datasaving requirements to call on the correct functions.
    
    structure:
        - Sec 1:  Parameter Choice
                    Generate all sequence types and define parameter ranges for N, t_n, mu and Delta
        - Sec 2:  Data Save Path  
                    Ddefine data save filepath
        - Sec 3:  Tailoring
                    Answer list of options relating to exactly what values should be calculated and saved, in what precision and using which solver type.
    
    usage instructions:
        1) Ensure the local_machine environment is initialised by running simulations/data_collection/local_machine/__init__.jl ONCE per Julia session.
        2) Complete Sec 1 to define your desired parameter ranges
        3) Check Sec 2 to ensure the data save path is correct (it will typically self-generate a folder name based on the parameters chosen)
        4) Complete Sec 3 to choose what calculations to perform, in what precision and using which solver type (for more details on solver types see the lead comments in solvers.jl)
        5) Run this script to execute the simulations and save the data
"""

using .SeqGen
using .ParamCombGen
using .LocalSolv: UserOptions
using BSON: @save, @load


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
N_range = [20]
#collect(Int, range(55,55,1))
#floor.(Int, log_range(2.0, 3.0, 7.0, 10))


## t_n Range (combine any number of different hopping ranges)
t1_range = collect(range(1.0, 1.0, 1))
t2_range = collect(range(2.0, 2.0, 11))
t_ranges = [t1_range, t2_range]
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
mu_range = collect(range(0.0, 3.0, length=11))


## Delta Range (even single value must be a vector type)
Delta_range = [0.1] 
#collect(range(0.0, 1.0, 51)) 
#log_range(10.0, 0.0, 0.0, 1)


## Sequence Chunking (use this to generate sequence samples which can be compared)
# chunk_length = 100
# seq_start_indicies = [1]
# sequences = ParamCombGen.split_sequence(normal_sequence, chunk_length, seq_start_indicies)
# if maximum(N_range) > chunk_length
#     error("The largest number of sites (N = $(maximum(N_range))) is incompatible with the chosen sequence chunk sizes ($chunk_length)")
# end


## Sequence Definition (chosoe which standard hopping sequence to use)
sequences = [golden_sequence]
sequence_name = "GQC"


## Set Precision for hp calculations
BigFloat_precision = 512


## Set chunk size for data saving
chunk_size = 1000



###########################################################
################## Sec 2: Data Save Path ##################
###########################################################


root_path = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data"
folder_name = "test_folder_hp_all_calcs"
path = "$(root_path)/$(folder_name)/"

# Create the folder if it doesn't exist
isdir(path) || mkpath(path)



###########################################################
#################### Sec 2: Tailoring #####################
###########################################################


function get_user_options()
    return UserOptions(
        true,    # calc_mp
        true,   # calc_ipr
        true,   # calc_mbs_energy_gap
        true,   # calc_loc_len
        :hp,     # calc_precision: :hp, :np
        :maj_np, # save_evecs: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :all_np, # save_evals: :all_np, :all_hp, :maj_np, :maj_hp, :none
        :generic # solver_type: :generic, :mu_loop, :N_loop
    )
end

# --- Main execution ---
opts = get_user_options()
run_selected_solver(opts, N_range, t_combinations, mu_range, Delta_range, sequences, sequence_name, path; precision=BigFloat_precision, chunk_size=chunk_size)
