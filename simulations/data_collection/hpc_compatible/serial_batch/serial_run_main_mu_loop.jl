include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/hp_serial_solver.jl")
include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/sequence_gen.jl")

using .hpSerialSolver
using .SeqGen

using DelimitedFiles


row_index = parse(Int64,ARGS[1]);


# seqeunce generation
N_seq_seed = 1024
normal_sequence = SeqGen.normal_SeqGen(N_seq_seed)
golden_sequence = SeqGen.golden_SeqGen(N_seq_seed)
silver_sequence = SeqGen.silver_SeqGen(N_seq_seed)
thue_morse_sequence = SeqGen.thue_morse_SeqGen(N_seq_seed)
plastic_sequence = SeqGen.plastic_SeqGen(N_seq_seed)
trib_sequence = SeqGen.tribonacci_SeqGen(N_seq_seed)


function read_parameters(filename::String, row_index::Int)

    raw_data = readdlm(filename, '\t', String; header=true)
    data = raw_data[1][row_index, :]


    N = parse(Int, data[1])
    mu_start = parse(Float64, data[2])
    mu_end = parse(Float64, data[3])
    mu_length = parse(Int, data[4])
    Delta = parse(Float64, data[5])
    sequence = data[6]
    t_columns = [parse(Float64, data[i]) for i in 7:9] 

    return N, mu_start, mu_end, mu_length, Delta, sequence, t_columns
end

##############################################
# # check here # #
##############################################
params_filename = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/param_sets/mu_loop_sets/final_abundance_sets/params_PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0].dat"
##############################################

N, mu_start, mu_end, mu_length, Delta, sequence_name, t_columns = read_parameters(params_filename, row_index)

# Print parsed variables
println("N: ", N)
println("mu length: ", mu_length)
println("mu start: ", mu_start)
println("mu end: ", mu_end)
println("Delta: ", Delta)
println("sequence name: ", sequence_name)
println("t_columns (t1, t2, t3): ", t_columns)


##############################################
# # check here # #
##############################################
## change this for PQC!! to tn = t_columns[1:3], else [1:2]
tn = t_columns[1:3]

mu_range = collect(range(mu_start, mu_end, mu_length))
sequence = trib_sequence

precision = 512
data_save_filepath = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/results/mu_loop_results/final_abundance_results/TRB_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]/"
##############################################


hpSerialSolver.hp_serial_solver_mu_loop(
    N, 
    tn, 
    mu_range,
    Delta,
    sequence,
    sequence_name,
    precision,
    data_save_filepath
)

