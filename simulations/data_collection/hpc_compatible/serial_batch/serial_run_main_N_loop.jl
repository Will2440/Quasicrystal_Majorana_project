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

    N_start = parse(Int, data[1])
    N_end = parse(Int, data[2])
    N_length = parse(Int, data[3])
    mu_start = parse(Float64, data[4])
    mu_end = parse(Float64, data[5])
    mu_length = parse(Int, data[6])
    Delta = parse(Float64, data[7])
    sequence = data[8]
    t_columns = [parse(Float64, data[i]) for i in 9:11] 

    return N_start, N_end, N_length, mu_start, mu_end, mu_length, Delta, sequence, t_columns
end

##############################################
# # check here # #
##############################################
params_filename = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/param_sets/N_loop_sets/final_fractality/params_GQC_N13-610-9_mu0.0-3.0-201_rho1.5-1.5-1_Delta1.0.dat"
##############################################

N_start, N_end, N_length, mu_start, mu_end, mu_length, Delta, sequence_name, t_columns = read_parameters(params_filename, row_index)

# Print parsed variables
println("N length: ", N_length)
println("N start: ", N_start)
println("N end: ", N_end)
println("mu length: ", mu_length)
println("mu start: ", mu_start)
println("mu end: ", mu_end)
println("Delta: ", Delta)
println("sequence name: ", sequence_name)
println("t_columns (t1, t2, t3): ", t_columns)


##############################################
# # check here # #
##############################################
## change this for PQC!! to tn = t_columns[1:3]
tn = t_columns[1:3]

mu_range = collect(range(mu_start, mu_end, mu_length))
N_range = floor.(Int, collect(range(N_start, N_end, N_length)))
sequence = plastic_sequence

precision = 512
data_save_filepath = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/results/N_loop_results/final_fractality/GQC_N13-610-9_mu0.0-3.0-201_rho1.5-1.5-1_Delta1.0/"
##############################################


hpSerialSolver.hp_serial_solver_N_loop_threaded(
    N_range, 
    tn,
    mu_range,
    Delta,
    sequence,
    sequence_name,
    precision,
    data_save_filepath
)

