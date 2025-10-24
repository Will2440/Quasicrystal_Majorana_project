include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/hp_serial_solver.jl")
include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/sequence_gen.jl")

using .hpSerialSolver
using .SeqGen

using DelimitedFiles


file_index = parse(Int64,ARGS[1]);


###########################################################
############### Sec 1: seqeunce Generation ################
###########################################################
N_seq_seed = 1024
normal_sequence = SeqGen.normal_SeqGen(N_seq_seed)
golden_sequence = SeqGen.golden_SeqGen(N_seq_seed)
silver_sequence = SeqGen.silver_SeqGen(N_seq_seed)
thue_morse_sequence = SeqGen.thue_morse_SeqGen(N_seq_seed)
plastic_sequence = SeqGen.plastic_SeqGen(N_seq_seed)
trib_sequence = SeqGen.tribonacci_SeqGen(N_seq_seed)



# Read all parameter sets from the file
function read_all_parameters(filename::String)
    raw_data = readdlm(filename, '\t', String; header=true)
    data = raw_data[1]
    param_sets = Vector{Tuple}
    for row in 1:size(data, 1)
        N = parse(Int, data[row, 1])
        mu = parse(Float64, data[row, 2])
        rho = parse(Float64, data[row, 3])
        Delta = parse(Float64, data[row, 4])
        sequence = data[row, 5]
        t3 = parse(Float64, data[row, 6])
        push!(param_sets, (N, mu, rho, Delta, sequence, t3))
    end
    return param_sets
end

function check_param_sets(param_sets::Vector{Tuple})
    N_ref, _, _, Delta_ref, sequence_ref, t3_ref = param_sets[1]
    for (i, (N, _, _, Delta, sequence, t3)) in enumerate(param_sets)
        if N != N_ref
            println("Mismatch in N at row $i: $N != $N_ref")
        end
        if Delta != Delta_ref
            println("Mismatch in Delta at row $i: $Delta != $Delta_ref")
        end
        if t3 != t3_ref
            println("Mismatch in t3 at row $i: $t3 != $t3_ref")
        end
        if sequence != sequence_ref
            println("Mismatch in sequence at row $i: $sequence != $sequence_ref")
        end
    end

    return N_ref, Delta_ref, sequence_ref, t3_ref
end

function check_sequence(sequence_name::String)
    if sequence_name == "GQC"
        return golden_sequence
    elseif sequence_name == "SQSC"
        return silver_sequence
    elseif sequence_name == "TMQC"
        return thue_morse_sequence
    elseif sequence_name == "PQC"
        return plastic_sequence
    else
        error("Unknown sequence name: $sequence_name")
    end
end

##############################################
# # check here # #
##############################################
params_filename = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/param_sets/mu_loop_sets/final_abundance_sets/params_PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0].dat"
##############################################

param_sets = read_all_parameters(params_filename)
N, Delta, sequence_name, t3 = check_param_sets(param_sets)
sequence = check_sequence(sequence_name)

mu_rho_points = [(mu, rho) for (_, mu, rho, _, _, _) in param_sets]

# Print parsed variables
println("N: ", N)
println("Delta: ", Delta)
println("sequence name: ", sequence_name)
println("t3: ", t3)


##############################################
# # check here # #
##############################################
precision = 512
data_save_filepath = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/results/mu_loop_results/final_abundance_results/TRB_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]/"
##############################################


hpSerialSolver.hp_mu_rho_restricted_solver(
    N, 
    mu_rho_points,
    Delta,
    sequence,
    sequence_name,
    precision,
    data_save_filepath
)

