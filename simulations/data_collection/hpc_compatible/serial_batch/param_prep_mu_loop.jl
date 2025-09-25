include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/param_comb_gen.jl")

using .ParamCombGen
using DelimitedFiles



# # Parameter ranges
# golden_exponents = collect(3:10)
# N_range = SeqGen.golden_LengthCalc(golden_exponents)
N_range = [50] #floor.(Int, collect(range(2,100,25)))

t1s=1.0
t1e=1.0
t1l=1

t2s=0.0
t2e=10.0
t2l=201

t3s=3.14
t3e=3.14
t3l=1

mus=0.0
mue=10.0
mul=201

t1_range = collect(range(t1s, t1e, t1l))
t2_range = collect(range(t2s, t2e, t2l))
t3_range = collect(range(t3s, t3e, t3l))
t_ranges = [t1_range, t2_range, t3_range]
t_combinations_1 = ParamCombGen.t_ranges_combs(t_ranges)
# println(t_combinations_1)
# println(length(t_combinations_1))

# t1_range = collect(range(0.0, 1.0, 3))
# t2_range = collect(range(1.0, 1.0, 1))
# t_ranges = [t1_range, t2_range]
# t_combinations_2 = ParamCombGen.t_ranges_combs(t_ranges)
# # println(t_combinations_1)
# # println(length(t_combinations_1))

# t_combinations = vcat(t_combinations_1, t_combinations_2)
t_combinations = t_combinations_1
# println(length(t_combinations))

mu_range =  collect(range(mus, mue, mul))

Delta_range_1 = collect(range(0.25, 2.0, 8)) #ParamCombGen.log_range(2.0, -1.0, 3.0, 5)
# Delta_range_2 = collect(range(1.0, 10.0, 10)) #collect(range(0.05, 0.49, 45))
Delta_range_3 = [0.1]
Delta_range = vcat(Delta_range_3, Delta_range_1)
println(Delta_range)
# println(length(Delta_range))

sequence_name = "GQC"
N = N_range[1]

mu_start = mu_range[1]
mu_end = mu_range[end]
mu_length = length(mu_range)

# Generate combinations
combinations = [(N, mus, mue, mul, Delta, sequence_name, t_n[1], t_n[2], t_n[3]) for N in N_range for t_n in t_combinations for Delta in Delta_range]
println(length(combinations))

# Write to .dat file
filepath = "/user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/param_sets/mu_loop_sets/final_abundance_sets/"
filename = "$(filepath)params_$(sequence_name)_N$(N)_mu$(mus)-$(mue)-$(mul)_rho$(t2s)-$(t2e)-$(t2l)_Delta$([Delta for Delta in Delta_range]).dat"
open(filename, "w") do io
    writedlm(io, ["N mu_start mu_end mu_length Delta sequence t1 t2 t3"])  # Write the header row
    for combo in combinations
        writedlm(io, [combo])  # Ensure each combination is written as a row
    end
end

