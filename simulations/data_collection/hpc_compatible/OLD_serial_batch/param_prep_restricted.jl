# include("/user/home/hb21877/Majorana_solver_for_BlueCrystal/param_comb_gen.jl")
include("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/param_comb_gen.jl")

using .ParamCombGen
using DelimitedFiles



###########################################################
################# Sec 1: Parameter Choice #################
###########################################################

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

mu_range =  collect(range(mus, mue, mul))

Delta_range = collect(range(0.1, 2.0, 2))

sequence_name = "GQC"

N = N_range[1]



###########################################################
################# Sec 2: Cut Parameters ###################
###########################################################

# This is specific to a 2-phopping system, adjust as needed for other cases
rho_min = t2s / t1e
rho_max = t2e / t1s
rho_step = (rho_max - rho_min) / t2l
rho_range = collect(range(rho_min, rho_max, t2l))

xs = mu_range
ys = rho_range
include("/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/auxilliary/param_restriction_cuts.jl")

cuts = GQC_D01_cuts # see auxilliary/param_restriction_cuts.jl for predefined cutting rules

mu_rho_restricted = ParamCombGen.restrict_range(xs, ys, cuts)



###########################################################
################ Sec 3: Generate Batches ##################
###########################################################

points_per_job = 1000 # how many (mu, rho) points to simulate per job

root = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/serial_batch/param_sets/restricted_sets/"
Delta_range_string = join([string(Delta) for Delta in Delta_range], ",")
# Delta_range_string_safe = replace(string(Delta_range_string), "." => "p")
run_folder = "test_abundance_sets/$(sequence_name)_N$(N)_mu$(mus)-$(mue)-$(mul)_rho$(t2s)-$(t2e)-$(t2l)_Delta$(Delta_range_string)/"
run_folder_safe = replace(string(run_folder), "." => "p")
folder_path = "$(root)$(run_folder_safe)"

# isdir(folder_path) || mkpath(folder_path)
# println("writing to folder: $(folder_path)")

function write_mu_rho_jobs(mu_rho_restricted::Vector{Tuple{Float64, Float64}}, N::Int, Delta::Float64, t3::Float64, sequence_name::String, points_per_job::Int, folder_path::String)
    
    num_jobs = ceil(Int, length(mu_rho_restricted) / points_per_job)

    for job_idx in 1:num_jobs

        start_idx = (job_idx - 1) * points_per_job + 1
        end_idx = min(job_idx * points_per_job, length(mu_rho_restricted))
        chunk = mu_rho_restricted[start_idx:end_idx]

        filename = "N$(N)_Delta$(Delta)_t3$(t3)_job$(job_idx)_npoints$(length(chunk)).dat"
        path = "$(folder_path)/$(filename)"
        isdir(path) || mkpath(path)

        open(filename, "w") do io
            writedlm(io, ["N mu rho Delta sequence t3"])  # Write the header row
            for (mu, rho) in chunk
                println(io, "$(N) $(mu) $(rho) $(Delta) $(sequence_name) $(t3)")
            end
        end
    end
end

for Delta in Delta_range
    for N in N_range
        for t3 in t3_range
            write_mu_rho_jobs(mu_rho_restricted, N, Delta, t3, sequence_name, points_per_job, folder_path)
        end
    end
end

println("finished")