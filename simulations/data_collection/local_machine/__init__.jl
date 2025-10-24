# Only run this ONCE per Julia session to start up the local_machine environment. Restart the Julia REPL to reset.
project_root = @__DIR__

include(joinpath(project_root, "../../data_collection/local_machine/solvers.jl"))
include(joinpath(project_root, "../../data_collection/local_machine/compute.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/sequence_gen.jl"))
include(joinpath(project_root, "../../data_collection/auxilliary/param_comb_gen.jl"))

using .LocalSolv
using .SeqGen
using .ParamCombGen
using BSON: @save, @load

print("local_machine environment initialised. Proceed to local_machine/main.jl to set parameters and run simulations.\n")