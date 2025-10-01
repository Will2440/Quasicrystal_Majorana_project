# Only run this ONCE per Julia session to start up the local_machine environment. Restart the Julia REPL to reset.
include("../../../simulations/data_collection/local_machine/solvers.jl")
include("../../../simulations/data_collection/local_machine/compute.jl")
include("../../../simulations/data_collection/auxilliary/sequence_gen.jl")
include("../../../simulations/data_collection/auxilliary/param_comb_gen.jl")

using .LocalSolv
using .SeqGen
using .ParamCombGen
using BSON: @save, @load

print("local_machine environment initialised. Proceed to local_machine/main.jl to set parameters and run simulations.\n")