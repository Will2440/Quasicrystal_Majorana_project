base_path = ENV["BOUNDARY_ARRAY_BASE_CONFIG"]
parameter = Symbol(ENV["BOUNDARY_ARRAY_PARAMETER"])
value = parse(Float64, ENV["BOUNDARY_ARRAY_VALUE"])
task_id = ENV["SLURM_ARRAY_TASK_ID"]
array_id = ENV["SLURM_ARRAY_JOB_ID"]

cfg = include(base_path)
cfg isa Dict{Symbol,Any} || error("Base configuration must return a Dict{Symbol,Any}")
isfinite(value) && value > 0 || error("Array fixed parameter must be finite and positive")
parameter in (:Delta, :t1_fixed, :t2_fixed) || error("Unsupported array fixed parameter: $parameter")
if cfg[:sweep_parameter] == :rho
    parameter == :Delta ||
        (parameter == :t1_fixed && cfg[:sweep_mode] == :vary_t2) ||
        (parameter == :t2_fixed && cfg[:sweep_mode] == :vary_t1) ||
        error("$parameter is not fixed in this rho sweep")
elseif cfg[:sweep_parameter] == :Delta
    parameter in (:t1_fixed, :t2_fixed) || error("$parameter is not fixed in a Delta sweep")
else
    error("Unsupported sweep parameter: $(cfg[:sweep_parameter])")
end
cfg[parameter] = value
cfg[:run_name] = "$(cfg[:run_name])_array$(array_id)_task$(task_id)"
cfg[:run_root_pointer_file] = "$(cfg[:run_root_pointer_file]).array$(array_id)_task$(task_id)"
cfg