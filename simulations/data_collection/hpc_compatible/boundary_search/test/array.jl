using Test

root = normpath(joinpath(@__DIR__, ".."))
base_config = joinpath(root, "config", "default_config.jl")
task_config = joinpath(root, "config", "array_task_config.jl")
launcher = joinpath(root, "runscripts", "boundary_array.sbatch")

function load_task(base, parameter, value, task="0")
    withenv("BOUNDARY_ARRAY_BASE_CONFIG" => base,
            "BOUNDARY_ARRAY_PARAMETER" => parameter,
            "BOUNDARY_ARRAY_VALUE" => value,
            "SLURM_ARRAY_JOB_ID" => "123",
            "SLURM_ARRAY_TASK_ID" => task) do
        include(task_config)
    end
end

@testset "array configuration" begin
    base = include(base_config)
    first_task = load_task(base_config, "Delta", "0.01")
    second_task = load_task(base_config, "Delta", "0.05", "2")
    @test first_task[:Delta] == 0.01
    @test second_task[:Delta] == 0.05
    @test first_task[:rho_start] == base[:rho_start]
    @test first_task[:rho_max] == base[:rho_max]
    @test first_task[:run_name] != second_task[:run_name]
    @test first_task[:run_root_pointer_file] != second_task[:run_root_pointer_file]
    @test first_task[:global_stop_file] == base[:global_stop_file]
    @test load_task(base_config, "t1_fixed", "1.5")[:t1_fixed] == 1.5
    @test_throws LoadError load_task(base_config, "t2_fixed", "1.5")
    @test_throws LoadError load_task(base_config, "Delta", "0")
    @test_throws LoadError load_task(base_config, "Delta", "NaN")
    @test_throws LoadError load_task(base_config, "rho_max", "2")

    mktempdir() do directory
        delta_config = joinpath(directory, "delta_config.jl")
        write(delta_config, "cfg = include($(repr(base_config)))\ncfg[:sweep_parameter] = :Delta\ncfg\n")
        @test load_task(delta_config, "t2_fixed", "3.0")[:t2_fixed] == 3.0
        @test_throws LoadError load_task(delta_config, "Delta", "0.1")
        vary_t1_config = joinpath(directory, "vary_t1_config.jl")
        write(vary_t1_config, "cfg = include($(repr(base_config)))\ncfg[:sweep_mode] = :vary_t1\ncfg\n")
        @test load_task(vary_t1_config, "t2_fixed", "3.0")[:t2_fixed] == 3.0
        @test_throws LoadError load_task(vary_t1_config, "t1_fixed", "3.0")
    end
end

@testset "Slurm array selection" begin
    mktempdir() do directory
        mkpath(joinpath(directory, "runscripts"))
        master = joinpath(directory, "runscripts", "boundary_master.sbatch")
        write(master, "printf '%s %s %s\\n' \"\$BOUNDARY_ARRAY_PARAMETER\" \"\$BOUNDARY_ARRAY_VALUE\" \"\$BOUNDARY_ARRAY_BASE_CONFIG\"\n")
        withenv("BOUNDARY_SEARCH_DIR" => directory, "SLURM_SUBMIT_DIR" => root,
                "SLURM_ARRAY_JOB_ID" => "123") do
            for (index, value) in (("0", "0.01"), ("2", "0.05"), ("3", "0.1"))
                withenv("SLURM_ARRAY_TASK_ID" => index) do
                    output = read(`bash $launcher $base_config`, String)
                    @test occursin("Delta $value $base_config", output)
                end
            end
            withenv("SLURM_ARRAY_TASK_ID" => "4") do
                @test !success(pipeline(`bash $launcher $base_config`; stdout=devnull, stderr=devnull))
            end
        end
    end
end