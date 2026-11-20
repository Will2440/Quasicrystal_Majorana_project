using LinearAlgebra
using Test
using JLD2

include(joinpath(@__DIR__, "solvers.jl"))
using .HPCSolv

@testset "High-precision real generic solver" begin
    N = 4
    hoppings = [1.0, 1.5]
    sequence = [1, 2, 1, 2]
    mu = 0.7
    Delta = 0.2
    precision = 128
    complex_matrix = HPCSolv.hp_create_bdg_hamiltonian(N, hoppings, mu, Delta, sequence, precision)
    real_matrix = HPCSolv.hp_create_bdg_hamiltonian(N, hoppings, mu, Delta, sequence, precision; real_matrix=true)

    @test eltype(real_matrix) == BigFloat
    @test real_matrix == complex_matrix
    complex_values, complex_vectors = eigen(Hermitian(complex_matrix))
    real_values, real_vectors = eigen(Hermitian(real_matrix))
    @test isapprox(real_values, complex_values; atol=big"1e-30")
    @test isapprox(HPCSolv.hp_calc_maj_mp(real_vectors), HPCSolv.hp_calc_maj_mp(complex_vectors); atol=1e-12)
    for calculate in (HPCSolv.hp_calc_maj_ipr, HPCSolv.hp_calc_maj_loc_len)
        @test isapprox(calculate(real_vectors), calculate(complex_vectors); atol=1e-12)
    end
    @test isapprox(
        HPCSolv.hp_calc_spec_func(real_values, real_vectors, 0.1, [0.0]),
        HPCSolv.hp_calc_spec_func(complex_values, complex_vectors, 0.1, [0.0]);
        atol=1e-12
    )
    real_scales = HPCSolv.hp_calc_local_scales_mp_allscales(real_vectors)
    complex_scales = HPCSolv.hp_calc_local_scales_mp_allscales(complex_vectors)
    @test length(real_scales) == length(complex_scales) == 10
    @test all(isapprox(real_scale.overall_MP, complex_scale.overall_MP; atol=1e-12)
              for (real_scale, complex_scale) in zip(real_scales, complex_scales))

    mktempdir() do directory
        for mode in (:none, :all_hp, :all_np, :maj_hp, :maj_np)
            options = UserOptions(true, true, true, true, true, true, 0.1, [0.0], :hp, mode, :none, :generic, false, 0)
            prefix = joinpath(directory, string(mode))
            HPCSolv.hp_generic_solver([N], [hoppings], [mu], [Delta], [sequence], "test", [(0.0, 0.0)], precision, 10, prefix, options, 1)
            results = JLD2.load(prefix * "_row1_thread_1_chunk_1.jld2", "results_df")
            @test isapprox(results.mp[1], HPCSolv.hp_calc_maj_mp(complex_vectors); atol=1e-12)
            if mode == :none
                @test ismissing(results.eigenvectors[1])
            else
                expected_type = mode in (:all_hp, :maj_hp) ? Complex{BigFloat} : ComplexF64
                expected_columns = mode in (:all_hp, :all_np) ? 2N : 2
                @test eltype(results.eigenvectors[1]) == expected_type
                @test size(results.eigenvectors[1]) == (2N, expected_columns)
            end
        end
    end

    mktempdir() do directory
        mus = [0.3, 0.5, 0.7, 0.9]
        options = UserOptions(true, false, false, false, false, false, 0.1, [0.0], :hp, :none, :none, :generic, false, 0)
        prefix = joinpath(directory, "threaded")
        HPCSolv.hp_generic_solver([N], [hoppings], mus, [Delta], [sequence], "test", [(0.0, 0.0)], precision, 1, prefix, options, 1)
        files = filter(path -> endswith(path, ".jld2"), readdir(directory; join=true))
        @test length(files) == length(mus)
        results = [JLD2.load(file, "results_df") for file in files]
        @test sort([result.mu[1] for result in results]) == mus
        for result in results
            matrix = HPCSolv.hp_create_bdg_hamiltonian(N, hoppings, result.mu[1], Delta, sequence, precision; real_matrix=true)
            values, vectors = eigen(Hermitian(matrix))
            @test isapprox(result.mp[1], HPCSolv.hp_calc_maj_mp(vectors); atol=1e-12)
        end
    end
end