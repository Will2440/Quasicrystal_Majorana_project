"""
    file name:   hpc_compatible/serial_batch/compute.jl
    created:     04/11/2025
    last edited: 05/11/2025

    overview:
        This file contains the dispatch function which calls on the appropriate solver function based on user input from main.jl.
	
	Latest edits:
		- Added support for :restricted solver type to allow parameter space restrictions during simulations.
		- Added sequence_id argument to track and save sequence identifiers alongside simulation results. (Used for Sturmian sequences.)
"""


using .HPCSolv

function run_selected_solver(
    opts::UserOptions,
    N_range::Vector{Int},
    t_n_range::Vector{Vector{Float64}},
    mu_range::Vector{Float64},
    Delta_range::Vector{Float64},
    sequences::Vector{Vector{Int}},
    sequence_name::String,
    filepath::String;
    precision::Int=512,
    chunk_size::Int=1000,
    param_restrictions::Union{Nothing, Vector{Tuple{Float64,Float64}}}=nothing,
    sequence_ids::Union{Nothing, Vector{Tuple{Float64,Float64}}}=nothing,
    rho_target::Float64=1.5,
    tbar::Float64=1.5,
    preserve_symbol::Symbol=:t1,
    row_index::Union{Int,Nothing}=nothing
)

    # helper: pick the first element for mu_loop and warn if multiple
    pick1(name, vec) = begin
        isempty(vec) && error("$name is empty; expected at least one value")
        if length(vec) > 1
            println("Note: :mu_loop uses only the first $name; ignoring $(length(vec)-1) others")
        end
        first(vec)
    end

    if opts.solver_type == :generic
        if opts.calc_precision == :np
            HPCSolv.np_generic_solver(N_range, t_n_range, mu_range, Delta_range, sequences, sequence_name, sequence_ids, chunk_size, filepath, opts, row_index)
        elseif opts.calc_precision == :hp
            HPCSolv.hp_generic_solver(N_range, t_n_range, mu_range, Delta_range, sequences, sequence_name, sequence_ids, precision, chunk_size, filepath, opts, row_index)
        elseif opts.calc_precision == :Arpack
            error("ARPACK-based generic solver not yet implemented.")
        else
            error("Unknown precision mode: $(opts.calc_precision)")
        end

    elseif opts.solver_type == :mu_loop
        # Use first of each range for a 1D mu loop
        N = pick1("N_range", N_range)
        t_n = pick1("t_n_range", t_n_range)
        Delta = pick1("Delta_range", Delta_range)
        seq = pick1("sequences", sequences)
        # sequence_id required for consistency (Sturmian usage); use first if provided
        seq_id = begin
            if sequence_ids === nothing || isempty(sequence_ids)
                println("Note: :mu_loop called without sequence_ids; using (NaN, NaN) as placeholder.")
                (NaN, NaN)
            else
                if length(sequence_ids) > 1
                    println("Note: :mu_loop uses only the first sequence_id; ignoring $(length(sequence_ids)-1) others")
                end
                first(sequence_ids)
            end
        end

        if opts.calc_precision == :np
            HPCSolv.np_mu_loop_solver(N, t_n, mu_range, Delta, seq, sequence_name, seq_id, chunk_size, filepath, opts, row_index)
        elseif opts.calc_precision == :hp
            HPCSolv.hp_mu_loop_solver(N, t_n, mu_range, Delta, seq, sequence_name, seq_id, precision, chunk_size, filepath, opts, row_index)
        elseif opts.calc_precision == :Arpack
            error("ARPACK-based mu_loop solver not yet implemented.")
        else
            error("Unknown precision mode: $(opts.calc_precision)")
        end

    elseif opts.solver_type == :N_loop
        if opts.calc_precision == :np
            error(":N_loop solver for normal precision not yet implemented.")
        elseif opts.calc_precision == :hp
            HPCSolv.hp_serial_solver_N_loop(N_range, t_n_range[1], mu_range, Delta_range[1], sequences[1], sequence_name, precision, filepath, opts, row_index)
        elseif opts.calc_precision == :Arpack
            error("ARPACK-based N_loop solver not yet implemented.")
        else
            error("Unknown precision mode: $(opts.calc_precision)")
        end

    elseif opts.solver_type == :restricted
        if opts.calc_precision != :np
            println(":restricted solver is only implemented for :np (got $(opts.calc_precision))")
            error("Unsupported precision for :restricted: $(opts.calc_precision)")
        end
        if param_restrictions === nothing
            error("param_restrictions must be provided for :restricted solver.")
        end
        HPCSolv.np_mu_rho_restricted_solver(N_range, t_n_range, mu_range, param_restrictions, Delta_range, sequences, sequence_name, chunk_size, filepath, opts, row_index)

    elseif opts.solver_type == :seq_scaled
        if opts.calc_precision != :np
            println(":seq_scaled solver is only implemented for :np (got $(opts.calc_precision))")
            error("Unsupported precision for :seq_scaled: $(opts.calc_precision)")
        end
        HPCSolv.np_seq_scaled_solver(N_range, rho_target, tbar, mu_range, Delta_range, sequences, sequence_name, sequence_ids, chunk_size, filepath, opts; preserve=preserve_symbol, row_index=row_index)

    else
        error("Unknown solver type: $(opts.solver_type)")
    end
end