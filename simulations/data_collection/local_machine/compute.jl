"""
    file name:   compute.jl
    created:     25/09/2025
    last edited: 25/09/2025

    overview:
        This file contains the dispatch function which calls on the appropriate solver function based on user input from main.jl.
"""


using .LocalSolv


function run_selected_solver(opts::UserOptions, N_range, t_n_range, mu_range, Delta_range, sequences, sequence_name, filepath; precision=512, chunk_size=1000, param_restrictions)

	if opts.solver_type == :generic
		if opts.calc_precision == :np
			LocalSolv.np_generic_solver(N_range, t_n_range, mu_range, Delta_range, sequences, sequence_name, chunk_size, filepath, opts)
		elseif opts.calc_precision == :hp
			LocalSolv.hp_generic_solver(N_range, t_n_range, mu_range, Delta_range, sequences, sequence_name, precision, chunk_size, filepath, opts)
		elseif opts.calc_precision == :Arpack
			error("ARPACK-based generic solver not yet implemented.")
		else
			error("Unknown precision mode: $(opts.calc_precision)")
		end
	elseif opts.solver_type == :mu_loop
		if opts.calc_precision == :np
			# Example: call np_mu_loop_solver (implement as needed)
			error(":mu_loop solver for normal precision not yet implemented.")
		elseif opts.calc_precision == :hp
			LocalSolv.hp_mu_loop_solver(
				N_range[1], t_n_range[1], mu_range, Delta_range[1], sequences[1], sequence_name, precision, filepath, opts)
		elseif opts.calc_precision == :Arpack
			error("ARPACK-based mu_loop solver not yet implemented.")
		else
			error("Unknown precision mode: $(opts.calc_precision)")
		end
	elseif opts.solver_type == :N_loop
		if opts.calc_precision == :np
			# Example: call np_serial_solver_N_loop (implement as needed)
			error(":N_loop solver for normal precision not yet implemented.")
		elseif opts.calc_precision == :hp
			LocalSolv.hp_serial_solver_N_loop(
				N_range, t_n_range[1], mu_range, Delta_range[1], sequences[1], sequence_name, precision, filepath, opts)
		elseif opts.calc_precision == :Arpack
			error("ARPACK-based N_loop solver not yet implemented.")
		else
			error("Unknown precision mode: $(opts.calc_precision)")
		end
    elseif opts.solver_type == :restricted
        if opts.calc_precision == :np
            if mu_rho_restricted === nothing
                error("mu_rho_restricted parameter must be provided for :np_restricted solver.")
            end
            LocalSolv.np_mu_rho_restricted_solver(N_range, t_n_range, mu_range, mu_rho_restricted, Delta_range, sequences, sequence_name, chunk_size, filepath, opts)
        else
            error(":restricted solver only supports normal precision (:np).")
        end
    else
        error("Unknown solver type: $(opts.solver_type)")
    end

end