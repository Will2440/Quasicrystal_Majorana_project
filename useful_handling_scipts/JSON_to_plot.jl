using JSON, Plots

function plot_eigenvalues_from_json(filename::String)
    # Read JSON file
    results = JSON.parsefile(filename)

    # Extract eigenvalues
    all_eigenvalues = [res["eigenvalues"] for res in results]

    # Determine the number of parameter sets
    num_sets = length(all_eigenvalues)

    # Create a plot
    plot(title="Eigenvalues of BdG Hamiltonian", xlabel="Index", ylabel="Eigenvalue", legend=false)

    # Plot eigenvalues for each parameter set
    for i in 1:num_sets
        scatter!(1:length(all_eigenvalues[i]), all_eigenvalues[i], marker=:circle, alpha=0.5)
    end

    display(plot!())  # Show plot
end

# Example usage:
filename = "simulations/system_eigensolver/results.json"
plot_eigenvalues_from_json(filename)
