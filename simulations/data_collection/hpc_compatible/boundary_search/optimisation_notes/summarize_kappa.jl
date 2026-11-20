using CSV, DataFrames
include("kappa_analysis.jl")
using .KappaAnalysis

length(ARGS) == 1 || error("Usage: summarize_kappa.jl kappa_map_directory")
directory = abspath(only(ARGS))
for path in sort(filter(path -> endswith(path, ".csv") && !endswith(path, "_summary.csv") &&
    !endswith(path, "_windows.csv"), readdir(directory; join=true)))
    data = CSV.read(path, DataFrame)
    sort!(data, [:u, :kappa])
    kappas = sort(unique(data.kappa))
    reference_kappa = only(filter(kappa -> isapprox(kappa, 1e-4; rtol=1e-10), kappas))
    reference = sort(data[data.kappa .== reference_kappa, :], :u)
    point_rows = NamedTuple[]
    for group in groupby(data, :u)
        group = sort(group, :kappa)
        nrow(group) == length(kappas) || error("Incomplete kappa sweep at u=$(group.u[1])")
        windows = positive_windows(group.phase_label)
        for (index, window) in enumerate(windows)
            push!(point_rows, (u=group.u[1], mu=group.mu[1], reference_phase=only(group[group.kappa .== reference_kappa, :phase_label]),
                window_index=index, first_kappa=group.kappa[window.first], last_kappa=group.kappa[window.last],
                lower_censored=window.lower_censored, upper_censored=window.upper_censored,
                previous_kappa=window.first == 1 ? missing : group.kappa[window.first - 1],
                next_kappa=window.last == length(kappas) ? missing : group.kappa[window.last + 1]))
        end
    end
    summaries = NamedTuple[]
    for kappa in kappas
        current = sort(data[data.kappa .== kappa, :], :u)
        current.u == reference.u || error("Kappa grids differ")
        metrics = agreement_counts(reference.phase_label, current.phase_label)
        push!(summaries, merge((kappa=kappa, min_localizer_gap=minimum(current.localizer_gap)), metrics))
    end
    CSV.write(replace(path, ".csv" => "_summary.csv"), DataFrame(summaries))
    if !isempty(point_rows)
        CSV.write(replace(path, ".csv" => "_windows.csv"), DataFrame(point_rows))
    end
    println(basename(path), " reference topological=", count(==(1), reference.phase_label),
        " points=", nrow(reference), " agreement kappas=", [row.kappa for row in summaries if row.topological_lost == 0 && row.trivial_changed == 0],
        " topological-preserving kappas=", [row.kappa for row in summaries if row.topological_lost == 0])
end