module PhasonInvariantHeatmap

using DataFrames
using Glob
using JLD2
using Plots
using Colors: RGB

const project_root = normpath(joinpath(@__DIR__, "../../.."))

const data_path = joinpath(
    project_root,
    "simulations/raw_data/specloc_phason_winding_recalc/hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-501-1.0_nbins1000_npb1_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-501)_Delta(0.05-0.05-1)_kappa(0.00025-0.00025-1)_x0s(250.0-250.0-1)"
)
const output_path = joinpath(
    project_root,
    "simulations/results/specloc_phason_winding_recalc",
    basename(data_path)
)

const target_N = 500
const target_t_n = [1.0, 1.5]
const target_Delta = 0.05
const target_kappa = 2.5e-4
const target_x0 = 250.0
const invariant_col = :specloc_winding

function load_jld2_results(folder_path::String)
    isdir(folder_path) || error("Data folder does not exist: $folder_path")
    file_paths = sort(glob("*.jld2", folder_path))
    isempty(file_paths) && error("No JLD2 files found in: $folder_path")

    dataframes = DataFrame[]
    for file_path in file_paths
        data = load(file_path)
        haskey(data, "results_df") || error("Missing results_df in $file_path")
        push!(dataframes, DataFrame(data["results_df"]))
    end

    println("Loaded $(length(dataframes)) JLD2 chunks from $folder_path")
    return vcat(dataframes...; cols=:union)
end

function same_hoppings(value, target; atol::Float64=1e-10)
    return value isa AbstractVector && length(value) == length(target) &&
           all(isapprox(a, b; atol=atol, rtol=0.0) for (a, b) in zip(value, target))
end

function log_min_abs_specloc_eigenvalue(spectrum, low_evals)
    eigenvalues = spectrum isa AbstractVector && !isempty(spectrum) ? spectrum : low_evals
    eigenvalues isa AbstractVector && !isempty(eigenvalues) || return NaN
    finite_abs_eigenvalues = filter(isfinite, abs.(eigenvalues))
    isempty(finite_abs_eigenvalues) && return NaN
    minimum_abs_eigenvalue = minimum(finite_abs_eigenvalues)
    return log10(minimum_abs_eigenvalue + 1e-15)
end

function build_invariant_grid(rows::DataFrame; value_col::Symbol)
    required_columns = (:mu, :sequence_id, value_col)
    missing_columns = filter(!in(propertynames(rows)), required_columns)
    isempty(missing_columns) || error("Missing required columns: $missing_columns")

    mu_values = sort(unique(Float64.(rows.mu)))
    phason_values = sort(unique(Float64(row.sequence_id[2]) for row in eachrow(rows)))
    mu_index = Dict(value => index for (index, value) in enumerate(mu_values))
    phason_index = Dict(value => index for (index, value) in enumerate(phason_values))
    values = fill(NaN, length(phason_values), length(mu_values))

    for row in eachrow(rows)
        sequence_id = row.sequence_id
        sequence_id isa Tuple && length(sequence_id) >= 2 ||
            error("Expected sequence_id to be a (phi, phason) tuple; got $sequence_id")

        mu = Float64(row.mu)
        phason = Float64(sequence_id[2])
        value = row[value_col]
        (ismissing(value) || !isfinite(value)) && continue

        row_index = phason_index[phason]
        col_index = mu_index[mu]
        previous = values[row_index, col_index]
        if isfinite(previous)
            isapprox(previous, value; atol=1e-10, rtol=1e-8) ||
                error("Conflicting invariant values at mu=$mu, phason=$phason: $previous and $value")
        else
            values[row_index, col_index] = Float64(value)
        end
    end

    return mu_values, phason_values, values
end

function main()
    mkpath(output_path)
    println("Loading spectral-localiser data from: $data_path")
    data = load_jld2_results(data_path)

    required_columns = (:N, :t_n, :Delta, :kappa, :x0, :sequence_id, invariant_col)
    all(column -> column in propertynames(data), required_columns) ||
        error("Input data is missing columns required to filter and plot: $required_columns")
    any(column -> column in propertynames(data), (:specloc_spectrum, :specloc_low_evals)) ||
        error("Input data must contain :specloc_spectrum or :specloc_low_evals to plot the spectral-localiser eigenvalue gap")

    rows = filter(data) do row
        row.N == target_N &&
        same_hoppings(row.t_n, target_t_n) &&
        isapprox(row.Delta, target_Delta; atol=1e-10, rtol=0.0) &&
        isapprox(row.kappa, target_kappa; atol=1e-12, rtol=1e-8) &&
        isapprox(row.x0, target_x0; atol=1e-8, rtol=0.0)
    end
    isempty(rows) && error("No rows matched N=$target_N, t_n=$target_t_n, Delta=$target_Delta, kappa=$target_kappa, x0=$target_x0")

    mu_values, phason_values, invariant_values = build_invariant_grid(rows; value_col=invariant_col)
    spectra = :specloc_spectrum in propertynames(rows) ? rows.specloc_spectrum : fill(missing, nrow(rows))
    low_evals = :specloc_low_evals in propertynames(rows) ? rows.specloc_low_evals : fill(missing, nrow(rows))
    rows.log_min_abs_specloc_eigenvalue = [
        log_min_abs_specloc_eigenvalue(spectra[index], low_evals[index])
        for index in 1:nrow(rows)
    ]
    _, _, log_min_abs_eigenvalue_values = build_invariant_grid(rows; value_col=:log_min_abs_specloc_eigenvalue)
    finite_invariants = sort(unique(filter(isfinite, vec(invariant_values))))
    isempty(finite_invariants) && error("No finite values found in $invariant_col for the selected rows")
    finite_log_min_abs_eigenvalues = filter(isfinite, vec(log_min_abs_eigenvalue_values))
    isempty(finite_log_min_abs_eigenvalues) && error("No finite spectral-localiser eigenvalues found for the selected rows")
    missing_cells = count(!isfinite, invariant_values)
    if missing_cells > 0
        @warn "Selected data do not fill the full mu-phason grid" missing_cells total_cells=length(invariant_values)
    end

    invariant_colormap = cgrad([RGB(1.0, 1.0, 1.0), RGB(181 / 255, 193 / 255, 211 / 255)])
    dark_blue_colormap = cgrad([RGB(1.0, 1.0, 1.0), RGB(7 / 255, 47 / 255, 107 / 255)])
    plot_title = "Spectral-localiser invariant (N=$target_N, t_n=$target_t_n, Delta=$target_Delta)"
    figure = heatmap(
        mu_values,
        phason_values,
        invariant_values;
        xlabel="Chemical potential, mu",
        ylabel="Phason",
        title=plot_title,
        color=invariant_colormap,
        clims=(0.0, 1.0),
        colorbar_title=string(invariant_col),
        colorbar_ticks=([0.0, 1.0], ["0", "1"]),
        framestyle=:box,
        grid=false,
        size=(1100, 650)
    )

    filename = joinpath(output_path, "specloc_winding_mu_vs_phason.svg")
    savefig(figure, filename)

    dark_blue_figure = heatmap(
        mu_values,
        phason_values,
        invariant_values;
        xlabel="Chemical potential, mu",
        ylabel="Phason",
        title=plot_title,
        color=dark_blue_colormap,
        clims=(0.0, 1.0),
        colorbar_title=string(invariant_col),
        colorbar_ticks=([0.0, 1.0], ["0", "1"]),
        framestyle=:box,
        grid=false,
        size=(1100, 650)
    )
    dark_blue_filename = joinpath(output_path, "specloc_winding_mu_vs_phason_dark_blue.svg")
    savefig(dark_blue_figure, dark_blue_filename)

    log_gap_colormap = cgrad(:viridis)
    log_gap_limits = extrema(finite_log_min_abs_eigenvalues)
    if log_gap_limits[1] == log_gap_limits[2]
        padding = max(abs(log_gap_limits[1]) * 0.05, 0.1)
        log_gap_limits = (log_gap_limits[1] - padding, log_gap_limits[2] + padding)
    end
    log_gap_figure = heatmap(
        mu_values,
        phason_values,
        log_min_abs_eigenvalue_values;
        xlabel="Chemical potential, mu",
        ylabel="Phason",
        title="Log minimum absolute spectral-localiser eigenvalue (N=$target_N, t_n=$target_t_n, Delta=$target_Delta)",
        color=log_gap_colormap,
        clims=log_gap_limits,
        colorbar_title="log10(min |eigenvalue|)",
        framestyle=:box,
        grid=false,
        size=(1100, 650)
    )
    log_gap_filename = joinpath(output_path, "specloc_log_min_abs_eigenvalue_mu_vs_phason.svg")
    savefig(log_gap_figure, log_gap_filename)

    println("Selected rows: $(nrow(rows)); grid: $(length(mu_values)) mu values x $(length(phason_values)) phasons")
    println("Invariant values: $finite_invariants")
    println("Log minimum absolute eigenvalue range: $(extrema(finite_log_min_abs_eigenvalues))")
    println("Saved heatmap to: $filename")
    println("Saved dark-blue heatmap to: $dark_blue_filename")
    println("Saved log minimum absolute eigenvalue heatmap to: $log_gap_filename")
end

main()

end