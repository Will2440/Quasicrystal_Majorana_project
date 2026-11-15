using DataFrames
using Glob
using JLD2
using Plots

const project_root = normpath(joinpath(@__DIR__, "../../.."))

const mp_data_path = joinpath(
    project_root,
    "simulations/raw_data/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N(5-10000-8)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-601)_Delta(0p05-0p05-1)"
)
const specloc_data_path = joinpath(
    project_root,
    "simulations/raw_data/spec_loc_qc-sc_comp/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-6001)_Delta(0.05-0.05-1)_kappa(0.00025-0.00025-1)_x0s(250.0-250.0-1)"
)
const output_path = joinpath(
    project_root,
    "simulations/results/spec_loc_qc-sc_comp/fig_2_test_plotting/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-6001)_Delta(0.05-0.05-1)_kappa(0.00025-0.00025-1)_x0s(250.0-250.0-1)"
)

const min_mp_N = 200
const max_mp_N = 1000
const target_t_n = [1.0, 1.5]
const target_Delta = 0.05
const specloc_N = 500
const specloc_x0 = 250.0

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

    return vcat(dataframes...; cols=:union)
end

function same_hoppings(value, target; atol::Float64=1e-10)
    return value isa AbstractVector && length(value) == length(target) &&
           all(isapprox(a, b; atol=atol, rtol=0.0) for (a, b) in zip(value, target))
end

function padded_limits(values)
    finite_values = filter(isfinite, Float64.(values))
    isempty(finite_values) && error("Cannot determine plot limits from empty or non-finite data.")

    lower, upper = extrema(finite_values)
    padding = lower == upper ? max(abs(lower) * 0.05, 0.1) : 0.05 * (upper - lower)
    return (lower - padding, upper + padding)
end

function sequence_ids_match(left_ids, right_ids; atol::Float64=1e-8)
    isempty(left_ids) && return false
    isempty(right_ids) && return false

    left = first(left_ids)
    right = first(right_ids)
    return left isa Tuple && right isa Tuple && length(left) == length(right) &&
           all(isapprox(Float64(a), Float64(b); atol=atol, rtol=0.0) for (a, b) in zip(left, right))
end

function main()
    mkpath(output_path)
    println("Loading MP data from: $mp_data_path")
    mp_data = load_jld2_results(mp_data_path)
    println("Loading spectral-localiser data from: $specloc_data_path")
    specloc_data = load_jld2_results(specloc_data_path)

    mp_rows = filter(mp_data) do row
        min_mp_N <= row.N <= max_mp_N &&
        same_hoppings(row.t_n, target_t_n) &&
        isapprox(row.Delta, target_Delta; atol=1e-10, rtol=0.0)
    end
    isempty(mp_rows) && error("No MP rows matched N in [$min_mp_N, $max_mp_N] and the selected hopping/Delta values.")

    specloc_rows = filter(specloc_data) do row
        row.N == specloc_N &&
        same_hoppings(row.t_n, target_t_n) &&
        isapprox(row.Delta, target_Delta; atol=1e-10, rtol=0.0) &&
        isapprox(row.x0, specloc_x0; atol=1e-8, rtol=0.0)
    end
    isempty(specloc_rows) && error("No spectral-localiser rows matched N=$specloc_N, x0=$specloc_x0, t_n=$target_t_n, Delta=$target_Delta.")

    specloc_sequence_ids = unique(specloc_rows.sequence_id)
    length(specloc_sequence_ids) == 1 || error("Expected one spectral-localiser sequence ID, found $specloc_sequence_ids.")
    target_sequence_id = only(specloc_sequence_ids)
    mp_sequence_ids = unique(mp_rows.sequence_id)
    if !sequence_ids_match(mp_sequence_ids, [target_sequence_id])
        @warn "MP and spectral-localiser sequence IDs differ; continuing with all selected MP rows." mp_sequence_ids target_sequence_id
    end

    sort!(mp_rows, [:N, :mu])
    sort!(specloc_rows, :mu)
    mu_limits = (
        max(minimum(mp_rows.mu), minimum(specloc_rows.mu)),
        min(maximum(mp_rows.mu), maximum(specloc_rows.mu))
    )
    mu_limits[1] < mu_limits[2] || error("MP and spectral-localiser data have no overlapping chemical-potential range.")
    mp_system_sizes = sort(unique(mp_rows.N))
    kappas = sort(unique(specloc_rows.kappa))
    println("MP system sizes included: $mp_system_sizes")
    println("MP rows: $(nrow(mp_rows)); spectral-localiser rows: $(nrow(specloc_rows))")
    println("Shared chemical-potential range: $mu_limits")
    println("Spectral-localiser sequence ID: $target_sequence_id")
    println("Generating one comparison plot for each of $(length(kappas)) kappa values.")

    mp_limits = padded_limits(mp_rows.mp)
    spectrum_panel = plot(
        xlabel="",
        ylabel="Spectral-localiser eigenvalue",
        title="Near-zero spectral-localiser spectrum across κ",
        ylims=(-0.5, 0.5),
        legend=:outerright,
        grid=false
    )
    log_gap_panel = plot(
        xlabel="Chemical potential, μ",
        ylabel="log10(minimum |eigenvalue|)",
        title="Spectral-localiser minimum-eigenvalue gap across κ",
        ylims=(-8.0, 0.0),
        legend=:outerright,
        grid=false
    )
    kappa_colors = [get(cgrad(:viridis), (i - 1) / max(length(kappas) - 1, 1)) for i in eachindex(kappas)]

    for (kappa_index, kappa) in enumerate(kappas)
        kappa_rows = filter(row -> isapprox(row.kappa, kappa; atol=1e-12, rtol=1e-8), specloc_rows)
        isempty(kappa_rows) && continue
        sort!(kappa_rows, :mu)
        winding_limits = padded_limits(kappa_rows.specloc_winding)
        plt = plot(
            xlabel="Chemical potential, μ",
            ylabel="Raw Majorana polarization",
            title="Majorana polarization and spectral-localiser winding (κ=$kappa)",
            ylims=mp_limits,
            xlims=mu_limits,
            legend=false,
            linewidth=2,
            grid=false
        )

        for system_size in mp_system_sizes
            rows = filter(:N => ==(system_size), mp_rows)
            sort!(rows, :mu)
            plot!(
                plt,
                Float64.(rows.mu),
                Float64.(rows.mp);
                label="",
                linewidth=1.7
            )
        end

        winding_axis = twinx(plt)
        plot!(
            winding_axis,
            Float64.(kappa_rows.mu),
            Float64.(kappa_rows.specloc_winding);
            ylabel="Spectral-localiser winding",
            ylims=winding_limits,
            xlims=mu_limits,
            label="",
            color=:black,
            linewidth=2.0,
            linestyle=:solid,
            legend=false,
            grid=false
        )
        xlims!(plt, mu_limits)

        plot!(
            spectrum_panel,
            Float64[],
            Float64[];
            label="κ=$kappa",
            color=kappa_colors[kappa_index],
            markersize=1
        )
        spectrum_mu = Float64[]
        spectrum_values = Float64[]
        for row in eachrow(kappa_rows)
            if row.specloc_low_evals isa AbstractVector
                for eval in row.specloc_low_evals
                    if -0.5 <= eval <= 0.5
                        push!(spectrum_mu, row.mu)
                        push!(spectrum_values, eval)
                    end
                end
            end
        end
        plot!(
            spectrum_panel,
            spectrum_mu,
            spectrum_values;
            seriestype=:scatter,
            label="",
            color=kappa_colors[kappa_index],
            markersize=0.7,
            markerstrokewidth=0,
            markeralpha=0.25
        )

        log_gap = log10.(abs.(kappa_rows.specloc_gap) .+ 1e-15)
        plot!(
            log_gap_panel,
            kappa_rows.mu,
            log_gap;
            label="κ=$kappa",
            color=kappa_colors[kappa_index],
            linewidth=1.5
        )

        kappa_spectrum_panel = plot(
            xlabel="",
            ylabel="Spectral-localiser eigenvalue",
            title="Near-zero spectral-localiser spectrum (κ=$kappa)",
            xlims=mu_limits,
            ylims=(-0.5, 0.5),
            legend=false,
            grid=false
        )
        plot!(
            kappa_spectrum_panel,
            spectrum_mu,
            spectrum_values;
            seriestype=:scatter,
            label="",
            color=kappa_colors[kappa_index],
            markersize=0.7,
            markerstrokewidth=0,
            markeralpha=0.25
        )
        kappa_log_gap_panel = plot(
            Float64.(kappa_rows.mu),
            log_gap;
            xlabel="Chemical potential, μ",
            ylabel="log10(minimum |eigenvalue|)",
            title="Spectral-localiser minimum-eigenvalue gap (κ=$kappa)",
            xlims=mu_limits,
            ylims=(-8.0, 0.0),
            label="",
            color=kappa_colors[kappa_index],
            linewidth=1.5,
            grid=false
        )
        kappa_spectral_plot = plot(
            kappa_spectrum_panel,
            kappa_log_gap_panel;
            layout=(2, 1),
            link=:x,
            size=(1100, 850),
            left_margin=8Plots.mm,
            right_margin=10Plots.mm
        )

        plot!(plt; size=(1100, 650), left_margin=8Plots.mm, right_margin=8Plots.mm)

        kappa_tag = replace(string(kappa), "." => "p")
        output_file = joinpath(output_path, "mp_and_specloc_winding_vs_mu_kappa_$(kappa_tag).svg")
        savefig(plt, output_file)
        println("Saved comparison plot for kappa=$kappa: $output_file")
        kappa_spectrum_output = joinpath(output_path, "specloc_spectrum_and_log_gap_kappa_$(kappa_tag).png")
        savefig(kappa_spectral_plot, kappa_spectrum_output)
        println("Saved individual spectrum/log-gap plot for kappa=$kappa: $kappa_spectrum_output")
    end

    spectral_superplot = plot(
        spectrum_panel,
        log_gap_panel;
        layout=(2, 1),
        link=:x,
        size=(1300, 1000),
        left_margin=8Plots.mm,
        right_margin=10Plots.mm
    )
    spectrum_output = joinpath(output_path, "specloc_spectrum_and_log_gap_all_kappas.png")
    savefig(spectral_superplot, spectrum_output)
    println("Saved all-kappa spectrum/log-gap superplot: $spectrum_output")
end

main()