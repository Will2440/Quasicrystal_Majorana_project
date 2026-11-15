using CSV
using DataFrames
using Glob
using JLD2
using Statistics

const project_root = normpath(joinpath(@__DIR__, "../../.."))
const specloc_input = joinpath(
    project_root,
    "simulations/raw_data/spec_loc_qc-sc_comp/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-6001)_Delta(0.05-0.05-1)_kappa(0.00025-0.00025-1)_x0s(250.0-250.0-1)"
)
const mp_diagnostics = joinpath(
    project_root,
    "simulations/results/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N(5-10000-8)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-601)_Delta(0p05-0p05-1)/mp_vs_mu_N100_1000_theilsen_diagnostics.csv"
)
const output_dir = joinpath(
    project_root,
    "simulations/results/spec_loc_qc-sc_comp/slopes_N10000_literal_0.61803_phason_0.61803-1-0.61803_mapto(2, 1)_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-6001)_Delta(0.05-0.05-1)_kappa(0.00025-0.00025-1)_x0s(250.0-250.0-1)"
)

function load_specloc_rows(folder_path::String)
    isdir(folder_path) || error("Spectral-localiser input folder not found: $folder_path")
    files = sort(glob("*.jld2", folder_path))
    isempty(files) && error("No JLD2 files found in: $folder_path")
    chunks = [DataFrame(load(file)["results_df"]) for file in files]
    return vcat(chunks...; cols=:union)
end

function same_hoppings(value, target; atol::Float64=1e-10)
    return value isa AbstractVector && length(value) == length(target) &&
           all(isapprox(a, b; atol=atol, rtol=0.0) for (a, b) in zip(value, target))
end

function transition_brackets(mu_values, phases)
    length(mu_values) == length(phases) || error("mu and phase columns have different lengths")
    any(ismissing, phases) && error("Cannot compare transitions with missing phase values")
    order = sortperm(mu_values)
    mu = Float64.(mu_values[order])
    phase = Int.(phases[order])
    length(unique(mu)) == length(mu) || error("Duplicate mu values found in a phase series")
    all(value -> value in (0, 1), phase) || error("Expected phase values to be 0 or 1")

    brackets = NamedTuple[]
    for index in 2:length(mu)
        if phase[index] != phase[index - 1]
            push!(brackets, (
                lower=mu[index - 1],
                upper=mu[index],
                midpoint=(mu[index - 1] + mu[index]) / 2,
                phase_before=phase[index - 1],
                phase_after=phase[index]
            ))
        end
    end
    return mu, brackets
end

function main()
    isfile(mp_diagnostics) || error("MP diagnostics CSV not found: $mp_diagnostics")
    specloc_data = load_specloc_rows(specloc_input)
    selected = filter(specloc_data) do row
        row.N == 500 &&
        isapprox(row.x0, 250.0; atol=1e-8, rtol=0.0) &&
        isapprox(row.kappa, 0.00025; atol=1e-12, rtol=1e-8) &&
        isapprox(row.Delta, 0.05; atol=1e-10, rtol=0.0) &&
        same_hoppings(row.t_n, [1.0, 1.5]) &&
        row.sequence_id isa Tuple && length(row.sequence_id) >= 2 &&
        isapprox(row.sequence_id[1], 0.6180339887498948; atol=1e-10, rtol=0.0) &&
        isapprox(row.sequence_id[2], 0.6180339887498948; atol=1e-10, rtol=0.0)
    end
    isempty(selected) && error("No rows matched the requested spectral-localiser parameters")

    specloc_mu, specloc_transitions = transition_brackets(selected.mu, selected.specloc_winding)
    mp_data = CSV.read(mp_diagnostics, DataFrame)
    mp_mu, mp_transitions = transition_brackets(mp_data.mu, mp_data.phase)
    length(specloc_transitions) == length(mp_transitions) ||
        error("Transition counts differ: specloc=$(length(specloc_transitions)), MP=$(length(mp_transitions)); refusing index-only pairing")

    specloc_step = median(diff(specloc_mu))
    mp_step = median(diff(mp_mu))
    resolution_bound = (specloc_step + mp_step) / 2
    comparison = DataFrame(
        transition_index=Int[],
        specloc_phase_before=Int[],
        specloc_phase_after=Int[],
        specloc_mu_lower=Float64[],
        specloc_mu_upper=Float64[],
        specloc_mu_midpoint=Float64[],
        mp_phase_before=Int[],
        mp_phase_after=Int[],
        mp_mu_lower=Float64[],
        mp_mu_upper=Float64[],
        mp_mu_midpoint=Float64[],
        mp_minus_specloc_mu=Float64[],
        absolute_midpoint_disagreement=Float64[],
        combined_half_step_resolution=Float64[],
        midpoint_exceeds_resolution=Bool[],
        bracket_gap=Float64[],
        brackets_overlap=Bool[]
    )

    for index in eachindex(specloc_transitions)
        specloc = specloc_transitions[index]
        mp = mp_transitions[index]
        (specloc.phase_before, specloc.phase_after) == (mp.phase_before, mp.phase_after) ||
            error("Transition directions differ at index $index; refusing index-only pairing")
        gap = max(0.0, max(specloc.lower, mp.lower) - min(specloc.upper, mp.upper))
        delta = mp.midpoint - specloc.midpoint
        push!(comparison, (
            index,
            specloc.phase_before,
            specloc.phase_after,
            specloc.lower,
            specloc.upper,
            specloc.midpoint,
            mp.phase_before,
            mp.phase_after,
            mp.lower,
            mp.upper,
            mp.midpoint,
            delta,
            abs(delta),
            resolution_bound,
            abs(delta) > resolution_bound + 1e-12,
            gap,
            gap <= 1e-12
        ))
    end

    mkpath(output_dir)
    output_csv = joinpath(output_dir, "mp_vs_specloc_transition_comparison.csv")
    CSV.write(output_csv, comparison)
    println("Spectral-localiser points: $(length(specloc_mu)), mu=$(extrema(specloc_mu)), step=$specloc_step")
    println("MP diagnostic points: $(length(mp_mu)), mu=$(extrema(mp_mu)), step=$mp_step")
    println("Transitions: $(nrow(comparison)); matching directions: $(count(comparison.specloc_phase_before .== comparison.mp_phase_before))")
    println("Mean absolute midpoint disagreement: $(mean(comparison.absolute_midpoint_disagreement))")
    println("Maximum absolute midpoint disagreement: $(maximum(comparison.absolute_midpoint_disagreement))")
    println("Transitions exceeding half-step resolution bound ($resolution_bound): $(count(comparison.midpoint_exceeds_resolution))")
    println("Transitions with non-overlapping sampled brackets: $(count(.!comparison.brackets_overlap))")
    println("Saved: $output_csv")
end

main()
