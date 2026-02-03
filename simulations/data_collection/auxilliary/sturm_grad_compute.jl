include(joinpath(@__DIR__, "sequence_gen.jl"))

using .SeqGen
using BSON: @save, @load
using Random
using Printf
using Plots
using DataFrames

## Sturmian slope sampling
K = 100
L = 4
println("Generating Sturmian slopes with K=$K and L=$L.")

# ## Standard tail repeats method (all 1's)
# tail_length = 50
# sturmian_slopes_df = SeqGen.sample_sturmian_slopes(K, L; tail_repeats=tail_length, measure=true) # for standard tail repeats method

# ## Custom tail option method
# tail_option = 0 # "random" or <nothing> or a particular non-negative integer
# tail_length = 50
# sturmian_slopes_df = SeqGen.sample_sturmian_slopes_custom(K, L; tail_option=tail_option, tail_length=tail_length, measure=true) # for custom tail repeats behaviour

## Option to use optimised balancing before sequence generation
tail_option = 1 # "random" or <nothing> or a particular non-negative integer
tail_length = 500
balance_option = false
bins = 1000
max_per_bin = 1
sturmian_slopes_df = SeqGen.sample_sturmian_slopes_custom_optBalance(K, L; tail_option=tail_option, tail_length=tail_length, measure=true, balanced=balance_option, bins=bins, max_per_bin=max_per_bin)


println("Sampled $(length(sturmian_slopes_df.slope)) Sturmian slopes with K=$K, L=$L and tail_length=$tail_length.")
println("Raw slope range: [$(minimum(sturmian_slopes_df.slope)), $(maximum(sturmian_slopes_df.slope))]")

## Add symmetric slopes (1 - phi)
comp = false
if comp == true
    slopes = collect(sturmian_slopes_df[!, :slope])
    symmetric_slopes = 1 .- slopes
    sturmian_slopes_df = vcat(
        sturmian_slopes_df,
        DataFrame(slope = symmetric_slopes);
        cols = :union
    )
    println("Appended $(length(symmetric_slopes)) complementary slopes; total slopes = $(nrow(sturmian_slopes_df))")
end

post_balance_opt = true
if post_balance_opt == true
    # Rebalance again after adding complementary slopes
    bins = 5000
    max_per_bin = 2
    rng = MersenneTwister(42)
    balanced_sturm_df = SeqGen.rebalance_slopes(sturmian_slopes_df; bins=bins, max_per_bin=max_per_bin, rng=rng)
    println("Rebalanced to $(length(balanced_sturm_df.slope)) slopes after adding complementary slopes.")
    println("Balanced slope range: [$(minimum(balanced_sturm_df.slope)), $(maximum(balanced_sturm_df.slope))]")
end


## Visualise
function plt_sturmian_slope_sampling_bins(
    slopes_or_df,
    savepath::String;
    slope_col::Symbol = :slope,
    K::Union{Nothing,Int}=nothing,
    L::Union{Nothing,Int}=nothing,
    size::Tuple{Int,Int} = (800,200),
    bin_width::Float64 = 0.01,
    normalize::Bool = true
)
    slopes = slopes_or_df isa DataFrame ? collect(slopes_or_df[!, slope_col]) : collect(slopes_or_df)
    n = length(slopes)
    if n == 0
        error("No slopes to plot")
    end

    @assert bin_width > 0 "bin_width must be positive"

    # bin edges covering [0,1]
    edges = collect(0.0:bin_width:1.0)
    if edges[end] < 1.0
        push!(edges, 1.0)
    end
    n_bins = length(edges) - 1

    # count slopes per bin
    counts = zeros(Int, n_bins)
    for s in slopes
        s_clamped = clamp(s, 0.0, 1.0)
        bin_idx = min(searchsortedlast(edges, s_clamped), n_bins)
        counts[bin_idx] += 1
    end

    # prepare bar plot values
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    heights = normalize ? Float64.(counts) ./ max(1, sum(counts)) : Float64.(counts)

    plt = Plots.bar(
        centers, heights;
        bar_width = bin_width * 0.9,
        xlabel = "slope (φ)",
        ylabel = normalize ? "relative frequency" : "count",
        title = "Binned Sturmian slope distribution",
        xlim = (0.0, 1.0),
        grid = true,
        size = size,
        legend = false
    )

    # annotation / key
    parts = String[]
    if K !== nothing
        push!(parts, "K=$(K)")
    end
    if L !== nothing
        push!(parts, "L=$(L)")
    end
    push!(parts, "n=$(n)")
    push!(parts, @sprintf("min=%.4f", minimum(slopes)))
    push!(parts, @sprintf("max=%.4f", maximum(slopes)))
    push!(parts, @sprintf("bin=%.4f", bin_width))
    push!(parts, normalize ? "freq=relative" : "freq=counts")
    key_text = join(parts, "\n")
    Plots.annotate!(plt, (0.98, 0.98, Plots.text(key_text, :right, 9)))

    Plots.savefig(plt, savepath)
    # Plots.display(plt)
    # return plt
end

outdir = joinpath(@__DIR__, "sturm_grad_sets")

# plt_sturmian_slope_sampling_bins(
#     sturmian_slopes_df,
#     joinpath(outdir, "sturmian_slope_sampling_bins_K$(K)_L$(L)_r$(tail_length)_comp-$(comp).png");
#     K=K,
#     L=L
# )

# plt_sturmian_slope_sampling_bins(
#     balanced_sturm_df,
#     joinpath(outdir, "balanced_sturmian_slope_sampling_bins_K$(K)_L$(L)_r$(tail_length)_comp-$(comp).png");
#     K=K,
#     L=L
# )

# balanced_sturm_df = sturmian_slopes_df

## save to BSON file
t0 = time_ns()
bal_name = "sturmian_slopes_K$(K)_L$(L)_balanced_bins$(bins)_mpb$(max_per_bin)_r$(tail_length)_comp-$(comp)_tailoption-$(tail_option).bson"
@save joinpath(outdir, bal_name) balanced_sturm_df
println("Saved balanced Sturmian slopes (took $((time_ns() - t0)/1e9) s) to $(joinpath(outdir, bal_name))")
# t0 = time_ns()
# raw_name = "sturmian_slopes_K$(K)_L$(L)_r$(tail_length)_comp-$(comp)_tailoption-$(tail_option).bson"
# @save joinpath(outdir, raw_name) sturmian_slopes_df
# println("Saved raw Sturmian slopes (took $((time_ns() - t0)/1e9) s) to $(joinpath(outdir, raw_name))")