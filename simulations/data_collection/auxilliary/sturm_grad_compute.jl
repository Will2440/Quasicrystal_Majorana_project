include(joinpath(@__DIR__, "sequence_gen.jl"))

using .SeqGen
using BSON: @save, @load
using Random
using Printf
using Plots
using DataFrames

## Sturmian slope sampling
K = 8
L = 3
tail_repeats = 50
sturmian_slopes_df = SeqGen.sample_sturmian_slopes(K, L; tail_repeats=tail_repeats, measure=true)
println("Sampled $(length(sturmian_slopes_df.slope)) Sturmian slopes with K=$K, L=$L and tail_repeats=$tail_repeats.")

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

bins = 500
max_per_bin = 1
rng = MersenneTwister(42)
balanced_sturm_df = SeqGen.rebalance_slopes(sturmian_slopes_df; bins=bins, max_per_bin=max_per_bin, rng=rng)
println("Rebalanced to $(length(balanced_sturm_df.slope)) slopes after removing overdense areas.")

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

plt_sturmian_slope_sampling_bins(
    sturmian_slopes_df,
    joinpath(outdir, "sturmian_slope_sampling_bins_K$(K)_L$(L)_r$(tail_repeats)_comp-$(comp).png");
    K=K,
    L=L
)

plt_sturmian_slope_sampling_bins(
    balanced_sturm_df,
    joinpath(outdir, "balanced_sturmian_slope_sampling_bins_K$(K)_L$(L)_r$(tail_repeats)_comp-$(comp).png");
    K=K,
    L=L
)


## save to BSON file
@save joinpath(outdir, "sturmian_slopes_K$(K)_L$(L)_balanced_bins$(bins)_mpb$(max_per_bin)_r$(tail_repeats)_comp-$(comp).bson") balanced_sturm_df
println("Saved balanced Sturmian slopes to $(joinpath(outdir, "sturmian_slopes_K$(K)_L$(L)_balanced.bson"))")
@save joinpath(outdir, "sturmian_slopes_K$(K)_L$(L)_r$(tail_repeats)_comp-$(comp).bson") sturmian_slopes_df
println("Saved raw Sturmian slopes to $(joinpath(outdir, "sturmian_slopes_K$(K)_L$(L).bson"))")