
include("functions.jl")

using .HofstadterButterfly
using CairoMakie
using Colors
using BSON

actual_butterfly_results_path = joinpath(@__DIR__, "..", "results")
data_path = joinpath(@__DIR__, "..", "raw_data")

function generate_custom_colormap(
    min_val::Int, max_val::Int
)
    qs = collect(min_val:max_val)
    
    lerp(c1::RGB, c2::RGB, t) = RGB((1-t)*c1.r + t*c2.r,
                                    (1-t)*c1.g + t*c2.g,
                                    (1-t)*c1.b + t*c2.b)

    cmap = fill(RGB(0.7,0.7,0.7), length(qs))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs)
    idx_pos1     = findfirst(==(  1), qs)
    idx_pos_high = length(qs)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1 (Light Blue)
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    # Set anchors
    cmap[idx_neg_high] = col_neg_high
    cmap[idx_pos_high] = col_pos_high
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end # Red for +1

    # Interpolate Negative Side (Min -> -1)
    if !isnothing(idx_neg1) && idx_neg1 > idx_neg_high
        span = idx_neg1 - idx_neg_high
        for (j,i) in enumerate(idx_neg_high+1:idx_neg1-1)
            t = j / span
            if t < 1/3
                cmap[i] = lerp(col_neg_high, col_neg_purp, t/(1/3))
            elseif t < 2/3
                cmap[i] = lerp(col_neg_purp, col_neg_dark, (t-1/3)/(1/3))
            else
                cmap[i] = lerp(col_neg_dark, col_neg_light, (t-2/3)/(1/3))
            end
        end
    end

    # Interpolate Positive Side (1 -> Max)
    if !isnothing(idx_pos1) && idx_pos_high > idx_pos1
        span = idx_pos_high - idx_pos1
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high-1)
            t = j / span
            if t < 0.5
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.5)
            else
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.5)/0.5)
            end
        end
    end

    # Interpolate Middle (-1 -> 1)
    if !isnothing(idx_neg1) && !isnothing(idx_pos1) && idx_pos1 - idx_neg1 > 1
        span = idx_pos1 - idx_neg1
        for (j,i) in enumerate(idx_neg1+1:idx_pos1-1)
            t = j / span
            cmap[i] = lerp(col_neg_light, RGB(1.0,0.0,0.0), t)
        end
    end

    return cmap
end



# # ####################################################################
# # ######## calculate and plot the Hofstadter butterfly spectrum ######
# # ####################################################################
# collect_data = true  # Set to false to load existing data

# if collect_data
#     qmax = 30
#     flux, energy = butterfly_spectrum(qmax)
#     data = Dict("flux" => flux, "energy" => energy)
#     BSON.bson(joinpath(data_path, "butterfly_qmax$(qmax).bson"), data)
# else
#     data = BSON.load(joinpath(data_path, "butterfly_qmax30.bson"))
#     flux = data["flux"]
#     energy = data["energy"]
# end

# f = CairoMakie.Figure(resolution=(800,600))
# ax = CairoMakie.Axis(f[1,1], xlabel="Energy", ylabel="Flux α", title="Hofstadter Butterfly")

# CairoMakie.scatter!(ax, energy, flux; markersize=1, color=:black)
# CairoMakie.save(joinpath(actual_butterfly_results_path, "butterfly.png"), f)




# ###########################################################################
# ### calculate and plot the Chern-coloured Hofstadter butterfly spectrum ###
# ###########################################################################
# collect_data = true  # Set to false to load existing data

# if collect_data
#     qmax = 50
#     nk = div(qmax, 2) 
#     flux, energy, colour = HofstadterButterfly.coloured_butterfly(qmax; nk=nk, tx=1.0, ty=1.0)
#     data = Dict("flux" => flux, "energy" => energy, "colour" => colour)
#     BSON.bson(joinpath(data_path, "chern_butterfly_qmax$(qmax)_nk$(nk).bson"), data)
# else
#     data = BSON.load(joinpath(data_path, "chern_butterfly_qmax50_nk25.bson"))
#     flux = data["flux"]
#     energy = data["energy"]
#     colour = data["colour"]
# end

# # FILTERING:
# # 1. Limit to small Chern numbers to remove numerical noise (artifacts from low nk grid)
# # 2. Most physically relevant gaps have C = 0, +/-1, +/-2
# max_chern = div(qmax, 2)  # Keep Chern numbers with |C| <= max_chern
# mask = abs.(colour) .<= max_chern

# flux_filt   = flux[mask]
# energy_filt = energy[mask]
# colour_filt = colour[mask]

# f = CairoMakie.Figure(resolution=(900,600))
# ax = CairoMakie.Axis(f[1,1], xlabel="Energy", ylabel="Flux α", title="Chern-Coloured Butterfly (|C| ≤ $max_chern)")

# # Use a discrete color map logic
# # We set the colorrange to ensure 0 is white/neutral and integers map to distinct colors
# cmap = :RdBu # Red-Blue diverging map
# clims = (-max_chern - 0.5, max_chern + 0.5)

# sc = CairoMakie.scatter!(ax, energy_filt, flux_filt; 
#     markersize=2, 
#     color=colour_filt, 
#     colormap=cmap, 
#     colorrange=clims
# )

# # Add a Colorbar with integer ticks (acting as a Legend)
# CairoMakie.Colorbar(f[1,2], sc, 
#     label="Chern Number", 
#     ticks = [-max_chern, 0, max_chern]
# )

# CairoMakie.save(joinpath(actual_butterfly_results_path, "chern_butterfly_qmax$(qmax)_nk$(nk)_maxchern$(max_chern).png"), f)


# #########################################################################
# ### calculate and plot the gap-coloured Hofstadter butterfly spectrum ###
# #########################################################################
# collect_data = true  # Set to false to load existing data

# if collect_data
#     # 1. Calculate with higher resolution (qmax) to find the gaps
#     qmax = 50
#     flux, gmin, gmax, glabel = HofstadterButterfly.coloured_butterfly_gaps_threaded(qmax; tx=1.0, ty=1.0)
#     data = Dict("flux" => flux, "gmin" => gmin, "gmax" => gmax, "glabel" => glabel)
#     BSON.bson(joinpath(data_path, "gap_coloured_butterfly_qmax$(qmax).bson"), data)
# else
#     data = BSON.load(joinpath(data_path, "gap_coloured_butterfly_qmax50.bson"))
#     flux = data["flux"]
#     gmin = data["gmin"]
#     gmax = data["gmax"]
#     glabel = data["glabel"]
# end

# # Recover q for filtering (since flux is p/q)
# # We use rationalize to recover the exact denominator used
# qs = [denominator(rationalize(f)) for f in flux]

# # 2. Filter: 
# #    a. "Exclude the rationals": Filter out simple rationals (small q) to see "irrational approximants"
# #    b. "Properly labelled gaps": Keep gap labels small relative to resolution
# min_q = 10       # Exclude simple fractions like 1/2, 1/3 ... 1/min_q
# max_label = qmax  # Only show gaps with Hall conductance <= max_label

# mask = (qs .> min_q) .& (abs.(glabel) .<= max_label)

# flux_filt = flux[mask]
# gmin_filt = gmin[mask]
# gmax_filt = gmax[mask]
# glabel_filt = glabel[mask]

# f = CairoMakie.Figure(resolution=(1000, 800))
# ax = CairoMakie.Axis(f[1,1], 
#     xlabel="Energy E", 
#     ylabel="Flux α = p/q", 
#     title="Hofstadter Butterfly (Irrational Approximants q > $min_q, |σ_xy| ≤ $max_label)",
#     backgroundcolor=:white
# )

# # 3. Calculate color limits explicitly
# clims = (-max_label, max_label)

# # Use rangebars to fill the vertical extent of the gaps
# rb = CairoMakie.rangebars!(ax, flux_filt, gmin_filt, gmax_filt; 
#     color=glabel_filt, 
#     colormap=:Spectral, 
#     colorrange=clims, 
#     linewidth=3.0 
# )

# CairoMakie.Colorbar(f[1, 2], colormap=:Spectral, colorrange=clims, label="Hall Conductance σ_xy")
# CairoMakie.save(joinpath(actual_butterfly_results_path, "gap_coloured_butterfly_qmax$(qmax)_ratqs$(min_q)_maxlabel$(max_label).png"), f)






#########################################################################
### calculate and plot the gap-coloured Hofstadter butterfly spectrum directly from diophantine solution (no berry  curvature integration) ###
#########################################################################

# ------------------------------------------
# Data Collection
# ------------------------------------------
collect_data = true  # Set to false to load existing data

qmax = 500

if collect_data
    # 1. Calculate with higher resolution (qmax) to find the gaps
    flux, gmin, gmax, glabel = HofstadterButterfly.coloured_butterfly_diophantine_gaps(qmax; tx=1.0, ty=1.0)
    data = Dict("flux" => flux, "gmin" => gmin, "gmax" => gmax, "glabel" => glabel)
    BSON.bson(joinpath(data_path, "gap_coloured_butterfly_diophantine_qmax$(qmax).bson"), data)
else
    data = BSON.load(joinpath(data_path, "gap_coloured_butterfly_diophantine_qmax$(qmax).bson"))
    flux = data["flux"]
    gmin = data["gmin"]
    gmax = data["gmax"]
    glabel = data["glabel"]
end
# ------------------------------------------


# ------------------------------------------
# Filtering
# ------------------------------------------
# Recover q for filtering (since flux is p/q)
# We use rationalize to recover the exact denominator used
qs = [denominator(rationalize(f)) for f in flux]

# 2. Filter: 
#    a. "Exclude the rationals": Filter out simple rationals (small q) to see "irrational approximants"
#    b. "Properly labelled gaps": Keep gap labels small relative to resolution
min_q = 10       # Exclude simple fractions like 1/2, 1/3 ... 1/min_q
max_label = 20 #qmax  # Only show gaps with Hall conductance <= max_label

mask = (qs .> min_q) #.& (abs.(glabel) .<= max_label)

flux_filt = flux[mask]
gmin_filt = gmin[mask]
gmax_filt = gmax[mask]
glabel_filt = glabel[mask]
# ------------------------------------------



# ------------------------------------------
# Plotting unormalised colormap
# ------------------------------------------
f = CairoMakie.Figure(resolution=(1000, 800))
ax = CairoMakie.Axis(f[1,1], 
    xlabel="Energy E", 
    ylabel="Flux α = p/q", 
    title="Hofstadter Butterfly (Irrational Approximants q > $min_q, |σ_xy| ≤ $max_label)",
    backgroundcolor=:white
)

# 3. Calculate color limits explicitly
clims = (-max_label, max_label)
custom_cmap = generate_custom_colormap(-max_label, max_label)

# Use linesegments to plot horizontal bars (swapped axes: energy on x, flux on y)
points = Vector{Point2f}()
colors = Vector{Int}()
for i in 1:length(flux_filt)
    push!(points, Point2f(gmin_filt[i], flux_filt[i]))
    push!(points, Point2f(gmax_filt[i], flux_filt[i]))
    push!(colors, glabel_filt[i])
    push!(colors, glabel_filt[i])  # Repeat color for each segment
end

CairoMakie.linesegments!(ax, points; 
    color=colors, 
    colormap=custom_cmap, 
    colorrange=clims, 
    linewidth=3.0 
)

CairoMakie.Colorbar(f[1, 2], colormap=custom_cmap, colorrange=clims, label="Hall Conductance σ_xy")
CairoMakie.save(joinpath(actual_butterfly_results_path, "gap_coloured_butterfly_diophantine_qmax$(qmax)_ratqs$(min_q)_maxlabel$(max_label).png"), f)



# ------------------------------------------
# Normalisation
# ------------------------------------------
# Normalise energy range per flux to [-1,1] for better visualisation
unique_fluxes = unique(flux_filt)
norm_gmin_filt = similar(gmin_filt)
norm_gmax_filt = similar(gmax_filt)

for uf in unique_fluxes
    indices = findall(flux_filt .== uf)
    local_gmins = gmin_filt[indices]
    local_gmaxs = gmax_filt[indices]
    local_min = minimum(local_gmins)
    local_max = maximum(local_gmaxs)
    
    # Avoid division by zero if local_min == local_max
    if local_min == local_max
        norm_gmin_filt[indices] .= 0.0
        norm_gmax_filt[indices] .= 0.0
    else
        norm_gmin_filt[indices] = (local_gmins .- local_min) ./ (local_max - local_min) .* 2 .- 1
        norm_gmax_filt[indices] = (local_gmaxs .- local_min) ./ (local_max - local_min) .* 2 .- 1
    end
end


# ------------------------------------------
# Plotting normalised colormap
# ------------------------------------------
f = CairoMakie.Figure(resolution=(1000, 800))
ax = CairoMakie.Axis(f[1,1], 
    xlabel="Normalised Energy E", 
    ylabel="Flux α = p/q", 
    title="Hofstadter Butterfly (Irrational Approximants q > $min_q, |σ_xy| ≤ $max_label)",
    backgroundcolor=:white
)

# Calculate color limits explicitly
clims = (-max_label, max_label)
custom_cmap = generate_custom_colormap(-max_label, max_label)

# Use linesegments with normalized energies
points_norm = Vector{Point2f}()
colors_norm = Vector{Int}()
for i in 1:length(flux_filt)
    push!(points_norm, Point2f(norm_gmin_filt[i], flux_filt[i]))
    push!(points_norm, Point2f(norm_gmax_filt[i], flux_filt[i]))
    push!(colors_norm, glabel_filt[i])
    push!(colors_norm, glabel_filt[i])  # Repeat color for each segment
end

CairoMakie.linesegments!(ax, points_norm; 
    color=colors_norm, 
    colormap=custom_cmap, 
    colorrange=clims, 
    linewidth=3.0 
)

CairoMakie.Colorbar(f[1, 2], colormap=custom_cmap, colorrange=clims, label="Hall Conductance σ_xy")
CairoMakie.save(joinpath(actual_butterfly_results_path, "gap_coloured_butterfly_diophantine_normalised_qmax$(qmax)_ratqs$(min_q)_maxlabel$(max_label).png"), f)