
include("bson_unpacker.jl")
include("standard_plotting.jl")

using DataFrames, Images, ImageMorphology, ImageFiltering, Statistics, Clustering, Plots

filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/mu_vs_rho_mp_heatmaps/GQC_N(50-50-1)_t1(1.0-1.0-101__t2(0.0-10.0-101)_mu(0.0-10.0-101)_Delta(0.1-0.1-1)/"
mp_tol = 0.9
df = unpack_bason_standard(filepath; mp_tol=mp_tol)

struct BandInfo
    id::Int
    slope_mean::Float64
    slope_std::Float64
    intercept_mean::Float64
    npoints::Int
    x_coords::Vector{Int}
    y_coords::Vector{Int}
end



"""
    extract_band_slopes(df; grid_res=(200,200), neigh=5, nclusters=nothing)

Takes dataframe with columns :mu_t, :rho, :mp_disc.
Returns cluster means of band slopes (dy/dx).

Arguments:
- grid_res :: tuple (nx, ny): resolution for gridding/interpolation
- neigh :: Int: half-window for local slope fit along skeleton
- nclusters :: Int or nothing: if known number of bands, force k-means clustering.
              Otherwise DBSCAN is used to auto-detect clusters.
"""
function extract_band_slopes(df; grid_res=(200,200), neigh=10, nclusters=4)
    # --- 1. make regular grid ---
    mu_vals = unique(sort(df.mu_t))
    rho_vals = unique(sort(df.rho))
    nx, ny = length(mu_vals), length(rho_vals)

    println("Data grid size: $nx x $ny")

    mu_index = Dict(mu_vals .=> 1:nx)
    rho_index = Dict(rho_vals .=> 1:ny)

    mp_grid = fill(0, nx, ny)
    for row in eachrow(df)
        i = mu_index[row.mu_t]
        j = rho_index[row.rho]
        mp_grid[i,j] = row.mp_disc == -1 ? 1 : 0   # binary mask
    end

    println("type of mp_grid: ", typeof(mp_grid))
    # println(mp_grid)

    # plt_df = DataFrame(mu_vals=mu_vals, rho_vals=rho_vals, mp_mask=mp_grid[:,:])
    # # println(plt_df)
    # filename = "heatmap_gradient_test.png"
    # plt_mp_heatmap(plt_df, :mp_mask, :mu_vals, :rho_cals, filename)

    plt = Plots.heatmap(
        mu_vals, rho_vals, mp_grid',
        xlabel = "mu_vals",
        ylabel = "rho_vals",
        title = "Binary Mask Heatmaps",
        # color = [:yellow :blue],  # 0 = white, 1 = black
        legend = false
    )
    Plots.savefig("binary_mask_heatmap.png")

    # --- 2. skeletonize ---
    # skel = thinning(mp_grid .> 0)   # boolean skeleton image
    mask = mp_grid .== 1   # Majorana pixels
    skel = thinning(mask; algo=GuoAlgo())

    println("type of mask: ", typeof(mask))
    
    # println(mask)

    println("Number of skeleton points: ", count(skel))

    coords = Tuple{Int,Int}[]
    for I in CartesianIndices(skel)
        if skel[I]
            push!(coords, (I[1], I[2]))
        end
    end

    # --- 3. local slope estimation ---
    slopes = Float64[]
    intercepts = Float64[]
    for (i,j) in coords
        # local window around (i,j)
        imin, imax = max(1,i-neigh), min(nx,i+neigh)
        jmin, jmax = max(1,j-neigh), min(ny,j+neigh)
        pts = [(x,y) for x in imin:imax, y in jmin:jmax if skel[x,y]]
        if length(pts) >= 3
            X = [p[1] for p in pts]
            Y = [p[2] for p in pts]
            m,c = linreg(X,Y)
            push!(slopes, m)
            push!(intercepts, c)
        end
    end

    println("Number of local fits: ", length(slopes))

    # --- 4. cluster slopes ---
    features = hcat(slopes, intercepts)
    X = Array(features')
    # if nclusters === nothing
    #     # DBSCAN auto-clustering in slope/intercept space
    #     labels = dbscan(features, 0.2, minpts=20).assignments
    # else
    #     res = kmeans(features', nclusters)
    #     labels = res.assignments
    # end
    res = kmeans(X, nclusters)
    labels = res.assignments

    bands = Dict()
    for (lab, m) in zip(labels, slopes)
        if lab == 0   # DBSCAN noise
            continue
        end
        if !haskey(bands, lab)
            bands[lab] = Float64[]
        end
        push!(bands[lab], m)
    end

    # band_slopes = Dict(lab => mean(ms) for (lab, ms) in bands)
    # return band_slopes


    bands_info = Dict{Int,BandInfo}()

    for band_id in unique(labels)

        # select indices belonging to this band
        idx = findall(labels .== band_id)

        # extract slopes, intercepts, and coordinates
        band_slopes = slopes[idx]
        band_intercepts = intercepts[idx]
        band_x = [coords[i][1] for i in idx]
        band_y = [coords[i][2] for i in idx]

        # create BandInfo
        bands_info[band_id] = BandInfo(
            band_id,
            mean(band_slopes),
            std(band_slopes),
            mean(band_intercepts),
            length(idx),
            band_x,
            band_y
        )
    end

    return mp_grid, bands_info
end



# --- simple linear regression helper ---
function linreg(X, Y)
    xm, ym = mean(X), mean(Y)
    sxx = sum((X .- xm).^2)
    if sxx == 0
        return (0.0, ym)
    end
    m = sum((X .- xm) .* (Y .- ym)) / sxx
    c = ym - m * xm
    return (m,c)
end




function plot_mp_with_bandss(mp_grid::AbstractMatrix, bands_info::Dict{Int64, BandInfo};
                            skeleton_color=:red,
                            line_colors=:auto,
                            line_width=2,
                            skeleton_size=2,
                            alpha_skel=0.6)
    """
    Plot Majorana phase map and overlay detected band skeletons and mean slope lines.
    
    Arguments:
    - mp_grid: 2D array of mp_disc values (0 = non-Majorana, -1 = Majorana)
    - bands_info: Dict of BandInfo structs for each detected band
    Keyword Arguments:
    - skeleton_color: color of skeleton points
    - line_colors: array of colors for each band line or :auto for automatic
    - line_width: width of the slope line
    - skeleton_size: marker size of skeleton points
    - alpha_skel: transparency of skeleton points
    """
    
    # Start with heatmap
    plt = heatmap(mp_grid';
                  c=:viridis,
                  xlabel="mu_t",
                  ylabel="rho",
                  title="Majorana Phase Map with Band Lines",
                  xlims=(0,100),
                  ylims=(0,100),
                  legend=false)
    
    # Overlay skeletons
    for (band_id, info) in bands_info
        scatter!(plt, info.x_coords, info.y_coords;
                 color=skeleton_color,
                 markersize=skeleton_size,
                 alpha=alpha_skel,
                 label="")
    end
    
    # Overlay mean slope lines
    band_ids = collect(keys(bands_info))
    n_bands = length(band_ids)
    
    # assign line colors
    if line_colors == :auto
        palette = distinguishable_colors(n_bands)
    else
        palette = line_colors
    end
    
    for (i, band_id) in enumerate(band_ids)
        info = bands_info[band_id]
        x0 = mean(info.x_coords)
        y0 = mean(info.y_coords)
        x_line = range(minimum(info.x_coords), stop=maximum(info.x_coords), length=100)
        y_line = y0 .+ info.slope_mean .* (x_line .- x0)
        
        plot!(plt, x_line, y_line;
              lw=line_width,
              color=palette[i],
              label="Band $band_id")
    end
    
    Plots.savefig(plt, "mp_with_bands.png")
end



              # slopes = extract_band_slopes(df)
# println("Detected band slopes (dy/dx):")
# for (lab, m) in slopes
#     println(" Band $lab: slope = $m")
# end


mp_grid, bands_info = extract_band_slopes(df; grid_res=(100,100), neigh=2, nclusters=3)

for (band_id, info) in bands_info
    println("Band $band_id:")
    println("  Mean slope: ", info.slope_mean)
    println("  Slope std: ", info.slope_std)
    println("  Mean intercept: ", info.intercept_mean)
    println("  Number of points: ", info.npoints)
    println("  x-range: ", minimum(info.x_coords), " -> ", maximum(info.x_coords))
    println("  y-range: ", minimum(info.y_coords), " -> ", maximum(info.y_coords))
    println()
end


println(typeof(bands_info))
plot_mp_with_bandss(mp_grid, bands_info)