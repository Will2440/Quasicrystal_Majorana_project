# extracts mp_tol ranges from .csv and plots them against slope data
using CSV
using DataFrames
using Plots


csv_filepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)/mp_tol_sweeps/mptarg-1.0_mptolrange1.0e-10-1.0-1001_N_rangenothing_Delta_rangenothing_tn_rangenothing_phi_rangenothing_phason_rangenothing/found_mp_tol_ranges.csv"
df = CSV.read(csv_filepath, DataFrame)

# plots log mp_tol ranges vs slope (phi)
function plot_mp_tol_vs_slope(csv_filepath)
    df = CSV.read(csv_filepath, DataFrame)
    
    # Parse mp_tol_range strings and extract min and max
    mins = Float64[]
    maxs = Float64[]
    for s in df.mp_tol_range
        parsed = eval(Meta.parse(s))
        if parsed isa Tuple
            push!(mins, parsed[1])
            push!(maxs, parsed[2])
        elseif parsed isa Vector && !isempty(parsed)
            push!(mins, parsed[1][1])
            push!(maxs, parsed[1][2])
        else
            push!(mins, NaN)
            push!(maxs, NaN)
        end
    end
    df.mp_tol_min = mins
    df.mp_tol_max = maxs
    
    # Assuming slope is the phi column
    p = plot()
    for row in eachrow(df)
        if !isnan(row.mp_tol_min)
            plot!([log10(row.mp_tol_min), log10(row.mp_tol_max)], [row.phi, row.phi], linewidth=2, color=:blue, label="")
        end
    end
    
    xlabel!("log(mp_tol)")
    ylabel!("slope (phi)")
    title!("mp_tol ranges vs slope")
    return p
end

p = plot_mp_tol_vs_slope(csv_filepath)
savefig(p, "mp_tol_vs_slope.png")