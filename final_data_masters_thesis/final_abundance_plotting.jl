using CSV
using DataFrames
using Glob
using Plots
using LaTeXStrings
using Measures
using Unzip

function combine_integrals(folder_path::String)
    file_paths = glob("*integrals.csv", folder_path)
    combined_rows = DataFrame(Delta=Float64[], mp_integral=Float64[], ipr_integral=Float64[], mbs_integral=Float64[])

    for file_path in file_paths
        filename = split(file_path, "/") |> last
        delta_match = match(r"(\d+p\d+)", filename)
        if delta_match === nothing
            println("Could not parse Delta from filename: $filename")
            continue
        end
        delta_str = replace(delta_match.match, "p" => ".")
        delta = parse(Float64, delta_str)

        df = CSV.read(file_path, DataFrame; delim=",", header=true)
        row = Dict(row.Integral_Type => row.Value for row in eachrow(df))

        push!(combined_rows, (
            Delta = delta,
            mp_integral = abs(row["mp_integral"]),
            ipr_integral = row["ipr_integral"],
            mbs_integral = row["mbs_integral"]
        ))
    end

    sort!(combined_rows, :Delta)

    return combined_rows
end



function plot_mp_integral_vs_delta(df::DataFrame)
    Plots.plot(
        df.Delta,
        df.mp_integral,
        label = L"\textnormal{GQC}",
        lw = 2,
        marker = :o,
        xlabel = L"\Delta",
        ylabel = L"\textnormal{Absolute \ Abundance}",
        legend = :topright,
        grid = false,
        size = (800, 600),
        labelfontsize = 20,
        tickfontsize = 12,
        legendfontsize = 20,
        left_margin = 5mm,
        bottom_margin = 5mm
        # color = :blue  # optional
    )
end



function plot_multiple_mp_integrals(dfs::Vector{DataFrame}, labels::Vector{LaTeXString};
                                    markers::Vector{Symbol} = repeat([:o], length(dfs)),
                                    colors::Vector{Symbol} = nothing)

    plt = Plots.plot(
        xlabel = L"\Delta / t_A",
        ylabel = L"\textnormal{Abundance}",
        legend = false, #:outerbottomleft,
        grid = false,
        size = (800, 600),
        labelfontsize = 20,
        tickfontsize = 12,
        legendfontsize = 20,
        left_margin = 5mm,
        bottom_margin = 5mm,
        # xlims = (0.0, 1.0)
    )

    for (i, df) in enumerate(dfs)
        filtered = [(Δ, mp) for (Δ, mp) in zip(df.Delta, df.mp_integral) if Δ <= 1.0]
        delta_filtered, mp_filtered = unzip(filtered)

        Plots.plot!(
            plt,
            delta_filtered, #df.Delta,
            mp_filtered, #df.mp_integral,
            # df.mbs_integral,
            label = labels[i],
            lw = 2,
            marker = markers[i],
            # color = colors === nothing ? nothing : colors[i]
        )
    end

    # display(plt)
    Plots.savefig("/Users/Will/Documents/FINAL_PROJECT/simulations/report_plots/abundance/final_abundance_plot_large_window_reduced.png")
end




# # Call the function
# GQC_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/GQC_N50_mu0.0-5.0-101_rho0.0-5.0-101_Delta0.1,0.2,0.3,0.4,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,15.0,20.0/"
# SQC_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/SQC_N50_mu0.0-5.0-101_rho0.0-5.0-101_Delta[0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]/"
# TMQC_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/TMQC_N50_mu0.0-5.0-101_rho0.0-5.0-101_Delta[0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0]/"
# GQC_df = combine_integrals(GQC_folder_path)
# SQC_df = combine_integrals(SQC_folder_path)
# TMQC_df = combine_integrals(TMQC_folder_path)
# # plot_mp_integral_vs_delta(combined_df)

# # plot_multiple_mp_integrals(
# #     [gqc_df, sqc_df, tmqc_df, pqc2_df, pqc3_df],
# #     [L"\textnormal{GQC}", L"\textnormal{SQC}", L"\textnormal{TMQC}", L"\textnormal{PQC} \ \sigma=2.0", L"\textnormal{PQC} \ \sigma=3.0"],
# #     markers = [:o, :s, :d, :x, :x],
# #     colors = [:blue, :red, :green, :orange, :purple]
# # )

# plot_multiple_mp_integrals([GQC_df, SQC_df, TMQC_df], 
#     [L"\textnormal{GQC}", L"\textnormal{SQC}", L"\textnormal{TMQC}"], 
#     markers = [:o, :s, :d], 
#     colors = [:blue, :red, :green]
# )



NCa_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/NC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
NCb_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/NC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]"
GQCa_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/GQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
GQCb_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/GQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]"
SQC_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/SQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
TMQC_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/TMQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
# PQC_1p5_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig1.5_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
PQC_2a_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
PQC_2b_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]"
PQC_3_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig3.0_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
PQC_4_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/PQC_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig4.0_Delta[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]/"
TRB_2_folder_path = "/Users/Will/Documents/FINAL_PROJECT/final_data_sharable/final_abundance_results/TRB_N50_mu0.0-10.0-201_rho0.0-10.0-201_sig2.0_Delta[2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]/"

# NC_df = vcat(combine_integrals(NCa_folder_path), combine_integrals(NCb_folder_path))
NC_df = combine_integrals(NCa_folder_path)
# GQC_df = vcat(combine_integrals(GQCa_folder_path), combine_integrals(GQCb_folder_path))
GQC_df = combine_integrals(GQCa_folder_path)
SQC_df = combine_integrals(SQC_folder_path)
TMQC_df = combine_integrals(TMQC_folder_path)
# PQC_1p5_df = combine_integrals(PQC_1p5_folder_path)
# PQC_2_df = vcat(combine_integrals(PQC_2a_folder_path), combine_integrals(PQC_2b_folder_path))
PQC_2_df = combine_integrals(PQC_2a_folder_path)
PQC_3_df = combine_integrals(PQC_3_folder_path)
PQC_4_df = combine_integrals(PQC_4_folder_path)
TRB_2_df = combine_integrals(TRB_2_folder_path)

plot_multiple_mp_integrals([GQC_df, SQC_df, TMQC_df, NC_df, PQC_3_df, PQC_4_df], 
    [L"\textnormal{GQC}", L"\textnormal{SQC}", L"\textnormal{TMQC}", L"\textnormal{NC}", L"\textnormal{PQC} \ \sigma=3.0", L"\textnormal{PQC} \ \sigma=4.0"], 
    markers = [:o, :s, :d, :o, :x, :x], 
    colors = [:blue, :red, :green]
)



println("finished")