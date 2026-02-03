using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter

include("bson_unpacker.jl")

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(110-110-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-301)_Delta(0.05-0.05-1)"
run_set_name = "sturm_range_full"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "raw_data", "np", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "results", "np", run_set_name, data_folder))
isdir(base_results_dir) || mkpath(base_results_dir)


######################################
############### Unpack? ##############
######################################

# data = process_bson_files_with_winding(folder_path_hof)
data = unpack_from_inspect(folder_path_hof)

println("Type of unpacked data: ", typeof(data))
println("Unpacked data keys: ", names(data))


function process_and_aggregate_winding(df::DataFrame)
    println("Processing raw winding data...")

    # 1. Parse t_n to get slope reliably if needed, or just use t_n directly if it's consistent.
    # We assume 't_n' is a vector [t, t_modulated]. 
    # Grouping by 't_n' directly works if vectors are identical.
    
    # 2. Extract phason from sequence_id for sorting
    # sequence_id is Tuple{Float64, Float64} -> (slope, phason)
    df.phason = [x[2] for x in df.sequence_id]
    
    # 3. Group by Physical Parameters
    # We group by everything EXCEPT the phason and the sequence specific data
    # Note: t_n is a Vector, which is groupable in recent DataFrames versions.
    grouped = groupby(df, [:N, :t_n, :mu, :Delta])
    
    # 4. Aggregate to find the final winding per group
    # Since 'winding' (or 'winding_phase') in the solver output was the *accumulated* phase
    # at that step, the value at the largest phason is the total winding for the cycle.
    
    processed_df = combine(grouped) do sdf
        # Sort by phason to ensure we follow the cycle 0 -> 2pi
        sort!(sdf, :phason)
        
        # The winding recorded in the rows is the running accumulator.
        # So we just take the LAST value.
        # Check column name: is it :winding_phase or :winding?
        # Based on solver: it was :winding. Based on your unpacker: it might be :winding_phase.
        # Let's assume :winding_phase based on your previous code.
        
        # If winding_phase is a scalar (Float64) per row:
        final_phase = last(sdf.winding_phase)
        
        # If winding_phase is a Vector per row (history):
        # final_phase = last(last(sdf.winding_phase)) 
        # But usually solver pushes the current scalar sum. Let's assume scalar per row.
        # If it's a vector per row, we need to know if that vector is the WHOLE history or just that step.
        # Assuming scalar per row based on "each row only has one phase winding value".
        
        return (
            final_winding_phase = final_phase,
            winding_number = final_phase / (2π),
            max_phason = last(sdf.phason),
            count = nrow(sdf)
        )
    end
    
    println("Aggregation complete. Reduced $(nrow(df)) raw rows to $(nrow(processed_df)) unique parameter sets.")
    return processed_df
end

aggregated_data = process_and_aggregate_winding(data)


# function add_winding_columns!(df::DataFrame)
#     # Extract the final accumulated winding phase from the list in :winding_phase
#     # :winding_phase is expected to be a Vector{Float64}
    
#     # helper to safely get last element
#     get_last(x) = isempty(x) ? NaN : last(x)

#     df.final_winding_phase = [get_last(wp) for wp in df.winding_phase]
    
#     # Compute the winding number as final_winding_phase / 2π
#     df.winding_number = df.final_winding_phase ./ (2π)
# end
# add_winding_columns!(data)
# println("Added winding columns to DataFrame.")

######################################
############## Plotting ##############
######################################

function plt_winding_vs_mu(
    df::DataFrame,
    save_path::String;
    fixed_values...
)
    # Filter data for the specified parameters
    df_filtered = filter(row -> all(getfield(row, k) == v for (k, v) in fixed_values), df)
    
    # Check if we have duplicates for mu
    # If we do, we might need to pick one (e.g. the first one) or average them
    # For now, let's just sort and plot. If there are multiples, it usually means 
    # multiple sequences or parameters were mixed in.
    
    # OPTIONAL: Deduplicate by averaging or picking unique
    # df_unique = combine(groupby(df_filtered, :mu), :winding_number => first => :winding_number)
    # sort!(df_unique, :mu)
    # data_to_plot = df_unique

    sort!(df_filtered, :mu)

    # improved scatter plot
    plt = scatter(
        df_filtered.mu, 
        df_filtered.winding_number, 
        xlabel="Chemical Potential (μ)", 
        ylabel="Winding Number", 
        title="Winding Number vs Chemical Potential", 
        legend=false,
        markersize=1.5,   # Increased size slightly for visibility
        markerstrokewidth=0,
        color=:blue,
        grid=true
    )
    
    # Save the plot
    png(plt, save_path)
    println("Saved plot to: $save_path")
end

savepath = joinpath(base_results_dir, "winding_vs_mu.png")
plt_winding_vs_mu(
    aggregated_data, 
    savepath;
)

# function plt_winding_accumulator_vs_mu(
#     df::DataFrame,
#     save_path::String;
#     fixed_values...
# )
#     # Filter data for the specified fixed values (e.g., N, Delta, etc.)
#     df_filtered = filter(row -> all(getfield(row, k) == v for (k, v) in fixed_values), df)
    
#     # Sort by mu for better plotting
#     sort!(df_filtered, :mu)
    
#     # Prepare data for plotting: each element of winding_phase as a point
#     x_vals = Float64[]
#     y_vals = Float64[]
#     color_vals = Int[]
    
#     for row in eachrow(df_filtered)
#         mu_val = row.mu
#         for (idx, phase) in enumerate(row.winding_phase)
#             push!(x_vals, mu_val)
#             push!(y_vals, phase)
#             push!(color_vals, idx)
#         end
#     end
    
#     # Create the plot
#     plt = scatter(
#         x_vals, 
#         y_vals, 
#         xlabel = "Chemical Potential (μ)", 
#         ylabel = "Winding Phase", 
#         title = "Winding Phase Elements vs Chemical Potential", 
#         legend = false,
#         markersize = 0.5,
#         color = color_vals,  # Color by index in the winding_phase vector
#         colormap = :viridis,  # Use a colormap for the indices
#         grid = true
#     )
    
#     # Save the plot
#     png(plt, save_path)
#     println("Saved plot to: $save_path")
# end

# # Example call (adjust fixed_values as needed, e.g., N=100, Delta=0.05)
# savepath_accum = joinpath(base_results_dir, "winding_accumulator_vs_mu.png")
# plt_winding_accumulator_vs_mu(
#     data, 
#     savepath_accum
# )
