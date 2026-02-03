using BSON
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter

# Set logger level from JULIA_LOG_LEVEL environment variable, default to Info if not set
level_str = uppercase(get(ENV, "JULIA_LOG_LEVEL", "info"))
level = if level_str == "DEBUG"
    Logging.Debug
elseif level_str == "INFO"
    Logging.Info
elseif level_str == "WARN"
    Logging.Warn
elseif level_str == "ERROR"
    Logging.Error
else
    Logging.Info  # fallback
end
global_logger(ConsoleLogger(stderr, level))

include("plotting.jl")
include("processing.jl")
include("RAM_opt_unpacker.jl")

using .IDOPUnpacking
using .IDOPPlotting
using .IDOPProcessing

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "hof_style_slopes_N400_target_0.61803_tol_0.001_phason_0.0-101-1.0_nbins1000_npb1_N(200-200-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1201)_Delta(0p05-0p05-1)"
run_set_name = "QC-SC_comp_small_L_check"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", "bp_results", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "bp_results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define the unpacking mode: :unpack_by_phason, :unpack_by_slope, :load_phason, or :load_slope
# - :unpack_by_phason: Unpack and slice data with slices per phason (varying phi within slices)
# - :unpack_by_slope: Unpack and slice data with slices per phi (slope), containing range of phason
# - :load_phason / :load_slope: Load existing combo_to_files.bson without unpacking
# -------- IMPORTANT --------
unpack_mode = :load_phason
partitioned_unpack = false  # Use partitioned unpacking to reduce RAM usage
n_partitions = 2
# unify = true
# unify_folder = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/bp_results/phason_final_data/EXTRAS_hof_style_slopes_N1000_target_0.61803_tol_0.001_phason_0.0-501-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(1p5-1p5-1)_mu(0p0-3p0-1311)_Delta(0p05-0p05-1)"
# ---------------------------

eig_save_mode = :all        # Options: :all (full spectrum), :maj (central 2 only)
mirror_data_about_mu0 = false  # Whether to mirror data about mu=0 for unpacking

# Define ranges for repacking or plotting -- otherwise the full range will be used.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phason_range = nothing
repack_phi_range = nothing #[0.6180048661800487]

# MP discretisation parameters
mp_targ = -1.0
mp_tol = 1e-5

# IDOP plateau parameters
delta_idop_thresh = 1e-4
min_mu_span = 0.0
min_samples = 2

# Gap labelling parameters
n = 10
gap_label_p_range = collect(-n:n)
gap_label_q_max = n

# Determine base_slice_dir based on unpack_mode for consistent file locations
slice_suffix = if unpack_mode in [:unpack_by_phason, :load_phason]
    "_phason"
elseif unpack_mode in [:unpack_by_slope, :load_slope]
    "_slope"
else
    error("Invalid unpack_mode: $unpack_mode. Must be :unpack_by_phason, :unpack_by_slope, :load_phason, or :load_slope.")
end
base_slice_dir = normpath(joinpath(base_results_dir, "idop_base_slices", "mptarg$(mp_targ)_mptol$(mp_tol)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)$(slice_suffix)_muduplicated-$(mirror_data_about_mu0)"))

if unpack_mode == :unpack_by_phason
    if partitioned_unpack
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_phason_partitioned(folder_path_hof, base_slice_dir;
            n_partitions=n_partitions,
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range,
            eig_save_mode=eig_save_mode,
            mirror_data_about_mu0=mirror_data_about_mu0
        )
    else
        # Run unpacking and slicing per phason (fixed phi, varying phason)
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_phason(folder_path_hof, base_slice_dir;
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range,
            eig_save_mode=eig_save_mode,
            mirror_data_about_mu0=mirror_data_about_mu0
        )
    end
elseif unpack_mode == :unpack_by_slope
    # Run unpacking and slicing per phi (fixed phason, varying phi)
    if partitioned_unpack
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope_partitioned(folder_path_hof, base_slice_dir;
            n_partitions=n_partitions,
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phason_range_arg=repack_phason_range,
        mirror_data_about_mu0=mirror_data_about_mu0
    )
    else
        combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope(folder_path_hof, base_slice_dir;
            mp_targ=mp_targ,
            mp_tol=mp_tol,
            plt_N_range_arg=repack_N_range,
            plt_Delta_range_arg=repack_Delta_range,
            plt_tn_range_arg=repack_tn_range,
            plt_phason_range_arg=repack_phason_range,
            mirror_data_about_mu0=mirror_data_about_mu0
        )
    end
elseif unpack_mode in [:load_phason, :load_slope]
    # Load existing slices
    combo_data = BSON.load(joinpath(base_slice_dir, "combo_to_files.bson"))
    combo_to_files = combo_data[:combo_to_files]
    println("keys of comb_to_files.bson: ", keys(combo_data))
else
    error("Unexpected unpack_mode: $unpack_mode")
end


######################################
########### Extra variables ##########
######################################

## Defines outdir for saving plots in this script
outdir = normpath(joinpath(base_results_dir, "idop", "mptarg$(mp_targ)_mptol$(mp_tol)_IDOPminsamp$(min_samples)_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)_mpduplicated-$(mirror_data_about_mu0)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)

# Collect per-slice gap labelling metadata for export
label_records = Vector{Dict{Symbol, Any}}()

######################################
########## Plotting Execution ########
######################################

## Colour Colorscheme
function q_anchor_colormap(
    qs::Vector{Int}
)

    lerp(c1::RGB, c2::RGB, t) = RGB(
        (1-t)*c1.r + t*c2.r,
        (1-t)*c1.g + t*c2.g,
        (1-t)*c1.b + t*c2.b
    )
    qs_sorted = sort(qs)
    
    # --- SENTINEL HANDLING ---
    SC_GAP_Q = 9999
    has_sentinel = SC_GAP_Q in qs_sorted
    
    # Filter out sentinel for gradient calculation
    qs_phys = filter(!=(SC_GAP_Q), qs_sorted)
    # -------------------------

    cmap = fill(RGB(0.7,0.7,0.7), length(qs_phys))

    idx_neg_high = 1
    idx_neg1     = findfirst(==( -1), qs_phys)
    idx_pos1     = findfirst(==(  1), qs_phys)
    idx_pos_high = length(qs_phys)

    col_neg_light = RGB(0.35, 0.75, 1.00)   # q = -1
    col_neg_dark  = RGB(0.05, 0.10, 0.45)   # mid dark blue
    col_neg_purp  = RGB(0.25, 0.00, 0.50)   # dark purple
    col_neg_high  = RGB(0.78, 0.20, 1.00)   # bright purple (most negative)

    col_pos_light = RGB(1.00, 0.60, 0.10)   # light orange
    col_pos_yell  = RGB(1.00, 0.95, 0.40)   # light yellow
    col_pos_high  = RGB(0.00, 1.00, 0.00)   # green (most positive)

    cmap[idx_neg_high] = col_neg_high
    
    if !isnothing(idx_neg1); cmap[idx_neg1] = col_neg_light; end
    if !isnothing(idx_pos1); cmap[idx_pos1] = RGB(1.0,0.0,0.0); end

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

    if !isnothing(idx_pos1) && idx_pos_high > idx_pos1
        span = idx_pos_high - idx_pos1
        for (j,i) in enumerate(idx_pos1+1:idx_pos_high)
            t = j / span
            if t < 0.33
                cmap[i] = lerp(RGB(1.0,0.0,0.0), col_pos_light, t/0.33)
            elseif t < 0.66
                cmap[i] = lerp(col_pos_light, col_pos_yell, (t-0.33)/0.33)
            else
                cmap[i] = lerp(col_pos_yell, col_pos_high, (t-0.66)/0.34)
            end
        end
    end

    if !isnothing(idx_neg1) && !isnothing(idx_pos1) && idx_pos1 - idx_neg1 > 1
        span = idx_pos1 - idx_neg1
        for (j,i) in enumerate(idx_neg1+1:idx_pos1-1)
            t = j / span
            cmap[i] = lerp(col_neg_light, RGB(1.0,0.0,0.0), t)
        end
    end

    # --- RECONSTRUCT MAP WITH SENTINEL ---
    # Map q -> color
    color_dict = Dict(q => c for (q,c) in zip(qs_phys, cmap))
    
    if has_sentinel
        color_dict[SC_GAP_Q] = RGB(0.85, 0.85, 0.85) # Light Grey for SC gaps
    end
    
    # Return colors in the order of the original sorted input
    return [color_dict[q] for q in qs_sorted]
end

function generate_q_color_dict(q_range)
    # Ensure range includes 0 and typical winding numbers
    qs = sort(unique(collect(q_range)))
    colors = q_anchor_colormap(qs)
    return Dict(q => c for (q,c) in zip(qs, colors))
end

# Determine fixed param based on mode for consistent handling
if unpack_mode in [:unpack_by_phason, :load_phason]
    fixed_param = :phi
    fixed_param_name = "phi"
    varying_param = :phason
    varying_param_name = "phason"
    x_axis_param = varying_param  # For plotting: x = phason, y = mu_mp
else  # :unpack_by_slope or :load_slope
    fixed_param = :phason
    fixed_param_name = "phason"
    varying_param = :phi
    varying_param_name = "phi"
    x_axis_param = fixed_param  # For plotting: x = phason (fixed per slice), but accumulate across slices
end


@showprogress for ((N_val, Delta_val, t_n_tuple, fixed_val), f) in combo_to_files

    slice_suffix_plotname = "$(fixed_param_name)$(fixed_val)"

    # # One-off filter for specific phi
    # if fixed_param == :phi && abs(fixed_val - 0.6180048661800487) > 1e-1
    #     continue
    # end

    # Apply filtering based on fixed_param
    if !isnothing(repack_N_range) && !(N_val in repack_N_range)
        continue
    end
    if !isnothing(repack_Delta_range) && !(Delta_val in repack_Delta_range)
        continue
    end
    if !isnothing(repack_tn_range) && !(t_n_tuple in repack_tn_range)
        continue
    end
    # Filter based on fixed_param
    range_to_check = if fixed_param == :phason
        repack_phason_range
    elseif fixed_param == :phi
        repack_phi_range
    else
        nothing
    end
    if !isnothing(range_to_check) && !(fixed_val in range_to_check)
        continue
    end

    data = BSON.load(f)
    df_slice = data[:df_slice]
    t_n = collect(t_n_tuple)
    mu_values = Float64[]
    phason_vals = Float64[]
    phason_range = missing
    target_phason_actual = missing
    idop_labels_formatted = Vector{Tuple{Float64, Float64, Int}}()
    phi_vals = Float64[]
    phi_value_for_record = fixed_param == :phi ? fixed_val : missing
    phi_range = missing

    # --- FIX: Ensure fixed parameter column exists in DataFrame ---
    # When loading slices by phason, 'phi' is fixed and might not be in the DF columns.
    # When loading slices by slope, 'phason' is fixed and might not be in the DF columns.
    if !(fixed_param in propertynames(df_slice))
        df_slice[!, fixed_param] .= fixed_val
    end

    # --- DEBUG: summarise slice contents ---
    @info "Slice summary" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val n_rows=nrow(df_slice)
    if nrow(df_slice) > 0
        # unique values actually present in this slice
        varying_present = sort(unique(df_slice[!, varying_param]))
        mu_mp_lengths = length.(df_slice.mu_mp)
        min_mu_mp_len = minimum(mu_mp_lengths)
        max_mu_mp_len = maximum(mu_mp_lengths)
        total_mu_points = sum(mu_mp_lengths)

        if :mu_range in propertynames(df_slice)
            # Unpacker provides 'mu_range'
            mu_values = df_slice.mu_range[1]

            @info "mu_values obtained (length, min, max): ", length(mu_values), ", ", minimum(mu_values), ", ", maximum(mu_values)
            @info "min_pos_eigs length: ", length(df_slice.min_pos_eigs[1])
            @info "mu_mp length: ", length(df_slice.mu_mp[1])

            lens_pos = :min_pos_eigs in propertynames(df_slice) ? length.(df_slice.min_pos_eigs) : Int[]
            lens_neg = :min_neg_eigs in propertynames(df_slice) ? length.(df_slice.min_neg_eigs) : Int[]
            
            # Combine lengths and find the safe minimum
            all_lens = vcat(lens_pos, lens_neg)
            if :mu_mp in propertynames(df_slice); append!(all_lens, length.(df_slice.mu_mp)); end
                        
            # Determine the safe length to use
            safe_len = isempty(all_lens) ? length(mu_values) : minimum(all_lens)
            @info "Safe length determined from data: " safe_len

            # Ensure 'mu_values' column exists for compatibility
            if !(:mu_values in propertynames(df_slice))
                df_slice[!, :mu_values] = df_slice.mu_range
            end
        else
            # Fallback for legacy data
            n_points = length(df_slice.mu_mp[1])
            mu_values = collect(range(-3.0, 3.0, length=n_points))
            df_slice[!, :mu_values] = [mu_values for _ in 1:nrow(df_slice)]
            @warn "'mu_range' column missing; generated 'mu_values' form hardcoded range(-3.0, 3.0)."
        end

        if :phason in propertynames(df_slice)
            phason_col = df_slice.phason
            valid_idx = findall(!ismissing, phason_col)
            if !isempty(valid_idx)
                phason_vals = Float64[phason_col[i] for i in valid_idx]
                phason_range = (minimum(phason_vals), maximum(phason_vals))
            end
        end

        if :phi in propertynames(df_slice)
            phi_col = df_slice.phi
            phi_valid_idx = findall(!ismissing, phi_col)
            if !isempty(phi_valid_idx)
                phi_vals = Float64[phi_col[i] for i in phi_valid_idx]
                if ismissing(phi_value_for_record)
                    phi_value_for_record = first(phi_vals)
                end
                phi_range = (minimum(phi_vals), maximum(phi_vals))
            end
        end

        @info "Slice value ranges" $varying_param_name=varying_present min_mu_mp_len=min_mu_mp_len max_mu_mp_len=max_mu_mp_len total_mu_points=total_mu_points
    else
        @warn "df_slice is empty" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val
    end
    # --- END DEBUG ---



    # ## 1a) IDOP and gap labelling processing (used with unpack by SLOPE)
    # df_slice = IDOPProcessing.compute_idop_df!(df_slice; disc_variable=:mu_mp, rescale_idop=false)
    
    # # copy mu_mp column but with missing values removed
    # df_slice[!, :mu_mp_nomiss] = [filter(!ismissing, mu_mp) for mu_mp in df_slice.mu_mp]

    # ## currently manually add mu_values column
    # # mu_values = collect(range(-3.0, 3.0, 1201))
    # # df_slice[!, :mu_values] = [mu_values for _ in 1:nrow(df_slice)]
    
    # ## Dynamically create mu_values based on the actual data length
    # if nrow(df_slice) > 0
    #     n_points = length(df_slice.mu_mp[1]) # Get length from the first row
    #     mu_values = collect(range(-3.0, 3.0, length=n_points))
    #     df_slice[!, :mu_values] = [mu_values for _ in 1:nrow(df_slice)]
    # end
    # # println("Added mu_values column with length ", length(df_slice.mu_values[1]), " based on data.")

    # df_slice = IDOPProcessing.compute_idop_plateaus_all!(
    #     df_slice; 
    #     delta_idop_thresh=delta_idop_thresh, 
    #     min_mu_span=min_mu_span, 
    #     min_samples=min_samples
    # )

    # IDOPProcessing.filter_idop_plateaus_at_1!(df_slice; atol=1e-6)
    # IDOPProcessing.filter_idop_plateaus_at_0!(df_slice; atol=1e-6)

    # df_slice = IDOPProcessing.compute_gap_labels_qlim_new!(
    #     df_slice;
    #     p_range=gap_label_p_range,
    #     q_max=gap_label_q_max
    # )



    # ## 1b) IDOP and gap labelling processing
    # df_slice = IDOPProcessing.compute_idop_df!(df_slice; disc_variable=:mu_mp, rescale_idop=false)

    # df_slice = IDOPProcessing.compute_idop_plateaus_all!(
    #     df_slice; 
    #     delta_idop_thresh=delta_idop_thresh, 
    #     min_mu_span=min_mu_span, 
    #     min_samples=min_samples
    # )

    # IDOPProcessing.filter_idop_plateaus_at_1!(df_slice; atol=1e-6)
    # IDOPProcessing.filter_idop_plateaus_at_0!(df_slice; atol=1e-6)

    # # df_slice = IDOPProcessing.compute_gap_labels_qlim_new!(
    # #     df_slice;
    # #     p_range=gap_label_p_range,
    # #     q_max=gap_label_q_max
    # # )

    # df_slice = IDOPProcessing.compute_gap_labels_parsimonious!(
    #     df_slice;
    #     p_range=gap_label_p_range,
    #     q_max=gap_label_q_max,
    #     tol=1e-6,                 # Tolerance for parsimony
    #     merge_threshold=nothing   # Set to e.g. 1e-4 if you want to merge close gaps
    # )
    
    #############################################
    ##### Choose IDOP labels to record ########
    #############################################

    # --- Extract IDOP Labels for Phason ≈ 0.0 ---
    # Find row closest to phason = 0.0
    # target_phason = 0.0
    # if :phason in propertynames(df_slice)
    #     phason_col = df_slice.phason
    #     valid_idx = findall(!ismissing, phason_col)
    #     if !isempty(valid_idx)
    #         diffs = [abs(phason_col[i] - target_phason) for i in valid_idx]
    #         closest_idx = valid_idx[argmin(diffs)]
    #         target_phason_actual = phason_col[closest_idx]
    #         row_phason0 = df_slice[closest_idx, :]

    #         if !ismissing(row_phason0.gap_labels)
    #             for gl in row_phason0.gap_labels
    #                 if !ismissing(gl.q)
    #                     push!(idop_labels_formatted, (gl.E_low, gl.E_high, Int(gl.q)))
    #                 end
    #             end
    #         end
    #     end
    # end
    # # -------------------------------------------------

    # # --- Collect IDOP Labels for ALL Phasons ---
    # if :phason in propertynames(df_slice) && :gap_labels in propertynames(df_slice)
    #     # Iterate over all rows with valid phason and gap_labels
    #     for row in eachrow(df_slice)
    #         if !ismissing(row.phason) && !ismissing(row.gap_labels)
    #             current_phason = row.phason
    #             # Determine phi: if fixed_param is :phi, use fixed_val, else look in column
    #             current_phi = if hasproperty(row, :phi) && !ismissing(row.phi)
    #                 row.phi
    #             elseif fixed_param == :phi
    #                 fixed_val
    #             else
    #                 missing
    #             end

    #             # Format gap labels
    #             current_labels = Vector{Tuple{Float64, Float64, Int}}()
    #             for gl in row.gap_labels
    #                 if !ismissing(gl.q)
    #                     push!(current_labels, (gl.E_low, gl.E_high, Int(gl.q)))
    #                 end
    #             end

    #             record = Dict{Symbol, Any}(
    #                 :N => N_val,
    #                 :Delta => Delta_val,
    #                 :t_n => t_n,
    #                 :phi => current_phi,
    #                 :phason => current_phason,
    #                 :idop_labels => current_labels,
    #                 :run_slice_path => f
    #             )
    #             push!(label_records, record)
    #         end
    #     end
    # end
    # # -------------------------------------------------
    #######################################################################

    # ## 1c) Plot IDOP plateau check (Diagnostic) - Plot only closest to target_phason
    # check_outdir = joinpath(outdir, "plateau_checks")
    # isdir(check_outdir) || mkpath(check_outdir)
    
    # if :plateaus in propertynames(df_slice) && :phason in propertynames(df_slice)
    #     # Find row closest to target_phason
    #     target_phason_check = 0.3
    #     phason_col = df_slice.phason
    #     valid_idx = findall(!ismissing, phason_col)
        
    #     if !isempty(valid_idx)
    #         diffs = [abs(phason_col[i] - target_phason_check) for i in valid_idx]
    #         closest_row_idx = valid_idx[argmin(diffs)]
            
    #         # Slice the DataFrame to keep just that row
    #         df_subset = df_slice[closest_row_idx:closest_row_idx, :]
    #         actual_phason = phason_col[closest_row_idx]

    #         # Extract plateaus for this single row
    #         plateaus_float = Vector{Float64}[Float64[p.plateau for p in row_data] for row_data in df_subset.plateaus]
            
    #         check_filename = joinpath(check_outdir, "phason$(actual_phason)", "plateau_check_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val)_phason$(actual_phason).png")
            
    #         IDOPPlotting.plt_idop_plateaus_check(
    #             df_subset,
    #             plateaus_float,
    #             check_filename
    #         )
    #     end
    # end


    # ## 1c) Plot IDOP plateau check (Diagnostic) - Plot IDOPs for all phasons
    # check_outdir = joinpath(outdir, "plateau_checks")
    
    # if :plateaus in propertynames(df_slice) && :phason in propertynames(df_slice)
    #     # Iterate over all valid phasons in this slice
    #     phason_col = df_slice.phason
    #     valid_idx = findall(!ismissing, phason_col)
        
    #     if !isempty(valid_idx)
    #         # Create a subdirectory for this slice's diagnostics
    #         slice_subdir = joinpath(check_outdir, "N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$(fixed_param_name)$(fixed_val)")
    #         isdir(slice_subdir) || mkpath(slice_subdir)

    #         for i in valid_idx
    #             # Slice the DataFrame to keep just that row
    #             df_subset = df_slice[i:i, :]
    #             actual_phason = phason_col[i]

    #             # Extract plateaus for this single row
    #             plateaus_float = Vector{Float64}[Float64[p.plateau for p in row_data] for row_data in df_subset.plateaus]
                
    #             check_filename = joinpath(slice_subdir, "plateau_check_phason$(actual_phason).png")
                
    #             IDOPPlotting.plt_idop_plateaus_check(
    #                 df_subset,
    #                 plateaus_float,
    #                 check_filename
    #             )
    #         end
    #     end
    # end


    # ## 2.0) Normalise mu_mp to [-1, 1] for plotting
    # df_slice = IDOPProcessing.normalise_scatter!(
    #     df_slice, 
    #     :mu_mp,
    #     :mu_mp_norm; 
    #     norm_range=(-1.0, 1.0)
    # )

    # coeff_filename = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/hof_style_slopes/hof_style_slopes_N500_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-5p0-7)_mu(0p0-0p0-1)_Delta(0p05-0p5-6)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_interp_norm/smooth_norm_coeffs_Delta0.05_t11.0_t22.0_phason0.0_thresh0.001,0.004_prange20_qmax20.txt"
    # df_slice = IDOPProcessing.apply_smooth_normalisation!(
    #     df_slice,
    #     coeff_filename;
    #     target_col=:mu_mp,
    #     out_col=:mu_mp_norm
    # )

    # ## 2a) Plot discrete phase projections
    # try
    #     phase_proj_outdir = joinpath(outdir, "phase_projections")
    #     isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)
        
    #     IDOPPlotting.plt_discrete_phase_projections(
    #         df_slice,
    #         :mu_mp,
    #         :phi,
    #         joinpath(phase_proj_outdir, "phase_proj_NN$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # plt_xlims=(-1.0,1.0),
    #         # plt_ylims=(0.0,1.0)
    #         # N=N_val,
    #         # Delta=Delta_val,
    #         # t_n=t_n,
    #         # phason=phason_val
    #     )

    # catch e
    #     @warn "plt_discrete_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    # ## 2b) Plot normalised discrete phase projections
    # try
    #     phase_proj_outdir = joinpath(outdir, "phase_projections")
    #     isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)
        
    #     IDOPPlotting.plt_discrete_phase_projections(
    #         df_slice,
    #         :mu_mp_norm,
    #         :phi,
    #         joinpath(phase_proj_outdir, "norm_phase_proj_NN$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         plt_xlims=(-1.0,1.0),
    #         plt_ylims=(0.0,1.0)
    #         # N=N_val,
    #         # Delta=Delta_val,
    #         # t_n=t_n,
    #         # phason=phason_val
    #     )

    # catch e
    #     @warn "plt_discrete_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    # ## N/A) Plot mu_mp vs phason for fixed N, Delta, t_n and slope (phi)
    # try
    #     mu_mp_phason_outdir = joinpath(outdir, "mu_mp_vs_phason")
    #     isdir(mu_mp_phason_outdir) || mkpath(mu_mp_phason_outdir)
        
    #     IDOPPlotting.plt_mu_mp_vs_phason(
    #         df_slice,
    #         joinpath(mu_mp_phason_outdir, "mu_mp_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # $fixed_param => fixed_val  # Filter by fixed param (phi for by_phason mode)
    #     )
    # catch e
    #     @warn "plt_mu_mp_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end


    # ## N/A) Plot gap eigs vs phason
    # # Process gap eigenvalues
    # df_slice = IDOPProcessing.process_gap_eigenvalues(df_slice)

    # try
    #     gap_eigs_phason_outdir = joinpath(outdir, "gap_eigs_vs_phason")
    #     isdir(gap_eigs_phason_outdir) || mkpath(gap_eigs_phason_outdir)
        
    #     IDOPPlotting.plt_gap_eigs_vs_phason(
    #         df_slice,
    #         joinpath(gap_eigs_phason_outdir, "gap_eigs_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png")
    #     )
    # catch e
    #     @warn "plt_gap_eigs_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    # # 2) Plot lowest eigs vs phason -- USE unpack_by_phason --
    # # Process lowest eigenvalues
    # try
    #     eigs_phason_outdir = joinpath(outdir, "lowest_eigs_vs_phason")
    #     isdir(eigs_phason_outdir) || mkpath(eigs_phason_outdir)
        
    #     IDOPPlotting.plt_lowest_eigs_vs_phason(
    #         df_slice,
    #         joinpath(eigs_phason_outdir, "lowest_eigs_vs_phason_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # $fixed_param => fixed_val
    #     )
    # catch e
    #     @warn "plt_lowest_eigs_vs_phason failed" N_val=N_val Delta_val=Delta_val t_n=t_n fixed_param=fixed_val error=e
    # end

    # ## 3a) Plot lowest eigs vs phason with winding rectangles (use with unpack by PHASON)
    # try
    #     eigs_phason_outdir = joinpath(outdir, "lowest_eigs_vs_phason_nodal_winding")
    #     isdir(eigs_phason_outdir) || mkpath(eigs_phason_outdir)

    #     # n_points = length(df_slice.min_pos_eigs[1])
    #     # mu_values = collect(range(-3.0, 3.0, length=n_points))

    #     qmax=gap_label_q_max
    #     qrange = -qmax:qmax
    #     q_color_dict = generate_q_color_dict(qrange)
        
    #     gap_windings_pos = IDOPProcessing.calculate_winding_gaps_nodal_lines(df_slice, mu_values, :min_pos_eigs)
    #     gap_windings_neg = IDOPProcessing.calculate_winding_gaps_nodal_lines(df_slice, mu_values, :min_neg_eigs)
        
    #     plot_filename = joinpath(eigs_phason_outdir, "lowest_eigs_nodal_winding_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png")
    #     IDOPPlotting.plt_lowest_eigs_vs_phason_with_winding_rects(
    #         df_slice, 
    #         plot_filename, 
    #         gap_windings_pos, 
    #         gap_windings_neg, 
    #         mu_values,
    #         q_color_dict
    #     )
    # catch e
    #     @warn "plot_eigs_with_winding_rects failed" N_val=N_val Delta_val=Delta_val t_n=t_n fixed_param=fixed_val error=e
    #     showerror(stderr, e, catch_backtrace())
    # end



    ## 3b) Plot lowest eigs vs phason with winding rectangles AND IDOP labels
    gap_windings = Vector{Tuple{Float64, Float64, Int}}()
    try
        eigs_phason_outdir = joinpath(outdir, "lowest_eigs_vs_phason_nodal_winding_with_IDOP")
        isdir(eigs_phason_outdir) || mkpath(eigs_phason_outdir)

        qmax = gap_label_q_max
        qrange = -qmax:qmax
        q_color_dict = generate_q_color_dict(qrange)
        
        # Calculate Winding (Nodal Lines)
        # gap_windings_pos = IDOPProcessing.calculate_winding_gaps_slice_counting(df_slice, mu_values, :min_pos_eigs)
        # gap_windings_neg = IDOPProcessing.calculate_winding_gaps_slice_counting(df_slice, mu_values, :min_neg_eigs)

        log_process = true

        # gap_windings, fit_results = IDOPProcessing.calculate_winding_gaps_sine_fit(
        #     df_slice, 
        #     mu_values, 
        #     :min_pos_eigs, 
        #     :min_neg_eigs;
        #     log_process=log_process,
        #     fit_strategy=:fft,
        # )

        gap_windings_pos = gap_windings
        # gap_windings_neg = gap_windings
        
        plot_filename = joinpath(eigs_phason_outdir, "lowest_eigs_winding_idop_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).svg")
        
        # IDOPPlotting.plt_lowest_eigs_winding_and_idop(
        #     df_slice, 
        #     plot_filename, 
        #     gap_windings_pos, 
        #     gap_windings_neg, 
        #     idop_labels_formatted, # Pass IDOP labels here
        #     mu_values,
        #     q_color_dict;
        #     log_process=log_process,
        #     fit_results=fit_results
        # )

        # IDOPPlotting.plt_lowest_eigs_winding_single(
        #     df_slice,
        #     plot_filename,
        #     gap_windings_pos,
        #     idop_labels_formatted, # Pass IDOP labels
        #     mu_values,
        #     q_color_dict;
        #     log_process=log_process,
        #     fit_results=nothing, #fit_results,
        #     side=:pos
        # )

        IDOPPlotting.plt_simple_lowest_eigs_heatmap(
            df_slice,
            plot_filename,
            mu_values;
            log_process=log_process,
            # threshold=-4.0,
            side=:pos
        )
    catch e
        @warn "plt_lowest_eigs_winding_and_idop failed" N_val=N_val Delta_val=Delta_val t_n=t_n fixed_param=fixed_val error=e
        showerror(stderr, e, catch_backtrace())
    end


    ## 4) plot raw mp vs phason heatmap
    # powers=collect(1:1:10)
    # mp_tol_range = 10.0 .^ (-powers)
    mp_tol_range=[2e-1]
    for mp_tol in mp_tol_range 
        try
            heatmap_outdir = joinpath(outdir, "mu_mp_heatmaps")
            isdir(heatmap_outdir) || mkpath(heatmap_outdir)
            # this function will plot mu_mp vs phason as a heatmap
            IDOPPlotting.plt_simple_mu_mp_heatmap(
                df_slice,
                joinpath(heatmap_outdir, "mu_mp_heatmap_mptol$(mp_tol)_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).svg"),
                mu_values,
                mp_tol=mp_tol,
                size=(1000,200)
            )
        catch e
            @warn "plt_simple_mu_mp_heatmap failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
        end
    end


    # if ismissing(target_phason_actual)
    #     target_phason_actual = target_phason
    # end

    # record = Dict{Symbol, Any}(
    #     :N => N_val,
    #     :Delta => Delta_val,
    #     :t_n => t_n,
    #     :phi => phi_value_for_record,
    #     :phi_range => phi_range,
    #     :target_phason_requested => target_phason,
    #     :target_phason_actual => target_phason_actual,
    #     :phason_range => phason_range,
    #     :idop_labels => idop_labels_formatted,
    #     :winding_labels => gap_windings,
    #     :run_slice_path => f
    # )
    # push!(label_records, record)


    # ## 6) Gap labelled plots
    # # println("example gap labels:) ", df_slice.gap_labels[1:5])
    # # Define colourscheme
    # unique_qs = sort(unique(g.q for g in Iterators.flatten(df_slice.gap_labels) if !ismissing(g.q)))
    # unique_qs = Int.(unique_qs)
    # # println("typeof(unique_qs) = ", typeof(unique_qs), ", unique_qs = ", unique_qs)
    # q_signed_custom = q_anchor_colormap(unique_qs)
    
    # try
    #     coloured_gaps_outdir = joinpath(outdir, "coloured_gaps")
    #     isdir(coloured_gaps_outdir) || mkpath(coloured_gaps_outdir)

    #     colour_mode = :qled_notp  # Options: :qled_p, :qled_notp, :absq_p, :absq_notp

    #     IDOPPlotting.plt_coloured_gaps_q_optionality(
    #         df_slice,
    #         joinpath(coloured_gaps_outdir, "qled_gaps_mu_vs_phi_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
    #         # p_range=nothing,
    #         # q_max=nothing,
    #         colour_mode=colour_mode,
    #         cmap=q_signed_custom, #:RdBu,
    #         atol=1e-8,
    #         rtol=1e-6,
    #         verbose=false,
    #         # Delta=Delta,
    #         # phason=phason,
    #         # t_n=t_n
    #     )
    # catch e
    #     @warn "plt_coloured_gaps_q_optionality failed"  N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    # end

    df_slice = nothing
    data = nothing
    GC.gc()
end

if !isempty(label_records)
    label_dataset_df = DataFrame(label_records)
    label_dataset_path = joinpath(outdir, "refined_gap_label_records.bson")
    BSON.@save label_dataset_path label_dataset_df
    preview_rows = min(3, nrow(label_dataset_df))
    println("Stored gap label dataset at: ", label_dataset_path)
    if preview_rows > 0
        println("Preview of gap label dataset:")
        # PrettyTables.pretty_table(label_dataset_df[1:preview_rows, [:N, :Delta, :t_n, :phi, :target_phason_actual, :phason_range]])
        PrettyTables.pretty_table(label_dataset_df[1:preview_rows, [:phi, :phason, :idop_labels]])
    end
else
    println("No gap label records collected; dataset export skipped.")
end

println("\nAll IDOP processing complete. Results in: ", outdir)