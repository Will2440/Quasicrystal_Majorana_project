using BSON
using CSV
using DataFrames
using PrettyTables
using Printf
using Logging
using Colors, ColorSchemes, Plots
using ProgressMeter

## data variable for found mp_tol ranges
struct MPTolResult
    N::Int64
    Delta::Float64
    t_n::Tuple{Float64, Float64}
    phi::Float64
    phason::Float64
    qc_region_count::Int64
    mp_tol_ranges::Union{Nothing, Tuple{Float64, Float64}}
end

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
include("../idos/processing.jl")

using .IDOPUnpacking
using .IDOPPlotting
using .IDOPProcessing
using .IDOSProcessing

######################################
########## Path Determination ########
######################################

## Determine the data folder for saving, should be consistent with raw data location
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)_FINALDATASHARED/"
run_set_name = "butterfly_data"

## Ensure data is being lifted from the correct raw data folder
folder_path_hof = normpath(joinpath(@__DIR__, "..", "..", "..", "raw_data", run_set_name, data_folder))

## Defines results directory for saving slices
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "..", "results", run_set_name, data_folder))

######################################
############### Unpack? ##############
######################################

# Define the unpacking mode: :unpack_by_slope or :load_slope
# - :unpack_by_slope: Unpack and slice data with slices per phi (slope), containing range of phason
# - :load_slope: Load existing combo_to_files.bson without unpacking
# -------- IMPORTANT --------
unpack_mode = :load_slope
# ---------------------------

eig_save_mode = :maj        # Options: :all (full spectrum), :maj (central 2 only)
mirror_data_about_mu0 = true  # Whether to mirror data about mu=0 for unpacking

# Define ranges for repacking or plotting -- otherwise the full range will be used.
repack_N_range = nothing
repack_Delta_range = nothing
repack_tn_range = nothing
repack_phason_range = nothing
repack_phi_range = nothing #[0.6180048661800487]

slice_suffix = "_slope"
base_slice_dir = normpath(joinpath(base_results_dir, "idop_base_slices", "N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)$(slice_suffix)_muduplicated-$(mirror_data_about_mu0)"))


if unpack_mode == :unpack_by_slope
    combo_to_files = IDOPUnpacking.unpack_and_slice_by_slope_raw_mp(
        folder_path_hof, 
        base_slice_dir;
        plt_N_range_arg=repack_N_range,
        plt_Delta_range_arg=repack_Delta_range,
        plt_tn_range_arg=repack_tn_range,
        plt_phason_range_arg=repack_phason_range,
        mirror_data_about_mu0=mirror_data_about_mu0
    )
elseif unpack_mode == :load_slope
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
outdir = normpath(joinpath(base_results_dir, "idop", "continuous_mp_N_range$(repack_N_range)_Delta_range$(repack_Delta_range)_tn_range$(repack_tn_range)_phason_range$(repack_phason_range)_muduplicated-$(mirror_data_about_mu0)"))
isdir(outdir) || mkpath(outdir)
println("Saving results to folder: ", outdir)

# Continuous colour-plot normalisation controls
# `plot_normalised_dataset=true` uses the IDOS bandwidth CSV (same approach as main_MBspecialised.jl).
plot_normalised_dataset = true
norm_filename = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/bp_results/MBData_final/hof_style_slopes_N1000_phason_0.0-101-1.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-0p0-1)_Delta(0p05-0p05-1)/idos/mu_rangenothing/plat_thresh0.001,0.004_maxGapLabel20_mptol0.1/norm_typeouter_plotNormeigs_outer_norm/eigs_bandwidth_norms_mu0.0_Delta0.05_t_n[1.0, 2.0]_phason=0.0.csv"
auto_build_norm_csv_if_missing = true

function build_local_norm_csv!(
    combo_to_files,
    csv_path::String,
    fixed_param::Symbol,
;
    eigs_col::Symbol=:eigenvalues,
    tol::Float64=1e-6,
    raw_data_path::Union{String, Nothing}=nothing,
)
    function _build_norm_csv_from_raw!(path::String, out_csv::String)
        println("Building normalization CSV from raw data path: ", path)
        df_raw = IDOPUnpacking.process_bson_files(path)
        isempty(df_raw) && error("Could not build local normalization CSV from raw data: no rows unpacked.")

        raw_cols = propertynames(df_raw)

        if !(:phi in raw_cols)
            df_raw[!, :phi] = getindex.(df_raw.sequence_id, 1)
        end
        if !(:phason in raw_cols)
            df_raw[!, :phason] = getindex.(df_raw.sequence_id, 2)
        end

        IDOSProcessing.compute_eigs_bandwidth_norm!(
            df_raw;
            eigs_col=:eigenvalues,
            out_col=:__eigs_outer_norm_tmp,
            tol=tol,
            min_out_col=:global_min,
            max_out_col=:global_max,
        )

        parameter_cols = [c for c in [:phi, :phason, :mu, :Delta, :t_n, :N] if c in propertynames(df_raw)]
        norm_df = unique(select(df_raw, parameter_cols..., :global_min, :global_max))

        if nrow(norm_df) == 0
            error("Could not build local normalization CSV from raw data: no valid eigenvalue rows found.")
        end

        isdir(dirname(out_csv)) || mkpath(dirname(out_csv))
        CSV.write(out_csv, norm_df)
        println("Built local normalization CSV: ", out_csv, " (rows=", nrow(norm_df), ")")
    end

    first_slice_cols = nothing
    if !isempty(combo_to_files)
        first_file = first(values(combo_to_files))
        first_data = BSON.load(first_file)
        first_df = first_data[:df_slice]
        first_slice_cols = propertynames(first_df)
    end

    if first_slice_cols !== nothing && !(eigs_col in first_slice_cols)
        has_min_pair = (:min_neg_eigs in first_slice_cols) && (:min_pos_eigs in first_slice_cols)
        if !has_min_pair
            if isnothing(raw_data_path)
                error("Normalization fallback requires eigenvalue data, but slices do not contain $(eigs_col), :min_neg_eigs, or :min_pos_eigs.")
            end
            _build_norm_csv_from_raw!(raw_data_path, csv_path)
            return
        end
    end

    norm_dfs = DataFrame[]

    for ((N_val, Delta_val, t_n_tuple, fixed_val), f) in combo_to_files
        data = BSON.load(f)
        df_slice = data[:df_slice]

        if !(fixed_param in propertynames(df_slice))
            df_slice[!, fixed_param] .= fixed_val
        end
        if !(:N in propertynames(df_slice))
            df_slice[!, :N] .= Int(N_val)
        end
        if !(:Delta in propertynames(df_slice))
            df_slice[!, :Delta] .= Float64(Delta_val)
        end
        if !(:t_n in propertynames(df_slice))
            df_slice[!, :t_n] .= string(collect(t_n_tuple))
        end

        eigs_col_eff = eigs_col
        if !(eigs_col in propertynames(df_slice)) && (:min_neg_eigs in propertynames(df_slice)) && (:min_pos_eigs in propertynames(df_slice))
            eigs_col_eff = :__eigenvalues_from_min
            reconstructed = Vector{Union{Vector{Float64}, Missing}}(undef, nrow(df_slice))
            for i in 1:nrow(df_slice)
                negs = df_slice[i, :min_neg_eigs]
                poss = df_slice[i, :min_pos_eigs]

                vals = Float64[]
                if !(negs === missing)
                    for x in negs
                        if !(x === missing)
                            push!(vals, Float64(x))
                        end
                    end
                end
                if !(poss === missing)
                    for x in poss
                        if !(x === missing)
                            push!(vals, Float64(x))
                        end
                    end
                end

                reconstructed[i] = isempty(vals) ? missing : sort!(vals)
            end
            df_slice[!, eigs_col_eff] = reconstructed
        end

        # Use the exact normalization function used in idos/main_bigData.jl
        IDOSProcessing.compute_eigs_bandwidth_norm!(
            df_slice;
            eigs_col=eigs_col_eff,
            out_col=:__eigs_outer_norm_tmp,
            tol=tol,
            min_out_col=:global_min,
            max_out_col=:global_max,
        )

        parameter_cols = [c for c in [:phi, :phason, :mu, :Delta, :t_n, :N] if c in propertynames(df_slice)]
        norm_save_df = select(df_slice, parameter_cols..., :global_min, :global_max)
        norm_save_df = unique(norm_save_df)
        if nrow(norm_save_df) > 0
            push!(norm_dfs, norm_save_df)
        end

        df_slice = nothing
        data = nothing
        GC.gc()
    end

    if isempty(norm_dfs)
        error("Could not build local normalization CSV: no valid eigenvalue rows found.")
    end

    norm_df = unique(vcat(norm_dfs...; cols=:union))
    isdir(dirname(csv_path)) || mkpath(dirname(csv_path))
    CSV.write(csv_path, norm_df)
    println("Built local normalization CSV: ", csv_path, " (rows=", nrow(norm_df), ")")
end

if plot_normalised_dataset && !isfile(norm_filename)
    if auto_build_norm_csv_if_missing
        fallback_norm_dir = joinpath(outdir, "norm_cache")
        generated_norm_filename = joinpath(fallback_norm_dir, "eigs_bandwidth_norms_generated.csv")
        println("Normalization CSV not found at configured path. Building local CSV at: ", generated_norm_filename)
        norm_filename = generated_norm_filename
    else
        error("Normalization CSV not found: $norm_filename")
    end
end


######################################
########## Plotting Execution ########
######################################


# Determine fixed param based on mode for consistent handling
fixed_param = :phason
fixed_param_name = "phason"
varying_param = :phi
varying_param_name = "phi"
x_axis_param = fixed_param  # For plotting: x = phason (fixed per slice), but accumulate across slices

function _norm_csv_has_param_cols(csv_path::String)
    if !isfile(csv_path)
        return false
    end
    norm_df_check = CSV.read(csv_path, DataFrame)
    cols = propertynames(norm_df_check)
    return any(c -> c in cols, [:phi, :phason, :mu, :Delta, :t_n, :N])
end

if plot_normalised_dataset && auto_build_norm_csv_if_missing
    needs_build = !isfile(norm_filename) || !_norm_csv_has_param_cols(norm_filename)
    if needs_build
        if isfile(norm_filename)
            @warn "Normalization CSV exists but has no parameter columns. Rebuilding." norm_filename=norm_filename
        end
        build_local_norm_csv!(combo_to_files, norm_filename, fixed_param; raw_data_path=folder_path_hof)
    end
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

    if !(fixed_param in propertynames(df_slice))
        df_slice[!, fixed_param] .= fixed_val
    end


    # --- FIX: Ensure fixed parameter column exists in DataFrame ---
    # When loading slices by phason, 'phi' is fixed and might not be in the DF columns.
    # When loading slices by slope, 'phason' is fixed and might not be in the DF columns.
    if !(fixed_param in propertynames(df_slice))
        df_slice[!, fixed_param] .= fixed_val
    end

    # --- DEBUG: summarise slice contents ---
    @info "Slice summary" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val n_rows=nrow(df_slice)
    if nrow(df_slice) > 0
        # unique phi values present in this slice
        varying_present = sort(unique(df_slice[!, varying_param]))
        @info "Unique $varying_param_name values in slice:" varying_present

        # mu_range summary
        if :mu_range in propertynames(df_slice)
            mu_lengths = length.(df_slice.mu_range)
            min_mu_len = minimum(mu_lengths)
            max_mu_len = maximum(mu_lengths)
            @info "mu_range vector lengths: min=$min_mu_len, max=$max_mu_len"
            @info "First mu_range: min=$(minimum(df_slice.mu_range[1])), max=$(maximum(df_slice.mu_range[1]))"
        end

        # mp_raw summary
        if :mp_raw in propertynames(df_slice)
            mp_lengths = length.(df_slice.mp_raw)
            min_mp_len = minimum(mp_lengths)
            max_mp_len = maximum(mp_lengths)
            @info "mp_raw vector lengths: min=$min_mp_len, max=$max_mp_len"
            # Print a sample of mp_raw values for the first phi
            @info "Sample mp_raw (first phi):" df_slice.mp_raw[1][1:min(5, length(df_slice.mp_raw[1]))]
        end

            # Print phi and phason ranges
        if :phi in propertynames(df_slice)
            phi_vals = [df_slice.phi[i] for i in 1:nrow(df_slice) if !ismissing(df_slice.phi[i])]
            if !isempty(phi_vals)
                @info "phi range: min=$(minimum(phi_vals)), max=$(maximum(phi_vals))"
            end
        end
        if :phason in propertynames(df_slice)
            phason_vals = [df_slice.phason[i] for i in 1:nrow(df_slice) if !ismissing(df_slice.phason[i])]
            if !isempty(phason_vals)
                @info "phason range: min=$(minimum(phason_vals)), max=$(maximum(phason_vals))"
            end
        end
    else
        @warn "df_slice is empty" N=N_val Delta=Delta_val t_n=t_n $fixed_param_name=fixed_val
    end
    # --- END DEBUG ---

    # filter out phi=0.0 and 1.0 if present (these are trivial cases)
    df_slice = filter(row -> !(row.phi == 0.0 || row.phi == 1.0), df_slice)

    # Optionally apply the same CSV-driven outer normalization used by main_MBspecialised.jl
    x_plot_col = :mu_range
    if plot_normalised_dataset
        IDOPProcessing.apply_outer_normalisation!(
            df_slice,
            norm_filename;
            target_col=:mu_range,
            out_col=:mu_range_norm,
            verbose_warnings=false,
        )
        x_plot_col = :mu_range_norm
    end

    ## 2a) Plot continuous (raw MP value) phase projections
    try
        phase_proj_outdir = joinpath(outdir, "phase_projections")
        isdir(phase_proj_outdir) || mkpath(phase_proj_outdir)

        plot_tag = plot_normalised_dataset ? "NORM" : "RAW"
        plot_xlims = plot_normalised_dataset ? (-1.0, 1.0) : nothing

        IDOPPlotting.plt_continuous_phase_projections(
            df_slice,
            x_plot_col,
            :phi,
            :mp_raw,
            joinpath(phase_proj_outdir, "HIGHRES_CONTINUOUS_$(plot_tag)_phase_proj_N$(N_val)_Delta$(Delta_val)_tn$(t_n)_$fixed_param_name$(fixed_val).png");
            log_colour_scale=true,
            plt_xlims=plot_xlims,
            clims=(-1.0, 0.0),
            cmap=:viridis,
        )

    catch e
        @warn "plt_continuous_phase_projections failed" N_val=N_val Delta_val=Delta_val t_n=t_n $fixed_param_name=fixed_val error=e
    end



    df_slice = nothing
    data = nothing
    GC.gc()
end


println("\nAll IDOP processing complete. Results in: ", outdir)