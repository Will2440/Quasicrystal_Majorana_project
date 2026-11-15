"""
    file name:   main.jl
    created:     2026-09-22
    last edited: 2026-09-22

    overview:
        Main script to process and plot spectral localizer butterfly data.
        Simplified version of the IDOP processing for pure specloc invariant visualization.
    
    structure:
        - Sec 1:  Path Determination
                  Define raw data and results directories
        - Sec 2:  Data Loading
                  Load all chunk files from raw data folder
        - Sec 3:  Data Extraction
                  Extract specloc gap and winding number values
        - Sec 4:  Plotting
                  Create butterfly plots and save to results directory
    
    usage instructions:
        1) Set the data_folder variable to match your raw data folder name
        2) Adjust base_results_dir if needed
        3) Run: julia main.jl
"""

using Logging
using ProgressMeter
using JLD2
using DataFrames

## Set logger level
level = Logging.Info
global_logger(ConsoleLogger(stderr, level))

include("plotting.jl")
include("processing.jl")

using .SpecLocPlotting
using .SpecLocProcessing


######################################
########## Processing Options ########
######################################

## Mirroring: Automatically mirror data about mu=0 if only positive mu detected
mirror_mu_enabled = true

## Plot selection: defaults preserve the full raw + normalized plot set
function parse_plot_options(args)
    selected_x0 = nothing
    plot_mode = :both
    plot_format = "png"
    plot_colour = "072f6b4d"
    plot_nxn = nothing
    plot_linewidth = 1.0
    plot_marker = :square

    for arg in args
        bare_format = lowercase(strip(arg))
        if bare_format in ("png", "svg", "pdf")
            plot_format = bare_format
            continue
        end

        parts = split(arg, '='; limit=2)
        length(parts) == 2 || error("Expected key=value argument, got: $arg")
        key = lowercase(strip(parts[1]))
        value = strip(parts[2])

        if key == "x0"
            selected_x0 = lowercase(value) == "all" ? nothing : parse(Float64, value)
        elseif key == "normalised"
            normalised_value = lowercase(value)
            normalised_value in ("true", "false") || error("normalised must be true or false")
            plot_mode = normalised_value == "true" ? :normalised : :raw
        elseif key == "format"
            plot_format = lowercase(value)
            plot_format in ("png", "svg", "pdf") || error("format must be png, svg, or pdf")
        elseif key in ("colour", "color")
            plot_colour = startswith(value, "#") ? value[2:end] : value
            occursin(r"^(?:[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$", plot_colour) ||
                error("colour must be a 6- or 8-digit hexadecimal RGB(A) value")
        elseif key == "nxn"
            plot_nxn = parse(Int, value)
            plot_nxn > 0 || error("nxn must be a positive pixel dimension")
        elseif key == "linewidth"
            plot_linewidth = parse(Float64, value)
            isfinite(plot_linewidth) && plot_linewidth > 0 || error("linewidth must be a positive finite number")
        elseif key == "marker"
            plot_marker = Symbol(lowercase(value))
            plot_marker in (:circle, :square, :diamond, :utriangle, :dtriangle, :hexagon, :cross, :x, :star5) ||
                error("marker must be circle, square, diamond, utriangle, dtriangle, hexagon, cross, x, or star5")
        else
            error("Unknown argument '$key'. Supported arguments: x0=<value|all>, normalised=true|false, format=png|svg|pdf, colour=<hex>, nxn=<pixels>, linewidth=<size>, marker=<shape>")
        end
    end

    return selected_x0, plot_mode, plot_format, plot_colour, plot_nxn, plot_linewidth, plot_marker
end

selected_x0, plot_mode, plot_format, plot_colour, plot_nxn, plot_linewidth, plot_marker = parse_plot_options(ARGS)
plot_raw_enabled = plot_mode != :normalised
plot_normalised_enabled = plot_mode != :raw
normalize_mu_enabled = plot_normalised_enabled

selected_x0_label = isnothing(selected_x0) ? "all" : string(selected_x0)
resolution_label = isnothing(plot_nxn) ? "default" : "$(plot_nxn)x$(plot_nxn) px"
@info "Plot selection: x0=$selected_x0_label, mode=$plot_mode, colour=#$plot_colour, resolution=$resolution_label, linewidth=$plot_linewidth, marker=$plot_marker"


######################################
########## Path Determination ########
######################################

## Determine the data folder for the raw specloc butterfly data
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p0-1201)_Delta(0p05-0p05-1)_x0(218p75-281p25-5)_kappa(0p0001-0p0001-1)"

## Ensure data is being lifted from the correct raw data folder
folder_path_specloc = normpath(joinpath(@__DIR__, "..", "..", "raw_data", "specloc_butterfly", data_folder))

## Define results directory for saving plots
base_results_dir = normpath(joinpath(@__DIR__, "..", "..", "results", "specloc_butterfly", data_folder))


######################################
############ Data Loading ############
######################################

@info "Loading spectral localizer data from: $folder_path_specloc"

# Load all chunk files
dataframes = SpecLocProcessing.load_chunk_files(folder_path_specloc)

if isempty(dataframes)
    @error "No data loaded from folder: $folder_path_specloc"
    exit(1)
end

@info "Loaded $(length(dataframes)) chunk files"


######################################
######## Data Extraction/Org #########
######################################

@info "Extracting spectral localizer data..."

# Extract specloc data from all chunks
specloc_df = SpecLocProcessing.extract_specloc_data(dataframes)

@info "Extracted $(nrow(specloc_df)) data points"

# Display summary statistics
if nrow(specloc_df) > 0
    n_mu = length(unique(specloc_df.mu))
    n_phi = length(unique(specloc_df.phi))
    n_delta = length(unique(specloc_df.Delta))
    n_x0 = length(unique(specloc_df.x0))
    nominal_grid = n_mu * n_phi * n_delta
    @info "Data summary:"
    @info "  mu range: $(minimum(specloc_df.mu)) to $(maximum(specloc_df.mu))"
    @info "  phi (slope) range: $(minimum(specloc_df.phi)) to $(maximum(specloc_df.phi))"
    @info "  phason range: $(minimum(specloc_df.phason)) to $(maximum(specloc_df.phason))"
    @info "  specloc_invariant range: $(minimum(specloc_df.specloc_invariant)) to $(maximum(specloc_df.specloc_invariant))"
    @info "  specloc_gap range: $(minimum(specloc_df.specloc_gap)) to $(maximum(specloc_df.specloc_gap))"
    @info "  Unique N values: $(unique(specloc_df.N))"
    @info "  Unique Delta values: $(unique(specloc_df.Delta))"
    @info "  Unique x0 values: $(sort(unique(specloc_df.x0)))"
    @info "  Nominal mu×phi×Delta grid: $n_mu × $n_phi × $n_delta = $nominal_grid"
    if nrow(specloc_df) > nominal_grid
        multiplicity = round(nrow(specloc_df) / nominal_grid, digits=2)
        @info "  Detected multiplicity per (mu,phi,Delta): ~$(multiplicity)x (likely from additional swept parameter(s), e.g. x0)"
    end
    
    # Check unique specloc_invariant values
    unique_invariants = sort(unique(skipmissing(specloc_df.specloc_invariant)))
    @info "  Unique specloc_invariant values: $unique_invariants"
    if all(inv -> inv in [0.0, 1.0], unique_invariants)
        @info "  ✓ Confirmed: specloc_invariant contains only 0's and 1's (as expected for a topological invariant)"
    else
        @warn "  ⚠ Unexpected: specloc_invariant contains values other than 0 and 1: $unique_invariants"
    end
else
    @warn "No specloc data was extracted!"
    exit(1)
end

# Aggregate data into two products:
# 1) per-x0 maps (one map per x0)
# 2) averaged-over-x0 maps
agg_by_x0_df = SpecLocProcessing.aggregate_specloc_data(specloc_df; by_mean=true, average_over_x0=false)
agg_avg_x0_df = SpecLocProcessing.aggregate_specloc_data(specloc_df; by_mean=true, average_over_x0=true)

@info "Aggregated per-x0 data into $(nrow(agg_by_x0_df)) records"
@info "Aggregated x0-averaged data into $(nrow(agg_avg_x0_df)) records"


######################################
########## Mirroring (Optional) ######
######################################

if mirror_mu_enabled
    @info "Applying mu mirroring..."
    agg_by_x0_df = SpecLocProcessing.mirror_mu_data(agg_by_x0_df)
    agg_avg_x0_df = SpecLocProcessing.mirror_mu_data(agg_avg_x0_df)
    @info "After mirroring: per-x0 data has $(nrow(agg_by_x0_df)) records, x0-avg has $(nrow(agg_avg_x0_df)) records"
end


######################################
####### Normalization (Optional) #####
##  Using the IDOP normalization cache ##
######################################

# Initialize bounds outside the if block to avoid soft scope ambiguity
eigs_bounds_df = nothing

if normalize_mu_enabled
    @info "="^50
    @info "Computing eigenvalue normalization from raw butterfly_data..."
    @info "="^50
    
    norm_start = time()
    
    norm_cache_path = normpath(joinpath(@__DIR__, "..", "..", "results", "butterfly_data",
        "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-5p0-146627)_Delta(0p05-0p05-1)_FINALDATASHARED",
        "idop", "continuous_mp_N_rangenothing_Delta_rangenothing_tn_rangenothing_phason_rangenothing_muduplicated-true",
        "norm_cache", "eigs_bandwidth_norms_generated.csv"))
    
    if !isfile(norm_cache_path)
        @warn "IDOP normalization cache not found at: $norm_cache_path"
        @warn "Normalization disabled. Continuing with raw mu axis..."
        normalize_mu_enabled = false
    else
        try
            global eigs_bounds_df
            @info "Loading IDOP mu=0 eigenvalue bounds from: $norm_cache_path"
            eigs_bounds_df = SpecLocProcessing.load_eigs_bandwidth_norms(norm_cache_path)
            @info "  Loaded bounds for $(nrow(eigs_bounds_df)) parameter combinations"
        catch e
            @error "Failed to load IDOP eigenvalue normalization cache" error=e
            global eigs_bounds_df
            eigs_bounds_df = nothing
        end
        
        # If extraction succeeded, apply normalization
        if !isnothing(eigs_bounds_df)
            @info "Applying eigenvalue-based normalization to per-x0 data..."
            agg_by_x0_df = SpecLocProcessing.apply_eigs_bandwidth_normalisation!(
                agg_by_x0_df,
                eigs_bounds_df
            )
            
            @info "Applying eigenvalue-based normalization to x0-averaged data..."
            agg_avg_x0_df = SpecLocProcessing.apply_eigs_bandwidth_normalisation!(
                agg_avg_x0_df,
                eigs_bounds_df
            )
            
            norm_elapsed = time() - norm_start
            @info "✓ Eigenvalue-based normalization complete in $(round(norm_elapsed, digits=1))s"
        else
            @warn "Normalization disabled. Continuing with raw mu axis..."
            normalize_mu_enabled = false
        end
    end
end


######################################
########## Plot Creation #############
######################################

# Show data statistics before plotting
@info "="^50
@info "Data Statistics:"
@info "  Original extracted points: $(nrow(specloc_df))"
@info "  Aggregated per-x0 points: $(nrow(agg_by_x0_df))"
@info "  Aggregated x0-averaged points: $(nrow(agg_avg_x0_df))"
@info "  Reduction factor (to per-x0): $(round(nrow(specloc_df) / nrow(agg_by_x0_df), digits=1))x"
@info "  Reduction factor (to x0-avg): $(round(nrow(specloc_df) / nrow(agg_avg_x0_df), digits=1))x"
@info "  Unique mu values: $(length(unique(agg_by_x0_df.mu)))"
@info "  Unique phi values: $(length(unique(agg_by_x0_df.phi)))"
@info "  Unique x0 values: $(length(unique(agg_by_x0_df.x0)))"
@info "="^50

# Determine fixed param based on the data
# Since we have one phason, one rho (t2/t1), one Delta, we slice by phi (slope)
fixed_param = :phi
varying_param = :mu

# Create output directory
outdir = normpath(joinpath(base_results_dir, "butterfly_plots"))
isdir(outdir) || mkpath(outdir)
outdir_per_x0 = joinpath(outdir, "per_x0")
outdir_avg_x0 = joinpath(outdir, "x0_averaged")
isdir(outdir_per_x0) || mkpath(outdir_per_x0)
isdir(outdir_avg_x0) || mkpath(outdir_avg_x0)

# Create subdirectories for normalized plots if needed
if normalize_mu_enabled
    outdir_per_x0_norm = joinpath(outdir_per_x0, "normalized")
    outdir_avg_x0_norm = joinpath(outdir_avg_x0, "normalized")
    isdir(outdir_per_x0_norm) || mkpath(outdir_per_x0_norm)
    isdir(outdir_avg_x0_norm) || mkpath(outdir_avg_x0_norm)
end

@info "Creating plots in: $outdir"


## PHASE 1: Per-x0 plots (full resolution, computed from fresh per-x0 data)
@info "="^50
@info "PHASE 1: Creating FULL RESOLUTION PER-X0 PLOTS..."
@info "="^50

per_x0_plot_df = agg_by_x0_df
if !isnothing(selected_x0)
    per_x0_plot_df = filter(row -> isapprox(row.x0, selected_x0; atol=1e-6, rtol=0), agg_by_x0_df)
    isempty(per_x0_plot_df) && error("No specloc data found for x0=$selected_x0")
end

@showprogress "Plotting per-x0 butterflies..." for (k, df_subset) in pairs(groupby(per_x0_plot_df, [:Delta, :x0]))
    Δ = df_subset.Delta[1]
    x0_val = df_subset.x0[1]
    N_val = df_subset.N[1]

    try
        x0_str = replace(string(round(x0_val, digits=6)), "." => "p")
        if plot_raw_enabled
            plot_start = time()
            plot_path = joinpath(outdir_per_x0, "specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), "." => "p"))_x0$(x0_str).$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path;
                title_suffix="(Δ=$(Δ), x0=$(round(x0_val, digits=6)))",
                value_col=:specloc_invariant,
                mu_col=:mu,
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            elapsed = time() - plot_start
            @info "✓ Per-x0 plot (raw) for Δ=$(Δ), x0=$(round(x0_val, digits=6)) saved in $(round(elapsed, digits=2))s"
        end
        
        # NORMALIZED MU AXIS PLOT (if enabled)
        if plot_normalised_enabled && normalize_mu_enabled && ("mu_specloc_invariant_norm" in names(df_subset))
            plot_start = time()
            plot_path_norm = joinpath(outdir_per_x0_norm, "specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), "." => "p"))_x0$(x0_str)_mu_normalized.$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path_norm;
                title_suffix="(Δ=$(Δ), x0=$(round(x0_val, digits=6)), μ NORMALIZED)",
                value_col=:specloc_invariant,
                mu_col=:mu_specloc_invariant_norm,
                plt_xlims=(-1.0, 1.0),
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            elapsed = time() - plot_start
            @info "✓ Per-x0 plot (normalized μ axis) for Δ=$(Δ), x0=$(round(x0_val, digits=6)) saved in $(round(elapsed, digits=2))s"
        end
    catch e
        @warn "Failed to plot per-x0 butterfly for Δ=$(Δ), x0=$(x0_val): $e"
    end
end

if isnothing(selected_x0)
@info ""
@info "Per-x0 plots complete. Now generating x0-averaged plots..."
@info ""


## PHASE 2: Fast diagnostic plot (x0-averaged, downsampled)
@info "="^50
@info "PHASE 2: Creating FAST DIAGNOSTIC PLOT (x0-averaged, downsampled 5x)..."
@info "="^50

fast_plot_start = time()

try
    downsampled_df = SpecLocProcessing.downsample_data(agg_avg_x0_df, 5)
    
    for (Delta_val, df_subset) in pairs(groupby(downsampled_df, :Delta))
        Δ = df_subset.Delta[1]
        N_val = df_subset.N[1]
        
        if plot_raw_enabled
            plot_path = joinpath(outdir_avg_x0, "FAST_specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), "." => "p"))_x0_averaged_downsampled.$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path;
                title_suffix="(Δ=$(Δ), x0-averaged, DOWNSAMPLED 5x)",
                value_col=:specloc_invariant,
                mu_col=:mu,
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            @info "✓ Fast diagnostic plot (raw) saved in $(round(time() - fast_plot_start, digits=2))s"
        end
        
        # NORMALIZED MU AXIS DIAGNOSTIC PLOT (if enabled)
        if plot_normalised_enabled && normalize_mu_enabled && ("mu_specloc_invariant_norm" in names(df_subset))
            plot_path_norm = joinpath(outdir_avg_x0_norm, "FAST_specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), "." => "p"))_x0_averaged_downsampled_mu_normalized.$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path_norm;
                title_suffix="(Δ=$(Δ), x0-averaged, μ NORMALIZED, DOWNSAMPLED 5x)",
                value_col=:specloc_invariant,
                mu_col=:mu_specloc_invariant_norm,
                plt_xlims=(-1.0, 1.0),
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            @info "✓ Fast diagnostic plot (normalized μ axis) saved in $(round(time() - fast_plot_start, digits=2))s"
        end
    end
catch e
    @warn "Fast diagnostic plot failed: $e"
end

@info ""
@info "If the fast diagnostic plot looks good, the full-resolution plot will render next."
@info "Full resolution plot may take 5-15 minutes depending on system..."
@info ""


## PHASE 3: Full resolution x0-averaged plots
@info "="^50
@info "PHASE 3: Creating FULL RESOLUTION x0-AVERAGED PLOTS..."
@info "="^50

@showprogress "Plotting x0-averaged butterflies..." for (Delta_val, df_subset) in pairs(groupby(agg_avg_x0_df, :Delta))
    Δ = df_subset.Delta[1]
    N_val = df_subset.N[1]

    try
        if plot_raw_enabled
            plot_start = time()
            plot_path = joinpath(outdir_avg_x0, "specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), ".", "p"))_x0_averaged.$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path;
                title_suffix="(Δ=$(Δ), x0-averaged)",
                value_col=:specloc_invariant,
                mu_col=:mu,
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            elapsed = time() - plot_start
            @info "✓ x0-averaged plot (raw) for Δ=$(Δ) saved in $(round(elapsed, digits=2))s"
        end
        
        # NORMALIZED MU AXIS PLOT (if enabled)
        if plot_normalised_enabled && normalize_mu_enabled && ("mu_specloc_invariant_norm" in names(df_subset))
            plot_start = time()
            plot_path_norm = joinpath(outdir_avg_x0_norm, "specloc_invariant_butterfly_N$(N_val)_Delta$(replace(string(Δ), "." => "p"))_x0_averaged_mu_normalized.$(plot_format)")
            SpecLocPlotting.plot_specloc_butterfly(
                df_subset,
                plot_path_norm;
                title_suffix="(Δ=$(Δ), x0-averaged, μ NORMALIZED)",
                value_col=:specloc_invariant,
                mu_col=:mu_specloc_invariant_norm,
                plt_xlims=(-1.0, 1.0),
                nonwhite_colour=plot_colour,
                nxn=plot_nxn,
                linewidth=plot_linewidth,
                marker=plot_marker
            )
            elapsed = time() - plot_start
            @info "✓ x0-averaged plot (normalized μ axis) for Δ=$(Δ) saved in $(round(elapsed, digits=2))s"
        end
    catch e
        @warn "Failed to plot x0-averaged butterfly for Δ=$(Δ): $e"
    end
end
end


# ## 2) Create combined plot with all data (if multiple Delta values exist)
# if length(unique(agg_df.Delta)) > 1
#     try
#         plot_path = joinpath(outdir, "specloc_invariant_combined_all_Delta.png")
#         SpecLocPlotting.plot_specloc_butterfly(
#             agg_df,
#             plot_path;
#             title_suffix="(All Δ)"
#         )
#         @info "Created combined plot for all Delta values"
#     catch e
#         @warn "Failed to create combined plot: $e"
#     end
# end


######################################
############ Summary #################
######################################

@info "Specloc butterfly processing complete!"
@info "Results saved to: $outdir"
@info "Total per-x0 data points processed: $(nrow(agg_by_x0_df))"
@info "Total x0-averaged data points processed: $(nrow(agg_avg_x0_df))"
@info "mu range: $(minimum(agg_by_x0_df.mu)) to $(maximum(agg_by_x0_df.mu))"
@info "phi (slope) range: $(minimum(agg_by_x0_df.phi)) to $(maximum(agg_by_x0_df.phi))"
@info "phason range: $(minimum(agg_by_x0_df.phason)) to $(maximum(agg_by_x0_df.phason))"
