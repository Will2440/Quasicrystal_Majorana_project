
# Slice Boundary Processing

Post-processing tools for the `dynamic_slice_search` output format. These scripts use the physical `mu` values and refined transition brackets already saved by the search; they do not reconstruct coordinates from a config file.

## Scripts

### Export transition data

```bash
julia export_slice_boundaries.jl <run_directory> [output.csv]
```

Reads `transition_boundaries.csv`, joins each transition to its slice gap count, and writes `boundary_data.csv` in the run directory by default. Rows include slice index, rho, tracked edge ID, transition position and bracket in mu, adjacent phase labels, localizer gap, and trivial-interval count.

### Plot the completed boundary and gap-count staircase

```bash
julia plot_slice_boundary.jl <run_directory> [--save output.png] [--points]
```

Plots refined transition positions and their mu bracket widths against rho, connects tracked boundaries, and shows the trivial-interval count in an aligned panel. `--points` overlays all evaluated phase labels.

### Plot algorithm progression

```bash
julia plot_slice_progression.jl <run_directory> [--save output_dir] [--interval N] [--max-slice N]
```

Creates a frame through each selected completed slice. Frames show accumulated evaluated phase points, detected transition tracks, and the latest gap count. If `--save` is omitted, plots are displayed interactively.

### Fit the trivial gap-count curve

```bash
julia fit_slice_gap_count.jl <run_directory> <max_gap_count> <lin|exp|poly> [poly_order] [--plot-excluded true|false] [normalise=true|half_norm=true]
```

Fits the initial point and the first point at each upward step, omitting plateau samples. Those fit points are highlighted in orange; other points under the cutoff are green, and excluded points are grey by default. The fitted equation and parameters appear in the fit-line legend label. Set `normalise=true` to scale the y-axis to 0-1, or `half_norm=true` to scale it to 0-0.5. These options scale plotted and fitted values by the maximum dataset count; the cutoff remains in raw count units. The modes are mutually exclusive and normalized plots receive distinct filename suffixes. `--plot-excluded false` hides points above the cutoff. Polynomial fits require an order, for example `poly 2`. Gap counts are boundary counts minus one and are checked to be integer-valued.

### Extract and plot gap plateaus

```bash
julia plot_gap_plateaus.jl <run_directory> [max_gap_index]
```

Defaults to the first 20 gap births and stops at the first decrease in adjusted trivial gap count. It writes mirrored boundary and gap-plateau CSV/plot outputs. Adjacent phase-transition boundaries identify bounded trivial intervals; intervals are matched between rho slices so births are indexed in rho order. Each gap's rho plateau width is the observed span of rho samples at that count level. Gap centers and their mirrored partners are ranked by μ and assigned normalized y values.

Each rho slice is normalized independently by its own outermost absolute boundary μ, so mirrored boundary rows span ±1 at every rho. The normalized transition CSV contains normalized positions, brackets, and the per-slice scale. The plateau CSV includes `mu_center_normalised` and the onset-slice `mu_normalisation_scale`; `gap_plateau_normalised_mu_overview.png` plots centers normalized using their onset slice's scale and ranks them by normalized μ.

## Input layout

The run directory should contain:

```text
run_directory/
├── slice_summary.csv
├── transition_boundaries.csv
└── slices/
    ├── slice_0001/point_results.csv
    ├── slice_0001/predictions.csv
    └── ...
```

Onset-refinement slices are included in progression in the order they were evaluated. Their transition `track_id` is zero because they were inserted after the initial sequential sweep and are not currently matched into the momentum tracks.

## Requirements

Julia with `CSV`, `DataFrames`, and `Plots`. The exporter and shared reader only require CSV and DataFrames; plotting scripts require Plots and a working plotting backend.
