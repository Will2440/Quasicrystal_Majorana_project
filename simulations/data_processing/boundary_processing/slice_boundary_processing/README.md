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
