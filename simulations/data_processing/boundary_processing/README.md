# Boundary Processing

Post-processing tools for visualizing topological phase boundaries from dynamic search results.

## Overview

These scripts extract and visualize the topological phase boundary from the adaptive refinement search pipeline. The boundary is plotted on a **mu vs rho** phase diagram, where:
- **ρ = t₂/t₁** (hopping ratio parameter)
- **μ** (chemical potential in meV or energy units)
- **Color**: Red = topological phase (ν=1), Blue = trivial phase (ν=0)

## Scripts

### 1. `plot_phase_boundary.jl` — Single Phase Diagram

Creates a complete phase diagram plot from all accumulated data in a run.

**Usage:**
```bash
julia plot_phase_boundary.jl <run_dir> [--iter N] [--save PATH] [--show]
```

**Arguments:**
- `<run_dir>` — Path to the dynamic search run output directory (contains `point_results_compact.csv`)
- `--iter N` — Use only data up to iteration N (default: all data)
- `--save PATH` — Save plot to file instead of displaying
- `--show` — Display plot interactively (default if no `--save`)

**Examples:**

View complete phase diagram:
```bash
julia plot_phase_boundary.jl ../../data_collection/hpc_compatible/dynamic_search_singlejob/runs/dynamic_mb_search_singlejob_20260920_171150 --show
```

Save after first 5 iterations:
```bash
julia plot_phase_boundary.jl ../../data_collection/hpc_compatible/dynamic_search_singlejob/runs/dynamic_mb_search_singlejob_20260920_171150 --iter 5 --save phase_iter5.png
```

### 2. `plot_phase_progression.jl` — Iteration Animation Frames

Creates individual plots for each iteration, showing how the boundary refines over time.

**Usage:**
```bash
julia plot_phase_progression.jl <run_dir> [--save OUTPUT_DIR] [--interval N]
```

**Arguments:**
- `<run_dir>` — Path to dynamic search run
- `--save OUTPUT_DIR` — Save all plots to directory (creates if needed)
- `--interval N` — Save every Nth iteration (default: 1, save all)

**Examples:**

Generate progression plots to directory:
```bash
julia plot_phase_progression.jl ../../data_collection/hpc_compatible/dynamic_search_singlejob/runs/dynamic_mb_search_singlejob_20260920_171150 --save ./phase_frames --interval 2
```

This creates `phase_iter_001.png`, `phase_iter_003.png`, etc., showing boundary refinement step-by-step.

## Data Format

### Input Files

The scripts read from the run directory structure:

```
run_directory/
├── point_results_compact.csv          # All evaluated points
├── iter_001/
│   ├── cells_input.csv               # Cells before refinement
│   ├── cells_next.csv                # Cells after refinement
│   └── cells_resolved.csv            # Fully resolved cells
├── iter_002/
│   ├── cells_resolved.csv
│   └── ...
└── ...
```

### `point_results_compact.csv` Format

| Column | Type | Description |
|--------|------|-------------|
| `point_key` | String | Format `rho\|u` (e.g., `1.5\|0.5`) |
| `phase_label` | Int | 0 = trivial, 1 = topological |
| `topo_fraction` | Float | Fraction of samples classified as topological |
| `invariant_mean` | Float | Mean spectral localiser invariant |
| `gap_min` | Float | Minimum gap across sampled positions |
| `...` | ... | Additional metrics |

## Coordinate Transformation

The scripts internally convert from (ρ, u) coordinates to (ρ, μ) by:

1. **From ρ to t₁, t₂:**
   - If `sweep_mode = :vary_t2`: `t₁ = 1.0`, `t₂ = ρ × t₁`
   - If `sweep_mode = :vary_t1`: `t₂ = 1.0`, `t₁ = t₂ / ρ`

2. **From u to μ:**
   ```
   μ_c(ρ) = 0.5 × 4 × t̄
   where t̄ = (1 - γ) × t₁ + γ × t₂
         γ = golden_ratio_conjugate ≈ 0.618
   μ = u × μ_c(ρ)
   ```

**Note:** To adjust these parameters (e.g., if the run used different settings), edit the constants at the top of each script:

```julia
const GAMMA = 2 / (1 + sqrt(5))  # From sequence
const T1_FIXED = 1.0              # From config
const T2_FIXED = 1.0
const SWEEP_MODE = :vary_t2
```

## Interpreting the Plots

### Boundary Location
- **Phase boundary** lies between red and blue point clouds
- Denser sampling near boundary = higher refinement resolution
- Sparse regions = low topological signal variability

### Point Density
- High density = algorithm identified complex region requiring refinement
- Low density = simple phase region, few boundary points needed

### Iteration Progression (with `plot_phase_progression.jl`)
- Early iterations: Coarse grid, large regions
- Later iterations: Dense boundary trace, coarse interiors
- Stable geometry = algorithm converged

## Requirements

Julia 1.7+ with packages:
- `CSV`
- `DataFrames`
- `Plots` (with compatible backend, e.g., GR or PyPlot)
- `Statistics` (stdlib)

Install if needed:
```julia
julia -e 'using Pkg; Pkg.add(["CSV", "DataFrames", "Plots"])'
```

## Customization

### Change Plot Style

Edit `create_plot()` function in either script to:
- Change colors: `color=:red` → any valid Plots.jl color
- Adjust transparency: `alpha=0.6` → `0.0` to `1.0`
- Modify marker: `markersize=3` → other values or `:circ`, `:square`, etc.
- Change plot size: `size=(900, 600)` → different dimensions

### Generate Contour Plot

To overlay a fitted surface instead of scatter:
```julia
# After loading data, fit a surface to the topological boundary
using Interpolations
# ... fit model
contour!(p, rhos_grid, mus_grid, boundary_surface)
```

(Example implementation available on request.)

## Troubleshooting

**"Results file not found"**
- Ensure `run_directory` is correct and contains `point_results_compact.csv`
- Run should be in progress or completed

**Plot appears empty**
- Check that `mu` values are being calculated correctly
- Verify constants match the actual run configuration

**Julia package not found**
- Install missing dependencies: `julia -e 'using Pkg; Pkg.add("PackageName")'`

## Contact

For modifications or new plotting needs, see the `standard_plotting.jl` utilities in the parent directory for common plotting helpers.
