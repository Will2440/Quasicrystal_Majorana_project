# MP Floor-Normalized Plots

The plotting and processing code is in [`plot_log10_mp_vs_normalized_mu.jl`](plot_log10_mp_vs_normalized_mu.jl). Run it from the repository root with:

```sh
julia --project=. simulations/data_processing/hofstadter_processing/idop/plot_log10_mp_vs_normalized_mu.jl
```

## Data and Axes

The script loads all `.bson` files in the configured raw-data directory. Each file contains a `results_df` table. Rows are grouped by system size `N` and chemical potential `mu`; MP values at duplicate `(N, mu)` points are averaged.

The plotted vertical quantity is

$$
y_N(\mu) = \log_{10}(MP_N(\mu) + 1).
$$

Values for which `MP + 1 <= 0` are undefined and omitted from the multi-size overlays. The horizontal coordinate is `mu / mu_c`.

## Reference Chemical Potential

A shared `mu_c` is found from the `N=1000` data. After sorting its MP curve by `mu`, the script finds the last point where `abs(MP) > MP_ZERO_TOL` (currently `1e-3`) and uses the next sampled `mu` as `mu_c`. Thus, all subsequent sampled MP values remain within the configured zero tolerance. This same `mu_c` is used for every N.

## Noise-Floor Estimate

For each N, the script considers finite `y_N` values in the displayed interval `0 <= mu / mu_c <= 1`. It estimates the curve's floor as the 5th percentile of those values (`MP_FLOOR_QUANTILE = 0.05`). Using a percentile instead of the absolute minimum makes the estimate less sensitive to an isolated numerical dip.

Two combined overlays are written:

- `log10_mp_floor_normalized_overlay.svg` subtracts the estimated floor: `y_N - F_N`. Each curve's estimated floor is at zero; its absolute log-scale range remains size-dependent.
- `log10_mp_range_normalized_overlay.svg` maps the estimated floor to zero and `log10(MP + 1) = 0` to one: `(y_N - F_N) / (-F_N)`. This gives every N a shared normalized vertical scale for comparing dome profiles relative to its estimated floor and the common upper endpoint. Values below the percentile floor remain visible below zero rather than being clipped.

The script also saves one unshifted `log10(MP + 1)` SVG per N. All figures are written under the configured `OUTPUT_DIR` in `simulations/results/`.
