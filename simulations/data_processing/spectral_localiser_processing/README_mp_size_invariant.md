# MP size-trend invariant

`new_data_mp_plotting.jl` estimates the sign of the MP size trend for each chemical potential, using only `100 <= N <= 1000` (here, `N` is the system length `L`). Duplicate observations at the same `(mu, N)` are averaged before fitting.

The slope is the Theil-Sen estimate: the median of all pairwise slopes `(MP[j] - MP[i]) / (N[j] - N[i])` across the available sizes. This is less sensitive to isolated numerical outliers than ordinary least-squares fitting and does not assume a linear MP-versus-size relationship. Positive slope is assigned phase 0; negative slope is assigned phase 1; zero or insufficient data is left unclassified.

Use the diagnostics CSV to inspect how well a single sign describes the data at each `mu`:

- `direction_agreement` is the fraction of non-tied size pairs whose MP ordering agrees with the Theil-Sen slope sign. Values near 1 indicate a consistent trend; lower values indicate non-monotonicity or potential outliers.
- `kendall_tau` is the signed concordant-minus-discordant fraction over non-tied pairs. Values near `+1` or `-1` indicate strongly increasing or decreasing trends; values near 0 indicate weak ordering.
- `robust_slope_dMP_dN`, `n_sizes`, and the contributing `min_N`/`max_N` are also recorded for each `mu`.

This sign is a robust trend summary, not proof of strict monotonicity or a linear scaling law. Check the agreement and rank statistic before interpreting a shaded interval strongly. In the generated data, all 1,201 `mu` values had 7–10 sizes; mean direction agreement was 0.916 (range 0.667–1.000), with 457 phase-0 and 744 phase-1 assignments.

## Outputs

The plot overlays the MP curves with green phase-0 and red phase-1 shading. The existing all-size plots are not overwritten.

[Open the invariant plot](../../results/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N%285-10000-8%29_t1%281p0-1p0-1%29_t2%281p5-1p5-1%29_mu%280p0-3p0-601%29_Delta%280p05-0p05-1%29/mp_vs_mu_N100_1000_theilsen_invariant.svg) · [Per-mu diagnostics CSV](../../results/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N%285-10000-8%29_t1%281p0-1p0-1%29_t2%281p5-1p5-1%29_mu%280p0-3p0-601%29_Delta%280p05-0p05-1%29/mp_vs_mu_N100_1000_theilsen_diagnostics.csv)

![MP versus chemical potential with size-trend phase shading](../../results/large_Ns_range/hof_style_slopes_N10000_target_0.61803_tol_0.0001_phason_0.0-1-0.0_nbins1000_npb1_N%285-10000-8%29_t1%281p0-1p0-1%29_t2%281p5-1p5-1%29_mu%280p0-3p0-601%29_Delta%280p05-0p05-1%29/mp_vs_mu_N100_1000_theilsen_invariant.svg)