# Boundary search investigation

## Systematic study, 2026-11-19

The new study treats resolution and kappa as separate axes. The resolution
scan measures 6401 full-range u points at each parameter point with no
prediction, exclusions or transition bisections. Nested grids with 401,
801, 1601 and 3201 nominal points and half-step offsets are reconstructed
from the same fine grid. `resolution_comparison.csv` records observed counts,
physical mu spacing and fine-grid phase runs not hit by each coarser grid.
Offsets retain both endpoints, so shifted grids can have an extra interval
and point. These are comparisons to the finite 6401-point reference, not a
proof of convergence. Even equal counts can hide different intervals or
two unresolved transitions. The reported minimum sample width is the number
of fine-grid samples in a phase-1 run times fine spacing, not a rigorous
lower bound on its true width.

At N=500, Delta=0.05, kappa=1e-4 and u=[0,1.05], the first completed
resolution tables give:

| rho | nominal baseline | mu spacing | observed phase-0 runs, anchored / shifted | fine-grid runs missed, anchored / shifted |
| ---: | ---: | ---: | ---: | ---: |
| 2.0 | 401 | 0.008495 | 8 / 8 | 0 / 0 |
| 2.0 | 1601 | 0.002124 | 8 / 8 | 0 / 0 |
| 2.0 | 6401 | 0.000531 | 8 / reference | 0 / reference |
| 3.5 | 401 | 0.013362 | 14 / 14 | 4 / 4 |
| 3.5 | 801 | 0.006681 | 16 / 18 | 2 / 0 |
| 3.5 | 1601 | 0.003340 | 18 / 18 | 0 / 0 |
| 3.5 | 3201 | 0.001670 | 18 / 18 | 0 / 0 |
| 3.5 | 6401 | 0.000835 | 18 / reference | 0 / reference |
| 5.0 | 1601 | 0.004557 | 7 / 8 | 4 / 3 |
| 5.0 | 3201 | 0.002279 | 11 / 9 | 0 / 2 |
| 5.0 | 6401 | 0.001139 | 11 / reference | 0 / reference |
| 3.5, Delta=0.02 | 1601 | 0.003340 | 13 / 12 | 2 / 3 |
| 3.5, Delta=0.02 | 3201 | 0.001670 | 14 / 15 | 1 / 0 |
| 3.5, Delta=0.02 | 6401 | 0.000835 | 15 / reference | 0 / reference |

For rho=3.5 the narrowest interior phase-1 run contains five fine-grid
samples, a sample-width proxy of 0.004176 in physical mu. A mu spacing
near 0.00334 works for both tested offsets against this reference; spacing
0.00668 does not. This is a case-specific empirical threshold, not a law
relating count to resolution, nor evidence that all true intervals were found.
At rho=5 and at rho=3.5, Delta=0.02, the 3201-point offsets disagree and
the finest grid contains interior topological runs of only one sample. The
rho=5 reference has six interior runs with at most two samples; the low-Delta
reference has four. Neither 6401-point reference can by itself support a
coarser two-offset threshold. The table's last four rows use kappa=1e-4; the two
Delta=0.02 rows explicitly differ from the Delta=0.05 cases.

Job 20238936 completed in 8:03 on `test` using four CPUs, with batch peak
RSS 1180408 KiB. Job 20239243 extends the two unresolved cases to 12801
points, with four CPUs, 8 GiB requested and a 20-minute limit.
Its rho=5 case has since completed: both offsets of the 6401-point grid
find the same 11 runs as its 12801-point reference, with no reference runs
missed. The physical spacings are 0.0011393 and 0.00056965 respectively;
the narrowest interior phase-1 run has only three 12801-grid samples. This
is an empirical 6401-point threshold against that reference, not a proof.
For rho=3.5, Delta=0.02 and kappa=1e-4, both 6401-point offsets match
the 12801-point reference at 15 runs; the nominal physical mu spacing is
0.00083511 and the reference spacing is 0.00041755. The narrowest interior
phase-1 run contains only two fine samples, so this apparent plateau is
particularly fragile. Job 20239243 completed in 10:06 with four CPUs and
peak batch RSS 1230824 KiB.
At rho=2.5 the 6401-point reference finds 14 runs, and the coarsest tested
passing pair is 1601 points at mu spacing 0.0025293; the 801-point offsets
find 13 and 14 runs. An additional intermediate-rho job 20239383 uses the
same four-CPU, 8-GiB, 20-minute test allocation.

Completed rho survey at N=500, Delta=0.05, kappa=1e-4, u=[0,1.05].
"Passing" means both tested grid offsets hit every run in the finite finest
reference and give its count and ordered phases; it does not mean proof of
physical convergence. The minimum sample counts are for interior phase-0 and
phase-1 runs in that reference. A value of 1 or 2 is a warning that the
reference itself may be too coarse.

| rho | finest points | reference phase-0 runs | coarsest passing points | passing mu spacing | min phase-0 / phase-1 samples |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2.0 | 6401 | 8 | 401 | 0.008495 | 68 / 147 |
| 2.5 | 6401 | 14 | 1601 | 0.002529 | 5 / 16 |
| 3.0 | 6401 | 14 | 801 | 0.005870 | 17 / 5 |
| 3.5 | 6401 | 18 | 1601 | 0.003340 | 9 / 5 |
| 4.0 | 6401 | 16 | 3201 | 0.001873 | 19 / 2 |
| 4.5 | 6401 | 14 | 3201 | 0.002076 | 15 / 2 |
| 5.0 | 12801 | 11 | 6401 | 0.001139 | 25 / 3 |

At rho=3.5, Delta=0.02, the 12801-point reference has 15 runs; the
6401-point offsets pass at spacing 0.0008351, and its minimum interior
phase-0 / phase-1 sample counts are 10 / 2. The sampled threshold is not
monotone in rho or count; in particular rho=2.5 and 3.0 have the same
reference count but different thresholds. Do not extrapolate a universal
spacing from this table or treat its decreasing counts as physical evidence.
Job 20239383 completed on `test` in 6:09 with four CPUs and peak batch RSS
1470148 KiB.

`kappa_map.jl` independently evaluates 17 logarithmic kappa values from
1e-7 to 1e-3 at full-grid phase-run midpoints and 33 evenly spaced u points.
It records each individual localizer label, invariant and localizer gap;
no kappa averaging or voting is used. Kappa-window endpoints in the output
are grid-sampled and bracketed by neighboring kappa values, and windows
touching the kappa limits are censored. `summarize_kappa.jl` compares every
sampled label with the same system's kappa=1e-4 labels. That reference is a
measurement, not ground truth. N=300 uses the same mu positions as the
N=500 reference scan, but its agreement comparison uses its own N=300 labels.

Job 20238921 completed on `test` in 1:21 with four CPUs and batch peak RSS
1108484 KiB. At rho=5, Delta=0.05, N=500, 54 sampled points agree with their
1e-4 labels at kappas 3.162e-5, 5.623e-5 and 1e-4. At rho=3.5,
Delta=0.05, 69 points agree only at 1e-4, though all 21 reference-phase-1
points remain phase 1 at several smaller kappas while two reference-phase-0
points change. At rho=3.5, Delta=0.02, 58 points agree only at 1e-4; four
of its twelve reference-phase-1 points have multiple separated phase-1
windows over the sampled kappas. At rho=5, Delta=0.05, N=300, 54 points
agree at 1e-4 and 1.778e-4, but only one reference point is phase 1.

Eight of ten reference-phase-1 windows at rho=5, Delta=0.05, N=500 touch
the smallest tested kappa; their lower edges have not been established.
The exact heuristic Delta/N is 1e-4 for Delta=0.05 but 4e-5 for Delta=0.02;
the logarithmic kappa map brackets the latter rather than sampling it. An
additional full-grid comparison includes that exact value. Kappa stability
at representative points does not ensure the full set of boundaries agrees,
nor identify the physically correct phase without an independent benchmark.

The first grid submission 20238919 failed immediately because of a Julia
conditional-iterator parse error; it performed no simulation. Corrected job
20238936 and aligned kappa-grid job 20238986 use only the test partition.
The scripts and synthetic analysis tests are in this directory; generated
outputs are in the ignored `runs` directory and batch logs in `../logs`.

At rho=5, Delta=0.05, N=500, 3201 aligned baseline points give 14 phase-0
runs for kappa 3.162e-5 versus 11 at 1e-4. Only five of the 3201 sampled
labels differ, all from phase 0 at 1e-4 to phase 1 at 3.162e-5. These five
points were not among the 54 kappa-map points; representative-point agreement
therefore did not establish full-grid count agreement.
At kappa 5.623e-5 the same rho=5 grid has 13 runs and four changed labels
relative to 1e-4. At rho=3.5, Delta=0.02, the exact Delta/N value 4e-5
gives 21 runs compared to 14 at 1e-4 on the identical 3201-point grid; 16
reference-phase-1 samples become phase 0 and nine phase-0 samples become
phase 1. Neither count has been independently established as the correct
physical count. The aligned-grid job 20238986 completed on `test` in 5:04,
four CPUs, peak batch RSS 1232944 KiB.
At 6401 aligned rho=5 points, 3.162e-5 gives 15 runs versus 11 at 1e-4,
and ten measured labels differ (all phase 0 to phase 1). The alternative
kappa count changed from 14 at 3201 points to 15 at 6401, so it too is not
yet shown to be converged. These two fine grids compare kappa sensitivity,
not physical truth.
At rho=3.5, Delta=0.02, 6401 aligned points give 27 runs for Delta/N=4e-5
versus 15 at 1e-4; 44 labels differ, with 27 phase-1 to phase-0 and 17
phase-0 to phase-1 changes. The heuristic count increases from 21 to 27
when refining 3201 to 6401 points, so it remains unresolved. Job 20239303
completed these two fine-kappa scans in 6:40 on `test` using four CPUs,
with batch peak RSS 1222648 KiB. Job 20239404 samples kappa more densely
around 1e-4 for the rho=5 and low-Delta representative points.
Of the 44 changed aligned labels in this low-Delta comparison, 12 have
u>=1 and 32 lie below 1. Kappa sensitivity is therefore not confined to the
edge of the chosen u window; extending the window remains a separate check.

Refined-window job 20239404 completed on `test` in 1:56 using four CPUs,
with batch peak RSS 1114492 KiB. For the 58 sampled rho=3.5, Delta=0.02
points, the contiguous sampled values agreeing with all 1e-4 reference
labels run from 9.170e-5 to 1.24186e-4. The adjacent tested values
8.781e-5 and 1.29684e-4 have disagreements. For the 54 sampled rho=5,
Delta=0.05 points, agreement spans sampled values 3.3096e-5 to 1.2086e-4.
These are discrete agreement sets, not continuous certified windows, and
do not establish a correct reference. In particular full-grid comparisons
already disprove representative-point agreement as a guarantee for rho=5.

Dense cross-check job 20239013 completed on `test` in 15 seconds using four
CPUs, with peak batch RSS 504204 KiB. At six selected N=500 points, including
rho=5 kappa-sensitive points and rho=3.5, Delta=0.02 at kappas 1e-5, 4e-5,
and 1e-4, direct dense eigenvalue signatures agreed with sparse signatures.
The largest relative difference in the localizer gap was approximately
1.6e-11. These checks do not prove every sparse evaluation accurate.

The paper's gap-based kappa bounds rely on assumptions and quantities not
established here. In addition, the localizer-gap routine currently reports
approximate eigenvalues without exposing residual/error bounds. Consequently,
these tests provide numerical convergence evidence, not a rigorous proof of a
global kappa window or a mathematically certified absence of narrower gaps.

### Practical acceptance protocol

1. Select N, sequence, x0, rho, Delta, u window and a single kappa. Scan a
  full baseline, then refine by factors of two. Compare both half-grid
  offsets against a finer reference and require the same phase runs and
  boundary ordering, not merely the same total count. Inspect boundary
  brackets separately at a tolerance small compared with the narrowest
  reported interval; count agreement alone does not converge locations.
2. Record the physical mu step, narrowest sampled phase-0 and phase-1 runs,
  endpoint-censored runs, and missed reference runs. Repeat at representative
  parameter values spanning the intended rho/Delta domain. The largest
  spacing that passes is an empirical threshold at that point, not a proven
  predictor for other points; use a safety factor and explicitly mark cases
  whose finest scan still has one- or two-sample runs as unresolved.
3. Sweep kappas independently across bulk phase interiors and near boundaries,
  including Delta/N, and record every label and localizer gap. Recheck
  candidate kappas on aligned dense mu grids and at other chain sizes and x0
  positions. Compare identities and brackets of the transitions, not only
  representative points or total counts. Near transitions, kappa can move
  boundaries, so exact label agreement over all continuous mu would require
  coincident boundaries; compare well-separated cores while separately
  bounding boundary displacement to the scientific tolerance.
4. For a proof rather than empirical convergence, obtain independently
  validated phase labels or applicable bulk-gap estimates, plus certified
  lower bounds on the localizer gap and numerically bounded kappa/coordinate
  variation. A pointwise estimated gap or a few dense checks cannot rule out
  additional closures between points. If those prerequisites cannot be met,
  report the tested resolution/window and residual uncertainty rather than
  asserting that the heatmap is converged.

For this implementation at fixed rho, Delta, kappa and x0, the localizer
depends linearly on mu and its mu derivative has operator norm one. Weyl's
inequality therefore gives a potential certification route: if a *certified
lower bound* g at a measured mu is positive, the localizer cannot close or
change signature at any mu within distance less than g. At two samples with
certified lower bounds g_left and g_right, g_left + g_right greater than
their separation covers the whole interval and excludes an unseen transition.
Equality is insufficient. A numerical Ritz value for the smallest absolute
eigenvalue is not automatically a certified lower bound, and the current
solver does not expose a bound or check all eigenvalues for missed near-zero
values. Dense cross-checks at six points do not change this. Closing that
verification gap, handling the actual classifier, and checking the stated
assumptions are required before advertising a mathematical certificate.

For reproducible post-processing, after loading Julia 1.10.3 run
`julia --project=. optimisation_notes/summarize_resolution.jl RUN_DIRECTORY`
or `julia --project=. optimisation_notes/compare_aligned.jl REFERENCE_RUN
COMPARISON_RUN` from the boundary-search directory. The latter requires
identical baseline coordinates, physical parameters and relevant configs and
reports every changed sample. `survey_summary.jl RUN_DIRECTORY [...]` reports
the coarsest tested two-offset spacing matching the fine reference, or
`unresolved` if no pair passes. Run
`summarize_kappa.jl` on a completed kappa-map output directory to write
per-point sampled windows and per-kappa agreement tables.

## Constraints and objective

Submit investigation jobs only to `test`, with at most 32 cores and under one
hour. Begin with modest resources. Existing runs and other user jobs are not
modified. The goal is to concentrate or increase sampling in the remaining
topological regions as trivial intervals grow, not merely reduce point counts.

## Initial observations

- The baseline planner already redistributes its point target across retained
  intervals. Exclusion audits and prediction bands cost additional evaluations.
  A proposed change that removed excluded grid points without redistributing
  them was reversed because it reduced discovery resolution.
- At fixed normalized spacing, physical mu spacing grows with the bandwidth
  estimate. Small gaps can therefore become harder to resolve as rho increases.
- Narrow regions can be missed both by finite baseline grids and by finite
  audits inside excluded intervals. Refined boundaries certify observed label
  changes, not completeness of the gap count.
- Existing job 20221156 used 24 CPUs and completed in 47:49; its final rho=10
  slice reported one phase-0 run, no transitions and 1001 baseline points.
  This does not distinguish insufficient resolution, finite-size or kappa
  effects from genuine properties of the measured finite system.

## Calibration

`calibration_config.jl` runs one full-baseline rho=3 slice at N=500,
Delta=0.05 and kappa=1e-4, using 201 points and two evaluation workers.
It is a runtime calibration, not a convergence test. Results are isolated in
this directory's `runs` subdirectory.

Job 20238378 completed on `test` in 1:10, peak batch RSS 743468 KiB. It
observed 13 phase-0 runs with 351 total evaluations; the 201-point grid is
too coarse to establish completeness.

## Independent full-scan comparisons

Job 20238484 used 2 CPUs, 4 GiB requested and 3:55 wall time; peak batch
RSS was 802256 KiB. All cases used N=500, Delta=0.05, u=[0,1.05], a
single kappa, x0 at the midpoint and no exclusions or tracking predictions.
Each is an independent single slice; values count observed phase-0 runs,
including any run touching a search-window endpoint.

| rho | baseline points | kappa | observed runs | total points |
| --- | ---: | ---: | ---: | ---: |
| 3.5 | 401 | 1e-4 | 14 | 536 |
| 3.5 | 1601 | 1e-4 | 18 | 1706 |
| 5.0 | 1601 | 1e-4 | 7 | 1637 |
| 5.0 | 1601 | 3e-4 | 3 | 1613 |

At rho=3.5, increasing the sampling target changes the count; this is direct
evidence that the coarse result is not converged. At rho=5, changing kappa
changes the measured localizer classification, but these two values do not
identify a correct kappa. The window, discretization, finite chain length,
choice of position and phase classification all remain relevant. Comparing
counts alone does not establish that the same boundaries are present.

The earlier run 20221156 first decreased from 19 to 18 at rho=3.31394;
its largest reported count was 21 at rho=3.42828, with 47 observed decreases
over its ordered slices. Some decreases occurred without exclusions, so
skipping is not the sole explanation. These are observations, not proof
that the true physical count must be monotone on a finite chain or finite
mu window.

## Changes and remaining checks

The optional `mu_baseline_spacing_max` increases the full-range baseline
target when the physical mu span grows. It is unset by default, applies to
fallback grids as well, and does not change the localizer or the budget
limits. With exclusions the target is distributed over retained intervals;
separate audit samples do not obey this spacing. This reduces one source of
resolution loss, but does not certify gap counts or guarantee monotonicity.

Coverage, tracking, sweep and integration tests passed. The test fixture
now fixes linear rho spacing to actually exercise its synthetic audit event;
the sweep test no longer assumes a specific user-configured skip default.

The localizer tutorial (https://arxiv.org/abs/2411.03515, sections III.2
and IV.1) describes the localizer gap as a local robustness measure and
explains that too small or too large a kappa can lead to incorrect finite-size
classification. Its sufficient bounds depend on a bulk spectral gap and
system size; they do not justify selecting a universal kappa for these
quasicrystalline scans. Compare individual boundary locations and localizer
gaps across kappa and N, not just aggregate counts.

Additional test job 20238651 compared two smaller kappas at rho=5, a finer
rho=5 grid and two Delta values at rho=3.5. It completed in 6:59 using 2
CPUs, with peak batch RSS 856904 KiB. See `convergence_extra.jl` for the
exact settings.

Completed cases at rho=5 and Delta=0.05: at 1601 baseline points, kappa
3e-5 gave 10 runs and 1e-5 gave 12 runs; at kappa 1e-4, increasing from
1601 to 3201 points increased the count from 7 to 11. For the 1e-4 versus
3e-4 comparison, only 4 of 1601 identical baseline locations changed labels,
all from phase 1 to phase 0, yet the count changed from 7 to 3. This shows
the importance of narrow intervening topological intervals and the need for
joint kappa and resolution convergence. At rho=3.5, Delta=0.02 and kappa
1e-4, the 1601-point full scan observed 13 runs versus 18 at Delta=0.05
and 8 at Delta=0.10; this comparison is not converged and does not establish
a physical decrease.

In the rho=5, 3201-point scan at kappa 1e-4, the narrowest observed interior
phase-1 interval has center-to-center width 0.000246094 in u. Its two
measured transition brackets imply a width between 0.000164063 and
0.000328125 in u; the baseline spacing is 0.000328125 and the transition
tolerance is 0.0001. The interval was detected but is close to the resolution
limit; smaller intervals can be missed. These bounds describe the measured
brackets, not an uncertainty bound on the physical phase diagram.

Finite-size job 20238718 completed on `test` in 1:32 using 2 CPUs, with
peak batch RSS 1038844 KiB. At rho=5, Delta=0.05, kappa=1e-4 and 1601
baseline points, it observed 2 runs at N=300 and 2 at N=400, compared with
7 at N=500 in job 20238484. The sequence prefix and midpoint also change
with N, so this is a finite-system sensitivity check, not an extrapolation
to infinite size. No N>500 calculations were performed.

Matched sweep job 20238660 completed on `test` in 2:31 with peak batch RSS
767376 KiB. It compared five rho values from 2.00 to 2.04 at Delta=0.05,
N=500 and kappa=1e-4, using 401 baseline points and a deliberately loose
onset tolerance so both runs use the same nominal slices. The no-skip and
skip runs both reported 8 phase-0 runs on every slice. Skipping activated
seven exclusions on each of the last three slices, with no fallbacks. All
75 ordered transition pairs had matching orientations and overlapping
brackets; the maximum absolute difference in bracket centers was
6.789e-5 in u. The no-skip run used 2644 points and the skip run used 2694:
audits and changed grid positions cost extra, while the baseline is focused
into retained regions. This limited comparison neither establishes a speedup
nor validates skipping for higher rho or narrower features.

For a defensible heatmap, first compare matched full scans at multiple grid
resolutions, kappas and N<=500; include wider u windows and several x0 values
where feasible. Then compare tracking/skipping against full scans on the same
sweep grid, recording disagreements, fallbacks and point counts. Resolve the
rho/Delta axes independently and flag unconverged cells; do not enforce
monotonicity by replacing measured counts. Full compute-partition campaigns
should only be designed after these checks and launched by the user.

## Candidate user-launched longer runs

After reviewing these results, one diagnostic campaign would hold Delta=0.05
and compare rho=2.0 to 3.5 on 16 linear slices with N=500, 1601 versus 3201
baseline points, initially with skipping disabled and kappa=1e-4. Keep
`rho_onset_tol` larger than 0.1 to compare just the nominal slices; a separate
run with tighter onset tolerance can refine observed changes. Repeat selected
slices with smaller and larger single kappas, compare boundary brackets and
localizer gaps, then compare an otherwise identical skip-enabled sweep. For a
physical-spacing experiment, choose `mu_baseline_spacing_max` from the
desired baseline mu resolution, and budget for the larger number of points
as rho increases. A second direction would fix rho and independently sweep
Delta downwards with the same convergence checks. These are diagnostic
starting points, not validated heatmap production settings or runtime
estimates. Requested wall time, memory, CPU allocation and account for any
compute job must be selected and submitted by the user; no compute jobs were
submitted for this investigation.