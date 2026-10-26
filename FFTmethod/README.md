# FFT method

The workflow is: particle tree → pair tree → TH3D source → FFT → C(q) fit.
The pair tree stays compact; the TH3D is made only when transforming it.

## Build and run

Enable FFT in `header_for_all_emissionsource.h`, or build with `FFT=1`:

```bash
make -C drho_analyze_urqmd FFT=1 pairsource_urqmd.exe
make -C FFTmethod FFT=1
```

Requires ROOT and FFTW3. Run the analyser from the directory containing the
converted particle tree, for example `AuAu_27_tree.root`:

```bash
cd drho_analyze_urqmd
./pairsource_urqmd.exe 27 0
cd ..
FFTmethod/transform_fft.exe drho_analyze_urqmd/UrQMD_3d_source_0-10cent_all_27_FFT.root
FFTmethod/CqFit.exe drho_analyze_urqmd/UrQMD_3d_source_0-10cent_all_27_Cq_FFT.root
```

The final file ends in `_Cq_fitresults_FFT.root`. FFT outputs are overwritten
when rerun. All names use `_FFT` to keep them separate from D(rho) results;
custom `--output` names also get this suffix automatically. The fit uses the
per-kT limits `fftQmaxFit[NKT]` in `header_for_all_emissionsource.h` by default.

For a batch run, set the switches at the top of `do_all_urqmd.sh`:

```bash
FFTmethod=true
do_conversion=false # true only if particle trees still need conversion
do_analysis=true
do_transform=true
do_fitting=true
FFT_systematics=true
do_final_plots=true
do_comparison=true # compare with D(rho) when both final ROOT outputs are available
```

Then run `./do_all_urqmd.sh 2` for at most two concurrent jobs. It uses the
energies listed in `energies_low` and `energies_high`, and overrides the
header with `FFT=1` (`FFT=0` when `FFTmethod=false`). FFT averages all stored
events; the D(rho) averaging/rho-range scans and plotting are skipped.
Pair trees, C(q) and fit results end up in `analysed/`, with logs in `logfiles/`.
The script creates all output folders, including `figs/fitting/lcms/FFT`,
`figs/syserr/FFT`, `figs/comparison/Drho_FFT` and `FFTmethod/results`, so a fresh checkout needs no empty
directories copied from another machine. The plotting programs also create
their own output directories when run separately.
Enabled stages regenerate their outputs, including pair trees already in
`analysed/`. Set `do_analysis=false` to reuse pair trees in `analysed/` or
`drho_analyze_urqmd/`, and `do_transform=false` to fit existing C(q) files.
Set `do_final_plots=false` for a transform-only run. For final plotting only,
turn analysis, conversion, transformation and fitting off, leaving
`FFT_systematics=true` and `do_final_plots=true`.
Set `do_comparison=false` to skip the separate method comparison.
Input/processing failures exit with code 1. Invalid fits are saved, the other
variations and final calculator still run, and code 2 reports incomplete results.
No missing variation is silently counted as zero systematic uncertainty.

## Where to edit things

- `header_for_all_emissionsource.h`: FFT switch, grid size, spatial range and
  per-kT C(q) fit limits (`fftQmaxFit`), fit-range scales and grid variations.
  `ikt_indices_to_plot` selects kT bins for energy figures in both final plotters
  and the comparison (currently `{0, 4}`). Rebuild after editing.
- `pair_source_tree.h`: branch variables and `SetBranchAddress` calls only.
- `transform_fft.cc`: settings, tree-to-TH3D loop, normalization, FFT, output,
  in that order inside `main`.
- `CqFit.cpp`: `AltChi2` contains the model and chi-square; `main` selects
  the fit bins, sets Minuit parameters and writes results.
- `calc_and_plot_syserr_FFT.cpp`: read seven fits, calculate directional
  systematics, write graphs/table, and draw the final figures.
- `compare_methods.cpp`: read both final ROOT outputs and draw method comparisons.
- `fft_common.h`: filenames and copying file information.

Default source grid: **251 bins, -100 to +100 fm** on each axis.
One kT grid is processed at a time. The FFT plan is reused. Memory is about
0.9 GiB at this grid size; larger spatial ranges improve q resolution.

The transformer accepts pair trees or top-level TH3 count histograms. Options:

```text
--bins N --range R       change the spatial grid
--binning 0|1|2          baseline/coarse/fine source grid; adds a filename tag
--kt K                  process one kT bin (default all)
--charge positive       or negative / both (default both)
--first-event E --events N   select event indices [E,E+N)
--save-source           also write the unnormalised source TH3D
--cutoff X              suppress small |C-1| values (default off)
--error X               common C(q) bin error (default 1)
--output file_FFT.root  choose a separate output for a variation
```

The initial Au+Au fit limits in the common header are (GeV/c):

```cpp
// kT centers:                 0.20  0.25  0.30  0.35  0.40  0.45  0.50  0.55  0.60  0.65
const double fftQmaxFit[NKT] = {0.10, 0.10, 0.10, 0.12, 0.12, 0.12, 0.14, 0.14, 0.16, 0.16};
```

These are trial values, chosen somewhat below the light-ion examples in
`fromchat.txt`: a larger source suggests a narrower correlation peak, while
decreasing source radii with increasing kT motivates widening the fit range.
For context on the radius trend, see the
[STAR Au+Au study](https://arxiv.org/pdf/1403.4972). The numerical limits above
are our starting guesses, not values taken from that study. Check fit-range
stability for the actual samples; the same array initially applies at all energies
and for all analysis qLCMS cuts. Those analysis cuts remain independent settings.

Both the batch script and standalone fitter use this array automatically.
`--fit-range 1` or `--fit-range 2` applies the narrow/wide scale in the header
and adds `_narrowqFit` or `_wideqFit` to the default fit output name.
For manual variations, `--qmax 0.12` overrides it with a common limit, and
`--qmax 0.10,0.10,0.10,0.12,0.12,0.12,0.14,0.14,0.16,0.16` supplies a separate array.
The kT index is read from `fft_sources`, including when using
`--hist Cq_D_lcms_3d_kt2` to select one histogram. TH3 inputs without kT metadata
need an explicit scalar `--qmax`. Each result row records the actual `qmax` used.
The region is a sphere
`q_out² + q_side² + q_long² <= qmax²`, including q=0 and both signs.
The five parameters are **lambda, Rout, Rside, Rlong, alpha**.
The fit output contains five `TGraphAsymmErrors` ready for a final plotter:
`lambda_vs_KT_FFT`, `R_out_vs_KT_FFT`, `R_side_vs_KT_FFT`, `R_long_vs_KT_FFT`
and `alpha_vs_KT_FFT`. They use `kT_center` on the x axis, with horizontal errors
extending to the kT-bin edges. Vertical errors are the rescaled fit errors,
equal on both sides because the current fitter uses symmetric Hesse errors.
These are effective residual-based errors because Fourier bins are correlated.

Graphs contain valid fits with known kT indices, sorted by kT. Missing or invalid
bins are omitted, so graph point indices need not equal kT-bin indices. Fits at
parameter limits are included if marked valid; their boundary flags remain in
the diagnostic tree. Unlabelled TH3 fits remain available in the tree/TF3 models
but cannot be placed on a kT graph.

A final plotter can read a graph directly, for example:

```cpp
TFile file("analysed/UrQMD_3d_source_0-10cent_all_27_Cq_fitresults_FFT.root");
auto* alpha = file.Get<TGraphAsymmErrors>("alpha_vs_KT_FFT");
alpha->Draw("AP"); // keep file open while using the graph
```

The small `fit_results_FFT` tree is retained for fit status, raw/scaled errors
and the full covariance, alongside TF3 models.
Rerun `CqFit` on existing C(q) files to update the fit results and graphs. Exit code 2
means at least one fit was marked invalid; 1 means an input/processing error.

## Uniform errors and chi-square

Following the [Supplemental Material of arXiv:2610.29760](https://arxiv.org/pdf/2610.29760#page=6),
we fit with a common bin error and estimate its scale from the fit residuals.
Uniform weighting is a convention here: an FFT does not produce independent
bins with a known statistical error. Fourier-bin correlations remain omitted,
so these parameter uncertainties are estimates.

For the input bin error `sigma0` (the transformer defaults to 1), the procedure is:

```text
SSE = sum((Cq - model)^2) over the selected fit sphere
NDF = number of selected bins - number of fitted parameters
chi2_raw = SSE / sigma0^2
error_scale = sqrt(chi2_raw / NDF)
sigma_eff = sigma0 * error_scale = sqrt(SSE / NDF)
errors = errors_raw * error_scale
covariance = covariance_raw * error_scale^2
chi2 = SSE / sigma_eff^2                 # equals NDF by construction
```

A raw `chi2_raw/NDF = 0.00012` with `sigma0 = 1` simply means the residual
scale is about 0.011. It is not a calibrated goodness-of-fit statistic.
The rescaled ratio is 1 by definition and cannot independently demonstrate
that the model is a good description either. `residual_rms` and the data/model
plots retain information about the actual discrepancy.

`CqFit` already applied the parameter/covariance rescaling before this reporting
change; it must not be applied twice. New fit files use `cqfit_error_schema=2`:

| Field | Meaning |
| --- | --- |
| `chi2_raw`, `chi2_ndf_raw` | Before estimating the error scale, using `bin_error_raw` |
| `chi2`, `chi2_ndf` | After rescaling, using `bin_error`; ratio is 1 by construction |
| `bin_error_raw`, `bin_error` | Input placeholder error and effective residual-based error |
| `error_scale` | `sqrt(chi2_ndf_raw)`; applied once to parameter errors |
| `errors_raw`, `covariance_raw` | Hesse results converted to the input error scale |
| `errors`, `covariance` | Final residual-rescaled parameter errors and covariance |
| `residual_rms` | `sqrt(SSE / number_of_selected_bins)` |

The saved TF3 and fit plots use the rescaled chi-square, while raw values remain
available in the tree, plot annotations and log. Files without the schema marker
used `chi2` and `chi2_ndf` for the raw quantities. Rerunning only `CqFit` updates
these diagnostics; the existing C(q) can be reused. Central values and the
already rescaled parameter errors are unaffected. If SSE is exactly zero,
normalization to unity is undefined; chi-square is saved as NaN and the fit is
flagged invalid rather than inventing a finite uncertainty scale.

## Fit figures

`CqFit` also saves one three-panel PNG and PDF for every fitted kT bin, including
invalid fits labelled for inspection. The layout follows `example.png`: adjoining
out/side/long panels, black C(q) points, and a red curve from the single 3D fit.
The data are central-axis slices `(q_out,0,0)`, `(0,q_side,0)`, `(0,0,q_long)`.
The model is solid inside `|q| <= qmax` and dashed outside. The displayed C(q)
errors use the same `sqrt(chi2_raw/NDF)` scaling as the fit. Annotations show
the fitted parameters, rescaled and raw chi-square ratios, fit status and any
parameter-limit flag; the normalized ratio is labelled as such.
Energy, centrality and the analysis qLCMS-cut label come from the input metadata.

When run from the repository root (as the batch script does), figures go to:

```text
figs/fitting/lcms/FFT/<fit-output-filename-without-.root>/Cq_D_lcms_3d_kt0_FFT.png
figs/fitting/lcms/FFT/<fit-output-filename-without-.root>/Cq_D_lcms_3d_kt0_FFT.pdf
```

This keeps energies, analysis cuts and separately named fit variations apart.
`--plot-dir DIR` changes the parent folder (relative to the working directory);
`--no-plots` disables figures. Plotting runs in ROOT batch mode and needs no
display on the server. Each canvas is also saved in the fit ROOT file as
`fitplot_<Cq-histogram-name>_FFT` for inspection or later restyling.

## FFT systematics and final plots

`FFT_systematics=true` runs seven fits per configured energy, each with all
events and each kT bin kept separate. It automatically selects all three
analysis qLCMS cuts. The settings in the common header are:

```cpp
const double fftQmaxFitScale[3] = {1.0, 0.8, 1.2};
const int fftSourceBinsSyst[3] = {fftSourceBins, 151, 351}; // baseline is 251
```

The analysis qLCMS cut remains a source-selection systematic. It is independent
of the Fourier-space fitting limit. The fit-range variations multiply each
bin's own `fftQmaxFit` by 0.8 or 1.2, keeping the grid and analysis cut fixed.
The binning variations rebuild D(rho) from the signed pair tree and repeat the
FFT with 151 or 351 bins, keeping the spatial range at ±100 fm and using the
baseline fit range and analysis cut. Thus five transforms supply seven fits;
there is no Cartesian product of the variations and no event-count scan.

Changing source binning belongs in the transformer: the fitter cannot recover
a finer spatial grid from an already transformed C(q). At fixed spatial range,
this checks the spatial discretization; the q spacing stays
`pi*hbarc/fftSourceHalfRange`. `--binning` needs pairTree input and cannot be
combined with a separate `--bins` override for a coarse/fine run.
The 351-bin transform uses roughly 2.5 GiB per job; the script processes the
variations sequentially within each energy and limits concurrent energies by
`MAXJOBS`.

For `S = UrQMD_3d_source_0-10cent_all_<energy>`, these fit files are in `analysed/`:

| Run | Fit filename after `S` |
| --- | --- |
| Baseline | `_Cq_fitresults_FFT.root` |
| Strict analysis qLCMS | `_strictqLCMS_Cq_fitresults_FFT.root` |
| Loose analysis qLCMS | `_looseqLCMS_Cq_fitresults_FFT.root` |
| Narrow fit sphere | `_Cq_fitresults_narrowqFit_FFT.root` |
| Wide fit sphere | `_Cq_fitresults_wideqFit_FFT.root` |
| Coarse source grid | `_Cq_coarseFFT_fitresults_FFT.root` |
| Fine source grid | `_Cq_fineFFT_fitresults_FFT.root` |

Standalone examples (replace the energy with the sample you have):

```bash
FFTmethod/transform_fft.exe analysed/UrQMD_3d_source_0-10cent_all_27_FFT.root --binning 1
FFTmethod/CqFit.exe analysed/UrQMD_3d_source_0-10cent_all_27_Cq_FFT.root --fit-range 1
FFTmethod/calc_and_plot_syserr_FFT.exe --input-dir analysed --energies 27
```

`--fit-range` also scales an explicitly supplied `--qmax`. Custom `--output`
names override the automatic variation tags, so choose distinct names when
keeping manual variations. For baseline-only batch runs set
`FFT_systematics=false`, `FFT_qlcms_cuts=(0)` and `do_final_plots=false`.
The batch still uses only its `energies_low` and `energies_high` arrays; set
those arrays to the energies available on the machine. In particular, add
`5p2` there if you want the manual 5.2 GeV test included in a batch.

The final calculator uses the same default rule as the D(rho) calculator:
opposite-side deviations contribute separately upwards and downwards;
two deviations on the same side contribute their mean shift. The three
source contributions are then added in quadrature, separately for each side.
Set `use_larger_same_side_shift=true` near the top of the new calculator to
use the larger same-side shift instead. Statistical errors come from the
baseline fit; variations' statistical errors are not added to the systematics.

The outputs are separate from the D(rho) results:

- `FFTmethod/results/syserr_results_FFT.root`: statistical, systematic, total
  and individual-source `TGraphAsymmErrors` for lambda, alpha and the three radii.
  For example, `alpha_stat_27_FFT`, `alpha_syserr_27_FFT`, `alpha_total_27_FFT`,
  `alpha_qLCMS_27_FFT`, `alpha_qFit_27_FFT` and `alpha_binning_27_FFT` use mT
  on the x axis. `alpha_sNN_stat_kt0_FFT` and `alpha_sNN_syserr_kt0_FFT`
  give the collision-energy dependence at fixed kT.
- `FFTmethod/results/syserr_results_breakdown_FFT.csv`: one row per energy,
  kT and parameter, with absolute errors from each source, totals and status.
- `FFTmethod/results/syserr_results_graphs_FFT.txt`: ROOT graph names followed
  by `mT_GeV y err_up err_dn` or `sqrtSNN_GeV y err_up err_dn` and point rows,
  matching the D(rho) text-table convention. Statistical, systematic, total
  and individual-source graphs are exported separately, with 12 significant digits.
- `FFTmethod/results/syserr_results_percent_breakdown_FFT.csv` and
  `syserr_results_percent_table_FFT.tex`: mean relative systematic uncertainties
  per energy, source and parameter. These are arithmetic means over complete
  kT bins with nonzero central values; CSV bin counts make omissions explicit.
  The total is combined in quadrature per kT bin before averaging percentages.
- `figs/syserr/FFT/<energy>/`: mT plots and signed percentage contributions.
  `mT_radii_FFT` combines Rout, Rside and Rlong in three panels with a common
  vertical scale. Alpha, lambda and the contribution plots remain separate.
- `figs/syserr/FFT/`: `mT_overlay_radii_FFT` combines the three radius directions
  in three panels with a common vertical scale and matching energy colors.
  Alpha/lambda energy overlays versus mT and fixed-kT plots versus
  collision energy for `ikt_indices_to_plot`, all as PNG/PDF. Central values
  form connected lines; filled bands connect the systematic-error edges, as
  in the D(rho) plotter. Bars show statistical errors. All kT bins remain in
  the numerical ROOT/text/CSV outputs.

Each ROOT graph's object name agrees with its file key, and its axes carry
the quantity and units. The final calculator consumes the already scaled
baseline parameter errors without any further chi-square rescaling. Custom
`--output` names give the accompanying tables the same filename prefix.

The calculator matches graph points by kT, not graph point index. A point
missing from any required variation is omitted from the final graphs and
marked incomplete in the table. Exit code 2 and `fft_systematics_complete=0`
in the ROOT file report this. Metadata checks also reject mismatched energies,
centralities, event counts, grids, spatial extents, fit limits, or event subsets.
Regenerate old C(q)/fit files once with the updated programs to record these
settings; thereafter fit-only and final-only runs can reuse them.

The calculator also accepts `--plot-dir DIR`, `--output file_FFT.root` and
`--no-plots`. Without `--energies`, it selects energies with baseline fit files.
An optional `--compare-drho syserr_results.root` produces comparisons for alpha
and the three radii. The old file supplies systematic bands only; its N parameter
is not equated to the FFT lambda. Empty legacy points are skipped, and the
D(rho) file is opened for reading only.

## Signs, normalization and D(rho) compatibility

The original branches still store magnitudes. FFT analysis adds the byte
`rho_osl_signs`: bits 1, 2, 4 indicate negative out, side, long coordinates.
Three explicit `if` statements in the transformer restore the signs.
Old trees without this byte need re-analysis from the particle event tree.

Pair analysis discards kT outside `[ktbins[0], ktbins[NKT]]` (currently
0.175–0.675 GeV/c), keeping both endpoints. This applies to both methods.
New files record `source_kt_min` and `source_kt_max`; older files folded
out-of-range pairs into the edge bins. Those files must be regenerated from
the particle trees to apply this cut. Enable analysis, transformation and fitting
in the batch script to regenerate all stages with the new selection.

The convention is `C(q) = 1 + Re FFT(D)`, with `hbarc = 0.197327 GeV fm`.
The source counts are divided by their total **including all flow bins**.
There is no extra bin-volume factor or squared transform. Thus
`C(0) = 1 + retained_fraction`, recorded in `fft_sources` for each histogram.
Check this fraction and grid convergence when choosing a spatial range.
Events are combined by adding pair counts, as in the original D(rho) method.
Centrality still comes from the input sample; the analyser retains its
existing `0-10` label rather than introducing new centrality selection.

The D(rho) fitter can read old trees or the magnitude branches of FFT trees.
Its optional final argument supplies an explicit input filename. Its usual
output names and `calc_and_plot_syserr.cpp` remain as before.

## Comparing D(rho) and FFT

From the repository root, after running both final plotters:

```bash
make -C FFTmethod compare_methods.exe
FFTmethod/compare_methods.exe
```

This reads `syserr_results.root` and
`FFTmethod/results/syserr_results_FFT.root`. No analysis or fits are repeated.
To select files, energies or the kT bins used in energy plots:

```bash
FFTmethod/compare_methods.exe \
  --drho syserr_results.root \
  --fft FFTmethod/results/syserr_results_FFT.root \
  --energies 3p0,5p2,27 --kt 0,4
```

By default, all energies present in both outputs are used. Energy figures use
`ikt_indices_to_plot` from the common header, currently kT bins 0 and 4 (centers
0.20 and 0.40 GeV/c), matching the D(rho) final plotter. The kT plots always
include every valid bin; `--kt` overrides the selection for energy plots.
The plotter creates its directories automatically and saves PNG/PDF figures in
`figs/comparison/Drho_FFT/`:

- `<energy>/kT_alpha_<energy>_comparison_FFT`: single panel overlaying the
  methods' alpha values.
- `<energy>/kT_radii_<energy>_comparison_FFT`: Rout, Rside and Rlong in three
  panels with a common vertical scale, overlaying the methods in each panel.
- `sNN_radii_kt<K>_comparison_FFT`: Rout, Rside and Rlong in three panels
  with a common vertical scale, overlaying the methods.
- Single panels for alpha versus energy and for the derived quantities versus
  kT and energy: `R_RMS`, `Rout_over_Rside`, `Rout2_minus_Rside2`.
- Normalizations overlaid on a single panel versus kT and energy, labelled
  `lambda*` (λ*). The original source-fit N and correlation-fit lambda values
  and errors are plotted directly; their saved graph names remain unchanged.

All normalization-parameter axes use λ*, including the FFT final plots and
the parameter graphs written by `CqFit`.

D(rho) uses solid red lines, filled markers and `/` hatching; FFT uses dashed
blue lines, open markers and `\\` hatching. Lines guide the eye. Connected bands show
systematic uncertainties and bars show the baseline errors saved by each final
plotter. D(rho)'s bars are the existing spreads across accepted event-group fits;
FFT's bars are residual-rescaled fit errors. Their meanings differ.
Bands join the upper and lower systematic errors between adjacent valid points,
as in the D(rho) plotter. A lone valid point uses an uncertainty box.
The D(rho) final plotter now also saves these existing baseline graphs as
`<parameter>_stat_<energy>`. Older ROOT files still work: their baseline bars
are omitted, flagged on the figures, and recorded as `nan` in the CSV.

The averaged radius is defined consistently here as
`R_RMS = sqrt((Rout² + Rside² + Rlong²)/3)`, computed from each method's reported
radii. It does not average over kT or energies. In D(rho) this can differ from
the saved `R` graph, which averages per-fit RMS radii over event groups.
Derived errors use first-order propagation of the saved marginal errors,
separately for baseline and systematic uncertainties. Negative derivatives
swap the asymmetric error directions (for example, increasing Rside decreases
Rout/Rside). Radius correlations are unavailable in these final graphs and are
neglected, so the derived uncertainty bands are approximations.

The plotted graphs, canvases, input filenames and conventions are saved in
`FFTmethod/results/comparison_Drho_FFT.root`, with a tidy CSV beside it.
Graphs use kT coordinates even though the final plotters save mT; missing bins
are matched by coordinate and never replaced with zeros. Negative and zero
values of Rout²−Rside² are retained. Existing input files are read only.
`--output NAME`, `--plot-dir DIR` and `--no-plots` are also available.
The batch script runs the comparison after either method's final plotter when
`do_comparison=true` (the default). It first checks that both final ROOT files
exist and are nonempty; otherwise it prints a skip message. The comparison uses
the script's selected energies and the common header's kT selection for energy
figures. Its log is
`logfiles/comparison_Drho_FFT.log`. Invalid ROOT files or comparison/build
failures exit with code 1; incomplete FFT fits still retain exit code 2 when
the comparison succeeds. FFT processing failures and failures in either final
plotter stop the script before comparison.

To run only the comparison from saved final ROOT files, set these switches
at the beginning of `do_all_urqmd.sh` and then run the script normally:

```bash
FFTmethod=true
do_analysis=false
do_conversion=false
do_transform=false
do_fitting=false
do_final_plots=false
do_comparison=true
```

This builds only the comparison executable and leaves the analyses and final
plotters untouched. For a transform-only or baseline-only run, turn comparison
off too if saved final outputs should not be plotted.

Existing PNG/PDF files from previous layouts are retained; a fresh `--plot-dir`
shows only the figures produced with the current selection.

Run the numerical and compatibility checks with `make -C FFTmethod check`.
The test code lives under `tests/`, separately from the analysis workflow.

This work was supported by the NKFIH grant NKKP ADVANCED 152097.
