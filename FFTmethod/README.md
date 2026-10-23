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
FFT_qlcms_cuts=(0 1 2) # or (0) for just the baseline
```

Then run `./do_all_urqmd.sh 2` for at most two concurrent jobs. It uses the
energies listed in `energies_low` and `energies_high`, and overrides the
header with `FFT=1` (`FFT=0` when `FFTmethod=false`). FFT averages all stored
events; the D(rho) averaging/rho-range scans and plotting are skipped.
Pair trees, C(q) and fit results end up in `analysed/`, with logs in `logfiles/`.
Enabled stages regenerate their outputs, including pair trees already in
`analysed/`. Set `do_analysis=false` to reuse pair trees in `analysed/` or
`drho_analyze_urqmd/`, and `do_transform=false` to fit existing C(q) files.
Failed jobs or invalid fits make the batch exit nonzero; fit results marked
invalid are still saved for inspection. FFT final parameter plotting and
systematic aggregation remain to be developed.

## Where to edit things

- `header_for_all_emissionsource.h`: FFT switch, grid size, spatial range and
  per-kT C(q) fit limits (`fftQmaxFit`). Rebuild after editing.
- `pair_source_tree.h`: branch variables and `SetBranchAddress` calls only.
- `transform_fft.cc`: settings, tree-to-TH3D loop, normalization, FFT, output,
  in that order inside `main`.
- `CqFit.cpp`: `AltChi2` contains the model and chi-square; `main` selects
  the fit bins, sets Minuit parameters and writes results.
- `fft_common.h`: filenames and copying file information.

Default source grid: **251 bins, -100 to +100 fm** on each axis.
One kT grid is processed at a time. The FFT plan is reused. Memory is about
0.9 GiB at this grid size; larger spatial ranges improve q resolution.

The transformer accepts pair trees or top-level TH3 count histograms. Options:

```text
--bins N --range R       change the spatial grid
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
and the full covariance, alongside TF3 models. Errors scale by `sqrt(chi2/NDF)`.
Rerun `CqFit` on existing C(q) files to update the fit results and graphs. Exit code 2
means at least one fit was marked invalid; 1 means an input/processing error.

## Fit figures

`CqFit` also saves one three-panel PNG and PDF for every fitted kT bin, including
invalid fits labelled for inspection. The layout follows `example.png`: adjoining
out/side/long panels, black C(q) points, and a red curve from the single 3D fit.
The data are central-axis slices `(q_out,0,0)`, `(0,q_side,0)`, `(0,0,q_long)`.
The model is solid inside `|q| <= qmax` and dashed outside. The displayed C(q)
errors use the same `sqrt(chi2/NDF)` scaling as the fit, and the annotations show
the fitted parameters, raw chi-square, fit status and any parameter-limit flag.
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

Run the numerical and compatibility checks with `make -C FFTmethod check`.
The test code lives under `tests/`, separately from the analysis workflow.

This work was supported by the NKFIH grant NKKP ADVANCED 152097.
