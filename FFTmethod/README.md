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

The final file ends in `_Cq_fitresults_FFT.root`. All FFT files use `_FFT`
and refuse to overwrite existing files. The fit uses the per-kT limits
`fftQmaxFit[NKT]` in `header_for_all_emissionsource.h` by default.

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
Existing FFT pair trees are reused, including a manual test left in
`drho_analyze_urqmd/`. Set `do_analysis=false` to use only existing pair trees,
or `do_transform=false` to fit existing C(q) files. C(q) and fit outputs still
refuse overwrites. Failed jobs or invalid fits make the batch exit nonzero;
fit results marked invalid are still saved for inspection. FFT final plotting
and systematic aggregation remain to be developed.

## Where to edit things

- `header_for_all_emissionsource.h`: FFT switch, grid size, spatial range and
  per-kT C(q) fit limits (`fftQmaxFit`). Rebuild after editing.
- `pair_source_tree.h`: branch variables and `SetBranchAddress` calls only.
- `transform_fft.cc`: settings, tree-to-TH3D loop, normalization, FFT, output,
  in that order inside `main`.
- `CqFit.cpp`: `AltChi2` contains the model and chi-square; `main` selects
  the fit bins, sets Minuit parameters and writes results.
- `fft_common.h`: three short functions for filenames and copying file information.

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
Results, raw/scaled errors and covariance are saved in `fit_results_FFT`,
alongside TF3 models. Errors scale by `sqrt(chi2/NDF)`; these are effective
residual-based errors because Fourier bins are correlated. Exit code 2
means at least one fit was marked invalid; 1 means an input/processing error.

## Signs, normalization and D(rho) compatibility

The original branches still store magnitudes. FFT analysis adds the byte
`rho_osl_signs`: bits 1, 2, 4 indicate negative out, side, long coordinates.
Three explicit `if` statements in the transformer restore the signs.
Old trees without this byte need re-analysis from the particle event tree.

Pair analysis discards kT outside `[ktbins[0], ktbins[NKT]]` (currently
0.175–0.675 GeV/c), keeping both endpoints. This applies to both methods.
New files record `source_kt_min` and `source_kt_max`; older files folded
out-of-range pairs into the edge bins. Those files must be regenerated from
the particle trees to apply this cut. Move old FFT pair trees and downstream
C(q)/fit files aside first: the batch script reuses existing pair trees.

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
