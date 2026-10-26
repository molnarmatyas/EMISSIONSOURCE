#!/bin/bash

# Example, if ~40--60 GB RAM & >3 cores available: TODO optimize fitting code to be less memory hungry, while not sacrificing speed
# time ./do_all_urqmd.sh 3 &> logdoall.log &

# Run ./do_all_urqmd.sh [MAXJOBS]; select stages and method here.
FFTmethod=true # true: signed pair trees -> FFT -> C(q) fits; false: original D(rho) workflow
do_conversion=false # whether to convert from .f19 to root tree files
do_analysis=false
do_transform=true # FFT only; false reuses C(q) files already in analysed/
do_fitting=true # fit C(q) in FFT mode, otherwise run the original Levy fits
FFT_systematics=true # FFT: all three qLCMS cuts, two fit-range and two grid variations
do_final_plots=true # FFT: combine systematics and draw the final parameter plots
do_comparison=true # both methods: compare final outputs if both ROOT files exist

# FFT fit limits: edit fftQmaxFit[NKT] in header_for_all_emissionsource.h.
FFT_qlcms_cuts=(0 1 2) # used when FFT_systematics=false; use (0) for just the baseline

# Base directory of this script (absolute)
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$BASEDIR" || exit 1
# Sanity check
echo "Base directory: $BASEDIR"
sleep 10

# The splitting has to be set in the common header as well!!!
energies_low=("3p0" "3p2" "3p5" "3p9" "4p5") # these may be with higher statistics - 100k full instead of 10k
energies_high=("7p7" "9p2" "11p5" "14p5" "19p6" "27")
energies=("${energies_low[@]}" "${energies_high[@]}")
energy_list=$(IFS=,; echo "${energies[*]}")

# The script's method choice overrides the header for this build.
fft_build=0
source_suffix=""
qlcms_types=(0 1 2)
qlcms_suffixes=("" "_strictqLCMS" "_looseqLCMS")
if [ "$FFTmethod" = true ]; then
  fft_build=1
  source_suffix="_FFT"
  qlcms_types=("${FFT_qlcms_cuts[@]}")
  if [ "$FFT_systematics" = true ]; then qlcms_types=(0 1 2); fi
fi
if [ "${#qlcms_types[@]}" -eq 0 ]; then
  echo "Select at least one qLCMS cut."
  exit 1
fi
for cut in "${qlcms_types[@]}"; do
  case "$cut" in
    0|1|2) ;;
    *) echo "qLCMS cuts must be 0, 1 or 2."; exit 1 ;;
  esac
done

mkdir -p analysed logfiles figs/fitting/lcms/FFT figs/syserr/FFT figs/comparison/Drho_FFT FFTmethod/results \
         levyfit/exe levyfit/object/deps levyfit/results || exit 1

# Called after either method's final plotter; also works with saved results.
compare_methods() {
  if [ "$do_comparison" != true ]; then return 0; fi
  local drho="$BASEDIR/syserr_results.root"
  local fft="$BASEDIR/FFTmethod/results/syserr_results_FFT.root"
  local file missing=false
  for file in "$drho" "$fft"; do
    if [ ! -s "$file" ]; then
      echo "Skipping method comparison: missing or empty $file"
      missing=true
    fi
  done
  if [ "$missing" = true ]; then return 0; fi
  make -C "$BASEDIR/FFTmethod" compare_methods.exe || return 1
  echo "Comparing D(rho) and FFT for energies: $energy_list"
  "$BASEDIR/FFTmethod/compare_methods.exe" --drho "$drho" --fft "$fft" --energies "$energy_list" \
    > "$BASEDIR/logfiles/comparison_Drho_FFT.log" 2>&1 || {
      echo "Method comparison failed; see logfiles/comparison_Drho_FFT.log."
      return 1
    }
  echo "Comparison figures: figs/comparison/Drho_FFT/; ROOT/CSV: FFTmethod/results/comparison_Drho_FFT.*"
}

# JOB CONTROL ------
# Enable job control so `jobs` builtin works inside this script
set -m

# Ensure background children are killed if the script exits or is interrupted
trap 'echo "Script exiting - terminating background jobs..."; jobs -p | xargs -r kill; wait' EXIT SIGINT SIGTERM

# Simple job limiter: keep at most MAXJOBS background processes in this script
# Usage: ./do_all_urqmd.sh [MAXJOBS]  (optional). Default is 3 to be safe with memory.
if [ -n "$1" ]; then
  MAXJOBS="$1"
else
  MAXJOBS=3
fi
echo "MAXJOBS set to $MAXJOBS"
if ! [[ "$MAXJOBS" =~ ^[1-9][0-9]*$ ]]; then
  echo "MAXJOBS must be a positive integer."
  exit 1
fi
limit_jobs() {
  # Wait while number of background jobs is >= MAXJOBS
  while [ "$(jobs -p | wc -l)" -ge "$MAXJOBS" ]; do
    sleep 2
  done
}

# Wait for every job in a conversion/analysis/FFT stage and report failures.
job_pids=()
wait_for_stage() {
  local failed=0 status
  for pid in "${job_pids[@]}"; do
    if wait "$pid"; then :; else
      status=$?
      if [ "$status" -eq 2 ] && [ "$failed" -ne 1 ]; then failed=2; else failed=1; fi
    fi
  done
  job_pids=()
  return "$failed"
}

# Helper to run a command in background and redirect stdout/stderr to logfile
run_bg() {
  local cmd="$1"
  local logfile="$2"
  eval "$cmd" &> "$logfile" &
  # throttle
  limit_jobs
}

# Move files matching a glob pattern only if matches exist.
move_matching_png() {
  local pattern="$1"
  local dest="$2"
  shopt -s nullglob
  local files=( $pattern )
  shopt -u nullglob
  if [ "${#files[@]}" -gt 0 ]; then
    mv "${files[@]}" "$dest"/
  fi
}
# END JOB CONTROL ------

cd "$BASEDIR/drho_analyze_urqmd" || exit 1
if [ "$do_analysis" = true ]; then
  make FFT="$fft_build" pairsource_urqmd.exe || exit 1
fi
if [ "$do_conversion" = true ]; then
  make converter_f19.exe || exit 1
fi

# BEWARE - NOT SURE IF ENOUGH RESOURCES FOR CONVERSION OR ANALYSIS IN PARALLEL
if [ "$do_conversion" = true ]; then
  for ienergy in "${energies[@]}"; do
    echo "Converting .f19 to .root for energy: ${ienergy}"
    ./converter_f19.exe "${ienergy}"  &> "${ienergy}.log" &
    job_pids+=("$!")
    limit_jobs
  done
  wait_for_stage || { echo "Conversion failed; see converter logs."; exit 1; }
fi

analysedname="UrQMD_3d_source_0-10cent_all_"
if [ "$do_analysis" = true ]; then
  for ienergy in "${energies[@]}"; do
    echo "Creating source for energy: ${ienergy}"
    for cut in "${qlcms_types[@]}"; do
      ./pairsource_urqmd.exe "$ienergy" "$cut" > "$BASEDIR/logfiles/analysis_${ienergy}${qlcms_suffixes[$cut]}${source_suffix}.log" 2>&1 &
      job_pids+=("$!")
      limit_jobs
    done
  done
  wait_for_stage || { echo "Pair analysis failed; see logfiles/analysis_*.log."; exit 1; }
  # Move only the requested pair trees, after all writers have finished.
  for ienergy in "${energies[@]}"; do
    for cut in "${qlcms_types[@]}"; do
      source_file="${analysedname}${ienergy}${qlcms_suffixes[$cut]}${source_suffix}.root"
      if [ -f "$source_file" ]; then mv "$source_file" "$BASEDIR/analysed/" || exit 1; fi
    done
  done
fi
cd "$BASEDIR" || exit 1

# FFT workflow: one fully event-averaged source per energy/qLCMS cut.
# C(q) and fit results go to analysed/; no D(rho) plot folders are cleared.
if [ "$FFTmethod" = true ]; then
  if [ "$do_transform" = true ] || [ "$do_fitting" = true ] || [ "$do_final_plots" = true ]; then
    make -C "$BASEDIR/FFTmethod" FFT=1 || exit 1
  fi
  range_suffixes=("" "_narrowqFit" "_wideqFit")
  binning_suffixes=("" "_coarseFFT" "_fineFFT")
  # A failed Minuit fit is saved and reported, while the remaining variations run.
  fit_fft() {
    "$BASEDIR/FFTmethod/CqFit.exe" "$1" --fit-range "$2" > "$3" 2>&1
    case $? in
      0) ;;
      2) fft_fit_status=2 ;;
      *) return 1 ;;
    esac
  }
  for energy in "${energies[@]}"; do
    (
      fft_fit_status=0
      for cut in "${qlcms_types[@]}"; do
        stem="${analysedname}${energy}${qlcms_suffixes[$cut]}"
        source_file="$BASEDIR/analysed/${stem}_FFT.root"
        if [ ! -f "$source_file" ]; then source_file="$BASEDIR/drho_analyze_urqmd/${stem}_FFT.root"; fi
        cq_file="$BASEDIR/analysed/${stem}_Cq_FFT.root"
        if [ "$do_transform" = true ]; then
          echo "Transforming $energy, qLCMS=$cut"
          "$BASEDIR/FFTmethod/transform_fft.exe" "$source_file" --binning 0 --output "$cq_file" > "$BASEDIR/logfiles/transform_${energy}${qlcms_suffixes[$cut]}_FFT.log" 2>&1 || exit 1
        fi
        if [ "$do_fitting" = true ]; then
          echo "Fitting C(q) for $energy, qLCMS=$cut"
          ranges=(0)
          if [ "$cut" -eq 0 ] && [ "$FFT_systematics" = true ]; then ranges=(0 1 2); fi
          for range in "${ranges[@]}"; do
            fit_fft "$cq_file" "$range" "$BASEDIR/logfiles/fit_${energy}${qlcms_suffixes[$cut]}${range_suffixes[$range]}_FFT.log" || exit 1
          done
        fi
      done
      if [ "$FFT_systematics" = true ]; then
        stem="${analysedname}${energy}"
        source_file="$BASEDIR/analysed/${stem}_FFT.root"
        if [ ! -f "$source_file" ]; then source_file="$BASEDIR/drho_analyze_urqmd/${stem}_FFT.root"; fi
        for grid in 1 2; do
          cq_file="$BASEDIR/analysed/${stem}_Cq${binning_suffixes[$grid]}_FFT.root"
          if [ "$do_transform" = true ]; then
            echo "Transforming $energy, grid=$grid"
            "$BASEDIR/FFTmethod/transform_fft.exe" "$source_file" --binning "$grid" --output "$cq_file" > "$BASEDIR/logfiles/transform_${energy}${binning_suffixes[$grid]}_FFT.log" 2>&1 || exit 1
          fi
          if [ "$do_fitting" = true ]; then
            fit_fft "$cq_file" 0 "$BASEDIR/logfiles/fit_${energy}${binning_suffixes[$grid]}_FFT.log" || exit 1
          fi
        done
      fi
      exit "$fft_fit_status"
    ) &
    job_pids+=("$!")
    limit_jobs
  done
  wait_for_stage
  fft_status=$?
  if [ "$fft_status" -eq 1 ]; then echo "FFT processing failed; see logfiles/*_FFT.log."; exit 1; fi
  if [ "$FFT_systematics" = true ] && [ "$do_final_plots" = true ]; then
    "$BASEDIR/FFTmethod/calc_and_plot_syserr_FFT.exe" --input-dir "$BASEDIR/analysed" --energies "$energy_list" > "$BASEDIR/logfiles/syserr_FFT.log" 2>&1
    case $? in
      0) ;;
      2) fft_status=2 ;;
      *) echo "FFT systematics failed; see logfiles/syserr_FFT.log."; exit 1 ;;
    esac
  fi
  compare_methods || exit 1
  echo "FFT stages finished. Outputs are in analysed/; logs are in logfiles/*_FFT.log."
  if [ "$do_fitting" = true ]; then echo "Fit figures are in figs/fitting/lcms/FFT/."; fi
  if [ "$fft_status" -ne 0 ]; then echo "Some fits or systematic points are incomplete; see logs and the FFT results table."; fi
  exit "$fft_status"
fi

#root.exe -b -q 'Average_Drho.cpp(100)' # why doesn't this work?! # I guess this rage-comment was related to why it doesn't work except for the default values # Anyways, this is looong deprecated

cd levyfit
make clean
make exe/EbE_or_Eavg_1d3d_fit.exe

nevt_avg_default=0 # default value; if <1, using different for each energy, acc. to header 
nevt=10000 # number of events to use for fitting
nevt_highstat=100000

# list of NEVT_AVG values
nevt_avgs=(10 25 50 100 200 500 1000 5000 10000)
nevt_avgs_highstat=(20000 25000 50000 100000)

for avg in "${nevt_avgs[@]}"; do
  echo "Preparing folders for nevt_avg: ${avg}"
  rm -rf $BASEDIR/figs/fitting/lcms/AVG${avg}/
  mkdir -p $BASEDIR/figs/fitting/lcms/AVG${avg}/
done
for avg in "${nevt_avgs_highstat[@]}"; do
  echo "Preparing folders for high-statistics nevt_avg: ${avg}"
  rm -rf $BASEDIR/figs/fitting/lcms/AVG${avg}/
  mkdir -p $BASEDIR/figs/fitting/lcms/AVG${avg}/
done

if [ "$do_fitting" = true ]; then
  echo "Starting fitting for all nevt_avg values..."
  # Fitting for different nevt_avg for systematics
  for avg in "${nevt_avgs[@]}"; do
    echo "Fitting for nevt_avg: ${avg}"
    for energy in "${energies_low[@]}"; do
      echo "Fitting for high-statistics energy ${energy}"
      # Run the fit and then move produced AVG images into their folder only after the fit finishes
      run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${avg} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}.log"
    done
    for energy in "${energies_high[@]}"; do
      echo "Fitting for energy ${energy}"
      run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${avg} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}.log"
    done
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.png" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.pdf" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
  done

  for avg in "${nevt_avgs_highstat[@]}"; do
    echo "Fitting for >10k nevt_avg: ${avg}"
    for energy in "${energies_low[@]}"; do
      echo "Fitting for high-statistics energy ${energy} with nevt_avg: ${avg}"
      run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${avg} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}_highstat.log"
    done
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.png" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*AVG${avg}*.pdf" "$BASEDIR/figs/fitting/lcms/AVG${avg}"
  done

  # Wait for any remaining background jobs from nevt_avg sweep
  wait
fi # end of fitting block (nevt_avg systematics)

cd ..

root.exe -b -q plot_param_vs_nevt_avg.cpp\(-1\)
for ikt in {0..9}; do
  echo "Plotting parameter vs NEVT_AVG graphs for KT bin: ${ikt}"
  # -1 means all KT bins averaged, otherwise give the KT bin index (0..9)
  root.exe -b -q plot_param_vs_nevt_avg.cpp\(${ikt}\)
done


# qLCMS systematics
echo "Preparing folders for qLCMS systematics"
rm -rf $BASEDIR/figs/fitting/lcms/defaultQlcms/
rm -rf $BASEDIR/figs/fitting/lcms/strictQlcms/
rm -rf $BASEDIR/figs/fitting/lcms/looseQlcms/
mkdir -p $BASEDIR/figs/fitting/lcms/defaultQlcms
mkdir -p $BASEDIR/figs/fitting/lcms/strictQlcms
mkdir -p $BASEDIR/figs/fitting/lcms/looseQlcms
if [ "$do_fitting" = true ]; then
  echo "Starting fitting for qLCMS systematics..."
  for energy in "${energies_low[@]}"; do
    echo "Fitting for qLCMS systematics, energy ${energy}"
    # Parallel execution
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 1 0 1" "$BASEDIR/logfiles/fit_log_${energy}_strictqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 2 0 1" "$BASEDIR/logfiles/fit_log_${energy}_looseqLCMS.log"
    # We started up to 3 jobs here; throttle will keep the overall concurrency <= MAXJOBS
    # Wait here to ensure all qLCMS systematics for this energy finish before moving on
    wait
    #mv $BASEDIR/figs/fitting/lcms/*strictqLCMS*.png $BASEDIR/figs/fitting/lcms/defaultQlcms/
    #mv $BASEDIR/figs/fitting/lcms/*looseqLCMS*.png $BASEDIR/figs/fitting/lcms/defaultQlcms/
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictqLCMS*.png" "$BASEDIR/figs/fitting/lcms/strictQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/strictQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looseqLCMS*.png" "$BASEDIR/figs/fitting/lcms/looseQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looseqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/looseQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.png" "$BASEDIR/figs/fitting/lcms/defaultQlcms" # move remaining default
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultQlcms" # move remaining default
    echo "Fitting for qLCMS systematics, energy ${energy} done."
  done
  for energy in "${energies_high[@]}"; do
    echo "Fitting for qLCMS systematics, energy ${energy}"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 1 0 1" "$BASEDIR/logfiles/fit_log_${energy}_strictqLCMS.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 2 0 1" "$BASEDIR/logfiles/fit_log_${energy}_looseqLCMS.log"
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictqLCMS*.png" "$BASEDIR/figs/fitting/lcms/strictQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/strictQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looseqLCMS*.png" "$BASEDIR/figs/fitting/lcms/looseQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looseqLCMS*.pdf" "$BASEDIR/figs/fitting/lcms/looseQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.png" "$BASEDIR/figs/fitting/lcms/defaultQlcms"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultQlcms"
    echo "Fitting for qLCMS systematics, energy ${energy} done."
  done
fi # end of qLCMS systematics fitting block

# rhofitmax systematics
echo "Preparing folders for rhofitmax systematics"
rm -rf $BASEDIR/figs/fitting/lcms/defaultrhoFitMax/
rm -rf $BASEDIR/figs/fitting/lcms/strictrhoFitMax/
rm -rf $BASEDIR/figs/fitting/lcms/looserhoFitMax/
mkdir -p $BASEDIR/figs/fitting/lcms/defaultrhoFitMax
mkdir -p $BASEDIR/figs/fitting/lcms/strictrhoFitMax
mkdir -p $BASEDIR/figs/fitting/lcms/looserhoFitMax
if [ "$do_fitting" = true ]; then
  echo "Starting fitting for rhofitmax systematics..."
  for energy in "${energies_low[@]}"; do
    echo "Fitting for rhofitmax systematics, energy ${energy}"
    # Parallel execution
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 1 1" "$BASEDIR/logfiles/fit_log_${energy}_strictrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt_highstat} ${nevt_avg_default} 0 2 1" "$BASEDIR/logfiles/fit_log_${energy}_looserhoFitMax.log"
    wait
    # Move remaining default
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictrhoFitMax*.png" "$BASEDIR/figs/fitting/lcms/strictrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictrhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/strictrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looserhoFitMax*.png" "$BASEDIR/figs/fitting/lcms/looserhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looserhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/looserhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.png" "$BASEDIR/figs/fitting/lcms/defaultrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultrhoFitMax"
    echo "Fitting for rhofitmax systematics, energy ${energy} done."
  done
  for energy in "${energies_high[@]}"; do
    echo "Fitting for rhofitmax systematics, energy ${energy}"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 0 1 1" "$BASEDIR/logfiles/fit_log_${energy}_strictrhoFitMax.log"
    run_bg "cd \"$BASEDIR/levyfit\" && exe/EbE_or_Eavg_1d3d_fit.exe 11 \"${energy}\" 1 ${nevt} ${nevt_avg_default} 0 2 1" "$BASEDIR/logfiles/fit_log_${energy}_looserhoFitMax.log"
    wait
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictrhoFitMax*.png" "$BASEDIR/figs/fitting/lcms/strictrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*strictrhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/strictrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looserhoFitMax*.png" "$BASEDIR/figs/fitting/lcms/looserhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*looserhoFitMax*.pdf" "$BASEDIR/figs/fitting/lcms/looserhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.png" "$BASEDIR/figs/fitting/lcms/defaultrhoFitMax"
    move_matching_png "$BASEDIR/figs/fitting/lcms/*.pdf" "$BASEDIR/figs/fitting/lcms/defaultrhoFitMax"
    echo "Fitting for rhofitmax systematics, energy ${energy} done."
  done
fi # end of rhofitmax systematics fitting block

# These plots only for informative purposes, with error bars showing stddev or sg like that
# mkdir -p $BASEDIR/alphaNR_vs_kt
# for ienergy in "${energies_low[@]}"; do
#   echo "Plotting individual graphs for energy: ${ienergy}"
#   #root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,${nevt_avg_default}\)
#   root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,5000\)
# done
# for ienergy in "${energies_high[@]}"; do
#   echo "Plotting individual graphs for energy: ${ienergy}"
#   root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,5000\)
# done

# Final summary plots
# -------------------
echo "Plotting param vs sqrt(sNN)"
#root.exe -b -q plot_alphaNR_allcent.cpp\(${nevt_avg_default}\)
# root.exe -b -q plot_alphaNR_allcent.cpp\(5000\) # sort of deprecated with UrQMD
# Calculate & collect all systematics
root.exe -b -q calc_and_plot_syserr.cpp\(-1\) || { echo "D(rho) systematics failed; comparison skipped."; exit 1; }
compare_methods || exit 1
# Plot rhofitmax values vs KT, additionally
#~/.venvs/jupyter/bin/python3.12 rhofitmax_vs_kt_snn.py
# or:
python3 rhofitmax_vs_kt_snn.py

