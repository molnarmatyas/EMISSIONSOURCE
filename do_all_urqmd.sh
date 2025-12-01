#!/bin/bash

# Base directory of this script (absolute)
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Sanity check
echo "Base directory: $BASEDIR"
sleep 10

energies=("3p0" "3p2" "3p5" "3p9" "4p5" "7p7" "9p2" "11p5" "14p5" "19p6" "27")

cd urqmd-3.4/drho_analyze
make clean
make pairsource_urqmd.exe

analysedname="UrQMD_3d_source_0-10cent_all_"
for ienergy in "${energies[@]}"; do
  echo "Creating source for energy: ${ienergy}"
  ./pairsource_urqmd.exe "${ienergy}" 0 # 0 = default qLCMS cut
  # For systematics:
  mv ${analysedname}${ienergy}.root ${analysedname}${ienergy}_defaultqLCMS.root
  ./pairsource_urqmd.exe "${ienergy}" 1 # 1 = strict qLCMS cut
  mv ${analysedname}${ienergy}.root ${analysedname}${ienergy}_strictqLCMS.root
  ./pairsource_urqmd.exe "${ienergy}" 2 # 2 = loose qLCMS cut
  mv ${analysedname}${ienergy}.root ${analysedname}${ienergy}_looseqLCMS.root
  mv ${analysedname}${ienergy}_defaultqLCMS.root ${analysedname}${ienergy}.root # restore default-named file
done
mv ${analysedname}*.root ../../analysed/
cd ../..

#root.exe -b -q 'Average_Drho.cpp(100)' # why doesn't this work?! # I guess this rage-comment was related to why it doesn't work except for the default values

cd levyfit
make clean
make exe/onedim_EbE_or_Eavg_fit.exe

# Enable job control so `jobs` builtin works inside this script
set -m

# Ensure log directory exists
mkdir -p $BASEDIR/logfiles

# Ensure background children are killed if the script exits or is interrupted
trap 'echo "Script exiting — terminating background jobs..."; jobs -p | xargs -r kill; wait' EXIT SIGINT SIGTERM

# Simple job limiter: keep at most MAXJOBS background processes in this script
# Usage: ./do_all_urqmd.sh [MAXJOBS]  (optional). Default is 5 to be safe with memory.
if [ -n "$1" ]; then
  MAXJOBS="$1"
else
  MAXJOBS=5
fi
echo "MAXJOBS set to $MAXJOBS"
limit_jobs() {
  # Wait while number of background jobs is >= MAXJOBS
  while [ "$(jobs -p | wc -l)" -ge "$MAXJOBS" ]; do
    sleep 2
  done
}

# Helper to run a command in background and redirect stdout/stderr to logfile
run_bg() {
  local cmd="$1"
  local logfile="$2"
  eval "$cmd" &> "$logfile" &
  # throttle
  limit_jobs
}

nevt_avg_default=500

# list of NEVT_AVG values
nevt_avgs=(10 25 50 100 200 500 1000 5000 10000)

for avg in "${nevt_avgs[@]}"; do
  echo "Preparing folders for nevt_avg: ${avg}"
  rm -rf $BASEDIR/figs/fitting/lcms/AVG${avg}/
  mkdir -p $BASEDIR/figs/fitting/lcms/AVG${avg}/
done

# Default fitting
# for ienergy in "${energies[@]}"; do
#   echo "Fitting for energy: ${ienergy}"
#   exe/onedim_EbE_or_Eavg_fit.exe 11 ${ienergy} 1 10000 ${nevt_avg_default} 0 0 1 &> logfiles/fit_log_${ienergy}.log # avg. by 100--10000 events
#   mv ../figs/fitting/lcms/*AVG1000*.png ../figs/fitting/lcms/AVG1000/
# done

# Fitting for different nevt_avg for systematics
for avg in "${nevt_avgs[@]}"; do
  echo "Fitting for nevt_avg: ${avg}"
  for energy in "${energies[@]}"; do
    echo "Fitting for energy ${energy}"
    # Run the fit and then move produced AVG images into their folder only after the fit finishes
    run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${avg} 0 0 1 && mv $BASEDIR/figs/fitting/lcms/*AVG${avg}*.png $BASEDIR/figs/fitting/lcms/AVG${avg}/" "$BASEDIR/logfiles/fit_log_${energy}_nevtavg${avg}.log"
  done
done

# Wait for any remaining background jobs from nevt_avg sweep
wait

cd ..

root.exe -b -q plot_param_vs_nevt_avg.cpp\(-1\)
for ikt in {0..9}; do
  echo "Plotting individual graphs for KT bin: ${ikt}"
  # -1 means all KT bins averaged, otherwise give the KT bin index (0..9)
  root.exe -b -q plot_param_vs_nevt_avg.cpp\(${ikt}\)
done

# These plots only for informative purposes, with error bars showing stddev or sg like that
for ienergy in "${energies[@]}"; do
  echo "Plotting individual graphs for energy: ${ienergy}"
  root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,${nevt_avg_default}\)
done

# qLCMS systematics
for energy in "${energies[@]}"; do
  echo "Fitting for qLCMS systematics, energy ${energy}"
  # Parallel execution
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 0 0 1" "$BASEDIR/logfiles/fit_log_${energy}_defaultqLCMS.log"
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 1 0 1" "$BASEDIR/logfiles/fit_log_${energy}_strictqLCMS.log"
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 2 0 1" "$BASEDIR/logfiles/fit_log_${energy}_looseqLCMS.log"
  # We started up to 3 jobs here; throttle will keep the overall concurrency <= MAXJOBS
  # Optionally wait here to ensure all qLCMS systematics for this energy finish before moving on
  wait
  echo "Fitting for qLCMS systematics, energy ${energy} done."
done

# rhofitmax systematics
for energy in "${energies[@]}"; do
  echo "Fitting for rhofitmax systematics, energy ${energy}"
  # Parallel execution
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 0 0 1" "logfiles/fit_log_${energy}_defaultrhoFitMax.log"
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 0 1 1" "logfiles/fit_log_${energy}_strictrhoFitMax.log"
  run_bg "cd \"$BASEDIR/levyfit\" && exe/onedim_EbE_or_Eavg_fit.exe 11 \"${energy}\" 1 10000 ${nevt_avg_default} 0 2 1" "logfiles/fit_log_${energy}_looserhoFitMax.log"
  wait
  echo "Fitting for rhofitmax systematics, energy ${energy} done."
done

# TODO: calculate & collect all systematics

# Final summary plots
echo "Plotting param vs sqrt(sNN)"
root.exe -b -q plot_alphaNR_allcent.cpp\(${nevt_avg_default}\)

root.exe -b -q calc_and_plot_syserr.cpp