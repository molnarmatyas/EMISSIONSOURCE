#!/bin/bash

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

# Default fitting
cd levyfit
make clean
make exe/onedim_EbE_or_Eavg_fit.exe

nevt_avg_default=500
for avg in "${nevt_avgs[@]}"; do
  echo "Preparing folders for nevt_avg: ${avg}"
  rm -rf ../figs/fitting/lcms/AVG${avg}/
  mkdir -p ../figs/fitting/lcms/AVG${avg}/
done

# for ienergy in "${energies[@]}"; do
#   echo "Fitting for energy: ${ienergy}"
#   exe/onedim_EbE_or_Eavg_fit.exe 11 ${ienergy} 1 10000 ${nevt_avg_default} 0 0 1 &> logfiles/fit_log_${ienergy}.log # avg. by 100--10000 events
#   mv ../figs/fitting/lcms/*AVG1000*.png ../figs/fitting/lcms/AVG1000/
# done

# Fitting for different nevt_avg for systematics
nevt_avgs=(10 25 50 100 200 500 1000 5000 10000) # default to be removed
                                                 # also 1 as that does not work and lasts too long

for energy in "${energies[@]}"; do
  echo "Fitting for energy ${energy}"
  for avg in "${nevt_avgs[@]}"; do
    echo "Fitting for nevt_avg: ${avg}"
    exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${avg} 0 0 1 &> logfiles/fit_log_${energy}_nevtavg${avg}.log # avg. by nevt events
    mv ../figs/fitting/lcms/*AVG${avg}*.png ../figs/fitting/lcms/AVG${avg}/
  done
done

cd ..


for ienergy in "${energies[@]}"; do
  echo "Plotting individual graphs for energy: ${ienergy}"
  root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,${nevt_avg_default}\)
done


# TODO: after finding out the best value for nevt_avg_default, do qLCMS systematics too with that nevt_avg
root.exe -b -q plot_param_vs_nevt_avg.cpp\(-1\) # -1 means all KT bins averaged, otherwise give the KT bin index (0..9)

# qLCMS systematics
for energy in "${energies[@]}"; do
  echo "Fitting for qLCMS systematics, energy ${energy}"
  # Parallel execution
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 0 0 1 &> logfiles/fit_log_${energy}_defaultqLCMS.log & # default qLCMS cut
  pid1=$!
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 1 0 1 &> logfiles/fit_log_${energy}_strictqLCMS.log & # strict qLCMS cut
  pid2=$!
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 2 0 1 &> logfiles/fit_log_${energy}_looseqLCMS.log & # loose qLCMS cut
  pid3=$!
  # Wait for all 3 to finish
  wait $pid1
  wait $pid2
  wait $pid3
  echo "Fitting for qLCMS systematics, energy ${energy} done."
done

# rhofitmax systematics
for energy in "${energies[@]}"; do
  echo "Fitting for rhofitmax systematics, energy ${energy}"
  # Parallel execution
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 0 0 1 &> logfiles/fit_log_${energy}_defaultrhoFitMax.log & # default rhofitmax
  pid1=$!
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 0 1 1 &> logfiles/fit_log_${energy}_strictrhoFitMax.log & # strict rhofitmax
  pid2=$!
  exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${nevt_avg_default} 0 2 1 &> logfiles/fit_log_${energy}_looserhoFitMax.log & # loose rhofitmax
  pid3=$!
  # Wait for all 3 to finish
  wait $pid1
  wait $pid2
  wait $pid3
  echo "Fitting for rhofitmax systematics, energy ${energy} done."
done

# TODO: calculate & collect all systematics

# Final summary plots
echo "Plotting param vs sqrt(sNN)"
root.exe -b -q plot_alphaNR_allcent.cpp\(${nevt_avg_default}\)