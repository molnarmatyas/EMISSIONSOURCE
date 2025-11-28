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

for ienergy in "${energies[@]}"; do
  echo "Fitting for energy: ${ienergy}"
  exe/onedim_EbE_or_Eavg_fit.exe 11 ${ienergy} 1 10000 1000 1 &> fit_log_${ienergy}.log # avg. by 100--10000 events
  rm -rf ../figs/fitting/lcms/AVG1000/
  mkdir -p ../figs/fitting/lcms/AVG1000/
  mv ../figs/fitting/lcms/*AVG1000*.png ../figs/fitting/lcms/AVG1000/
done

# Fitting for different nevt_avg for systematics
nevt_avgs=(1 10 25 50 100 200 500 5000 10000) # removed 1000 as that is the default

for energy in "${energies[@]}"; do
  echo "Fitting for energy ${energy}"
  for avg in "${nevt_avgs[@]}"; do
    echo "Fitting for nevt_avg: ${avg}"
    exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${avg} 1 &> fit_log_${energy}_nevtavg${avg}.log # avg. by nevt events
    rm -rf ../figs/fitting/lcms/AVG${avg}/
    mkdir -p ../figs/fitting/lcms/AVG${avg}/
    mv ../figs/fitting/lcms/*AVG${avg}*.png ../figs/fitting/lcms/AVG${avg}/
  done
done

cd ..


# for ienergy in "${energies[@]}"; do
#   echo "Plotting individual graphs for energy: ${ienergy}"
#   root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,1000\)
# done

# root.exe -b -q plot_alphaNR_allcent.cpp
