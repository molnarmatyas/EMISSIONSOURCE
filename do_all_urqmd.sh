#!/bin/bash

# cd urqmd-3.4/drho_analyze
# ./pairsource_fit.exe 3p0 && ./pairsource_fit.exe 3p2 && ./pairsource_fit.exe 3p5 && ./pairsource_fit.exe 3p9 && ./pairsource_fit.exe 4p5 && ./pairsource_fit.exe 7p7
# cd ../..

#root.exe -b -q 'Average_Drho.cpp(100)' # why doesn't this work?! # I guess this rage-comment was related to why it doesn't work except for the default values

cd levyfit
make clean
make exe/onedim_EbE_or_Eavg_fit.exe

energies=("3p0" "3p2" "3p5" "3p9" "4p5" "7p7" "9p2" "11p5" "14p5" "19p6" "27")

for ienergy in "${energies[@]}"; do
  echo "Fitting for energy: ${ienergy}"
  exe/onedim_EbE_or_Eavg_fit.exe 11 ${ienergy} 1 10000 100 1 &> fit_log_${ienergy}.log # avg. by 100--10000 events
done

cd ..


# for ienergy in "${energies[@]}"; do
#   echo "Plotting individual graphs for energy: ${ienergy}"
#   root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"${ienergy}\",true,1000\)
# done

# root.exe -b -q plot_alphaNR_allcent.cpp
