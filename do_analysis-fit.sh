#!/bin/bash

# time ./do_analysis-fit.sh &> fullanalysis_fullfit.log &

#centleg=("0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10")
centleg=(0 1 2 3 4 5 6 7 8 9 10 11)
energies=("7.7GeV" "9.2GeV" "11.5GeV" "14.5GeV" "19.6GeV" "27GeV" "39GeV" "62.4GeV" "130GeV" "200GeV")

# Compile fitting code
cd levyfit
make clean
make exe/onedim_EbE_fit.exe
cd ..

# ANALYSIS AND FITTING
for ienergy in "${energies[@]}"; do
  for cent in "${centleg[@]}"; do
    # Analysis for centrality classes
    root.exe -b -q lcmsonly_epos_pair_source_3d.C\(40,100,\"${ienergy}\",${cent}\)
    # Fitting
    cd levyfit
    exe/onedim_EbE_fit.exe ${cent} "${ienergy}" 40 100
    cd ..
  done
done

