#!/bin/bash

cp ArSc_output_30.root urqmd-3.4/drho_analyze/AuAu_30A_tree.root
cd urqmd-3.4/drho_analyze
# Run the rho dist creation
make clean
make pairsource_urqmd.exe
./pairsource_urqmd.exe "30A"
mv UrQMD_3d_source_0-10cent_all_30A.root ../../analysed/
rm AuAu_30A_tree.root
cd ../..
# Now do the fitting
cd levyfit
make clean
make exe/onedim_EbE_or_Eavg_fit.exe
exe/onedim_EbE_or_Eavg_fit.exe 11 "30A" 1 10000 10000 1 &> fit_log_30A.log # avg. by 1000 events: set second 10k to 1000
# Plotting param vs mT
cd ..
root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"30A\",true,10000\) # do not forget to set to same NEVT_AVG