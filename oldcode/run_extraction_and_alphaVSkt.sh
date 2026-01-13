#!/bin/bash

./extract_params.sh 7.7GeV
root.exe -b -q 'plot_alpha_vs_kt.cpp("7.7GeV")'

./extract_params.sh 9.2GeV
root.exe -b -q 'plot_alpha_vs_kt.cpp("9.2GeV")'

./extract_params.sh 19.6GeV
root.exe -b -q 'plot_alpha_vs_kt.cpp("19.6GeV")'

./extract_params.sh 200GeV
root.exe -b -q 'plot_alpha_vs_kt.cpp("200GeV")'