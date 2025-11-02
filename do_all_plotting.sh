#!/bin/bash
root.exe -b -q collect_alphavsR.cpp

root.exe -b -q  Print_nphist.cpp\(\"7.7GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"9.2GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"11.5GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"14.5GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"19.6GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"27GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"39GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"62.4GeV\"\)
root.exe -b -q  Print_nphist.cpp\(\"200GeV\"\)

root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"7.7GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"9.2GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"11.5GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"14.5GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"19.6GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"27GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"39GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"62.4GeV\"\)
root.exe -b -q plot_alpha_vs_kt_EbE.cpp\(\"200GeV\"\)

root.exe -b -q allenergies_paramgraph_EbE.cpp

root.exe -b -q averaged_paramgraph_EbE.cpp

root.exe -b -q plot_alphaNR_allcent.cpp
