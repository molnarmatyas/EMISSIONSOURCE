#ifndef HEADER_FOR_ALL_EMISSIONSOURCE_H
#define HEADER_FOR_ALL_EMISSIONSOURCE_H

//#include <array>    // std::array
//#include <cstddef>  // std::size_t

// Mass^2 values
const double Mass2_pi = 0.019479835;
const double Mass2_ka = 0.24371698032;

// kT bins
const int NKT = 10;
const double ktbins[NKT+1] = {0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675};
const double kT_center[NKT] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65};

// Energies
const char* energies[] = {"3p0","3p2","3p5","3p9","4p5","7p7","9p2","11p5","14p5","19p6","27"};
const int NENERGIES = static_cast<int>(sizeof(energies) / sizeof(energies[0]));
const double energydouble[NENERGIES] = {3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0};

// from pairsource_urqmd.cc (avoid leading underscore in global scope)
const char* qLCMS_cut[3] = {"default", "strict", "loose"};
const double qLCMS_cut_values[3] = {0.15, 0.05, 0.25}; // GeV/c

// from onedim_EbE_or_Eavg_fit.cc
const double B[3] = {2500.0, 1600.0, 3600.0}; // rho_fitmax limits: default, strict, loose
const double rfitmax_systlimits[3] = {100.0, 50.0, 150.0}; // simpler limits

// bporfy kT bins and centers
// const int NKT = 8;
// const double kT_center[NKT] = {0.0452, 0.0904, 0.1299, 0.1649, 0.2044, 0.2633, 0.3820, 0.6221};
// const double ktbins[NKT + 1] = {0.0, 0.07, 0.11, 0.15, 0.18, 0.23, 0.3, 0.5, 1.0}; // bporfy
// const char* energies[] = {"30A"};
// const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
// const double energydouble[NENERGIES] = {30.0};

#endif // HEADER_FOR_ALL_EMISSIONSOURCE_H
