#ifndef HEADER_FOR_ALL_EMISSIONSOURCE_H
#define HEADER_FOR_ALL_EMISSIONSOURCE_H

#include <array>    // std::array
#include <cstddef>  // std::size_t

// Mass^2 values
inline constexpr double Mass2_pi = 0.019479835;
inline constexpr double Mass2_ka = 0.24371698032;

// kT bins
inline constexpr int NKT = 10;
inline constexpr std::array<double, NKT + 1> ktbins = {
    0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675
};

// Centers computed at compile time via lambda
inline constexpr std::array<double, NKT> kT_center =  {
    std::array<double, NKT> centers{};
    for (std::size_t ii = 0; ii < static_cast<std::size_t>(NKT); ++ii) {
        centers[ii] = 0.5 * (ktbins[ii] + ktbins[ii + 1]);
    }
    return centers;
}();

// Energies
inline constexpr std::array<const char*, 11> energies = {
    "3p0","3p2","3p5","3p9","4p5","7p7","9p2","11p5","14p5","19p6","27"
};
inline constexpr int NENERGIES = static_cast<int>(energies.size());
inline constexpr std::array<double, NENERGIES> energydouble = {
    3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0
};

// from pairsource_urqmd.cc (avoid leading underscore in global scope)
inline constexpr std::array<const char*, 3> qLCMS_cut = {"default", "strict", "loose"};
inline constexpr std::array<double, 3> qLCMS_cut_values = {0.15, 0.05, 0.25}; // GeV/c

// from onedim_EbE_or_Eavg_fit.cc
inline constexpr std::array<double, 3> B = {2500.0, 1600.0, 3600.0}; // rho_fitmax limits: default, strict, loose
inline constexpr std::array<double, 3> rfitmax_systlimits = {100.0, 50.0, 150.0}; // simpler limits


// bporfy kT bins and centers
// const int NKT = 8;
// const double kT_center[NKT] = {0.0452, 0.0904, 0.1299, 0.1649, 0.2044, 0.2633, 0.3820, 0.6221};
// const double ktbins[NKT + 1] = {0.0, 0.07, 0.11, 0.15, 0.18, 0.23, 0.3, 0.5, 1.0}; // bporfy
// const char* energies[] = {"30A"};
// const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
// const double energydouble[NENERGIES] = {30.0};

#endif // HEADER_FOR_ALL_EMISSIONSOURCE_H
