const double Mass2_pi = 0.019479835;
const double Mass2_ka = 0.24371698032;

const int NKT = 10;
const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675}; // Matyas
const char* energies[] = {"3p0","3p2","3p5","3p9","4p5","7p7","9p2","11p5","14p5","19p6","27"};
const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
const double energydouble[NENERGIES] = {3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0};

//const int NKT = 8;
//const double kT_center[NKT] = {0.0452, 0.0904, 0.1299, 0.1649, 0.2044, 0.2633, 0.3820, 0.6221};
//const double ktbins[NKT + 1] = {0.0, 0.07, 0.11, 0.15, 0.18, 0.23, 0.3, 0.5, 1.0}; // bporfy
// const char* energies[] = {"30A"};
// const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
// const double energydouble[NENERGIES] = {30.0};