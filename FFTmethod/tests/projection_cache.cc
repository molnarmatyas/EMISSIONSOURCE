// Exercise the existing fitter's actual cache/loading path without loading
// physics tables or running its minimizer/plotting main program.
#define main drho_cli_main
#include "../../levyfit/sources/EbE_or_Eavg_1d3d_fit.cc"
#undef main
#include <stdexcept>

void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
int main(int argc, char** argv) {
  try {
    for(int f = 1; f < argc; ++f) {
      TFile file(argv[f], "READ");
      check(!file.IsZombie(), "Cannot open projection fixture");
      TH1* loaded[NOSL + 1]{};
      loadHistograms(false, 1, 0, 0, 0, &file, true, loaded);
      TH1F expected("expected", "", g_Nbins, g_logbins);
      expected.SetDirectory(nullptr);
      for(int axis = 0; axis <= NOSL; ++axis) {
        expected.Reset();
        if(axis == 0) { expected.Fill(std::sqrt(6.f)); expected.Fill(std::sqrt(6.f)); }
        if(axis == 1) { expected.Fill(1.); expected.Fill(1.); }
        if(axis == 2) { expected.Fill(2.); expected.Fill(1.); }
        if(axis == 3) { expected.Fill(1.); expected.Fill(2.); }
        check(loaded[axis] && loaded[axis]->GetEntries() == 2, "Projection counts changed");
        for(int bin = 0; bin <= g_Nbins + 1; ++bin)
          check(loaded[axis]->GetBinContent(bin) == expected.GetBinContent(bin), "Projection bin content changed");
      }
      loadHistograms(true, 1, 2, 0, 0, &file, true, loaded);
      check(loaded[0]->GetEntries() == 3, "Averaging did not add both events");
      loadHistograms(false, 1, 0, 0, 0, &file, true, loaded);
      check(loaded[0]->GetEntries() == 2, "Clone modification corrupted cache");
      loadHistograms(false, 1, 1, 0, 0, &file, true, loaded);
      for(auto* h : loaded) check(h == nullptr, "Empty first event retained previous averaging block");
      loadHistograms(false, 1, 2, 0, 0, &file, true, loaded);
      check(loaded[0]->GetEntries() == 1, "New averaging block reused old pairs");
      for(auto* h : loaded) delete h;
      clearHistCache();
      g_pairTree = nullptr; g_cachedFile = nullptr;
    }
    std::cout << "PASS existing D(rho) cache, projection bins, averaging and empty-block regression\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
