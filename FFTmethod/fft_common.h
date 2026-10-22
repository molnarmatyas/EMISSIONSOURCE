#ifndef FFT_COMMON_H
#define FFT_COMMON_H
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include "../header_for_all_emissionsource.h"

inline TString fftOutputName(TString input, const char* suffix)
{
  if(input.EndsWith(".root")) input.Resize(input.Length() - 5);
  if(input.EndsWith("_FFT")) input.Resize(input.Length() - 4);
  return input + suffix + "_FFT.root";
}

inline bool checkFFTOutput(TString name)
{
  if(!name.EndsWith("_FFT.root") || !gSystem->AccessPathName(name)) {
    std::cerr << "Use a new output filename ending in _FFT.root: " << name << std::endl;
    return false;
  }
  return true;
}

inline void copyFFTInfo(TFile* input, TFile* output)
{
  output->cd();
  const char* names[] = {"source_events", "pair_source_schema", "source_charge_mode",
    "source_qlcms_cut", "source_energy", "source_centrality", "source_kt_min", "source_kt_max",
    "fft_input", "fft_convention", "fft_cutoff", "fft_bin_error", "fft_hbarc"};
  for(const char* name : names) {
    TObject* object = input->Get(name);
    if(object) object->Write(name);
  }
}
#endif
