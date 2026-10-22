// Build: make -C FFTmethod FFT=1 transform_fft.exe
#include "fft_common.h"
#include "../pair_source_tree.h"
#include <TH3.h>
#include <TKey.h>
#include <TClass.h>
#include <TNamed.h>
#include <TParameter.h>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
using namespace std;

int main(int argc, char** argv)
{
  if(argc < 2 || string(argv[1]) == "--help") {
    cout << "Usage: transform_fft.exe input.root [cutoff] [options]\n"
         << "  --bins N --range R    source grid (defaults in global header)\n"
         << "  --kt K --charge both|negative|positive\n"
         << "  --first-event E --events N    event interval [E,E+N)\n"
         << "  --cutoff X --error X --save-source --output file_FFT.root\n";
    return argc < 2 ? 1 : 0;
  }
  if(!useFFTmethod) { cerr << "FFT is disabled in the global header/build.\n"; return 1; }

  // Settings: all events, both charges, all kT bins by default.
  TString inputName = argv[1], outputName = fftOutputName(inputName, "_Cq");
  int nbins = fftSourceBins, onlyKT = -1, charge = -1;
  Long64_t firstEvent = 0, events = -1;
  double range = fftSourceHalfRange, cutoff = 0., binError = 1.;
  bool saveSource = false;
  for(int i = 2; i < argc; ++i) {
    string option = argv[i];
    if(option == "--save-source") { saveSource = true; continue; }
    if(i == 2 && option.find("--") != 0) { cutoff = atof(argv[i]); continue; }
    if(i + 1 == argc) { cerr << "Missing value for " << option << endl; return 1; }
    string value = argv[++i];
    if(option == "--output") outputName = value;
    else if(option == "--bins") nbins = atoi(value.c_str());
    else if(option == "--range") range = atof(value.c_str());
    else if(option == "--cutoff") cutoff = atof(value.c_str());
    else if(option == "--error") binError = atof(value.c_str());
    else if(option == "--kt") onlyKT = atoi(value.c_str());
    else if(option == "--first-event") firstEvent = atoll(value.c_str());
    else if(option == "--events") events = atoll(value.c_str());
    else if(option == "--charge") {
      if(value == "both") charge = -1;
      else if(value == "negative") charge = 0;
      else if(value == "positive") charge = 1;
      else { cerr << "Unknown charge: " << value << endl; return 1; }
    }
    else { cerr << "Unknown option: " << option << endl; return 1; }
  }
  if(nbins < 3 || nbins > 1287 || nbins % 2 == 0 || range <= 0 || !isfinite(range) ||
     cutoff < 0 || !isfinite(cutoff) || binError <= 0 || !isfinite(binError) ||
     onlyKT < -1 || onlyKT >= NKT || firstEvent < 0 || events == 0 || events < -1) {
    cerr << "Check grid size (odd), range/error (positive), cutoff and pair selection.\n"; return 1;
  }
  if(!checkFFTOutput(outputName)) return 1;
  TFile* input = TFile::Open(inputName, "READ");
  if(!input || input->IsZombie()) { cerr << "Cannot open " << inputName << endl; return 1; }
  TTree* tree = dynamic_cast<TTree*>(input->Get("pairTree"));
  PairSourceData pair;
  vector<TString> names;
  vector<int> ktBins;

  // Find populated kT bins, or take the TH3 source names from an older file.
  if(tree) {
    if(!setPairSourceBranches(tree, pair, true)) return 1;
    Long64_t counts[NKT] = {};
    for(Long64_t i = 0; i < tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if(pair.ievent < firstEvent || (events >= 0 && pair.ievent - firstEvent >= events)) continue;
      if(charge >= 0 && pair.ich != charge) continue;
      if(pair.iKT >= 0 && pair.iKT < NKT) ++counts[pair.iKT];
    }
    for(int kt = 0; kt < NKT; ++kt) {
      if(!counts[kt] || (onlyKT >= 0 && kt != onlyKT)) continue;
      names.push_back(Form("D_lcms_3d_kt%d", kt)); ktBins.push_back(kt);
    }
  } else {
    if(onlyKT >= 0 || charge >= 0 || firstEvent || events >= 0) {
      cerr << "Pair selections require pairTree input.\n"; return 1;
    }
    TIter next(input->GetListOfKeys());
    TKey* key;
    while((key = (TKey*)next())) {
      TClass* type = TClass::GetClass(key->GetClassName());
      if(type && type->InheritsFrom(TH3::Class()) && find(names.begin(), names.end(), key->GetName()) == names.end()) {
        names.push_back(key->GetName()); ktBins.push_back(-1);
      }
    }
  }
  if(names.empty()) { cerr << "No sources in the requested selection.\n"; return 1; }

  // Small table describing the source normalization and selection of each C(q).
  string sourceName, correlationName;
  int kt, bins;
  double total, inGrid, retained;
  TTree summary("fft_sources", "FFT source normalization and selection");
  summary.SetDirectory(nullptr);
  summary.Branch("source", &sourceName); summary.Branch("correlation", &correlationName);
  summary.Branch("kt", &kt); summary.Branch("bins", &bins); summary.Branch("charge", &charge);
  summary.Branch("first_event", &firstEvent); summary.Branch("events", &events);
  summary.Branch("total_with_flow", &total); summary.Branch("in_grid", &inGrid);
  summary.Branch("retained_fraction", &retained);
  TFile* output = nullptr;
  int N = 0, written = 0;
  fftw_complex *fftIn = nullptr, *fftOut = nullptr;
  fftw_plan plan = nullptr;

  for(size_t ihist = 0; ihist < names.size(); ++ihist) {
    kt = ktBins[ihist];
    TH3* hDrho;
    if(tree) {
      hDrho = new TH3D(names[ihist], "Pair counts;#rho_{out} [fm];#rho_{side} [fm];#rho_{long} [fm]",
                       nbins,-range,range, nbins,-range,range, nbins,-range,range);
      hDrho->SetDirectory(nullptr);
      for(Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if(pair.iKT != kt || (charge >= 0 && pair.ich != charge)) continue;
        if(pair.ievent < firstEvent || (events >= 0 && pair.ievent - firstEvent >= events)) continue;
        double rhoout = pair.rhoout, rhoside = pair.rhoside, rholong = pair.rholong;
        // Restore the signs stored by pairsource_urqmd (out=1, side=2, long=4).
        if(pair.signs & 1) rhoout = -rhoout;
        if(pair.signs & 2) rhoside = -rhoside;
        if(pair.signs & 4) rholong = -rholong;
        if(pair.signs > 7 || !isfinite(rhoout) || !isfinite(rhoside) || !isfinite(rholong)) {
          cerr << "Invalid separation/signs in pairTree.\n"; return 1;
        }
        hDrho->Fill(rhoout, rhoside, rholong);
      }
    } else hDrho = (TH3*)input->Get(names[ihist]);
    hDrho->SetDirectory(nullptr);

    // FFT needs equal odd bin counts, uniform axes and nonnegative source counts.
    bins = hDrho->GetNbinsX();
    if(bins < 3 || bins % 2 == 0 || hDrho->GetNbinsY() != bins || hDrho->GetNbinsZ() != bins) {
      cerr << "Use equal, odd bin counts in all three directions.\n"; return 1;
    }
    TAxis* axes[] = {hDrho->GetXaxis(), hDrho->GetYaxis(), hDrho->GetZaxis()};
    for(TAxis* axis : axes)
      for(int b = 1; b <= bins; ++b)
        if(!isfinite(axis->GetBinWidth(b)) || axis->GetBinWidth(b) <= 0 ||
           abs(axis->GetBinWidth(b) - axis->GetBinWidth(1)) > 1.e-9 * axis->GetBinWidth(1)) {
          cerr << "FFT needs uniformly spaced spatial bins.\n"; return 1;
        }
    for(int b = 0; b < hDrho->GetNcells(); ++b)
      if(!isfinite(hDrho->GetBinContent(b)) || hDrho->GetBinContent(b) < 0) {
        cerr << "Invalid source count.\n"; return 1;
      }
    total = hDrho->Integral(0,bins+1, 0,bins+1, 0,bins+1); // include ALL flow bins!
    inGrid = hDrho->Integral(1,bins, 1,bins, 1,bins);
    if(total == 0) { delete hDrho; continue; }
    if(!isfinite(total) || inGrid == 0) { cerr << "No usable counts inside the source grid.\n"; return 1; }
    retained = inGrid / total;

    // Allocate once; only rebuild the plan if the next histogram has another size.
    if(N != bins) {
      if(plan) { fftw_destroy_plan(plan); fftw_free(fftIn); fftw_free(fftOut); }
      N = bins;
      size_t cells = (size_t)N * N * N;
      fftIn = fftw_alloc_complex(cells); fftOut = fftw_alloc_complex(cells);
      if(!fftIn || !fftOut) { cerr << "Not enough memory for FFT.\n"; return 1; }
      plan = fftw_plan_dft_3d(N,N,N, fftIn,fftOut, FFTW_BACKWARD,FFTW_ESTIMATE);
      if(!plan) { cerr << "Cannot create FFTW plan.\n"; return 1; }
    }
    for(int x = 0; x < N; ++x) for(int y = 0; y < N; ++y) for(int z = 0; z < N; ++z) {
      size_t index = ((size_t)x * N + y) * N + z;
      fftIn[index][0] = hDrho->GetBinContent(x+1,y+1,z+1) / total;
      fftIn[index][1] = 0.;
    }
    fftw_execute(plan);

    double dq[3], origin[3];
    for(int axis = 0; axis < 3; ++axis) {
      dq[axis] = 2.*acos(-1.)*fftHbarC / (N*axes[axis]->GetBinWidth(1));
      origin[axis] = axes[axis]->GetBinCenter(1);
    }
    sourceName = hDrho->GetName(); correlationName = "Cq_" + sourceName;
    TH3D* hCq = new TH3D(correlationName.c_str(), "LCMS C(q);q_{out} [GeV/c];q_{side} [GeV/c];q_{long} [GeV/c]",
      N,-N*dq[0]/2,N*dq[0]/2, N,-N*dq[1]/2,N*dq[1]/2, N,-N*dq[2]/2,N*dq[2]/2);
    hCq->SetDirectory(nullptr);
    for(int x = 0; x < N; ++x) for(int y = 0; y < N; ++y) for(int z = 0; z < N; ++z) {
      int mx = x <= N/2 ? x : x-N, my = y <= N/2 ? y : y-N, mz = z <= N/2 ? z : z-N;
      size_t index = ((size_t)x * N + y) * N + z;
      double phase = (mx*dq[0]*origin[0] + my*dq[1]*origin[1] + mz*dq[2]*origin[2]) / fftHbarC;
      double real = fftOut[index][0]*cos(phase) - fftOut[index][1]*sin(phase);
      if(abs(real) < cutoff) real = 0.;
      // D is already a pair-separation distribution: use Re(FFT), not |FFT|^2.
      hCq->SetBinContent(mx+N/2+1, my+N/2+1, mz+N/2+1, 1.+real);
      hCq->SetBinError(mx+N/2+1, my+N/2+1, mz+N/2+1, binError);
    }
    if(!output) {
      output = TFile::Open(outputName, "CREATE");
      if(!output || output->IsZombie()) { cerr << "Cannot create " << outputName << endl; return 1; }
      TParameter<int>("fft_complete", 0).Write();
    }
    output->cd(); hCq->Write();
    if(saveSource) hDrho->Write();
    summary.Fill(); ++written;
    cout << sourceName << ": pairs=" << total << ", retained=" << retained << ", C(0)=" << 1.+retained << endl;
    delete hCq; delete hDrho;
  }
  if(plan) { fftw_destroy_plan(plan); fftw_free(fftIn); fftw_free(fftOut); }
  if(tree) tree->ResetBranchAddresses();
  if(!written) { cerr << "No nonempty sources to transform.\n"; return 1; }
  copyFFTInfo(input, output);
  TNamed("fft_input", inputName.Data()).Write();
  TNamed("fft_convention", "C(q)=1+Re FFT(D); normalize source counts including flow bins").Write();
  TParameter<double>("fft_cutoff", cutoff).Write();
  TParameter<double>("fft_bin_error", binError).Write();
  TParameter<double>("fft_hbarc", fftHbarC).Write();
  summary.Write(); TParameter<int>("fft_complete", 1).Write("fft_complete", TObject::kOverwrite);
  output->Close(); input->Close(); delete output; delete input;
  cout << "Wrote " << written << " correlations to " << outputName << endl;
  return 0;
}
