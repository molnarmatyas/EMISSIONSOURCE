// Build: make -C FFTmethod FFT=1 CqFit.exe
#include "fft_common.h"
#include <TH3.h>
#include <TF3.h>
#include <TKey.h>
#include <TClass.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const int NPARS = 5; // lambda, Rout, Rside, Rlong, alpha
const double fitWeight = 0.01; // fixed numerical weight for a stable Minuit Hessian
vector<double> qout2, qside2, qlong2, cqValues; // selected bins, q^2/hbarc^2

double AltChi2(const double* par)
{
  double chi2 = 0.;
  for(size_t i = 0; i < cqValues.size(); ++i) {
    double qR2 = qout2[i]*par[1]*par[1] + qside2[i]*par[2]*par[2] + qlong2[i]*par[3]*par[3];
    double theory = 1. + par[0]*exp(-pow(qR2, par[4]/2.));
    double residual = (cqValues[i] - theory) / fitWeight;
    chi2 += residual*residual;
  }
  return chi2;
}

int main(int argc, char** argv)
{
  if(argc < 2 || string(argv[1]) == "--help") {
    cout << "Usage: CqFit.exe input_Cq_FFT.root [--qmax Q] [--hist NAME] [--output file_FFT.root]\n"
         << "Default: per-kT spherical |q| limits fftQmaxFit in header_for_all_emissionsource.h.\n"
         << "--qmax overrides these with one limit in GeV/c, or NKT comma-separated limits.\n";
    return argc < 2 ? 1 : 0;
  }
  if(!useFFTmethod) { cerr << "FFT is disabled in the global header/build.\n"; return 1; }
  TString inputName = argv[1], outputName = fftOutputName(inputName, "_fitresults");
  TString onlyHist = "";
  vector<double> qmaxValues(fftQmaxFit, fftQmaxFit + NKT);
  for(int i = 2; i < argc; ++i) {
    string option = argv[i];
    if(i + 1 == argc) { cerr << "Missing value for " << option << endl; return 1; }
    string value = argv[++i];
    if(option == "--output") outputName = value;
    else if(option == "--hist") onlyHist = value;
    else if(option == "--qmax") {
      qmaxValues.clear();
      istringstream values(value);
      string entry;
      while(getline(values, entry, ',')) qmaxValues.push_back(atof(entry.c_str()));
    } else { cerr << "Unknown option: " << option << endl; return 1; }
  }
  if(qmaxValues.size() != 1 && qmaxValues.size() != NKT) {
    cerr << "Supply --qmax with one limit or NKT comma-separated limits.\n"; return 1;
  }
  for(double q : qmaxValues) if(q <= 0 || !isfinite(q)) { cerr << "Invalid qmax.\n"; return 1; }
  if(!checkFFTOutput(outputName)) return 1;
  TFile* input = TFile::Open(inputName, "READ");
  if(!input || input->IsZombie()) { cerr << "Cannot open " << inputName << endl; return 1; }
  TParameter<int>* complete = (TParameter<int>*)input->Get("fft_complete");
  if(complete && complete->GetVal() != 1) { cerr << "FFT input is incomplete.\n"; return 1; }

  // Read which kT bin belongs to each C(q) (needed for kT-dependent qmax).
  map<string,int> ktByName;
  TTree* sources = (TTree*)input->Get("fft_sources");
  if(sources) {
    string* name = nullptr;
    int kt;
    sources->SetBranchAddress("correlation", &name);
    sources->SetBranchAddress("kt", &kt);
    for(Long64_t i = 0; i < sources->GetEntries(); ++i) {
      sources->GetEntry(i); ktByName[*name] = kt;
    }
    sources->ResetBranchAddresses();
    delete name;
  }

  // Result branches. Arrays use the parameter order given above.
  string histogram;
  int kt, status, covarianceStatus;
  Long64_t npoints, NDF;
  bool valid, atLimit;
  double qmax, chi2, chi2ndf, scale, binError, edm;
  double par[NPARS], errorsRaw[NPARS], errors[NPARS], covRaw[25], cov[25];
  TTree results("fit_results_FFT", "3D Levy fits of C(q)");
  results.SetDirectory(nullptr);
  results.Branch("histogram", &histogram); results.Branch("kt", &kt);
  results.Branch("status", &status); results.Branch("covariance_status", &covarianceStatus);
  results.Branch("valid", &valid); results.Branch("at_parameter_limit", &atLimit);
  results.Branch("npoints", &npoints); results.Branch("ndf", &NDF); results.Branch("qmax", &qmax);
  results.Branch("chi2", &chi2); results.Branch("chi2_ndf", &chi2ndf);
  results.Branch("error_scale", &scale); results.Branch("bin_error_raw", &binError); results.Branch("edm", &edm);
  results.Branch("parameters", par, "parameters[5]/D");
  results.Branch("errors_raw", errorsRaw, "errors_raw[5]/D"); results.Branch("errors", errors, "errors[5]/D");
  results.Branch("covariance_raw", covRaw, "covariance_raw[25]/D"); results.Branch("covariance", cov, "covariance[25]/D");
  TFile* output = nullptr;
  vector<string> fitted;
  int failed = 0;
  TIter next(input->GetListOfKeys());
  TKey* key;
  while((key = (TKey*)next())) {
    TString name = key->GetName();
    TClass* type = TClass::GetClass(key->GetClassName());
    if(!name.BeginsWith("Cq_") || !type || !type->InheritsFrom(TH3::Class())) continue;
    if(onlyHist != "" && name != onlyHist) continue;
    histogram = name.Data();
    if(find(fitted.begin(), fitted.end(), histogram) != fitted.end()) continue;
    kt = ktByName.count(histogram) ? ktByName[histogram] : -1;
    if(qmaxValues.size() > 1 && (kt < 0 || kt >= NKT)) {
      cerr << "Per-kT qmax needs kT information in fft_sources; use --qmax Q for unlabelled TH3 inputs.\n"; return 1;
    }
    qmax = qmaxValues.size() == 1 ? qmaxValues[0] : qmaxValues[kt];
    TH3* hCq = (TH3*)input->Get(name);
    hCq->SetDirectory(nullptr);
    TAxis* axes[] = {hCq->GetXaxis(), hCq->GetYaxis(), hCq->GetZaxis()};
    for(TAxis* axis : axes)
      if(qmax > min(-axis->GetXmin(), axis->GetXmax())) {
        cerr << "Requested fit sphere extends beyond the C(q) grid.\n"; return 1;
      }

    // Select bins once, so Minuit only loops over points inside the fit sphere.
    qout2.clear(); qside2.clear(); qlong2.clear(); cqValues.clear(); binError = 0.;
    for(int x = 1; x <= hCq->GetNbinsX(); ++x) {
      double qo = hCq->GetXaxis()->GetBinCenter(x);
      if(abs(qo) > qmax) continue;
      for(int y = 1; y <= hCq->GetNbinsY(); ++y) {
        double qs = hCq->GetYaxis()->GetBinCenter(y);
        if(qo*qo + qs*qs > qmax*qmax) continue;
        for(int z = 1; z <= hCq->GetNbinsZ(); ++z) {
          double ql = hCq->GetZaxis()->GetBinCenter(z);
          if(qo*qo + qs*qs + ql*ql > qmax*qmax) continue;
          double value = hCq->GetBinContent(x,y,z), error = hCq->GetBinError(x,y,z);
          if(binError == 0) binError = error;
          if(!isfinite(value) || !isfinite(error) || error <= 0 || abs(error-binError) > 1.e-9*binError) {
            cerr << "C(q) needs finite values and positive, uniform bin errors.\n"; return 1;
          }
          qout2.push_back(qo*qo/(fftHbarC*fftHbarC));
          qside2.push_back(qs*qs/(fftHbarC*fftHbarC));
          qlong2.push_back(ql*ql/(fftHbarC*fftHbarC)); cqValues.push_back(value);
        }
      }
    }
    npoints = cqValues.size(); NDF = npoints - NPARS;
    if(NDF <= 0) { cerr << "Too few bins inside the fit sphere.\n"; return 1; }

    // Fit lambda, Rout, Rside, Rlong, alpha, as in the original Minuit sketch.
    ROOT::Math::Functor objective(&AltChi2, NPARS);
    ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
    minimizer.SetFunction(objective);
    minimizer.SetMaxFunctionCalls(100000); minimizer.SetTolerance(1.e-8); minimizer.SetErrorDef(1.);
    double lambdaStart = max(0.05, min(1.5, hCq->GetBinContent(hCq->FindBin(0.,0.,0.))-1.));
    minimizer.SetLimitedVariable(0, "lambda", lambdaStart, 0.01, 0., 2.);
    minimizer.SetLimitedVariable(1, "R_out", 4., 0.05, 0.05, 100.);
    minimizer.SetLimitedVariable(2, "R_side", 4., 0.05, 0.05, 100.);
    minimizer.SetLimitedVariable(3, "R_long", 6., 0.05, 0.05, 100.);
    minimizer.SetLimitedVariable(4, "alpha", 1.5, 0.01, 0.2, 2.);
    bool minimized = minimizer.Minimize(), hesse = minimizer.Hesse();
    status = minimizer.Status(); covarianceStatus = minimizer.CovMatrixStatus(); edm = minimizer.Edm();

    // Convert the fixed numerical weight back to the input error, then rescale
    // errors by sqrt(chi2/NDF). Final errors do not depend on that common weight.
    double weightRatio2 = pow(binError/fitWeight, 2);
    chi2 = AltChi2(minimizer.X()) / weightRatio2;
    chi2ndf = chi2/NDF; scale = sqrt(chi2ndf);
    valid = minimized && hesse && status == 0 && covarianceStatus == 3 && isfinite(chi2ndf);
    atLimit = false;
    double lower[NPARS] = {0.,0.05,0.05,0.05,0.2}, upper[NPARS] = {2.,100.,100.,100.,2.};
    for(int p = 0; p < NPARS; ++p) {
      par[p] = minimizer.X()[p];
      errorsRaw[p] = sqrt(max(0., minimizer.CovMatrix(p,p)*weightRatio2));
      errors[p] = errorsRaw[p]*scale;
      if(par[p]-lower[p] < 1.e-4 || upper[p]-par[p] < 1.e-4) atLimit = true;
      if(!isfinite(par[p]) || !isfinite(errors[p])) valid = false;
      for(int q = 0; q < NPARS; ++q) {
        covRaw[p*NPARS+q] = minimizer.CovMatrix(p,q)*weightRatio2;
        cov[p*NPARS+q] = covRaw[p*NPARS+q]*chi2ndf;
      }
    }
    if(!output) {
      output = TFile::Open(outputName, "CREATE");
      if(!output || output->IsZombie()) { cerr << "Cannot create " << outputName << endl; return 1; }
    }
    output->cd();
    TF3 fit(Form("fit_%s_FFT", name.Data()),
      Form("1+[0]*exp(-pow((x*x*[1]*[1]+y*y*[2]*[2]+z*z*[3]*[3])/(%.17g*%.17g),[4]/2))", fftHbarC,fftHbarC),
      -qmax,qmax, -qmax,qmax, -qmax,qmax);
    fit.SetParNames("lambda", "R_out", "R_side", "R_long", "alpha");
    fit.SetParameters(par); fit.SetParErrors(errors); fit.SetNDF(NDF); fit.SetChisquare(chi2); fit.Write();
    results.Fill(); fitted.push_back(histogram);
    if(!valid) ++failed;
    cout << name << ": kT bin=" << kt << ", qmax=" << qmax << " GeV/c, status=" << status << ", chi2/ndf=" << chi2ndf
         << ", alpha=" << par[4] << ", radii=" << par[1] << ',' << par[2] << ',' << par[3] << endl;
    delete hCq;
  }
  if(fitted.empty()) { cerr << "No selected Cq_ histograms.\n"; return 1; }
  copyFFTInfo(input, output);
  if(sources) { TTree* copy = sources->CloneTree(-1, "fast"); copy->Write(); delete copy; }
  TNamed("cqfit_input", inputName.Data()).Write();
  TNamed("cqfit_parameter_order", "lambda,R_out,R_side,R_long,alpha").Write();
  TNamed("cqfit_error_convention", "errors=errors_raw*sqrt(chi2/ndf); covariance=covariance_raw*chi2/ndf; Fourier bins are correlated").Write();
  TNamed("cqfit_selection", "Spherical |q|<=qmax, including q=0 and both signs").Write();
  results.Write(); TParameter<int>("fft_fit_complete", 1).Write();
  output->Close(); input->Close(); delete output; delete input;
  cout << "Wrote " << fitted.size() << " fits to " << outputName << " (" << failed << " invalid)\n";
  return failed ? 2 : 0;
}
