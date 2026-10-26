// Build: make -C FFTmethod FFT=1 CqFit.exe
#include "fft_common.h"
#include <TH3.h>
#include <TF3.h>
#include <TKey.h>
#include <TClass.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <limits>
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

// Central-axis slices of the same 3D data/model used in the fit.
bool plotCqFit(TH3* cq, TF3& fit, int kt, double errorScale, double rawChi2Ndf, bool valid, bool atLimit,
               TFile* input, const TString& plotName)
{
  TString folder = gSystem->DirName(plotName);
  if(gSystem->AccessPathName(folder) && gSystem->mkdir(folder, true) != 0) {
    cerr << "Cannot create plot folder " << folder << endl; return false;
  }
  const char* directions[] = {"out", "side", "long"};
  const char* axisTitles[] = {"(q_{out},0,0) [GeV/c]", "(0,q_{side},0) [GeV/c]", "(0,0,q_{long}) [GeV/c]"};
  TAxis* axes[] = {cq->GetXaxis(), cq->GetYaxis(), cq->GetZaxis()};
  int zero[] = {axes[0]->FindFixBin(0.), axes[1]->FindFixBin(0.), axes[2]->FindFixBin(0.)};
  double qmax = fit.GetXmax(), qPlot = max(0.30, 1.4*qmax);
  for(TAxis* axis : axes)
    qPlot = min(qPlot, min(-axis->GetBinCenter(1), axis->GetBinCenter(axis->GetNbins())));

  TGraphErrors data[3];
  TGraph fullCurve[3], fittedCurve[3];
  double ymin = 0.75, ymax = 2.25;
  for(int d = 0; d < 3; ++d) {
    data[d].SetName(Form("Cq_slice_%s_FFT", directions[d]));
    data[d].SetTitle(Form(";%s;%s", axisTitles[d], d == 0 ? "C(q)" : ""));
    data[d].SetMarkerStyle(20); data[d].SetMarkerSize(0.65); data[d].SetLineWidth(1);
    for(int b = 1; b <= axes[d]->GetNbins(); ++b) {
      double q = axes[d]->GetBinCenter(b);
      if(abs(q) > qPlot) continue;
      int bins[] = {zero[0], zero[1], zero[2]}; bins[d] = b;
      double value = cq->GetBinContent(bins[0],bins[1],bins[2]);
      double error = cq->GetBinError(bins[0],bins[1],bins[2])*errorScale;
      if(!isfinite(value) || !isfinite(error)) continue;
      int point = data[d].GetN();
      data[d].SetPoint(point, q, value); data[d].SetPointError(point, 0., error);
      ymin = min(ymin, value-error-0.08); ymax = max(ymax, value+error+0.12);
    }
    fullCurve[d].SetName(Form("Cq_fullfit_%s_FFT", directions[d]));
    fittedCurve[d].SetName(Form("Cq_fit_%s_FFT", directions[d]));
    for(int point = 0; point <= 600; ++point) {
      double q[] = {0.,0.,0.};
      q[d] = -qPlot + 2.*qPlot*point/600.;
      fullCurve[d].SetPoint(point, q[d], fit.Eval(q[0],q[1],q[2]));
      q[d] = -qmax + 2.*qmax*point/600.;
      fittedCurve[d].SetPoint(point, q[d], fit.Eval(q[0],q[1],q[2]));
    }
    fullCurve[d].SetLineColor(kRed+1); fullCurve[d].SetLineWidth(3); fullCurve[d].SetLineStyle(2);
    fittedCurve[d].SetLineColor(kRed+1); fittedCurve[d].SetLineWidth(3);
  }

  TString sample = "LCMS correlation, FFT";
  TNamed* energy = (TNamed*)input->Get("source_energy");
  TNamed* centrality = (TNamed*)input->Get("source_centrality");
  if(energy) {
    TString number = energy->GetTitle(); number.ReplaceAll("p", ".");
    sample = "#sqrt{s_{NN}} = " + number + " GeV";
  }
  if(centrality) sample += TString(", ") + centrality->GetTitle() + "%";
  TParameter<int>* cut = (TParameter<int>*)input->Get("source_qlcms_cut");

  // Wider first pad leaves room for the shared y axis, as in example.png.
  TCanvas canvas(Form("fitplot_%s_FFT", cq->GetName()), "3D C(q) fit", 1800, 620);
  double edges[] = {0.,0.37,0.685,1.};
  for(int d = 0; d < 3; ++d) {
    canvas.cd();
    TPad* pad = new TPad(Form("fitpad_%s_FFT", directions[d]), "", edges[d],0.,edges[d+1],1.);
    pad->SetLeftMargin(d == 0 ? 0.15 : 0.);
    pad->SetRightMargin(d == 2 ? 0.025 : 0.);
    pad->SetBottomMargin(0.18); pad->SetTopMargin(0.025); pad->SetTicks(1,1);
    pad->Draw(); pad->cd();
    data[d].SetMinimum(ymin); data[d].SetMaximum(ymax);
    data[d].Draw("AP");
    data[d].GetXaxis()->SetLimits(-1.20*qPlot,1.20*qPlot); // space for labels at panel joins
    for(TAxis* axis : {data[d].GetXaxis(), data[d].GetYaxis()}) {
      axis->SetTitleFont(43); axis->SetTitleSize(30);
      axis->SetLabelFont(43); axis->SetLabelSize(25);
      axis->SetNdivisions(505);
    }
    data[d].GetXaxis()->SetTitleOffset(1.15);
    data[d].GetYaxis()->SetTitleOffset(1.4);
    if(d > 0) data[d].GetYaxis()->SetLabelSize(0.);
    fullCurve[d].Draw("L SAME"); fittedCurve[d].Draw("L SAME"); data[d].Draw("P SAME");

    double left = d == 0 ? 0.18 : 0.04;
    TLatex text; text.SetNDC(); text.SetTextFont(43); text.SetTextSize(24);
    if(d == 0) {
      text.DrawLatex(left,0.925,"UrQMD, ^{197}Au+^{197}Au");
      text.DrawLatex(left,0.88,sample);
      text.DrawLatex(left,0.72,Form("#alpha = %.3f #pm %.2g", fit.GetParameter(4),fit.GetParError(4)));
      text.DrawLatex(left,0.66,Form("#lambda^{*} = %.3f #pm %.2g", fit.GetParameter(0),fit.GetParError(0)));
      text.DrawLatex(left,0.60,Form("#chi^{2}_{sc}/NDF = %.2f", fit.GetChisquare()/fit.GetNDF()));
    } else if(d == 1) {
      text.DrawLatex(left,0.935,"Like-sign #pi#pi, event averaged");
      if(kt >= 0 && kt < NKT)
        text.DrawLatex(left,0.88,Form("%.3f < k_{T} < %.3f GeV/c", ktbins[kt],ktbins[kt+1]));
      text.DrawLatex(left,0.72,Form("#sigma_{C} = %.3g (rescaled)", cq->GetBinError(zero[0],zero[1],zero[2])*errorScale));
      text.DrawLatex(left,0.66,Form("#chi^{2}_{raw}/NDF = %.2g", rawChi2Ndf));
      text.DrawLatex(left,0.60,"Scaled by construction");
    } else {
      text.DrawLatex(left,0.935,Form("Spherical fit: |q| #leq %.2f GeV/c", qmax));
      if(cut && cut->GetVal() >= 0 && cut->GetVal() < 3)
        text.DrawLatex(left,0.88,Form("Analysis q_{LCMS}: %s", qLCMS_cut[cut->GetVal()]));
      text.SetTextColor(valid ? kBlack : kRed+1);
      text.DrawLatex(left,0.72,valid ? "Fit valid" : "INVALID FIT (see fit_results_FFT)");
      if(atLimit) { text.SetTextColor(kOrange+7); text.DrawLatex(left,0.66,"Parameter at limit"); }
      text.SetTextColor(kBlack);
      text.DrawLatex(left,0.60,Form("NDF = %d", fit.GetNDF()));
    }
    text.DrawLatex(left,0.80,Form("R_{%s} = (%.2f #pm %.2g) fm", directions[d],fit.GetParameter(d+1),fit.GetParError(d+1)));
    TLegend legend(left,0.205,0.98,0.275);
    legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(43); legend.SetTextSize(25);
    if(d == 0) legend.AddEntry(&data[d],"C(q), UrQMD","pe");
    if(d == 1) legend.AddEntry(&fittedCurve[d],"3D L#acute{e}vy fit","l");
    if(d == 2) legend.AddEntry(&fullCurve[d],"Extrapolation outside fit range","l");
    legend.DrawClone();
  }
  canvas.SaveAs(plotName + ".png"); canvas.SaveAs(plotName + ".pdf");
  canvas.Write();
  return !gSystem->AccessPathName(plotName + ".png") && !gSystem->AccessPathName(plotName + ".pdf");
}

int main(int argc, char** argv)
{
  if(argc < 2 || string(argv[1]) == "--help") {
    cout << "Usage: CqFit.exe input_Cq_FFT.root [--qmax Q] [--hist NAME] [--output file_FFT.root]\n"
         << "Default: per-kT spherical |q| limits fftQmaxFit in header_for_all_emissionsource.h.\n"
         << "--qmax overrides these with one limit in GeV/c, or NKT comma-separated limits.\n"
         << "--fit-range 0|1|2 uses the baseline/narrow/wide qmax scale from the header.\n"
         << "Plots: --plot-dir DIR (default figs/fitting/lcms/FFT) or --no-plots.\n";
    return argc < 2 ? 1 : 0;
  }
  if(!useFFTmethod) { cerr << "FFT is disabled in the global header/build.\n"; return 1; }
  TString inputName = argv[1], outputName = fftOutputName(inputName, "_fitresults");
  TString onlyHist = "";
  TString plotDir = "figs/fitting/lcms/FFT";
  bool makePlots = true;
  bool customOutput = false;
  int fitRange = 0;
  gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  vector<double> qmaxValues(fftQmaxFit, fftQmaxFit + NKT);
  for(int i = 2; i < argc; ++i) {
    string option = argv[i];
    if(option == "--no-plots") { makePlots = false; continue; }
    if(i + 1 == argc) { cerr << "Missing value for " << option << endl; return 1; }
    string value = argv[++i];
    if(option == "--output") { outputName = fftOutputName(value.c_str(), ""); customOutput = true; }
    else if(option == "--hist") onlyHist = value;
    else if(option == "--fit-range") fitRange = atoi(value.c_str());
    else if(option == "--plot-dir") plotDir = value;
    else if(option == "--qmax") {
      qmaxValues.clear();
      istringstream values(value);
      string entry;
      while(getline(values, entry, ',')) qmaxValues.push_back(atof(entry.c_str()));
    } else { cerr << "Unknown option: " << option << endl; return 1; }
  }
  if(fitRange < 0 || fitRange > 2) { cerr << "Fit range must be 0, 1 or 2.\n"; return 1; }
  double qmaxScale = fftQmaxFitScale[fitRange];
  for(double& q : qmaxValues) q *= qmaxScale;
  if(!customOutput) outputName = fftOutputName(inputName, (TString("_fitresults") + fftFitRangeSuffix[fitRange]).Data());
  if(qmaxValues.size() != 1 && qmaxValues.size() != NKT) {
    cerr << "Supply --qmax with one limit or NKT comma-separated limits.\n"; return 1;
  }
  for(double q : qmaxValues) if(q <= 0 || !isfinite(q)) { cerr << "Invalid qmax.\n"; return 1; }
  if(outputName == inputName) { cerr << "Input and output must have different names.\n"; return 1; }
  TString plotFolder = gSystem->BaseName(outputName); plotFolder.Resize(plotFolder.Length()-5);
  plotFolder = plotDir + "/" + plotFolder;
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
  double qmax, chi2, chi2ndf, chi2Raw, chi2ndfRaw, scale, binError, binErrorScaled, residualRMS, edm;
  double par[NPARS], errorsRaw[NPARS], errors[NPARS], covRaw[25], cov[25];
  TTree results("fit_results_FFT", "3D Levy fits of C(q)");
  results.SetDirectory(nullptr);
  results.Branch("histogram", &histogram); results.Branch("kt", &kt);
  results.Branch("status", &status); results.Branch("covariance_status", &covarianceStatus);
  results.Branch("valid", &valid); results.Branch("at_parameter_limit", &atLimit);
  results.Branch("npoints", &npoints); results.Branch("ndf", &NDF); results.Branch("qmax", &qmax);
  results.Branch("chi2", &chi2); results.Branch("chi2_ndf", &chi2ndf);
  results.Branch("chi2_raw", &chi2Raw); results.Branch("chi2_ndf_raw", &chi2ndfRaw);
  results.Branch("bin_error", &binErrorScaled); results.Branch("residual_rms", &residualRMS);
  results.Branch("error_scale", &scale); results.Branch("bin_error_raw", &binError); results.Branch("edm", &edm);
  results.Branch("parameters", par, "parameters[5]/D");
  results.Branch("errors_raw", errorsRaw, "errors_raw[5]/D"); results.Branch("errors", errors, "errors[5]/D");
  results.Branch("covariance_raw", covRaw, "covariance_raw[25]/D"); results.Branch("covariance", cov, "covariance[25]/D");

  // Ready for a final plotter; keep the tree above for fit diagnostics/covariance.
  const char* graphNames[NPARS] = {"lambda_vs_KT_FFT", "R_out_vs_KT_FFT", "R_side_vs_KT_FFT", "R_long_vs_KT_FFT", "alpha_vs_KT_FFT"};
  const char* graphTitles[NPARS] = {";k_{T} [GeV/c];#lambda^{*}", ";k_{T} [GeV/c];R_{out} [fm]",
    ";k_{T} [GeV/c];R_{side} [fm]", ";k_{T} [GeV/c];R_{long} [fm]", ";k_{T} [GeV/c];#alpha"};
  TGraphAsymmErrors graphs[NPARS];
  for(int p = 0; p < NPARS; ++p) {
    graphs[p].SetNameTitle(graphNames[p], graphTitles[p]);
    graphs[p].SetMarkerStyle(20);
  }
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

    // Uniform input errors have an arbitrary scale. Use the residual variance
    // to estimate that scale, as in arXiv:2610.29760, Supplemental Material.
    // Keep the raw chi-square; the rescaled one is NDF by construction.
    double weightRatio2 = pow(binError/fitWeight, 2);
    double sumSquares = AltChi2(minimizer.X())*fitWeight*fitWeight;
    chi2Raw = sumSquares/(binError*binError);
    chi2ndfRaw = chi2Raw/NDF; scale = sqrt(chi2ndfRaw);
    binErrorScaled = binError*scale;
    residualRMS = sqrt(sumSquares/npoints);
    chi2 = binErrorScaled > 0 ? sumSquares/(binErrorScaled*binErrorScaled) : numeric_limits<double>::quiet_NaN();
    chi2ndf = chi2/NDF;
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
        cov[p*NPARS+q] = covRaw[p*NPARS+q]*chi2ndfRaw;
      }
    }
    if(!output) {
      output = TFile::Open(outputName, "RECREATE");
      if(!output || output->IsZombie()) { cerr << "Cannot create " << outputName << endl; return 1; }
    }
    output->cd();
    TF3 fit(Form("fit_%s_FFT", name.Data()),
      Form("1+[0]*exp(-pow((x*x*[1]*[1]+y*y*[2]*[2]+z*z*[3]*[3])/(%.17g*%.17g),[4]/2))", fftHbarC,fftHbarC),
      -qmax,qmax, -qmax,qmax, -qmax,qmax);
    fit.SetParNames("lambda", "R_out", "R_side", "R_long", "alpha");
    fit.SetParameters(par); fit.SetParErrors(errors); fit.SetNDF(NDF); fit.SetChisquare(chi2); fit.Write();
    if(makePlots && !plotCqFit(hCq,fit,kt,scale,chi2ndfRaw,valid,atLimit,input,plotFolder + "/" + name + "_FFT")) {
      cerr << "Could not save fit plots in " << plotFolder << endl; return 1;
    }
    results.Fill(); fitted.push_back(histogram);
    if(valid && kt >= 0 && kt < NKT) {
      for(int p = 0; p < NPARS; ++p) {
        int point = graphs[p].GetN();
        graphs[p].SetPoint(point, kT_center[kt], par[p]);
        graphs[p].SetPointError(point, kT_center[kt]-ktbins[kt], ktbins[kt+1]-kT_center[kt], errors[p], errors[p]);
      }
    }
    if(!valid) ++failed;
    cout << name << ": kT bin=" << kt << ", qmax=" << qmax << " GeV/c, status=" << status << ", chi2/ndf=" << chi2ndf
         << " (rescaled by construction), raw=" << chi2ndfRaw << ", sigma_C=" << binErrorScaled << ", residual_RMS=" << residualRMS
         << ", alpha=" << par[4] << ", radii=" << par[1] << ',' << par[2] << ',' << par[3] << endl;
    delete hCq;
  }
  if(fitted.empty()) { cerr << "No selected Cq_ histograms.\n"; return 1; }
  copyFFTInfo(input, output);
  if(sources) { TTree* copy = sources->CloneTree(-1, "fast"); copy->Write(); delete copy; }
  TNamed("cqfit_input", inputName.Data()).Write();
  TNamed("cqfit_parameter_order", "lambda,R_out,R_side,R_long,alpha").Write();
  TParameter<int>("cqfit_error_schema", 2).Write();
  TNamed("cqfit_error_convention", "error_scale=sqrt(chi2_ndf_raw); bin_error=bin_error_raw*error_scale; errors=errors_raw*error_scale; covariance=covariance_raw*error_scale^2; chi2 and chi2_ndf use bin_error and are rescaled by construction; Fourier-bin correlations are neglected").Write();
  TNamed("cqfit_uncertainty_reference", "https://arxiv.org/pdf/2610.29760#page=6").Write();
  TNamed("cqfit_selection", "Spherical |q|<=qmax, including q=0 and both signs").Write();
  TParameter<int>("cqfit_range_variation", fitRange).Write();
  TParameter<double>("cqfit_qmax_scale", qmaxScale).Write();
  for(int p = 0; p < NPARS; ++p) { graphs[p].Sort(); graphs[p].Write(); }
  results.Write(); TParameter<int>("fft_fit_complete", 1).Write();
  output->Close(); input->Close(); delete output; delete input;
  cout << "Wrote " << fitted.size() << " fits to " << outputName << " (" << failed << " invalid)\n";
  return failed ? 2 : 0;
}
