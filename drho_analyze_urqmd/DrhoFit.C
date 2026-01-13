#include <TFile.h>
#include <TH1D.h>
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <TString.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "levy_calc.h"

const int NKT = 10;

const int NEVTS = 100;

TH1D *histarray[NKT];
TH1D *hist;

double MyChi2(const double* par) {
  double N = par[0];
  double R = par[1];
  double alpha = par[2];
  double chi2 = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    double x = hist->GetBinCenter(bin);
    double y = hist->GetBinContent(bin);
    double err = hist->GetBinError(bin);
    if (err == 0) continue;
    double Rcc = R * pow(2., 1. / alpha);
    double model = N * levy_calc(x / Rcc, Rcc, alpha);
    chi2 += pow((y - model) / err, 2);
  }
  return chi2;
}


void DrhoFit(const char* filename) {
  TFile* file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cerr << "Could not open file: " << filename << std::endl;
    return;
  }

  for (int iev = 0; iev < 1000; iev++) {
    for (int ikt = 0; ikt < 10; ikt++) {
      TString suffix = Form("_ev%d_ch1_KT%d", iev, ikt);
      TH1D *histtemp = (TH1D*)file->Get("D_LCMS" + suffix);
      if(iev%NEVTS==0) histarray[ikt] = (TH1D*)histtemp->Clone();
      else histarray[ikt]->Add(histtemp);

      if(iev%NEVTS==NEVTS-1) {
        hist = histarray[ikt];
        hist->Print("base");
        hist->Scale(1.0/hist->GetEntries());
        ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined );
        min.SetMaxFunctionCalls(1000000);
        min.SetMaxIterations(100000);
        min.SetTolerance(0.001);
        ROOT::Math::Functor f(&MyChi2,3);
        min.SetFunction(f);
        min.SetVariable(0, "N", 1.0, 0.1);
        min.SetVariable(1, "R", 5.0, 0.1);
        min.SetVariable(2, "alpha", 1.5, 0.1);
        min.SetVariableLimits(1, 3, 15);
        min.SetVariableLimits(2, 0.6, 1.99);
        min.Minimize();
        //min.PrintResults();
        cout << min.X()[0] << "\t" << min.X()[1] << "\t" << min.X()[2] << endl;
        histarray[ikt]->Reset();
      }
      
    }
  }
  file->Close();
}
