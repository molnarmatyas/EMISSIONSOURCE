// Standard Library
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ROOT Core
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TVector3.h>

// ROOT Graphics and UI
#include <TLatex.h>
#include <TLegend.h>
#include <TEllipse.h>

// ROOT Math
#include <Math/Factory.h>
#include <Math/Functor.h>

// Levy reader
#include "Levy_proj_reader.h"

using namespace std;

Levy_reader* myLevy_reader;

double fit_min = 1.;
double fit_max = 100.;

const int NPAR = 3; //alpha, R, N
int NDF; // For tracking the number of degrees of freedom in the fit

TH1* histogram;

double LevyAvg1DFunc(const double *x, const double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (2.*N/Rcc)*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

double myChi2(const double *params)
{
  NDF = 0;
  double chi2 = 0.0;
  double integral = histogram->Integral();
  for (int ibin = 1; ibin <= histogram->GetNbinsX(); ++ibin)
  {
    double x = histogram->GetXaxis()->GetBinCenter(ibin);
    if(x > fit_max || x < fit_min) continue;
    double binwidth = histogram->GetXaxis()->GetBinWidth(ibin);
    double data = histogram->GetBinContent(ibin);
    double error = histogram->GetBinError(ibin);
    if(error == 0.) continue; // Avoid division by zero
    double center = histogram->GetBinCenter(ibin);
    double phasespace = center*center*4*M_PI;
    double theor = LevyAvg1DFunc(&x, params)*binwidth*integral*phasespace;
    double chi = (data-theor)/error;
    chi2 += chi*chi;
    NDF++;
  }
  NDF -= 5;
  return chi2;
}

int main(int argc, char** argv)
{
	// Initialize Levy reader object
  myLevy_reader = new Levy_reader("levy_proj3D_values.dat"); // Download this file from https://csanad.web.elte.hu/phys/levy_proj3D_values.dat
	
	// Obtain histograms from data file
	//TFile* infile = new TFile("AuAu_7p7_drho_merged.root"); // An example file
	//histogram = (TH1*)infile->Get("D_LCMS_KT3"); // An example KT bin
	TFile* infile = new TFile("D_rho_0.root"); // An example file
	histogram = (TH1*)infile->Get("hDrhoEvt0Cent0kt3"); // An example KT bin
  for (int ibin = 1; ibin <= histogram->GetNbinsX(); ibin++)
  {
    double content = histogram->GetBinContent(ibin);
    double error = histogram->GetBinError(ibin);
    double binwidth = histogram->GetXaxis()->GetBinWidth(ibin);
    double center = histogram->GetBinCenter(ibin);
    double phasespace = center*center*4*M_PI;
    histogram->SetBinContent(ibin, content * binwidth * phasespace);
    histogram->SetBinError(ibin, error * binwidth * phasespace);
  }

  // Create the minimizer
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  // Set the minimizer properties
  minimizer->SetMaxFunctionCalls(10000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  
  // Create the function to be minimized
  ROOT::Math::Functor functor(&myChi2, NPAR);
  minimizer->SetFunction(functor);
	
  // Set the initial parameters and their steps
  minimizer->SetLimitedVariable(0, "alpha", 1.5, 0.01, 0.5, 2.0);
  minimizer->SetLimitedVariable(1, "R", 5.0, 0.01, 2.0, 12.0);
  minimizer->SetVariable(2, "N", 1., 0.01);
	
  // Minimize the function
  minimizer->Minimize();
  minimizer->PrintResults();
	
  // Save the results
  const double *params = minimizer->X();
  const double *errors = minimizer->Errors();
  const double chi2val = minimizer->MinValue();
  double CL = TMath::Prob(chi2val, NDF);
	cerr << "chi^2/NDF = " << chi2val << "/" << NDF << " -> C.L. = " << CL << endl;
	
	// Create canvas for plotting
  TCanvas *c1 = new TCanvas("c1", "", 1000, 1000);
  c1->SetLogy();
  c1->SetLogx();
  c1->SetLeftMargin(0.13);
	
	// ROOT functions for drawing the fit
  TF1* f_levy_func_fitted;
  TF1* f_levy_func_fullrange;
	// Rescale histogram for plotting
  double integral = histogram->Integral();
  for (int ibin = 1; ibin <= histogram->GetNbinsX(); ibin++)
  {
    double content = histogram->GetBinContent(ibin);
    double error = histogram->GetBinError(ibin);
    double binwidth = histogram->GetXaxis()->GetBinWidth(ibin);
    double center = histogram->GetBinCenter(ibin);
    double phasespace = center*center*4*M_PI;
    histogram->SetBinContent(ibin, content / binwidth / phasespace);
    histogram->SetBinError(ibin, error / binwidth / phasespace);
  }
  histogram->Scale(1.0 / integral);

  // Set histogram plotting properties and draw it
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitle("|#rho_{LCMS}|");
  histogram->GetYaxis()->SetTitle("D(#rho_{LCMS})");
  histogram->GetYaxis()->SetTitleOffset(1.5);
  histogram->SetStats(0);
  histogram->GetYaxis()->SetRangeUser(1e-10, 1e-3);
  histogram->GetXaxis()->SetRangeUser(0.1, 200.);
  histogram->SetMarkerColor(kBlue);
  histogram->SetMarkerSize(1);
  histogram->SetMarkerStyle(8);
  histogram->SetLineColor(kBlue);
  histogram->GetYaxis()->SetTitleSize(0.04);
  histogram->GetXaxis()->SetTitleSize(0.04);
  histogram->Draw("pe");
  
  // Create a new function with the fitted parameters
	f_levy_func_fitted = new TF1("f_levy_func_fitted", LevyAvg1DFunc, fit_min, fit_max, NPAR);
  f_levy_func_fitted->SetParameters(params[0], params[1], params[2]); // alpha, R, N
  f_levy_func_fitted->SetLineStyle(1); 
  f_levy_func_fitted->SetLineWidth(4);
  f_levy_func_fitted->SetLineColor(kRed);
  f_levy_func_fitted->Draw("same");
  
	// Another function to be drawin in the full range
	f_levy_func_fullrange = new TF1("f_levy_func_fullrange", LevyAvg1DFunc, 0.1, 200., NPAR);
  f_levy_func_fullrange->SetParameters(params[0], params[1], params[2]); // alpha, R, N
  f_levy_func_fullrange->SetLineStyle(2);
  f_levy_func_fullrange->SetLineWidth(4);
  f_levy_func_fullrange->SetLineColor(kRed);
  f_levy_func_fullrange->Draw("same");
  
	// Draw parameters on the plot
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.035);
  l.DrawLatex(0.2,0.55,  Form("#alpha = %.2f #pm %.2f",params[0],errors[0]));
  l.DrawLatex(0.2,0.50, Form("R = (%.2f #pm %.2f) fm",params[1], errors[1]));
  l.DrawLatex(0.2,0.45, Form("N = (%.4f #pm %.4f) fm",params[2], errors[2]));
  l.DrawLatex(0.2,0.40, Form("#chi^{2} / NDF = %.2f / %d",chi2val,NDF));
  l.DrawLatex(0.2,0.35, Form("C.L.: %.2f%%",CL * 100));
  
  TLegend *legend = new TLegend(0.65, 0.75, 0.9, 0.9); 
  legend->SetTextSize(0.04);
  legend->AddEntry(histogram, "D(#rho_{LCMS})", "ep");
  legend->AddEntry(f_levy_func_fitted, "Levy fit", "l");
  legend->Draw("same");
  c1->Print("rhist.pdf");
  
	delete c1;
  return 0;
}



