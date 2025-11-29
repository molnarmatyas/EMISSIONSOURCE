#include <Levy_proj_reader.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
#include <list>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <TLegend.h>
#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#include <TLatex.h>

using namespace std;

// Global variables
const int NKT = 16;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 4;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
//const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575};
const double ktbins[NKT + 1] = {0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,0.575};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.,0.,0.05,0.,0.,0.05,0.,0.,0.05};
const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.2,0.2,0.2,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* frames[3] = {"lab","lcms","pcms"};
const char* osl_labels[3] = {"out","side","long"};

const char* path = "~/work/figs/EPOS_project/3d_proj/";

std::vector<TH1*> histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

const double rfitmax = 80.;
const double rfitmin = 3.;

int thiskt = 0;
int thisframe = 0;
int NDF = 0;
// Define the fit function
double fitFunction(double *x, double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (N/Rcc)*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
//  return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  NDF = 0;
  double alpha = params[0];
  double R_out = params[1];
  double N_out = params[2];
  double R_side = params[3];
  double N_side = params[4];
  double R_long = params[5];
  double N_long = params[6];

  double chi2 = 0.0;

  for (int i = 0; i < 3; ++i)
  {
    TH1* hist = histograms[thiskt][thisframe][i];
    double R = (i == 0) ? R_out : (i == 1) ? R_side : R_long;
    double N = (i == 0) ? N_out : (i == 1) ? N_side : N_long;

    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
    {
      double x = hist->GetBinCenter(bin);
      double ex = hist->GetBinError(bin);
      if(x > rfitmax) continue;
      if(x < rfitmin) continue;
      double observed = hist->GetBinContent(bin);
      double expected = fitFunction(&x, new double[3]{alpha, R, N});
      if (ex > 0)
      {
        chi2 += pow((observed - expected)/ex, 2.);
      }
      NDF++;
    }
  }
  NDF -= 7;
  return chi2;
}

int main()
{
  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/star/u/kincses/work/UrQMD/levyfit/tables/levy_values_complex2.dat");
  cout << "levy reader created." << endl;

  // Open the file and get the histograms
  TFile *file = TFile::Open("../EPOS_3d_source.root");
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

//  TH3F* sourcehists[NKT][NFRAME];

  for (int ikt = 0; ikt < NKT; ++ikt)
  {
    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(45);
    Tl.SetNDC(kTRUE);
  
    TCanvas* canvas = new TCanvas("canvas", "", 1200, 1000);
    TPad *pad[9];
    pad[0] = new TPad("pad1","",0.,0.65,0.36,1.);
    pad[1] = new TPad("pad2","",0.36,0.65,0.68,1.);
    pad[2] = new TPad("pad3","",0.68,0.65,1.,1.);
    pad[3] = new TPad("pad4","",0.,0.36,0.36,0.65);
    pad[4] = new TPad("pad5","",0.36,0.36,0.68,0.65);
    pad[5] = new TPad("pad6","",0.68,0.36,1.,0.65);
    pad[6] = new TPad("pad7","",0.,0.,0.36,0.36);
    pad[7] = new TPad("pad8","",0.36,0.,0.68,0.36);
    pad[8] = new TPad("pad9","",0.68,0.,1.,0.36);
    for(int ipad = 0; ipad < 9; ipad++) pad[ipad]->Draw();

    for(int iframe = 0; iframe < 3; iframe++)
    {
      thiskt = ikt;
      thisframe = iframe;
      for(int iosl = 0; iosl < 3; iosl++)
      {
        // Form the histogram name
        TString histName = Form("pion_pair_source_avg_%s_ikt%i_ich0_%s", frames[iframe], ikt, osl_labels[iosl]);
        // Read the histogram
        TH1* sourcehist = dynamic_cast<TH1F*>(file->Get(histName));
        if (!sourcehist)
        {
          std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
          continue;
        }
        double integral = sourcehist->Integral();
        if (integral == 0)
        {
          std::cerr << "Warning: Histogram " << histName << " has zero integral." << std::endl;
          continue;
        }
        // Divide the histogram by its integral
        sourcehist->Scale(1.0 / integral);
        // Loop over all bins to divide by bin volume
        for (int x = 1; x <= sourcehist->GetNbinsX(); ++x)
        {
          // Get the bin content and error
          double content = sourcehist->GetBinContent(x);
          double error = sourcehist->GetBinError(x);
          // Calculate the bin volume
          double binVolume = sourcehist->GetXaxis()->GetBinWidth(x);
          // Update the bin content and error
          sourcehist->SetBinContent(x, content / binVolume);
          sourcehist->SetBinError(x, error / binVolume);
        }
        // Print some information about the histogram
        std::cout << "Processed Histogram " << histName << " with " << sourcehist->GetEntries() << " entries." << std::endl;
        histograms[ikt][iframe].push_back((TH1F*)sourcehist);
      }

      if (histograms[ikt][iframe].size() != 3)
      {
        std::cerr << "Error: Expected 3 histograms" << std::endl;
        return 1;
      }

      // Create the minimizer
      ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    
      // Set the minimizer properties
      minimizer->SetMaxFunctionCalls(10000);
      minimizer->SetMaxIterations(10000);
      minimizer->SetTolerance(0.001);
    
      // Create the function to be minimized
      ROOT::Math::Functor functor(&chiSquare, 7);
      minimizer->SetFunction(functor);
    
      // Set the initial parameters and their steps
      minimizer->SetVariable(0, "alpha", 1.5, 0.01);
      minimizer->SetVariable(1, "R_out", 10.0, 0.01);
      minimizer->SetVariable(2, "N_out", 10.0, 0.1);
      minimizer->SetVariable(3, "R_side", 10.0, 0.01);
      minimizer->SetVariable(4, "N_side", 10.0, 0.1);
      minimizer->SetVariable(5, "R_long", 10.0, 0.01);
      minimizer->SetVariable(6, "N_long", 10.0, 0.1);
    
      // Minimize the function
      minimizer->Minimize();
    
      // Print the results
      const double *results = minimizer->X();
      const double *errors = minimizer->Errors();
      const double chi2val = minimizer->MinValue();
      std::cout << "alpha: " << results[0] << std::endl;
      std::cout << "R_out: " << results[1] << std::endl;
      std::cout << "N_out: " << results[2] << std::endl;
      std::cout << "R_side: " << results[3] << std::endl;
      std::cout << "N_side: " << results[4] << std::endl;
      std::cout << "R_long: " << results[5] << std::endl;
      std::cout << "N_long: " << results[6] << std::endl;
      std::cout << "chi2: " << chi2val << std::endl;
      
      for(int iosl = 0; iosl < 3; iosl++)
      {
        int ipad = 3*iframe+iosl;
        gStyle->SetLabelSize(0.06,"Y");
        canvas->SetLogx(1);
        canvas->SetLogy(1);
        pad[ipad]->cd();
        gPad->SetRightMargin(rightmargins[ipad]);
        gPad->SetLeftMargin(leftmargins[ipad]);
        gPad->SetTopMargin(topmargins[ipad]);
        pad[ipad]->SetBottomMargin(bottommargins[ipad]);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        histograms[ikt][iframe][iosl]->SetStats(0);
        histograms[ikt][iframe][iosl]->SetTitle("");
        histograms[ikt][iframe][iosl]->GetXaxis()->SetRangeUser(0.9,199.);
        histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleSize(0.1);
        histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleSize(0.1);
        histograms[ikt][iframe][iosl]->GetYaxis()->SetLabelSize(0.08);
        if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitle("D(#rho)");
        if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleOffset(0.9);
        if(iframe == 2 && iosl == 2) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitle("#rho [fm]");
        if(iframe == 2 && iosl == 2) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleOffset(0.9);
        histograms[ikt][iframe][iosl]->SetMinimum(5.e-6);
        histograms[ikt][iframe][iosl]->SetMaximum(30.);
        TF1* f_levyfunc = new TF1(Form("levyfunc%i%i",iframe,iosl), fitFunction, 1., 5000., 3);
        f_levyfunc->SetParNames("alpha","R","norm");
        f_levyfunc->SetLineStyle(1);
        double alpha = results[0];
        double dalpha = errors[0];
        double Rosl = results[1+2*iosl];
        double dRosl = errors[1+2*iosl];
        double Nosl = results[2+2*iosl];
        double dNosl = errors[2+2*iosl];
        f_levyfunc->SetParameters(results[0],Rosl,Nosl);
        histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelSize(0.08);
        histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelOffset(-0.0);
        histograms[ikt][iframe][iosl]->Draw("pe");
        f_levyfunc->SetLineStyle(2);
        f_levyfunc->DrawCopy("same");
        f_levyfunc->SetLineStyle(1);
        f_levyfunc->SetRange(rfitmin,rfitmax);
        f_levyfunc->DrawCopy("same");
        double conflev = TMath::Prob(chi2val, NDF);
        Tl.SetTextSize(30);
        if(iframe == 0 && iosl == 0) Tl.DrawLatex(0.45, 0.82, "OUT");
        if(iframe == 0 && iosl == 1) Tl.DrawLatex(0.45, 0.82, "SIDE");
        if(iframe == 0 && iosl == 2) Tl.DrawLatex(0.45, 0.82, "LONG");
        if(iframe == 0 && iosl == 0) Tl.DrawLatex(0.25, 0.50, "LAB");
        if(iframe == 1 && iosl == 0) Tl.DrawLatex(0.25, 0.50, "LCMS");
        if(iframe == 2 && iosl == 0) Tl.DrawLatex(0.25, 0.45, "PCMS");
        Tl.SetTextSize(20);
        Tl.DrawLatex(0.60-0.05*iosl, 0.75+iframe*0.1, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
        Tl.DrawLatex(0.60-0.05*iosl, 0.69+iframe*0.1, Form("conf.lev. = %1.5f", conflev));
        Tl.DrawLatex(0.60-0.05*iosl, 0.63+iframe*0.1, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
        Tl.DrawLatex(0.60-0.05*iosl, 0.57+iframe*0.1, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));

      }
      // Clean up
      delete minimizer;
    }
    canvas->cd();
    Tl.SetTextSize(30);
    Tl.DrawLatex(0.05, 0.97, Form("EPOS 200 GeV 0#minus10%% AuAu PION PAIR SOURCE PROJECTIONS, #LTm_{T}#GT = %1.3f",sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
    canvas->SaveAs(Form("%sprojections_ikt%i_ich0.pdf", path, ikt));
  }
  file->Close();
  delete file;
  return 0;
}

