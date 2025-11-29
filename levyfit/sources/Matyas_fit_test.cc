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
const int NFILE = 400;
const int NEVT = 30;
const int NKT = 5;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 4;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const double ktbins[NKT + 1] = {0.,0.1,0.2,0.3,0.4,0.5};
//const double ktbins[NKT + 1] = {0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,0.575};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.,0.,0.05,0.,0.,0.05,0.,0.,0.05};
const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* frames[3] = {"lcms","pcms","lab"};
const char* osl_labels[3] = {"out","side","long"};

const char* path = "~/work/figs/EPOS_project/matyat_test/";

TH1* histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

const double rfitmax = 100.;
const double rfitmin = 5.;

int thiskt = 0;
int thisframe = 0;
int NDF = 0;
// Define the fit function
double fitFunction(const double *x, const double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (N/Rcc/Rcc/Rcc)*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
//  return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  NDF = 0;
  double chi2 = 0.0;

  TH1* hist = (TH1D*)histograms[thiskt][thisframe]->Clone();
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
  {
    double x = hist->GetBinCenter(bin);
    double ex = hist->GetBinError(bin);
    if(x > rfitmax) continue;
    if(x < rfitmin) continue;
    double observed = hist->GetBinContent(bin);
    double expected = fitFunction(&x, params);
    if (ex > 0)
    {
      chi2 += pow((observed - expected)/ex, 2.);
    }
    NDF++;
  }
  NDF -= 3;
  return chi2;
}
// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  NDF = 0;
  double logL = 0.0;
  double integral = histograms[thiskt][thisframe]->Integral();
  for (int ibin = 1; ibin <= histograms[thiskt][thisframe]->GetNbinsX(); ++ibin)
  {
    double x = histograms[thiskt][thisframe]->GetXaxis()->GetBinCenter(ibin);
    if(x > rfitmax) continue;
    if(x < rfitmin) continue;
    double binVolume = ( histograms[thiskt][thisframe]->GetXaxis()->GetBinWidth(ibin) );
      
    double observed = histograms[thiskt][thisframe]->GetBinContent(ibin);
    double expected = fitFunction(&x, params)*binVolume*integral*(x*x*4.*M_PI);
    if(expected <= 0) continue; // Avoid log(0) or negative expected values
    if(observed != 0)
      logL += expected + observed*log(observed/expected) - observed;
    NDF++;
  }
  NDF -= 3;
  return 2.0 * logL;
}

int main()
{
  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/star/u/kincses/work/EPOS_project/3D_analysis/levyfit/tables/levy_values_moreprecise.dat");
  cout << "levy reader created." << endl;

  TH1* alphahist[NKT];
  TH1* Rhist[NKT];
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    alphahist[ikt] = new TH1F(Form("alphahist_ikt%i",ikt),"",100,0.,2.);
    Rhist[ikt] = new TH1F(Form("Rhist_ikt%i",ikt),"",100,5.,15.);
  }
  // Open the file and get the histograms
  TFile *file = TFile::Open("./matyas_dr.root");
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

//  TH3F* sourcehists[NKT][NFRAME];
  TCanvas* canvas = new TCanvas("c1", "", 600, 600);

  for(int ifile = 1; ifile < 2; ifile++)
  {
    for(int ievt = 0; ievt < 1; ievt++)
      for (int ikt = 1; ikt < NKT; ++ikt)
      {
        cout << "ifile,ievt,ikt: " << ifile << "," << ievt << "," << ikt << "," << endl;
        TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(45);
        Tl.SetNDC(kTRUE);
    
        for(int iframe = 0; iframe < 3; iframe++)
        {
          thiskt = ikt;
          if(iframe != thisframe) continue;
          if(histograms[ikt][iframe]) histograms[ikt][iframe]->Reset();
          // Form the histogram name
          TString histName = Form("spatdistLCMS%i", ikt);
          // Read the histogram
          histograms[ikt][iframe] = dynamic_cast<TH1D*>(file->Get(histName));
          if (!histograms[ikt][iframe])
          {
            std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
            continue;
          }
          // Print some information about the histogram
          std::cout << "Processed Histogram " << histName << " with " << histograms[ikt][iframe]->GetEntries() << " entries." << std::endl;
//            sourcehist->Scale(1.0 / sourcehist->Integral(1,sourcehist->GetNbinsX()));
          cout << "Integral w overflow after normalization: " << histograms[ikt][iframe]->Integral(0,histograms[ikt][iframe]->GetNbinsX()+1) << endl;
    
          // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        
          // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
        
          // Create the function to be minimized
//          ROOT::Math::Functor functor(&logLikelihood, 3);
          ROOT::Math::Functor functor(&chiSquare, 3);
          minimizer->SetFunction(functor);
        
          // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.6, 0.01, 0.5, 2.0);
//          minimizer->SetFixedVariable(0, "alpha", 1.3);
          minimizer->SetLimitedVariable(1, "R", 8.5, 0.01, 0., 15.);
//          minimizer->SetFixedVariable(3, "R_long", 20.0);
//          minimizer->SetFixedVariable(4, "N_out", 1.);
//          minimizer->SetFixedVariable(5, "N_side", 1.);
//          minimizer->SetFixedVariable(6, "N_long", 1.);
          minimizer->SetVariable(2, "N", 1., 0.01);
        
            double integral = histograms[ikt][iframe]->Integral();//0,histograms[ikt][iframe][iosl]->GetNbinsX()+1);
            for (int x = 1; x <= histograms[ikt][iframe]->GetNbinsX(); ++x)
            {
              double content = histograms[ikt][iframe]->GetBinContent(x);
              double error = histograms[ikt][iframe]->GetBinError(x);
              double binVolume = histograms[ikt][iframe]->GetXaxis()->GetBinWidth(x);
              histograms[ikt][iframe]->SetBinContent(x, content / binVolume);
              histograms[ikt][iframe]->SetBinError(x, error / binVolume);
            }
            TF1* f_r2 = new TF1("f_r2","1./x/x/4./pi",0.,1.e8);
            histograms[ikt][iframe]->Multiply(f_r2,1.);
            histograms[ikt][iframe]->Scale(1.0 / integral);
            delete f_r2;

          // Minimize the function
          minimizer->Minimize();
        
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          std::cout << "alpha: "  << results[0] << std::endl;
          std::cout << "R: "  << results[1] << std::endl;
          std::cout << "N "       << results[2] << std::endl;
          std::cout << "chi2: "   << chi2val << std::endl;
          
            gStyle->SetLabelSize(0.06,"Y");
            canvas->SetLogx(1);
            canvas->SetLogy(1);
            canvas->SetRightMargin(0.01);
            canvas->SetLeftMargin(0.15);
            canvas->SetTopMargin(0.01);
            canvas->SetBottomMargin(0.15);

            histograms[ikt][iframe]->SetStats(0);
            histograms[ikt][iframe]->SetTitle("");
            histograms[ikt][iframe]->GetXaxis()->SetRangeUser(0.2,5000.);
            histograms[ikt][iframe]->GetXaxis()->SetTitleSize(0.1);
            histograms[ikt][iframe]->GetYaxis()->SetTitleSize(0.1);
            histograms[ikt][iframe]->GetYaxis()->SetLabelSize(0.08);
            histograms[ikt][iframe]->GetYaxis()->SetTitle("D(#rho)");
            histograms[ikt][iframe]->GetYaxis()->SetTitleOffset(0.9);
            histograms[ikt][iframe]->GetXaxis()->SetTitle("#rho [fm]");
            histograms[ikt][iframe]->GetXaxis()->SetTitleOffset(0.7);
            histograms[ikt][iframe]->SetMinimum(5.e-13);
            histograms[ikt][iframe]->SetMaximum(1.e-2);
            TF1* f_levyfunc = new TF1(Form("levyfunc%i",iframe), fitFunction, 0.1, 5000., 3);
            f_levyfunc->SetParNames("alpha","R","norm");
            f_levyfunc->SetLineStyle(1);
            double alpha = results[0];
            double dalpha = errors[0];
            double R = results[1];
            double dR = errors[1];
            double N = results[2];
            double dN = errors[2];

            alphahist[ikt]->Fill(alpha);
            Rhist[ikt]->Fill(R);

            f_levyfunc->SetParameters(alpha,R,N);
            histograms[ikt][iframe]->GetXaxis()->SetLabelSize(0.08);
            histograms[ikt][iframe]->GetXaxis()->SetLabelOffset(-0.03);
            histograms[ikt][iframe]->Draw("pe");
            f_levyfunc->SetLineStyle(2);
            f_levyfunc->DrawCopy("same");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->SetRange(rfitmin,rfitmax);
            f_levyfunc->DrawCopy("same");
            double conflev = TMath::Prob(chi2val, NDF);
            Tl.SetTextSize(20);
            Tl.DrawLatex(0.58-0.05, 0.78, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
            Tl.DrawLatex(0.58-0.05, 0.73, Form("conf.lev. = %1.5f", conflev));
            Tl.DrawLatex(0.58-0.05, 0.68, Form("R = (%1.2f #pm %1.2f) fm",R, dR));
            Tl.DrawLatex(0.58-0.05, 0.63, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
            Tl.DrawLatex(0.58-0.05, 0.58, Form("N = %1.2f #pm %1.2f",N, dN));
            delete minimizer;
            delete f_levyfunc;
          }
          // Clean up
          Tl.SetTextSize(30);
          Tl.DrawLatex(0.05, 0.94, Form("EPOS 200 GeV 0#minus10%% AuAu PION PAIR SOURCE, #LTm_{T}#GT = %1.3f, LCMS",sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
          canvas->SaveAs(Form("%s%s/onedsource_ifile%i_ievt%i_ikt%i_ich0.pdf", path, frames[thisframe], ifile, ievt, ikt));
          canvas->Clear();
        }
  }
  delete canvas;
  delete myLevy_reader;
  cout << "about to close input file." << endl;
  file->Close();
  delete file;
  cout << "input file closed, delete done." << endl;
  return 0;
}

