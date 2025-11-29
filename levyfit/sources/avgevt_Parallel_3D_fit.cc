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
const int NFILE = 400.;//9986;
const int NEVT = 30;
const int NKT = 10;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 4;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
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

const char* path = "~/work/figs/EPOS_project/3d_proj/";

std::vector<TH1*> histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

const double rfitmax = 100.;
const double rfitmin = 1.;

const bool donotdraw = false;

int thiskt = 0;
const int thisframe = 0;
int NDF = 0;
// Define the fit function
double fitFunction(double *x, double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (2.*N/Rcc)*(myLevy_reader->getValue_1d(alpha, x[0]/Rcc));
//  return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  NDF = 0;
  double alpha  = params[0];
  double R_out  = params[1];
  double R_side = params[2];
  double R_long = params[3];
  double N  = params[4];

  double chi2 = 0.0;

  for (int i = 0; i < 3; ++i)
  {
    double R = (i == 0) ? R_out : (i == 1) ? R_side : R_long;

    for (int bin = 1; bin <= histograms[thiskt][thisframe][i]->GetNbinsX(); ++bin)
    {
      double x = histograms[thiskt][thisframe][i]->GetBinCenter(bin);
      double ex = histograms[thiskt][thisframe][i]->GetBinError(bin);
      if(x > rfitmax) continue;
      if(x < rfitmin) continue;
      double observed = histograms[thiskt][thisframe][i]->GetBinContent(bin);
      double expected = fitFunction(&x, new double[3]{alpha, R, N});
      if (ex > 0)
      {
        chi2 += pow((observed - expected)/ex, 2.);
      }
      NDF++;
    }
  }
  NDF -= 5;
  return chi2;
}
// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  NDF = 0;
  double logL = 0.0;
  double alpha  = params[0];
  double R_out  = params[1];
  double R_side = params[2];
  double R_long = params[3];
  double N  = params[4];
  for (int i = 0; i < 3; ++i)
  {
    double R = (i == 0) ? R_out : (i == 1) ? R_side : R_long;
    double *par = new double[3];
    par[0] = alpha;
    par[1] = R;
    par[2] = N;
    double integral = histograms[thiskt][thisframe][i]->Integral();
    for (int ibin = 1; ibin <= histograms[thiskt][thisframe][i]->GetNbinsX(); ++ibin)
    {
      double x = histograms[thiskt][thisframe][i]->GetXaxis()->GetBinCenter(ibin);
      if(x > rfitmax) continue;
      if(x < rfitmin) continue;
      double binVolume = ( histograms[thiskt][thisframe][i]->GetXaxis()->GetBinWidth(ibin) );
        
      double observed = histograms[thiskt][thisframe][i]->GetBinContent(ibin);
      double expected = fitFunction(&x, par)*binVolume*integral;
      if(expected <= 0.) continue; // Avoid log(0) or negative expected values
      if(observed != 0.)
        logL += expected + observed*log(observed/expected) - observed;
      NDF++;
    }
    delete par;
  }
  NDF -= 5;
  return 2.0 * logL;
}

int main()
{
  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/star/u/kincses/work/EPOS_project/3D_analysis/levyfit/tables/levy_values_mostprecise_e24.dat");
//  myLevy_reader = new Levy_reader("/star/u/kincses/work/UrQMD/levyfit/tables/levy_values_complex2.dat");
  cout << "levy reader created." << endl;

  TH1* alphahist[NKT];
  TH1* Rohist[NKT];
  TH1* Rshist[NKT];
  TH1* Rlhist[NKT];
  TH1* Ravghist[NKT];
  TH1* conflevhist[NKT];
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    alphahist[ikt] = new TH1F(Form("alphahist_ikt%i",ikt),"",200,0.,2.);
    Rohist[ikt] = new TH1F(Form("Rohist_ikt%i",ikt),"",200,0.,20.);
    Rshist[ikt] = new TH1F(Form("Rshist_ikt%i",ikt),"",200,0.,20.);
    Rlhist[ikt] = new TH1F(Form("Rlhist_ikt%i",ikt),"",200,0.,20.);
    Ravghist[ikt] = new TH1F(Form("Ravghist_ikt%i",ikt),"",200,5.,15.);
    conflevhist[ikt] = new TH1F(Form("clhist_ikt%i",ikt),"",1000,0.,1.);
  }
  // Open the file and get the histograms
//  TFile *file = TFile::Open("../EPOS_3d_source_longrange_evtavg.root");
  TFile *file = TFile::Open("../EPOS_3d_source_010cent_summed_strictQcut.root");
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

//  TH3F* sourcehists[NKT][NFRAME];

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(45);
  Tl.SetNDC(kTRUE);
  TCanvas* canvas = new TCanvas("c1", "", 1200, 600);
  TPad *pad[3];
  pad[0] = new TPad("pad1","",0.,0.,0.36,1.);
  pad[1] = new TPad("pad2","",0.36,0.,0.68,1.);
  pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
  canvas->cd();
  for(int ipad = 0; ipad < 3; ipad++)
    pad[ipad]->Draw();
  for(int ifile = 1; ifile < 2; ifile++)
  {
//    if(!(ifile%20 == 1 || ifile%20 == 2)) continue;
    cout << "ifile: " << ifile << endl;
    for(int ievt = 0; ievt < 1; ievt++)
    {
      for (int ikt = 0; ikt < NKT; ++ikt)
      {
        thiskt = ikt;
//        cout << "ikt = " << thiskt << endl;
        for(int iframe = 0; iframe < 3; iframe++)
        {
          if(iframe != thisframe) continue;
//          cout << "iframe = " << thisframe << endl;
          if (histograms[ikt][iframe].size() == 3) 
          {
//            cout << "histograms size is 3, deleting." << endl;
            for(int i = 0; i < 3; i++) delete histograms[ikt][iframe][i];
            histograms[ikt][iframe].clear();
          }
          for(int iosl = 0; iosl < 3; iosl++)
          {
            // Form the histogram name
            TString histName = Form("pion_pair_source_avg_%s_ikt%i_%s", frames[iframe], ikt, osl_labels[iosl]);
            // Read the histogram
            if (!(file->Get(histName)))
              continue;
            // Print some information about the histogram
//            std::cout << "Processed Histogram " << histName << " with " << sourcehist->GetEntries() << " entries." << std::endl;
//            sourcehist->Scale(1.0 / sourcehist->Integral(1,sourcehist->GetNbinsX()));
//            cout << "Integral w overflow after normalization: " << sourcehist->Integral(0,sourcehist->GetNbinsX()+1) << endl;
            histograms[ikt][iframe].push_back(dynamic_cast<TH1F*>(file->Get(histName)));
          }
//          cout << "histograms read. " << endl;
    
          if (histograms[ikt][iframe].size() != 3)
          {
//            std::cerr << "Error: Expected 3 histograms" << std::endl;
            continue;
          }
    
          // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        
          // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
        
          // Create the function to be minimized
          ROOT::Math::Functor functor(&logLikelihood, 5);
//          ROOT::Math::Functor functor(&chiSquare, 5);
          minimizer->SetFunction(functor);
        
          // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.6, 0.01, 0.5, 2.0);
//          minimizer->SetFixedVariable(0, "alpha", 1.3);
          minimizer->SetLimitedVariable(1, "R_out", 8.5, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(2, "R_side", 8.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(3, "R_long", 9.0, 0.01, 0., 20.);
//          minimizer->SetFixedVariable(3, "R_long", 20.0);
//          minimizer->SetFixedVariable(4, "N_out", 1.);
//          minimizer->SetFixedVariable(5, "N_side", 1.);
//          minimizer->SetFixedVariable(6, "N_long", 1.);
          minimizer->SetVariable(4, "N", 1., 0.01);
        
          // Minimize the function
          minimizer->Minimize();
        
//          cout << "minimizing done. " << endl;
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          double conflev = TMath::Prob(chi2val, NDF);
//          std::cout << "alpha: "  << results[0] << std::endl;
//          std::cout << "R_out: "  << results[1] << std::endl;
//          std::cout << "R_side: " << results[2] << std::endl;
//          std::cout << "R_long: " << results[3] << std::endl;
//          std::cout << "N "       << results[4] << std::endl;
//          std::cout << "chi2: "   << chi2val << std::endl;
          conflevhist[ikt]->Fill(conflev);
          if(conflev > 0.0001 && results[1] > 1. && results[2] > 1. && results[3] > 1.)
          {
            alphahist[ikt]->Fill(results[0]);
            Rohist[ikt]->Fill(results[1]);
            Rshist[ikt]->Fill(results[2]);
            Rlhist[ikt]->Fill(results[3]);
            Ravghist[ikt]->Fill(sqrt( (results[1]*results[1] + results[2]*results[2] + results[3]*results[3]) / 3. ));
          }
//          cout << "histograms filled." << endl;             
          if(donotdraw) {delete minimizer; continue;}
          for(int iosl = 0; iosl < 3; iosl++)
          {
            canvas->cd();
//            cout << "iosl: " << iosl << endl;
            int ipad = iosl;
            gStyle->SetLabelSize(0.06,"Y");
            canvas->SetLogx(1);
            canvas->SetLogy(1);
//            cout << "about to cd into ipad" << endl;
            pad[ipad]->cd();
            gPad->SetRightMargin(rightmargins[ipad]);
            gPad->SetLeftMargin(leftmargins[ipad]);
            gPad->SetTopMargin(topmargins[ipad]);
            pad[ipad]->SetBottomMargin(bottommargins[ipad]);
            gPad->SetLogx(1);
            gPad->SetLogy(1);

//            cout << "about to scale histograms." << endl;
            double integral = histograms[ikt][iframe][iosl]->Integral();//0,histograms[ikt][iframe][iosl]->GetNbinsX()+1);
            histograms[ikt][iframe][iosl]->Scale(1.0 / integral);
            for (int x = 1; x <= histograms[ikt][iframe][iosl]->GetNbinsX(); ++x)
            {
              double content = histograms[ikt][iframe][iosl]->GetBinContent(x);
              double error = histograms[ikt][iframe][iosl]->GetBinError(x);
              double binVolume = histograms[ikt][iframe][iosl]->GetXaxis()->GetBinWidth(x);
              histograms[ikt][iframe][iosl]->SetBinContent(x, content / binVolume);
              histograms[ikt][iframe][iosl]->SetBinError(x, error / binVolume);
            }
//            cout << "histograms scaled." << endl;
            histograms[ikt][iframe][iosl]->SetStats(0);
            histograms[ikt][iframe][iosl]->SetTitle("");
            histograms[ikt][iframe][iosl]->GetXaxis()->SetRangeUser(0.2,9.e5);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleSize(0.1);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleSize(0.1);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetLabelSize(0.08);
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitle("D(#rho)");
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleOffset(0.9);
            if(iframe == 0 && iosl == 2) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitle("#rho [fm]");
            if(iframe == 0 && iosl == 2) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleOffset(0.7);
            histograms[ikt][iframe][iosl]->SetMinimum(5.e-11);
            histograms[ikt][iframe][iosl]->SetMaximum(3.);
            TF1* f_levyfunc = new TF1(Form("levyfunc%i%i",iframe,iosl), fitFunction, 0.1, 1.e6, 3);
            f_levyfunc->SetParNames("alpha","R","norm");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->SetNpx(10000);
            double alpha = results[0];
            double dalpha = errors[0];
            double Rosl = results[1+iosl];
            double dRosl = errors[1+iosl];
            double N = results[4];
            double dN = errors[4];
            f_levyfunc->SetParameters(results[0],Rosl,N);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelSize(0.08);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelOffset(-0.03);
            histograms[ikt][iframe][iosl]->Draw("pe");
            f_levyfunc->SetLineStyle(2);
            f_levyfunc->DrawCopy("same");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->SetRange(rfitmin,rfitmax);
            f_levyfunc->DrawCopy("same");
            Tl.SetTextSize(30);
            if(iframe == 0 && iosl == 0) Tl.DrawLatex(0.45, 0.87, "OUT");
            if(iframe == 0 && iosl == 1) Tl.DrawLatex(0.45, 0.87, "SIDE");
            if(iframe == 0 && iosl == 2) Tl.DrawLatex(0.45, 0.87, "LONG");
//            if(iframe == 0 && iosl == 0) Tl.DrawLatex(0.25, 0.50, "LAB");
//            if(iframe == 1 && iosl == 0) Tl.DrawLatex(0.25, 0.50, "LCMS");
//            if(iframe == 2 && iosl == 0) Tl.DrawLatex(0.25, 0.45, "PCMS");
            Tl.SetTextSize(20);
//            Tl.DrawLatex(0.58-0.05*iosl, 0.78+iframe*0.1, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.73+iframe*0.1, Form("conf.lev. = %1.5f", conflev));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.68+iframe*0.1, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.63+iframe*0.1, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.58+iframe*0.1, Form("N = %1.2f #pm %1.2f",N, dN));
            Tl.DrawLatex(0.58-0.05*iosl, 0.78, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
            Tl.DrawLatex(0.58-0.05*iosl, 0.73, Form("conf.lev. = %1.5f", conflev));
            Tl.DrawLatex(0.58-0.05*iosl, 0.68, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
            Tl.DrawLatex(0.58-0.05*iosl, 0.63, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
            Tl.DrawLatex(0.58-0.05*iosl, 0.58, Form("N = %1.2f #pm %1.2f",N, dN));
          }
        // Clean up
        delete minimizer;
        }
        if(!donotdraw)
        {
          canvas->cd();
          Tl.SetTextSize(30);
          Tl.DrawLatex(0.05, 0.94, Form("EPOS 200 GeV 0#minus10%% Au+Au PION PAIR SOURCE PROJECTIONS, #LTm_{T}#GT = %1.3f",sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
          canvas->SaveAs(Form("%s%s/projections_strictQcut_ikt%i_ich0.pdf", path, frames[thisframe], ikt));
//          cout << "canvas saved." << endl;
//          for(int ipad = 0; ipad < 3; ipad++) pad[ipad]->Clear();
//          cout << "pads cleared." << endl;
          canvas->Clear();
          pad[0] = new TPad("pad1","",0.,0.,0.36,1.);
          pad[1] = new TPad("pad2","",0.36,0.,0.68,1.);
          pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
          for(int ipad = 0; ipad < 3; ipad++) pad[ipad]->Draw();
//          cout << "canvas cleared." << endl;
        }
      }
    }
  }
  cout << "File loop ended." << endl;
  if(canvas) delete canvas;
  file->Close();
  if(file) delete file;
  TFile* file_output = new TFile(Form("fitresults_%s_010cent_avgevt_strictQcut.root",frames[thisframe]), "RECREATE");  
  file_output->cd();
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    cout << ikt << endl;
    alphahist[ikt]->Write();
    Rohist[ikt]->Write();
    Rshist[ikt]->Write();
    Rlhist[ikt]->Write();
    Ravghist[ikt]->Write();
    conflevhist[ikt]->Write();
  }
  file_output->Close();
  delete myLevy_reader;
  return 0;
}

