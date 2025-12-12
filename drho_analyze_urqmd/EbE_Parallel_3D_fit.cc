#include <include/Levy_proj_reader.h>
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

const int NKT = 10;
const double kt_limits[NKT+1]= {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};

// Global variables
const int NFILE = 2;//9986;
const int NEVT = 10000;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 4;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

//const float rightmargins[9] = {0.,0.,0.04,0.,0.,0.04,0.,0.,0.04};
//const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
//const float topmargins[9]  = {0.03,0.03,0.03,0.,0.,0.,0.,0.,0.};
//const float bottommargins[9]  = {0.16,0.16,0.16,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const float topmargins[4] = {0.12, 0.12, 0.05, 0.05};    
const float bottommargins[4] = {0.05, 0.05, 0.15, 0.15}; 
const float leftmargins[4] = {0.15, 0.05, 0.15, 0.05};   
const float rightmargins[4] = {0.05, 0.15, 0.05, 0.15};
const char* frames[3] = {"LCMS","pcms","lab"};
const char* osl_labels[4] = {"out","side","long","outlong"};

const char* path = "plots/";

std::vector<TH1D*> histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

const double rfitmax = 100.;
const double rfitmin = 1.;

const bool donotdraw = false;
const bool writetofile = true;

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
  double R_outlong_2  = params[4];
  double R_outlong_proj = std::sqrt(0.5 * R_out * R_out + 0.5 * R_long * R_long + R_outlong_2);
  double N  = params[5];
  for (int i = 0; i < 4; ++i)
  {
    double R = (i == 0) ? R_out : ((i == 1) ? R_side : ((i == 2) ? R_long : R_outlong_proj));
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
  NDF -= 6;
  return 2.0 * logL;
}

int main(int argc, char** argv)
{
  int ikt;
  if(argc < 2) {
    std::cerr << "Invalid input levy_fit <kt num> <evt avg> <beammom>" << std::endl; 
    return 0;
  }
  ikt = atoi(argv[1]); 
  int evt_avg = atoi(argv[2]); 
	int beamMomentum = std::atoi(argv[3]);

  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/afs/cern.ch/work/b/bporfy/UrQMD/LevyFit/bin_table/levy_values_mostprecise_e24.dat");
  cout << "levy reader created." << endl;

  TH1* alphahist[NKT];
  TH1* Rohist[NKT];
  TH1* Rshist[NKT];
  TH1* Rlhist[NKT];
  TH1* Ravghist[NKT];
  TH1* Rolhist[NKT];
  TH2* alphaRavghist[NKT];
  TH1* conflevhist[NKT];
  /*
  for(int ikt = 0; ikt < NKT; ikt++)
  {
  */
    alphahist[ikt] = new TH1D(Form("alphahist_ikt%i",ikt),"",200,0.,2.);
    alphaRavghist[ikt] = new TH2F(Form("alphaRavghist_ikt%i",ikt),"",200,0.,2.,200,2.,18.);
    Rohist[ikt] = new TH1D(Form("Rohist_ikt%i",ikt),"",200,0.,20.);
    Rshist[ikt] = new TH1D(Form("Rshist_ikt%i",ikt),"",200,0.,20.);
    Rlhist[ikt] = new TH1D(Form("Rlhist_ikt%i",ikt),"",200,0.,20.);
    Ravghist[ikt] = new TH1D(Form("Ravghist_ikt%i",ikt),"",200,5.,15.);
    Rolhist[ikt] = new TH1D(Form("Rolhist_ikt%i",ikt),"",400,-20.,20.);
    conflevhist[ikt] = new TH1D(Form("clhist_ikt%i",ikt),"",1000,0.,1.);
  /*
  }
  */
  // Open the file and get the histograms
  TFile *file = TFile::Open(Form("/afs/cern.ch/work/b/bporfy/UrQMD/out/ArSc%i_qLCMScut_qmaxlim_test.root", beamMomentum));
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

//  TH3F* sourcehists[NKT][NFRAME];

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(45);
  Tl.SetNDC(kTRUE);
  TCanvas* canvas = new TCanvas("c1", "", 1200, 600);
  TPad *pad[4];
  pad[0] = new TPad("pad1","",0.,0.5,0.5,1.);
  pad[1] = new TPad("pad2","",0.5,0.5,1.0,1.);
  pad[2] = new TPad("pad3","",0.0,0.,0.5,0.5);
  pad[3] = new TPad("pad4","",0.5,0.,1.0,0.5);
  canvas->cd();
  TLegend *leg = new TLegend(0.5,0.3,0.7,0.5);
  for(int ipad = 0; ipad < 4; ipad++)
    pad[ipad]->Draw();
  for(int ifile = 1; ifile < NFILE; ifile++)
  {
    for(int ievt = 0; ievt < NEVT; ievt += evt_avg)
    {
      /*
      for (int ikt = 0; ikt < NKT; ++ikt)
      {
      */
        thiskt = ikt;
//        cout << "ikt = " << thiskt << endl;
        for(int iframe = 0; iframe < 3; iframe++)
        {
          if(iframe != thisframe) continue;
//          cout << "iframe = " << thisframe << endl;
          if (histograms[ikt][iframe].size() == 4) 
          {
//            cout << "histograms size is 4, deleting." << endl;
            for(int i = 0; i < 4; i++) delete histograms[ikt][iframe][i];
            histograms[ikt][iframe].clear();
          }
          for(int iosl = 0; iosl < 4; iosl++)
          {
            // Form the histogram name
            TH1D* tempHist_evt_avg;
            for(int ievTemp = ievt; ievTemp < ievt+evt_avg; ievTemp++)
            {
              TString histName = Form("D_%s_%s_ev%i_ch0_KT%i", osl_labels[iosl], frames[iframe], ievTemp, ikt);
              TString histName_pos = Form("D_%s_%s_ev%i_ch1_KT%i", osl_labels[iosl], frames[iframe], ievTemp, ikt);
              // Read the histogram
              if (!(file->Get(histName)))
                continue;
              // Print some information about the histogram
//              std::cout << "Processed Histogram " << histName << " with " << sourcehist->GetEntries() << " entries." << std::endl;
//              sourcehist->Scale(1.0 / sourcehist->Integral(1,sourcehist->GetNbinsX()));
//              cout << "Integral w overflow after normalization: " << sourcehist->Integral(0,sourcehist->GetNbinsX()+1) << endl;
              TH1D* tempHist_neg = (TH1D*)file->Get(histName_pos);
              tempHist_neg->Add((TH1D*)file->Get(histName));
              if(ievTemp%evt_avg == 0) 
              {
                tempHist_evt_avg = (TH1D*)tempHist_neg->Clone();
              }
              else {
                tempHist_evt_avg->Add(tempHist_neg);
              }
            }
            histograms[ikt][iframe].push_back(tempHist_evt_avg);
          }
//          cout << "histograms read. " << endl;
    
          if (histograms[ikt][iframe].size() != 4)
          {
            std::cerr << "Error: Expected 4 histograms" << std::endl;
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
          minimizer->SetLimitedVariable(0, "alpha", 1.3, 0.01, 0.5, 2.0);
//          minimizer->SetFixedVariable(0, "alpha", 1.3);
          minimizer->SetLimitedVariable(1, "R_out", 4.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(2, "R_side", 5.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(3, "R_long", 9.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(4, "R_outlong_2", 5.0, 0.01, -400, 400.);
//          minimizer->SetFixedVariable(3, "R_long", 20.0);
//          minimizer->SetFixedVariable(4, "N_out", 1.);
//          minimizer->SetFixedVariable(5, "N_side", 1.);
//          minimizer->SetFixedVariable(6, "N_long", 1.);
//          minimizer->FixVariable(...);
          minimizer->SetVariable(5, "N", 1., 0.01);
        
          // Minimize the function
          minimizer->Minimize();
        
//          cout << "minimizing done. " << endl;
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          double conflev = TMath::Prob(chi2val, NDF);
          //std::cout << "alpha: "  << results[0] << std::endl;
          //std::cout << "R_out: "  << results[1] << std::endl;
          //std::cout << "R_side: " << results[2] << std::endl;
          //std::cout << "R_long: " << results[3] << std::endl;
          //std::cout << "N "       << results[4] << std::endl;
          //std::cout << "chi2: "   << chi2val << std::endl;
          conflevhist[ikt]->Fill(conflev);
          if(conflev > 0.0001 && results[1] > 1. && results[2] > 1. && results[3] > 1.)
          {
            alphahist[ikt]->Fill(results[0]);
            Rohist[ikt]->Fill(results[1]);
            Rshist[ikt]->Fill(results[2]);
            Rlhist[ikt]->Fill(results[3]);
            double Ravg = sqrt( (results[1]*results[1] + results[2]*results[2] + results[3]*results[3]) / 3. );
            Ravghist[ikt]->Fill(Ravg);
            Rolhist[ikt]->Fill(results[4]);
            alphaRavghist[ikt]->Fill(results[0],Ravg);
          }
//          cout << "histograms filled." << endl;             
          if(donotdraw) {delete minimizer; continue;}
          for(int iosl = 0; iosl < 4; iosl++)
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
            double padwdth = pad[ipad]->GetAbsWNDC();
            double padht = pad[ipad]->GetAbsHNDC();
            double xfactor = 2./(canvas->GetAbsWNDC()/padwdth);
            double yfactor = 2./(canvas->GetAbsHNDC()/padht);
//            double padwdth = pad[ipad]->GetAbsWNDC();
//            double xfactor =  3./(canvas->GetAbsWNDC()/padwdth);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetRangeUser(0.2,200.);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleSize(0.08/xfactor);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleSize(0.1/yfactor);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetLabelSize(0.08);
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitle("D(#rho)");
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleOffset(0.9);
            if(iframe == 0) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitle(Form("#rho_{%s} [fm]",osl_labels[iosl]));
            if(iframe == 0) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleOffset(0.9*xfactor);
            histograms[ikt][iframe][iosl]->SetMinimum(5.e-7);
            histograms[ikt][iframe][iosl]->SetMaximum(3.);
            TF1* f_levyfunc = new TF1(Form("levyfunc%i%i",iframe,iosl), fitFunction, 0.1, 5000., 3);
            f_levyfunc->SetParNames("alpha","R","norm");
            f_levyfunc->SetLineStyle(1);
            double alpha = results[0];
            double dalpha = errors[0];
            double Rosl = results[1+iosl];
            double dRosl = errors[1+iosl];
            double N = results[5];
            double dN = errors[5];
            f_levyfunc->SetParameters(results[0],Rosl,N);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelSize(0.08/xfactor);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelOffset(-0.033/xfactor/xfactor/xfactor);
            histograms[ikt][iframe][iosl]->Draw("pe");


            f_levyfunc->SetLineStyle(2);
            f_levyfunc->DrawCopy("same");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->SetRange(rfitmin,rfitmax);
            f_levyfunc->DrawCopy("same");

            if(iosl == 0) leg->AddEntry(histograms[ikt][iframe][iosl],"D(#rho)","PE");

            Tl.SetTextSize(28);
//            Tl.DrawLatex(0.58-0.05*iosl, 0.78+iframe*0.1, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.73+iframe*0.1, Form("conf.lev. = %1.5f", conflev));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.68+iframe*0.1, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.63+iframe*0.1, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.58+iframe*0.1, Form("N = %1.2f #pm %1.2f",N, dN));
            if(iosl == 0) Tl.DrawLatex((iosl == 0 ? 0.28 : 0.06), 0.44, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
            if(iosl == 0) Tl.DrawLatex((iosl == 0 ? 0.28 : 0.06), 0.38, Form("conf.lev. = %1.5f", conflev));
            if(iosl == 0) Tl.DrawLatex((iosl == 0 ? 0.28 : 0.06), 0.32, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
            Tl.DrawLatex((iosl == 0 ? 0.28 : 0.06), 0.26, Form("R_{%s} = (%1.2f #pm %1.2f) fm", osl_labels[iosl], Rosl, dRosl));
//            Tl.DrawLatex(0.28-0.07*iosl, 0.23, Form("N = %1.2f #pm %1.2f",N, dN));
          }
        // Clean up
        delete minimizer;
        }
        if(!donotdraw)
        {
          canvas->cd();
          Tl.SetTextSize(30);
          Tl.DrawLatex(0.10, 0.91, Form("UrQMD %i event",evt_avg) );
          Tl.DrawLatex(0.10, 0.85, "pair-source projections" );
          Tl.DrawLatex(0.385, 0.91, "0#minus10% ^{40}Ar#plus^{45}Sc" );
          Tl.DrawLatex(0.385, 0.85, "#sqrt{s_{NN}} #approx 17 GeV" );
//          Tl.DrawLatex(0.71, 0.91, "p_{T} [GeV/c] = 0.15#minus1.0");
//          Tl.DrawLatex(0.71, 0.91, "#pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
          Tl.DrawLatex(0.57, 0.91, "#pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
          Tl.DrawLatex(0.71, 0.85, Form("k_{T} [GeV/c] = %1.2f #minus %1.2f",ktbins[ikt],ktbins[ikt+1]) );
//          Tl.DrawLatex(0.71, 0.65, "|#eta| < 1.0"); 
//          Tl.DrawLatex(0.71, 0.59, "Q_{LCMS} < 0.5k_{T}"); 
//          Tl.DrawLatex(0.71, 0.59, "Q_{LCMS} < k_{T}"); 
//          Tl.DrawLatex(0.71, 0.85, Form("#LTm_{T}#GT = %1.3f",sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
          canvas->SaveAs(Form("%s%s/fits/projections_ifile%i_ievt%i-%i_ikt%i_ich0.pdf", path, frames[thisframe], ifile, ievt, ievt+evt_avg, ikt));
//          cout << "canvas saved." << endl;
//          for(int ipad = 0; ipad < 3; ipad++) pad[ipad]->Clear();
//          cout << "pads cleared." << endl;
          canvas->Clear();
          leg->Clear();
          pad[0] = new TPad("pad1","",0.,0.5,0.5,1.);
          pad[1] = new TPad("pad2","",0.5,0.5,1.0,1.);
          pad[2] = new TPad("pad3","",0.0,0.,0.5,0.5);
          pad[3] = new TPad("pad4","",0.5,0.,1.0,0.5);
          for(int ipad = 0; ipad < 4; ipad++) pad[ipad]->Draw();
//          cout << "canvas cleared." << endl;
        }
      /*
      }
      */
    }
  }
  cout << "File loop ended." << endl;
  if(canvas) delete canvas;
  file->Close();
  if(file) delete file;
  if(writetofile)
  {
    TFile* file_output = new TFile(Form("out/fitresults_%s_qLCMS_cut_010cent_all.root",frames[thisframe]), "RECREATE");  
    file_output->cd();
    /*
    for(int ikt = 0; ikt < NKT; ikt++)
    {
    */
      cout << ikt << endl;
      alphahist[ikt]->Write();
      alphaRavghist[ikt]->Write();
      Rohist[ikt]->Write();
      Rshist[ikt]->Write();
      Rlhist[ikt]->Write();
      Ravghist[ikt]->Write();
      Rolhist[ikt]->Write();
      conflevhist[ikt]->Write();
    /*
    }
    */
    file_output->Close();
  }
  delete myLevy_reader;
  return 0;
}

