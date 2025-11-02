// compile & run, e.g. for KT bin index 2: make clean && make all && ./levy_fit.exe 2

#include <iostream>
//#include <string>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "levy_calc.h" // on-the-fly Levy calculation

using namespace std;

TGraphErrors* gr;

double rmin;
double rmax;
const int NPARS = 3;
int NDF;


const int NKTBIN = 10;
//const char* bounds[6] = {"0.0","0.1","0.2","0.3","0.4","0.5"}; // K_T bin bounds
const double kt_limits[NKTBIN+1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};

const char *statuses[6] = {
                         "converged",
                         "cov. made pos.def.",
                         "Hesse invalid",
                         "Edm above max",
                         "call lim. reached",
                         "other failure"
                       };
const char *covstatuses[5] = {
                            "not available",
                            "not pos.def.",
                            "approximate",
                            "forced pos.def.",
                            "accurate"
                          };

// Implemented from eposSpatialDist.cpp
const Double_t Pi = TMath::Pi();
void rdivision(TH1F *h)
{
    Int_t nbins = h->GetXaxis()->GetNbins(); // 168
    h->Sumw2();
    //h->ComputeIntegral();
    //cerr << "Integral = " << h->Integral() << endl;
    double Integral = h->Integral();
    //cerr << "GetEntries pre: = " << h->GetEntries() << endl;
    for (Int_t i = 0; i <= nbins; i++)
    {
        auto BinValue = h->GetBinContent(i);
        auto BinCenter = h->GetBinCenter(i);
        auto BinWidth = h->GetBinWidth(i);
        auto BinError = h->GetBinError(i);
        BinError = sqrt(BinValue); //h->GetBinError(i);
        
        auto Scale = 4. * Pi * BinCenter * BinCenter * BinWidth;
        h->SetBinContent(i, BinValue / Scale);
        h->SetBinError(i, BinError / Scale);
    }
    //h->ComputeIntegral();
    //double Integral = h->Integral();

    //h->Print("all");
    h->Scale(1.0/Integral);
    cerr << "Integral = " << Integral << endl;
}

// Function to convert a histogram to a TGraphErrors
TGraphErrors* ConvertHistToGraph(TH1F* hist) {
  if (!hist) return nullptr;

  int nBins = hist->GetNbinsX();
  std::vector<double> x, y, ex, ey;

  for (int i = 1; i <= nBins; ++i) {
    double binCenter = hist->GetBinCenter(i);
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);
    x.push_back(binCenter);
    y.push_back(binContent);
    ex.push_back(0);       // No x error
    ey.push_back(binError); // y error from histogram
  }

  TGraphErrors* graph = new TGraphErrors(nBins, x.data(), y.data(), ex.data(), ey.data());
  graph->SetName(Form("%s_graph", hist->GetName()));
  graph->SetTitle(hist->GetTitle());

  return graph;
}

// Function to get a TGraphErrors for a specific KT bin
TGraphErrors* GetGraphForKT(const char* filename, int iKT) {
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open file " << filename << "!" << std::endl;
    return nullptr;
  }

  /*
  TIter next(file->GetListOfKeys());
  TKey *key;
  TGraphErrors* graph = nullptr;

  while ((key = (TKey*)next())) {
    TObject *obj = key->ReadObj();
    TH1F *hist = dynamic_cast<TH1F*>(obj);
    if (!hist) continue;

    std::string histname = hist->GetName();
    // Look for the _sqrtrho2 histogram with the correct KT bin
    if (histname.find("_sqrtrho2") != std::string::npos &&
        histname.find(Form("_ikt%d_", iKT)) != std::string::npos) {

      std::cout << "Found histogram: " << histname << std::endl;
      graph = ConvertHistToGraph(hist);
      break;  // Stop after finding the first match
    }
  }
  */

  // Construct the exact histogram name
  std::string histname = Form("pion_pair_source_avg_lcms_ifile1201_ievt0_ikt%d_sqrtrho2", iKT);
  // Get the histogram directly by name
  TH1F* hist = dynamic_cast<TH1F*>(file->Get(histname.c_str()));
  if (!hist) {
    std::cerr << "Error: Histogram " << histname << " not found in file!" << std::endl;
    file->Close();
    return nullptr;
  }
  std::cout << "Found histogram: " << histname << std::endl;

  rdivision(hist); // TODO: normalization on the TGraphErrors or similar (for logistic instead of chi^2 etc.)

  TGraphErrors* graph = ConvertHistToGraph(hist);
  
  file->Close();
  return graph;
}
// ---------------


double FitFunction(const double *x, const double *par)
{
	double r = x[0];
	double N = par[0];
  double R = par[1];
  double alpha = par[2];
  double Rcc = R*pow(2.,1. / alpha);
  return N*levy_calc(r/Rcc,Rcc,alpha);
}

double MyChi2(const double *par)
{
  double chi2 = 0;
  NDF = 0;
  for(int ix=1;ix<gr->GetN();ix++)
//for(int ix=1;ix <= 35; ix++)
  {
    double r = gr->GetX()[ix];
    if(r<rmin || r>rmax) continue;
    double exp = gr->GetY()[ix];
    double theor = FitFunction(&r,par);
    double err = gr->GetEY()[ix];
    if(err==0) continue;
    double chi = (exp-theor)/err;
    chi2 += chi*chi;
    NDF++;
  }
  NDF -= NPARS;
	cerr << ".";
  return chi2;
}

int main(const int argc, char **argv)
{
  cout << "Starting fitting D(r)..." << endl;
  // FIT BOUNDS ÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷÷
  rmin = 1.0;//0.5;
	rmax = 300.0;//7000.0;
	double rminplot = 0.1;
	double rmaxplot = 2000.;

  int iKT = 2;
  const char* energy = "200GeV";
  const char* inputfile = "EPOS_3d_source_010cent_all.root";
  if(argc>1)
  {
    iKT = (int)atoi(argv[1]);
    energy = argv[2];
    if(argc>3)
    {
      inputfile = argv[3];
    }
    if(argc==6)
    {
      rmin = (double)atof(argv[4]);
      rmax = (double)atof(argv[5]);
    }
  }
  
//  gr = new TGraphErrors("Drdata.txt","%lg %lg %lg");
  gr = new TGraphErrors(Form("distance_dists_EPOS_KT%d.txt",iKT),"%lg %lg %lg");
  //gr = GetGraphForKT(inputfile, iKT);
  if(!gr){
    std::cerr << "Error: No matching histogram found for KT bin " << iKT << std::endl;
    return 1;
  }
  std::cout << "Successfully converted histogram to TGraphErrors for KT bin " << iKT << std::endl;

  // Choose method upon creation between:
  // kMigrad, kSimplex, kCombined, kScan, kFumili
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kCombined );

  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.001);

  ROOT::Math::Functor f(&MyChi2,NPARS);

  min.SetFunction(f);

  // Set the free variables to be minimized!
  min.SetVariable(0,"N",      1.0 ,0.01);
  min.SetVariable(1,"R",      8.  ,0.01);
  min.SetVariable(2,"alpha",  1.1 ,0.01);
  //min.FixVariable(0);
	//min.SetVariableLimits(0,0.5,1.1);
	min.SetVariableLimits(1,3,15);
	min.SetVariableLimits(2,0.8,1.8);
  //min.FixVariable(1);
  //double chi2before = MyChi2(min.X());
	//cout << "Probability before: " << chi2before << "/" << NDF << "->" << TMath::Prob(chi2before,NDF) << endl;

  min.Minimize();
  min.PrintResults();
  min.ProvidesError();
  min.Hesse();
  const double *par = min.X();
  const double *err = min.Errors();
  double chi2 = MyChi2(par);
  cout << "Probability: " << chi2 << "/" << NDF << "->" << TMath::Prob(chi2,NDF) << endl;
  cout << "Parameters:" << endl;
  for(int ipar=0;ipar<NPARS;ipar++) cout << "par" << ipar << "=" << par[ipar] << "+-" << err[ipar] << endl;

  cout << "Fit status:" << endl;
  int fitstatus = min.Status();
  int fitcovstatus = min.CovMatrixStatus();
  if(fitstatus<0 || fitstatus>5) fitstatus=5;
  if(fitcovstatus<-1 || fitcovstatus>3) fitcovstatus=-1;
  if(fitstatus == 0 && fitcovstatus == 3) cout << "Fit converged, full accurate cov. matrix";
  else 
  {
    cout << "fit status: " << statuses[fitstatus] << endl;
    cout << "cov. matrix " << covstatuses[fitcovstatus+1] << endl;
  }
  cout << endl;
  cout << "(fitstatus=" << fitstatus << ",covstatus=" << fitcovstatus << ")" << endl;

// innen kiszedni az aszimm hibát, meg a többi paramétert is írja ki coutba

  double errup[NPARS] = {0};
  double errdn[NPARS] = {0};

  const char* paramnames[NPARS] = {"N", "R", "alpha"};
  cout << "Minos errors: " << endl;
  for(unsigned int ipar=0; ipar<NPARS; ipar++)
  {
    min.GetMinosError(ipar,errdn[ipar],errup[ipar]);
    cout << "Asym. err " << paramnames[ipar] << ": +" << errup[ipar] << " -" << errdn[ipar] << endl;
  }

  TF1* fitfunc0 = new TF1("fitfunc0",FitFunction,rminplot,rmaxplot,NPARS);
	fitfunc0->SetLineColor(kRed);
	fitfunc0->SetLineStyle(2);
	fitfunc0->SetLineWidth(1);
	fitfunc0->SetParameters(min.X());
	fitfunc0->SetMinimum(1e-20);
	fitfunc0->SetMaximum(5e3);
	fitfunc0->SetNpx(100);
	fitfunc0->SetTitle(Form("D(r) fit from EPOS4, %s", energy));
  fitfunc0->GetXaxis()->SetTitle("r [fm]");
  fitfunc0->GetYaxis()->SetTitle("D(r)");
  TF1* fitfunc = new TF1("fitfunc",FitFunction,rmin,rmax,NPARS);
	fitfunc->SetLineColor(kRed);
	fitfunc->SetParameters(min.X());
	fitfunc->SetNpx(100);

  TCanvas *c = new TCanvas();
	c->SetLogx(1);
	c->SetLogy(1);
	fitfunc0->Draw();
	fitfunc->Draw("SAME");
	gr->SetMarkerStyle(20);
	gr->Draw("EP");
	TLatex lt;
	lt.SetTextSize(0.05);
	lt.DrawLatexNDC(0.58,0.83,Form(       "R = (%.2f #pm %.2f) fm",par[1],err[1]));
	lt.DrawLatexNDC(0.58,0.78,Form(  "#alpha = %.2f #pm %.2f",     par[2],err[2]));
	lt.DrawLatexNDC(0.58,0.73,Form(       "N = %.2f #pm %.2f",     par[0],err[0]));
	lt.DrawLatexNDC(0.58,0.68,Form( "#chi^{2}/NDF = %.1f / %d",   chi2, NDF));
	lt.DrawLatexNDC(0.58,0.63,Form( "C.L. = %.2f%%",                  100*TMath::Prob(chi2,NDF)));
  lt.DrawLatexNDC(0.58,0.58,Form( "%.3f < K_{T} < %.3f GeV", kt_limits[iKT], kt_limits[iKT+1]));
	c->Print(Form("figs/Dr_levyfit-AuAu%s_KT%d.png",energy,iKT));
  //c->SaveAs(Form("levy_fit-KT_%d.pdf",iKT));

  return 0;
}

