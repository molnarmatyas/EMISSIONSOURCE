// compile & run, e.g. for KT bin index 2: make clean && make all && ./levy_fit.exe 2

#include <iostream>
//#include <string>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "levy_calc.h"

using namespace std;

TGraphErrors* gr;

double rmin;
double rmax;
const int NPARS = 3;
int NDF;

//const char* bounds[6] = {"0.0","0.1","0.2","0.3","0.4","0.5"}; // K_T bin bounds
const int nkT = 10;
const char* bounds[nkT+1] = {"0.175","0.225","0.275","0.325","0.375","0.425","0.475","0.525","0.575","0.625","0.675"};

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
  
  const char* energy = "9.2GeV";
  if(argc>1)
  {
    iKT = (int)atoi(argv[1]);
    energy = argv[2];
    if(argc==4)
    {
      rmin = (double)atof(argv[3]);
      rmax = (double)atof(argv[4]);
    }
  }
//  gr = new TGraphErrors("Drdata.txt","%lg %lg %lg");
  gr = new TGraphErrors(Form("distance_dists_EPOS_KT%d.txt",iKT),"%lg %lg %lg");

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
  lt.DrawLatexNDC(0.58,0.58,Form( "%s < K_{T} < %s GeV", bounds[iKT], bounds[iKT+1]));
	c->Print(Form("figs/Dr_levyfit-AuAu%s_KT%d.png",energy,iKT));
  //c->SaveAs(Form("levy_fit-KT_%d.pdf",iKT));

  return 0;
}

