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
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const int linestyles[4] = {1,9,7,1};
const int markerstyles[4] = {20,21,34,24};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* frames[3] = {"lcms","pcms","lab"};
const char* osl_labels[3] = {"out","side","long"};
const char* path = "~/work/figs/corehalotest/";

Levy_reader* myLevy_reader;

double fitFunction(double *x, double *par)
{
  double alphacc = par[0];
  double alphach = par[1];
  double Rcc = par[2];
  double Rch = par[3];
  double lambda = par[4];
  return (2.*lambda/Rcc)*(myLevy_reader->getValue_1d(alphacc, x[0]/Rcc)) + 2.*(2.*sqrt(lambda)*(1.-sqrt(lambda)))*(alphach/2./Rch/TMath::Gamma(1./alphach))*(exp(-pow(x[0]/Rch,alphach)));
//  return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

int main()
{
  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/star/u/kincses/work/EPOS_project/3D_analysis/levyfit/tables/levy_values_mostprecise_e24.dat");
//  myLevy_reader = new Levy_reader("/star/u/kincses/work/UrQMD/levyfit/tables/levy_values_complex2.dat");
  cout << "levy reader created." << endl;

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(45);
  Tl.SetNDC(kTRUE);
  TCanvas* canvas = new TCanvas("c1", "", 800, 800);

  gStyle->SetLabelSize(0.06,"Y");
  canvas->SetLogx(1);
  canvas->SetLogy(1);
  canvas->SetRightMargin(0.05);
  canvas->SetLeftMargin(0.15);
  canvas->SetTopMargin(0.05);
  canvas->SetBottomMargin(0.15);

  TF1* f_levyfunc = new TF1("levyfunc", fitFunction, 0.1, 100000., 5);
  f_levyfunc->SetParNames("alphacc","alphach","Rcc","Rch","lambda");
  f_levyfunc->SetLineStyle(1);
  double alphacc = 1.60;
  double alphach = 0.5;
  double Rcc = 7.*pow(2.,1./alphacc);
  double Rch = 5.e3;
  double lambda = 0.93;
  f_levyfunc->SetParameters(alphacc,alphach,Rcc,Rch,lambda);
  f_levyfunc->SetNpx(10000);
  f_levyfunc->SetLineStyle(1);
  f_levyfunc->SetRange(0.2,9.e5);
  f_levyfunc->GetYaxis()->SetRangeUser(5.e-11,3.);
  f_levyfunc->SetTitle("");
  f_levyfunc->Draw("");
//    Tl.DrawLatex(0.58-0.05*iosl, 0.68, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
//    Tl.DrawLatex(0.58-0.05*iosl, 0.63, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
//    Tl.DrawLatex(0.58-0.05*iosl, 0.58, Form("N = %1.2f #pm %1.2f",N, dN));
  Tl.SetTextSize(30);
  canvas->SaveAs(Form("%s/corehalotest.pdf", path));
  canvas->Clear();
  if(canvas) delete canvas;
  delete myLevy_reader;
  return 0;
}

