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
const int NFILE = 1;
const int NEVT = 1;
const int NKT = 9;
const int NFRAME = 1;
const int NCH = 2;
const int NCENT = 4;
const int NPAR = 5;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const double ktbins[NKT + 1] = {0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600};
//const double ktbins[NKT + 1] = {0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,0.575};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.,0.,0.05,0.,0.,0.05,0.,0.,0.05};
const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* osl_labels[3] = {"out","side","long"};

const char* parnames[NPAR] = {"alpha", "Ro", "Rs", "Rl", "N"};

const char* path = "~/work/figs/EPOS_project/3d_proj/lcms/";

TH3* sourcehist[NKT];
Levy_reader* myLevy_reader;

const double rfitmax_o = 100.;
const double rfitmax_s = 100.;
const double rfitmax_l = 100.;
const double rfitmin = 1.;

int thiskt = 0;
int NDF = 0;
// Define the fit function
double fitFunction(double *x, double *par)
{
  double xo = x[0];
  double xs = x[1];
  double xl = x[2];
  double alpha = par[0];
  double Rcc_o = par[1];
  double Rcc_s = par[2];
  double Rcc_l = par[3];
  double N = par[4];
  return (N / Rcc_o / Rcc_s / Rcc_l) * (myLevy_reader->getValue_3d(alpha, sqrt(xo*xo/Rcc_o/Rcc_o + xs*xs/Rcc_s/Rcc_s + xl*xl/Rcc_l/Rcc_l) ));
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  NDF = 0;
  double chi2 = 0.0;

  double integral = sourcehist[thiskt]->Integral();
  for (int ibin = 1; ibin <= sourcehist[thiskt]->GetNbinsX(); ++ibin)
    for (int jbin = 1; jbin <= sourcehist[thiskt]->GetNbinsY(); ++jbin)
      for (int kbin = 1; kbin <= sourcehist[thiskt]->GetNbinsZ(); ++kbin)
      {
        double ro = sourcehist[thiskt]->GetXaxis()->GetBinCenter(ibin);
        double rs = sourcehist[thiskt]->GetYaxis()->GetBinCenter(jbin);
        double rl = sourcehist[thiskt]->GetZaxis()->GetBinCenter(kbin);
        double x[3] = {ro,rs,rl};
        if(ro < rfitmin) continue;
        if(rs < rfitmin) continue;
        if(rl < rfitmin) continue;
        if(sqrt(ro*ro+rs*rs+rl*rl) < 3.*rfitmin) continue;
        double romax = rfitmax_o-thiskt;
        double rsmax = rfitmax_s-thiskt;
        double rlmax = rfitmax_l-thiskt;
        double r_ellipsoid = sqrt(ro*ro/romax/romax+rs*rs/rsmax/rsmax+rl*rl/rlmax/rlmax);
        if(r_ellipsoid > 1.) continue;
        double binVolume = ( sourcehist[thiskt]->GetXaxis()->GetBinWidth(ibin) * sourcehist[thiskt]->GetYaxis()->GetBinWidth(jbin) * sourcehist[thiskt]->GetZaxis()->GetBinWidth(kbin) );
        double error = sourcehist[thiskt]->GetBinError(ibin,jbin,kbin);
        if(error == 0) error = 1.0;
        double observed = sourcehist[thiskt]->GetBinContent(ibin,jbin,kbin);
        double expected = fitFunction(x, const_cast<double*>(params))*binVolume*integral;
        chi2 += pow((observed - expected)/error, 2.);
        NDF++;
      }
  NDF -= 5;
  return chi2;
}

// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  double logL = 0.0;
  double integral = sourcehist[thiskt]->Integral();
  for (int ibin = 1; ibin <= sourcehist[thiskt]->GetNbinsX(); ++ibin)
    for (int jbin = 1; jbin <= sourcehist[thiskt]->GetNbinsY(); ++jbin)
      for (int kbin = 1; kbin <= sourcehist[thiskt]->GetNbinsZ(); ++kbin)
      {
        double ro = sourcehist[thiskt]->GetXaxis()->GetBinCenter(ibin);
        double rs = sourcehist[thiskt]->GetYaxis()->GetBinCenter(jbin);
        double rl = sourcehist[thiskt]->GetZaxis()->GetBinCenter(kbin);
        double x[3] = {ro, rs, rl};
        if (ro < rfitmin) continue;
        if (rs < rfitmin) continue;
        if (rl < rfitmin) continue;
        if(sqrt(ro*ro+rs*rs+rl*rl) < 3.*rfitmin) continue;
        double romax = rfitmax_o-thiskt;
        double rsmax = rfitmax_s-thiskt;
        double rlmax = rfitmax_l-thiskt;
        double r_ellipsoid = sqrt(ro*ro/romax/romax+rs*rs/rsmax/rsmax+rl*rl/rlmax/rlmax);
        if (r_ellipsoid > 1.) continue;
        double binVolume = ( sourcehist[thiskt]->GetXaxis()->GetBinWidth(ibin) * sourcehist[thiskt]->GetYaxis()->GetBinWidth(jbin) * sourcehist[thiskt]->GetZaxis()->GetBinWidth(kbin) );
        
        double observed = sourcehist[thiskt]->GetBinContent(ibin, jbin, kbin);
        double expected = fitFunction(x, const_cast<double*>(params))*binVolume*integral;
        if (expected <= 0) continue; // Avoid log(0) or negative expected values

        if(observed != 0)
          logL += expected + observed*log(observed/expected) - observed;
        
      }
  return 2.0 * logL;
}

int main()
{
  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader("/star/u/kincses/work/UrQMD/levyfit/tables/levy_values_complex2.dat");
  cout << "levy reader created." << endl;

  TGraphErrors* params[NPAR];
  for(int ipar = 0; ipar < NPAR; ipar++) {params[ipar] = new TGraphErrors(NKT); params[ipar]->SetName(parnames[ipar]);}

  // Open the file and get the histograms
  TFile *file = TFile::Open("../EPOS_TH3_source_linbin_ebe.root");
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

  for (int ikt = 0; ikt < NKT; ++ikt)
  {
    thiskt = ikt;
    // Form the histogram name
    TString histName = Form("pion_pair_source_avg_lcms_ifile1_ikt%i", ikt);
    // Read the histogram
    sourcehist[ikt] = dynamic_cast<TH3F*>(file->Get(histName));
    if (!sourcehist[ikt])
    {
      std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
      continue;
    }
    sourcehist[ikt]->Sumw2();
//    std::cout << "Processing Histogram " << histName << " with " << sourcehist[ikt]->GetEntries() << " entries." << std::endl;
//    // Loop over all bins to divide by bin volume
//    for (int xbin = 1; xbin <= sourcehist[ikt]->GetNbinsX(); ++xbin)
//      for (int ybin = 1; ybin <= sourcehist[ikt]->GetNbinsY(); ++ybin)
//        for (int zbin = 1; zbin <= sourcehist[ikt]->GetNbinsZ(); ++zbin)
//        {
//          // Get the bin content and error
//          double content = sourcehist[ikt]->GetBinContent(xbin,ybin,zbin);
//          double error = sourcehist[ikt]->GetBinError(xbin,ybin,zbin);
//          // Calculate the bin volume
//          double binVolume = ( sourcehist[ikt]->GetXaxis()->GetBinWidth(xbin) * sourcehist[ikt]->GetYaxis()->GetBinWidth(ybin) * sourcehist[ikt]->GetZaxis()->GetBinWidth(zbin) );
//          // Update the bin content and error
//          sourcehist[ikt]->SetBinContent(xbin, ybin, zbin, content / binVolume);
//          sourcehist[ikt]->SetBinError(xbin, ybin, zbin, error / binVolume);
//        }
//    double integral = sourcehist[ikt]->Integral();//0, sourcehist[ikt]->GetNbinsX() + 1,
////                                                0, sourcehist[ikt]->GetNbinsY() + 1,
////                                                0, sourcehist[ikt]->GetNbinsZ() + 1);;
//    if (integral == 0)
//    {
//      std::cerr << "Warning: Histogram " << histName << " has zero integral." << std::endl;
//      continue;
//    }
//    sourcehist[ikt]->Scale(1.0 / integral);
    // Print some information about the histogram
    std::cout << "Processed Histogram " << histName << " with " << sourcehist[ikt]->GetEntries() << " entries." << std::endl;
    
    // Create the minimizer
//    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
       
    // Set the minimizer properties
    minimizer->SetMaxFunctionCalls(20000);
    minimizer->SetMaxIterations(20000);
    minimizer->SetTolerance(0.01);
//    minimizer->SetStrategy(2);
    minimizer->SetPrintLevel(1);
    minimizer->SetPrecision(1e-12);        
    // Create the function to be minimized
    ROOT::Math::Functor functor(&logLikelihood, 5);
    minimizer->SetFunction(functor);
        
    // Set the initial parameters and their steps
//    minimizer->SetFixedVariable(0, "alpha", 1.6);
    minimizer->SetVariable(0, "alpha", 1.55, 0.001);
    minimizer->SetVariable(1, "Rcc_out", 15.0-0.5*ikt, 0.001);
    minimizer->SetVariable(2, "Rcc_side", 15.0-0.5*ikt, 0.001);
    minimizer->SetVariable(3, "Rcc_long", 15.0-0.5*ikt, 0.001);
    minimizer->SetVariable(4, "N", 8.2, 0.001);
        
    // Minimize the function
    minimizer->Minimize();
        
    // Print the results
    const double *results = minimizer->X();
    const double *errors = minimizer->Errors();
    const double minval = minimizer->MinValue();
    const double chi2val = chiSquare(results);

//    std::cout << "alpha: "  << results[0] << std::endl;
    std::cout << "R_out: "  << results[1]/pow(2.,1./results[0]) << std::endl;
    std::cout << "R_side: " << results[2]/pow(2.,1./results[0]) << std::endl;
    std::cout << "R_long: " << results[3]/pow(2.,1./results[0]) << std::endl;
//    std::cout << "Rcc_out: "  << results[1] << std::endl;
//    std::cout << "Rcc_side: " << results[2] << std::endl;
//    std::cout << "Rcc_long: " << results[3] << std::endl;
//    std::cout << "N: "      << results[4] << std::endl;
    std::cout << "minval: "   << minval << std::endl;
    std::cout << "chi2val: "  << chi2val << std::endl;
    std::cout << "NDF: "    << NDF << std::endl;
    std::cout << "conflev_logl: " << TMath::Prob(minval, NDF) << endl;
    std::cout << "conflev_chi2: " << TMath::Prob(chi2val, NDF) << endl;

    double mt = sqrt( Mass2_pi + SQR((ktbins[ikt]+ktbins[ikt+1])/2.) );
    for(int ipar = 0; ipar < NPAR; ipar++)
    {
      params[ipar]->SetPoint(ikt, mt, results[ipar]);
      params[ipar]->SetPointError(ikt, 0., errors[ipar]);
    }
    // Clean up
    delete minimizer;
  }
  file->Close();
  delete file;
  TFile* file_output = new TFile("fitresults.root", "RECREATE");  
  file_output->cd();
  for(int ipar = 0; ipar < NPAR; ipar++)
    params[ipar]->Write();
  file_output->Close();
  return 0;
}

