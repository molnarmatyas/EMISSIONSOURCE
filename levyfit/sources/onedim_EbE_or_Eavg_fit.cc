// Compile: make exe/onedim_EbE_or_Eavg_fit.exe
// Run also from path/levyfit: exe/onedim_Ebe_or_Eavg_fit.exe **args, e.g. exe/onedim_EbE_or_Eavg_fit.exe 11 "27" 1 10000 1 1
// "path": parent directory of levyfit, analysed, figs directories
// make sure there is a results directory within levyfit dir!
// make sure to move .root input files to "analysed" directory (if not already there)!

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
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <TLegend.h>
#include <TString.h>
#include <TPad.h>
#include <TMath.h>
#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#include <TLatex.h>
#include "../../header_for_all_emissionsource.h"

// Note that the header also includes is3Dfit bool flag that decides 
// if this code will perform 1D or 3D fit.
// Several variables, commented out here, are defined in the header.

using namespace std;

// Global variables
const bool donotdraw = false;
//bool is3Dfit = false;
//const double Mass2_pi = 0.019479835;
//const int NKT = 10;
//const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
// const int NFRAME = 3;
// const int NOSL = 3;
// const int NCH = 2;
// const int NCENT = 10; // number of centrality classes
// const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};
//const int colors[NKT] = {632,416,600,432,616,400,8,9,42,46};
int NPARS = 3; // alpha, R, N, default for 1D fit
bool ikt_plotted[NKT] = {false}; // Track which ikt values have been plotted
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.,0.,0.05,0.,0.,0.05,0.,0.,0.05};
const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* frames[3] = {"lcms","pcms","lab"};
const char* osl_labels[3] = {"out","side","long"};

// !!! change to actual "path"
const char* path = "..";

//TH1* histograms[NKT][NFRAME][NOSL+1]; // +1 for rho proj histograms
TH1* histograms[NKT][NOSL+1]; // +1 for rho proj histograms
Levy_reader* myLevy_reader;

double rfitmax = rfitmax_systlimits[0];//60.; // TODO make adjustable for kT/cent/... & incr/decr if fit does not converge
double rfitmin = 3.;//1.;//5.;
// const double B[3] = {2500, 1600, 3600}; // for rho_fitmax limits: default, strict, loose
// const double rfitmax_systlimits[3] = {100., 50., 150.}; // for simpler rho_fitmax limits: default, strict, loose

int thiskt = 0;
int thisframe = 0;
int NDF = 0;

const char* statuses[6] = {"converged",
                           "cov._made_pos.def.",
                           "Hesse_invalid",
                           "Edm_above_max",
                           "call_lim._reached",
                           "other_failure"};
const char* covstatuses[4] = {"not_calculated",
                              "approximated",
                              "forced_pos.def.",
                              "accurate"};

// Define the fit function
double fitFunction(const double *x, const double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  if(is3Dfit)
    return (2.*N/Rcc)*(myLevy_reader->getValue_1d(alpha, x[0]/Rcc));
  else
    // counter-intuitively, the getValue_3d is used for 1D projections - and for 3D projections the getValue_1d
    //return (N/Rcc/Rcc/Rcc)*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
    return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc)); // should be faster with 1 division + 2 multiplications than 1 multiplication + 3 divisions
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  // params are alpha, R, N (1D) or alpha, Rout, Rside, Rlong, N (3D)
  NDF = 0;
  double chi2 = 0.0;

  int begin_i=0; // rho
  if(is3Dfit) begin_i = 1;; // rho_out, then side, long

  for(int i=0; i<NPARS-NOSL+1; i++) // NPARS-NOSL+1 = 3 for 3D fit, 1 for 1D fit
  {
    double *_params = new double[NOSL]; // alpha, R_i, N
    _params[0] = params[0]; // alpha
    _params[1] = params[i+1]; // R_i
    _params[2] = params[NPARS-1]; // N

    double integral = histograms[thiskt][i+begin_i]->Integral();
    for (int ibin = 1; ibin <= histograms[thiskt][i+begin_i]->GetNbinsX(); ++ibin)
    {
      double x = histograms[thiskt][i+begin_i]->GetXaxis()->GetBinCenter(ibin);
      double ex = histograms[thiskt][i+begin_i]->GetBinError(ibin);
      if(x > rfitmax) continue;
      if(x < rfitmin) continue;
      double binVolume = ( histograms[thiskt][i+begin_i]->GetXaxis()->GetBinWidth(ibin) );
      
      double observed = histograms[thiskt][i+begin_i]->GetBinContent(ibin);
      double expected = fitFunction(&x, _params)*binVolume*integral;
      if(!is3Dfit)
        expected *= (x*x*4.*M_PI); // for 1D fit, include the Jacobian here
      if (ex > 0)
      {
        chi2 += pow((observed - expected)/ex, 2.);
      }
      NDF++;
    }
    delete[] _params;
  }
  NDF -= NPARS;
  return chi2;
}
// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  // params are alpha, R, N (1D) or alpha, Rout, Rside, Rlong, N (3D)
  NDF = 0;
  double logL = 0.0;

  int begin_i=0; // rho
  if(is3Dfit)
  {
    begin_i = 1; // rho_out, then side, long
  }
  
  for(int i=0; i<NPARS-NOSL+1; i++) // NPARS-NOSL+1 = 3 for 3D fit, 1 for 1D fit
  {
    double *_params = new double[NOSL]; // alpha, R_i, N
    _params[0] = params[0]; // alpha
    _params[1] = params[i+1]; // R_i
    _params[2] = params[NPARS-1]; // N

    double integral = histograms[thiskt][i+begin_i]->Integral();
    for (int ibin = 1; ibin <= histograms[thiskt][i+begin_i]->GetNbinsX(); ++ibin)
    {
      double x = histograms[thiskt][i+begin_i]->GetXaxis()->GetBinCenter(ibin);
      if(x > rfitmax) continue;
      if(x < rfitmin) continue;
      double binVolume = ( histograms[thiskt][i+begin_i]->GetXaxis()->GetBinWidth(ibin) );
        
      double observed = histograms[thiskt][i+begin_i]->GetBinContent(ibin);
      double expected = fitFunction(&x, _params)*binVolume*integral;
      if(!is3Dfit)
        expected *= (x*x*4.*M_PI); // for 1D fit, include the Jacobian here
      if(expected <= 0.) continue; // Avoid log(0) or negative expected values
      if(observed != 0.)
        logL += expected + observed*log(observed/expected) - observed;
      NDF++;
    }
    delete[] _params;
  }
  NDF -= NPARS;
  return 2.0 * logL;
}

void Drho_from_rhohist(TH1* hist, int iosl)
{
  // Create D(rho) from rho histogram
  double integral = hist->Integral();//0,histograms[ikt][iosl]->GetNbinsX()+1);
  for (int x = 1; x <= hist->GetNbinsX(); ++x)
  {
    double content = hist->GetBinContent(x);
    double error = hist->GetBinError(x);
    double binVolume = hist->GetXaxis()->GetBinWidth(x);
    hist->SetBinContent(x, content / binVolume);
    hist->SetBinError(x, error / binVolume);
  }
  TF1* f_r2 = nullptr;
  if(!is3Dfit)
  {
    f_r2 = new TF1("f_r2","1./(x*x*4.*pi)",0.,1.e8); // Undo Jacobian from 1D fit
    hist->Multiply(f_r2,1.);
    delete f_r2;
  }
  hist->Scale(1.0 / integral);

  hist->SetStats(0);
  hist->SetTitle("");
  //hist->SetTitle(Form("EPOS4 pion-pion pair source, k_{T} #in [%.2f, %.2f] GeV/c, #sqrt{s_{NN}} = %s, %s%%", ktbins[ikt], ktbins[ikt + 1], energy, centleg[ICENT]));
  hist->GetXaxis()->SetRangeUser(0.2,5000.);//0.1,1000.); // Compare with Yan
  //hist->GetXaxis()->SetTitleSize(0.1);
  //hist->GetYaxis()->SetTitleSize(0.1);
  //hist->GetYaxis()->SetLabelSize(0.08);
  hist->GetYaxis()->SetTitle("D(#rho)");
  //hist->GetYaxis()->SetTitleOffset(0.9);
  const char* rho_label = is3Dfit ? Form("%s",osl_labels[iosl-1]) : "1D proj";
  hist->GetXaxis()->SetTitle(Form("#rho_{%s} [fm]", rho_label));
  //hist->GetXaxis()->SetTitleOffset(0.7);
  hist->SetMinimum(5.e-13); // -12 FIXED
  hist->SetMaximum(1.e-2);
  
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitleOffset(1.2);
}

int main(int argc, char *argv[])
{
  // Changes only for 3D fit
  if(is3Dfit) NPARS = 5; // alpha, Rout, Rside, Rlong, N; else alpha, R, N

  // Default args
  int ICENT=-1;
  const char* energy = "200GeV";
  int NFILEMAX = 10;//2;
  int NEVT = 25;//1000;
  int NEVT_AVG = 100; // number of events to be averaged over for each fit
  bool IsUrQMD = false;
  int qlcms_syst = 0; // 0: default, 1: strict, 2: loose
  int rho_fitmax_syst = 0; // 0: default, 1: strict, 2: loose
  // Args overridden
  if(argc < 2 || argc > 9)
  {
    cerr << "Correct usage: " << endl;
    cerr << "exe/programname.exe ICENT energy NFILEMAX NEVT NEVT_AVG qlcms rho_fitmax IsUrQMD[bool:true/false]" << endl;
    cerr << "(NEVT is number of events per 'file', NEVT_AVG is number of events over to be averaged)" << endl;
    cerr << "qlcms and rho_fitmax of systematics, default: 0, strict: 1, loose: 2" << endl;
    return 1;
  }
  if(argc > 1)
  {
    ICENT = (int)atoi(argv[1]);
  }
  if(argc > 2)
  {
    energy = argv[2];
  }
  if(argc > 3)
  {
    NFILEMAX = (int)atoi(argv[3]);
  }
  if(argc > 4)
  {
    NEVT = (int)atoi(argv[4]);
  }
  if(argc > 5)
  {
    NEVT_AVG = (int)atoi(argv[5]);
    if(NEVT_AVG < 1)
    {
      cerr << "NEVT_AVG must be at least 1!" << endl;
      return 1;
    }
    qlcms_syst = (int)atoi(argv[6]);
    rho_fitmax_syst = (int)atoi(argv[7]);
  }
  if(argc == 9)
  {
    IsUrQMD = static_cast<bool>(atoi(argv[8]));
  }

  // Correct indexing of centrality classes' array - implemented from lcmsonly_epos_pair_source.C macro 
  if(ICENT==NCENT || ICENT < 0)
  {
    ICENT = NCENT; // centleg element: "all"
    cout << "WARNING: all centrality classes together." << endl;
  }else if(ICENT > NCENT) // 0-10% together
  {
    ICENT = NCENT + 1;
    cout << "WARNING: 0-10 percent centrality together." << endl;
  }else{
    cout << centleg[ICENT] << " percent centrality class to be fitted." << endl;
  }

  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader(Form("%s/levyfit/tables/levy_values_moreprecise.dat",path));
  cout << "levy reader created." << endl;

  TH1* alphahist[NKT];
  TH1* Rhist[NKT];
  TH1* R_out_hist[NKT];
  TH1* R_side_hist[NKT];
  TH1* R_long_hist[NKT];
  TH1* Nhist[NKT];
  TH2* alpha_vs_R[NKT];
  // widened ranges to avoid under/overflow for low-energy fits
  TH2F* alpha_vs_R_all = new TH2F("alpha_vs_R_all","",100,0.0,2.0,100,2.5,15.);
  // Per-ikt confidence histograms and one collecting all kT bins
  TH1* confidencehist[NKT];
  TH1* confidencehist_all = new TH1F("confidencehist_all","",100,0.,1.);
  // Temporary storages
  // For writing the vectors into TGraphAsymmErrors
  Double_t ktbin_centers[NKT], xLow[NKT], xHigh[NKT];//, binWidths[NKT];
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    ktbin_centers[ikt] = (ktbins[ikt] + ktbins[ikt + 1]) / 2.;
    //binWidths[ikt] = ktbins[ikt + 1] - ktbins[ikt];
    double binWidth = ktbins[ikt + 1] - ktbins[ikt];
    xLow[ikt] = binWidth / 2.;
    xHigh[ikt] = binWidth / 2.;
  }
  // For collecting from temporary storages
  int Ngoodfits = 0;
  std::vector<TGraphAsymmErrors*> alpha_vs_KT_all;
  std::vector<TGraphAsymmErrors*> N_vs_KT_all;
  std::vector<TGraphAsymmErrors*> R_vs_KT_all;
  std::vector<TGraphAsymmErrors*> R_out_vs_KT_all;
  std::vector<TGraphAsymmErrors*> R_side_vs_KT_all;
  std::vector<TGraphAsymmErrors*> R_long_vs_KT_all;
  //vector<double> chi2_vec;
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    alphahist[ikt] = new TH1F(Form("alphahist_ikt%i",ikt),"",100,0.,2.);
    Rhist[ikt] = new TH1F(Form("Rhist_ikt%i",ikt),"",100,2.5,15.);
    R_out_hist[ikt] = new TH1F(Form("R_out_hist_ikt%i",ikt),"",100,2.5,15.);
    R_side_hist[ikt] = new TH1F(Form("R_side_hist_ikt%i",ikt),"",100,2.5,15.);
    R_long_hist[ikt] = new TH1F(Form("R_long_hist_ikt%i",ikt),"",100,2.5,15.);
    Nhist[ikt] = new TH1F(Form("Nhist_ikt%i",ikt),"",100,0.85,1.15);
    // per-ikt src range: increase Y max to cover larger R at low energies
    alpha_vs_R[ikt] = new TH2F(Form("alpha_vs_R_ikt%i",ikt),"",100,0.,2.,100,2.5,15.);
    confidencehist[ikt] = new TH1F(Form("confidencehist_ikt%i",ikt),"",100,0.,1.);
  }
  // Initialize histogram pointers to nullptr to avoid accidental dereference
  for(int iosl = 0; iosl < NOSL + 1; iosl++) // +1 for rho proj histograms
  {
    for (int ikt = 0; ikt < NKT; ++ikt) {
      for (int iframe = 0; iframe < NFRAME; ++iframe) {
        histograms[ikt][iosl] = nullptr;
      }
    }
  }

  const char* isPathUrqmd = IsUrQMD ? "UrQMD" : "EPOS";
  const char* qlcms_syst_label = (qlcms_syst == 0) ? "" : (qlcms_syst == 1) ? "_strictqLCMS" : "_looseqLCMS"; // in sh script there is _defaultqLCMS created, but renamed back to no suffix
 
  //TH3F* sourcehists[NKT][NFRAME];
  TCanvas* canvas; TLatex Tl; TPad *pad[NOSL]; TLegend * leg;
  Tl.SetNDC(kTRUE);

  if(!is3Dfit)
  {
    canvas = new TCanvas("c1", "", 600, 600);
    Tl.SetTextSize(0.04);//Tl.SetTextFont(43); Tl.SetTextSize(15);
    // Canvas setup
    gStyle->SetLabelSize(0.04, "XY"); // Set label size for both axes
    gStyle->SetTitleSize(0.05, "XY"); // Set title size for both axes
    gStyle->SetTitleOffset(1.2, "X"); // Adjust X-axis title offset
    gStyle->SetTitleOffset(1.5, "Y"); // Adjust Y-axis title offset
    
    canvas->SetLogx(1);
    canvas->SetLogy(1);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.15);
    canvas->SetTopMargin(0.10);
    canvas->SetBottomMargin(0.15);
  }
  else
  {
    canvas = new TCanvas("c1", "", 1200, 600);
    Tl.SetTextFont(43); Tl.SetTextSize(35);
    leg = new TLegend(0.5,0.3,0.7,0.5);
    pad[0] = new TPad("pad1","",0.,0.,0.37,1.);
  	pad[1] = new TPad("pad2","",0.37,0.,0.68,1.);
  	pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
    canvas->cd();
    for(int ipad = 0; ipad < NOSL; ipad++)
      pad[ipad]->Draw();
  }

  // Initialize parameter arrays OUTSIDE ievt loop so they accumulate across all events within an averaging block
  Double_t alpha_vec[NKT]={0};
  Double_t alpha_errup_vec[NKT]={0};
  Double_t alpha_errdn_vec[NKT]={0};
  Double_t R_vec[NKT]={0};
  Double_t R_errup_vec[NKT]={0};
  Double_t R_errdn_vec[NKT]={0};
  Double_t R_out_vec[NKT]={0};
  Double_t R_out_errup_vec[NKT]={0};
  Double_t R_out_errdn_vec[NKT]={0};
  Double_t R_side_vec[NKT]={0};
  Double_t R_side_errup_vec[NKT]={0};
  Double_t R_side_errdn_vec[NKT]={0};
  Double_t R_long_vec[NKT]={0};
  Double_t R_long_errup_vec[NKT]={0};
  Double_t R_long_errdn_vec[NKT]={0};
  Double_t N_vec[NKT]={0};
  Double_t N_errup_vec[NKT]={0};
  Double_t N_errdn_vec[NKT]={0};

  for(int iframe = 0; iframe < 3; iframe++)
  {
    if(iframe != thisframe) continue;

    for(int ifile = 1; ifile < NFILEMAX + 1; ifile++) // files were indexed from 1
    {
      // Open the file and get the histograms
      // FIXME no file numbering, needed for EPOS (soon for UrQMD as well)
      TFile *file = TFile::Open(Form("%s/analysed/%s_3d_source_%scent_all_%s%s.root", 
                                      path,isPathUrqmd,centleg[ICENT],energy,qlcms_syst_label));
      if (!file)
      {
        std::cerr << "Error opening file" << std::endl;
        return 1;
      }
      for(int ievt = 0; ievt < NEVT; ievt++)
      {
        bool averagingComplete = false; // NEW: track if this is the last event of averaging block
        // global (across files) event index used for averaging blocks that span files
        int globalEventIndex = (ifile - 1) * NEVT + ievt;
        
        for (int ikt = 0; ikt < NKT; ++ikt)
        {
          cout << "ifile,ievt,ikt: " << ifile << "," << ievt << "," << ikt << "," << endl;
          thiskt = ikt;
          
          TH1* temp_rhohist[4] = {nullptr, nullptr, nullptr, nullptr}; // sqrtrho2, Rout, Rside, Rlong

          if(IsUrQMD)
          {
            // Add both charges together...
            // pi- pi-
            // -------
            // Form the histogram name
            TString histName = Form("D_%s_ev%d_ch%d_KT%d", frames[iframe], ievt, 0, ikt);
            // Read the histogram
            temp_rhohist[0] = dynamic_cast<TH1D*>(file->Get(histName));
            // Check if histogram was successfully retrieved
            if (!temp_rhohist[0])
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }
            if(is3Dfit)
            {
              for(int iosl = 0; iosl < NOSL; iosl++)
              {
                histName = Form("D_%s_%s_ev%d_ch%d_KT%d", osl_labels[iosl], frames[iframe], ievt, 0, ikt);
                temp_rhohist[iosl+1] = dynamic_cast<TH1D*>(file->Get(histName));
              }
            }

            // pi+ pi+
            // -------
            // Form the histogram name
            histName = Form("D_%s_ev%d_ch%d_KT%d", frames[iframe], ievt, 1, ikt);
            // Read the histogram
            temp_rhohist[0]->Add(dynamic_cast<TH1D*>(file->Get(histName)));
            if (file->Get(histName) == nullptr)
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }
            if(is3Dfit)
            {
              for(int iosl = 0; iosl < NOSL; iosl++)
              {
                histName = Form("D_%s_%s_ev%d_ch%d_KT%d", osl_labels[iosl], frames[iframe], ievt, 1, ikt);
                temp_rhohist[iosl+1]->Add(dynamic_cast<TH1D*>(file->Get(histName)));
              }
            }
          }
          else
          {
            // With EPOS analysis code, both charges were already added together...
            // Form the histogram name
            TString histName = Form("pion_pair_source_avg_%s_ifile%i_ievt%i_ikt%i_sqrtrho2", frames[iframe], ifile, ievt, ikt);
            // Read the histogram
            temp_rhohist[0] = dynamic_cast<TH1F*>(file->Get(histName));
          
            // Check if histogram was successfully retrieved
            if (!temp_rhohist[0])
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }
            if(is3Dfit)
            {
              for(int iosl = 0; iosl < NOSL; iosl++)
              {
                histName = Form("pion_pair_source_avg_%s_ifile%i_ievt%i_ikt%i_%s", frames[iframe], ifile, ievt, ikt, osl_labels[iosl]);
                temp_rhohist[iosl+1] = dynamic_cast<TH1F*>(file->Get(histName));
              }
            }
          }

          //  --- AVERAGING logic ---
          // Only do resetting if first one of (NEVT_AVG) averaged. Use a global event
          // index so averaging can span across multiple files when NEVT_AVG > NEVT.
          if((globalEventIndex % NEVT_AVG == 0) || NEVT_AVG == 1)
          {
            for(int iosl = 0; iosl < NOSL + 1; iosl++) // +1 for rho proj
            {
              if(!temp_rhohist[iosl]) continue; // Skip if nullptr (e.g., 3D slices in 1D fit)
              if(histograms[ikt][iosl]) histograms[ikt][iosl]->Reset();
              histograms[ikt][iosl] = (TH1F*)temp_rhohist[iosl]->Clone();
              if(histograms[ikt][iosl]) histograms[ikt][iosl]->SetDirectory(nullptr); // avoid ROOT ownership issues
            }
            if(NEVT_AVG != 1) continue; // if averaging and only first event in avg. block, do not proceed to fitting yet
          }
          else
          {
            for(int iosl = 0; iosl < NOSL + 1; iosl++) // +1 for rho proj
            {
              if(!temp_rhohist[iosl]) continue; // Skip if nullptr (e.g., 3D slices in 1D fit)
              if(!histograms[ikt][iosl]) continue; // Safety: should have been cloned earlier
              histograms[ikt][iosl]->Add(temp_rhohist[iosl]); // add for averaging
            }
            if(globalEventIndex % NEVT_AVG != NEVT_AVG - 1)
            {
              continue; // do until last one of averaging
            }
            averagingComplete = true; // Set when last event of block reached
          }
          
          // Print some information about the histogram
          int entries = histograms[ikt][0]->GetEntries();
          std::cout << "If 1D, processed sqrtrho2 histogram with " << entries << " entries (= pairs, I guess)." << std::endl;
          // TODO based on that, with a global variable, the averaging could stop at a given number of pairs as well (though there will be some fluctuation, as we already can only loop over evts)
          //sourcehist->Scale(1.0 / sourcehist->Integral(1,sourcehist->GetNbinsX()));
          cout << "Its integral w overflow after normalization: " << histograms[ikt][0]->Integral(0,histograms[ikt][0]->GetNbinsX()+1) << endl;
    
          // FITTING PROCEDURE STARTING HERE
          // k_T- (or m_T-) dependent fit rannge
          //rfitmax = sqrt(ktbin_centers[ikt]*ktbin_centers[ikt] + Mass2_pi) * B[qlcms_syst];
          rfitmax = rfitmax_systlimits[rho_fitmax_syst]; // simpler limits, "by-look" better
          
          // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        
          // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
        
          // Create the function to be minimized
          ROOT::Math::Functor functor(&logLikelihood, NPARS);
          //ROOT::Math::Functor functor(&chiSquare, NPARS);
          minimizer->SetFunction(functor);
        
          // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.6, 0.01, 0.5, 2.0);
          //minimizer->SetFixedVariable(0, "alpha", 1.3);
          if(is3Dfit)
          {
            minimizer->SetLimitedVariable(1, "Rout", 4.0, 0.01, 0., 20.);
            minimizer->SetLimitedVariable(2, "Rside", 5.0, 0.01, 0., 20.);
            minimizer->SetLimitedVariable(3, "Rlong", 9.0, 0.01, 0., 20.);
          }
          else
          {
            minimizer->SetLimitedVariable(1, "R", 8.5, 0.01, 0., 15.);
          }
          minimizer->SetVariable(NPARS-1, "N", 1., 0.01); // this likes to become negative, not good
          minimizer->SetLimitedVariable(NPARS-1, "N", 1.0, 0.001, 0.0, 2.0); // FIXME is this even needed?
        
          // Minimize the function
          minimizer->Minimize();
          //minimizer->ProvidesError(); // !!!
          //minimizer->Hesse(); // !!!
                  
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          std::cout << "alpha: "  << results[0] << std::endl;
          if(is3Dfit)
          {
            std::cout << "Rout: "  << results[1] << std::endl;
            std::cout << "Rside: " << results[2] << std::endl;
            std::cout << "Rlong: " << results[3] << std::endl;
          }
          else
          {
            std::cout << "R: "  << results[1] << std::endl;
          }
          std::cout << "N: "      << results[NPARS-1] << std::endl;
          std::cout << "chi2: "   << chi2val << std::endl;

          // FIND OUT IF: BAD FITS
          // Find out hits in a given range vs. number of bins ONLY FOR rho HISTOGRAM
          int nBins = histograms[ikt][0]->GetNbinsX();
          Int_t binsInRange = 0;
          Double_t hitsInRange = 0;
          Double_t xMin=1.; Double_t xMax=100.;
          for(Int_t i=1; i<=nBins; i++)
          {
            Double_t binCenter = histograms[ikt][0]->GetBinCenter(i);
            if(binCenter >= xMin && binCenter <= xMax)
            {
              binsInRange++;
              hitsInRange += histograms[ikt][0]->GetBinContent(i);
            }
          }
          // Fit statuses
          int fitCovStatus = minimizer->CovMatrixStatus();
          int fitStatus = minimizer->Status();
          cout << "CovMatrixStatus: " << fitCovStatus << " (" << covstatuses[fitCovStatus] << ")"<< endl;
          cout << "Status: " << fitStatus << " (" << statuses[fitStatus] << ")"<< endl;
          double confidence = chi2val / NDF;
          double conflev = TMath::Prob(chi2val, NDF); // this gives the p-value TODO check if correct
          cout << "Confidence: " << confidence << endl;
          
          const char* fitQuality = "GOODFIT";
          if(fitCovStatus != 3) fitQuality = covstatuses[fitCovStatus];
          if(fitStatus != 0 && fitStatus != 3) fitQuality = statuses[fitStatus]; // only mark bad when status is neither 0 nor 3
          if(results[0] < 0.55 || results[0] > 1.95) fitQuality = "alpha_out_of_bounds";
          if(results[1] < 0.05 || results[1] > 14.95) fitQuality = "R_out_of_bounds"; // TODO check but probably enough for one rho
          if(results[NPARS-1] < 0.5 || results[NPARS-1] > 1.5) fitQuality = "N_out_of_bounds";
          if(NEVT_AVG==1 && confidence < 0.01) fitQuality = "conf_too_low"; // confidence < 0.01 should not be checked when averaged
          if(static_cast<Double_t>(binsInRange) * 0.5 > hitsInRange) fitQuality = "too_few_hits";
        
          // Here, among others, creating D(rho) from rho hists (normalising to 1)
          for(int iosl = 0; iosl < NOSL+1; iosl++) // +1 for rho proj
          {
            if(!is3Dfit && iosl>0) continue; // for 1D fit, only do the rho proj
            if(is3Dfit && iosl==0) continue; // for 3D fit, skip the rho proj
            Drho_from_rhohist(histograms[ikt][iosl], iosl);
          }

          // COLLECTING fit parameters into arrays for later averaging/saving
          double alpha = results[0];
          double dalpha = errors[0];
          double R=0., dR=0.;
          double R_out=0., dR_out=0., R_side=0., dR_side=0., R_long=0., dR_long=0.;
          if(!is3Dfit)
          {
            R = results[1];
            dR = errors[1];
          }
          else
          {
            R_out = results[1];
            dR_out = errors[1];
            R_side = results[2];
            dR_side = errors[2];
            R_long = results[3];
            dR_long = errors[3];
            R = TMath::Sqrt((R_out*R_out + R_side*R_side + R_long*R_long) / 3.); // sorry, not RMS, "/3" really needed
            if(R > 0.)
              dR = TMath::Sqrt( (R_out*R_out*dR_out*dR_out) + (R_side*R_side*dR_side*dR_side) + (R_long*R_long*dR_long*dR_long) ) / (3.*R);
            else
              dR = 0.;
          }
          double N = results[NPARS-1];
          double dN = errors[NPARS-1];
          // TODO minos errors?
          //double alpha_dn=0., alpha_up=0.;
          //minimizer->GetMinosError(0, alpha_dn, alpha_up);
          // etc. for R and N
          delete minimizer;
          
          // PLOTTING FIT RESULTS
          // --------------------
          if(!donotdraw)
          {
            if(!is3Dfit)
            {
            TF1* f_levyfunc = new TF1(Form("levyfunc%i",iframe), fitFunction, 0.1, 5000., 3);
            f_levyfunc->SetParNames("alpha","R","norm");
            f_levyfunc->SetLineStyle(1);
            
            f_levyfunc->SetParameters(alpha,R,N);
            histograms[ikt][0]->Draw("pe");
            f_levyfunc->SetLineStyle(2);
            f_levyfunc->DrawCopy("same");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->SetRange(rfitmin,rfitmax);
            f_levyfunc->DrawCopy("same");
            delete f_levyfunc;
            
            Tl.SetTextSize(0.04);
            Tl.DrawLatexNDC(0.58, 0.83, Form("R = (%.2f #pm %.2f) fm", R, dR));
            Tl.DrawLatexNDC(0.58, 0.78, Form("#alpha = %.2f #pm %.2f", alpha, dalpha));
            Tl.DrawLatexNDC(0.58, 0.73, Form("N = %.2f #pm %.2f", N, dN));
            Tl.DrawLatexNDC(0.58, 0.68, Form("#chi^{2}/NDF = %.1f / %d", chi2val, NDF));
            Tl.DrawLatexNDC(0.58, 0.63, Form("C.L. = %.2f%%", 100 * conflev));
            
            // Add a title using TLatex for better customization
            const char* isPathUrqmdTitle = IsUrQMD ? "UrQMD" : "EPOS4";
            TLatex title;
            title.SetTextAlign(12);  //centered
            title.SetTextSize(0.03);
            title.SetTextFont(42);
            title.SetNDC(true);
            title.DrawLatex(0.04, 0.95, Form("%s pion-pion pair source, k_{T}#in[%.3f, %.3f] GeV/c,#sqrt{s_{NN}} = %s, %s%%",
                                  isPathUrqmdTitle, ktbins[ikt], ktbins[ikt + 1], energy, centleg[ICENT]));
            } // end if 1D
            else
            {
              // 3D fit plotting
              for(int iosl = 0; iosl < NOSL; iosl++)
              {
                canvas->cd();
                int ipad = iosl;
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

                histograms[ikt][iosl+1]->SetStats(0);
                histograms[ikt][iosl+1]->SetTitle("");
                double padwdth = pad[ipad]->GetAbsWNDC();
                double xfactor =  3./(canvas->GetAbsWNDC()/padwdth);
                histograms[ikt][iosl+1]->GetXaxis()->SetRangeUser(0.2,200.);
                histograms[ikt][iosl+1]->GetXaxis()->SetTitleSize(0.08/xfactor);
                histograms[ikt][iosl+1]->GetYaxis()->SetTitleSize(0.1);
                histograms[ikt][iosl+1]->GetYaxis()->SetLabelSize(0.06);
                if(iframe == 0 && iosl == 0) histograms[ikt][iosl+1]->GetYaxis()->SetTitle("D(#rho)");
                if(iframe == 0 && iosl == 0) histograms[ikt][iosl+1]->GetYaxis()->SetTitleOffset(0.9);
                if(iframe == 0) histograms[ikt][iosl+1]->GetXaxis()->SetTitle(Form("#rho_{%s} [fm]",osl_labels[iosl]));
                if(iframe == 0) histograms[ikt][iosl+1]->GetXaxis()->SetTitleOffset(0.8*xfactor);
                histograms[ikt][iosl+1]->SetMinimum(5.e-7);
                histograms[ikt][iosl+1]->SetMaximum(3.);
                TF1* f_levyfunc = new TF1(Form("levyfunc%i%i",iframe,iosl), fitFunction, rfitmin, rfitmax, 3);
                TF1* f_levyfunc_full = new TF1(Form("levyfunc_full%i%i",iframe,iosl), fitFunction, 0.1, 5000., 3);
                f_levyfunc->SetParNames("alpha","R","norm");
                f_levyfunc->SetLineStyle(1);
                double Rosl = results[1+iosl];
                double dRosl = errors[1+iosl];
                //if(iosl == 3) Rosl = TMath::Sqrt(std::abs(results[1+iosl]));
                f_levyfunc_full->SetParameters(alpha, Rosl, N); // alpha, (Rout or Rside or Rlong), N
                f_levyfunc->SetParameters(alpha, Rosl, N);
                histograms[ikt][iosl+1]->GetXaxis()->SetLabelSize(0.06/xfactor);
                histograms[ikt][iosl+1]->GetXaxis()->SetLabelOffset(-0.025/xfactor/xfactor/xfactor);
                histograms[ikt][iosl+1]->Draw("pe");


                f_levyfunc_full->SetLineStyle(2);
                f_levyfunc_full->Draw("same");
                f_levyfunc->SetLineStyle(1);
                f_levyfunc->Draw("same");

                if(iosl == 0) leg->AddEntry(histograms[ikt][iosl+1],"D(#rho)","PE");

                Tl.SetTextSize(30);
                if(iosl == 0) Tl.DrawLatex(0.24, 0.40, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
                //if(iosl == 0) Tl.DrawLatex(0.24, 0.56, Form("conf.lev. = %1.5f", conflev));
                if(iosl == 0) Tl.DrawLatex(0.25, 0.33, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
                if(iosl == 0) Tl.DrawLatex(0.24, 0.26, Form("^{*}#lambda = %1.2f #pm %1.2f",N, dN));
                //if(iosl == 0) Tl.DrawLatex(0.24, 0.38, Form("Fit status: %s", statuses[fitstatus]));
                //if(iosl == 0) Tl.DrawLatex(0.24, 0.32, Form("Cov. matrix: %s", covstatuses[fitcovstatus+1]));
                //if(iosl == 0) Tl.DrawLatex(0.24, 0.26, Form("Edm %1.3f", minimizer->Edm()));
                Tl.DrawLatex((iosl == 0 ? 0.24 : 0.18), 0.20, Form("R_{%s} = (%1.2f #pm %1.2f) fm^{2}", osl_labels[iosl], Rosl, dRosl));
              } // end iosl loop

              canvas->cd();
              Tl.SetTextSize(30);
              const char* isPathUrqmdTitle = IsUrQMD ? "UrQMD" : "EPOS4";
              Tl.DrawLatex(0.08, 0.92, Form("%s %i events", isPathUrqmdTitle, NEVT_AVG) );
              Tl.DrawLatex(0.08, 0.86, "#pi#kern[-0.3]{{}^{#pm}}#pi#kern[-0.3]{{}^{#pm}} pair-source projections" );
              Tl.DrawLatex(0.4, 0.92, "0#minus10% ^{197}Au#plus^{197}Au" );
              Tl.DrawLatex(0.4, 0.86, Form("#sqrt{s_{NN}} = %s GeV", energy));
              Tl.DrawLatex(0.70, 0.90, Form("K_{T} [GeV/c] = %1.3f #minus %1.3f",ktbins[ikt],ktbins[ikt+1]) );
            } // end if 3D

            // Both for 1D and 3D: save the canvas as a PNG file
            if(!ikt_plotted[ikt])
            {
              const char* fitQualityTag = (strcmp(fitQuality,"GOODFIT") == 0) ? "" : Form("_BADFIT_%s",fitQuality);
              canvas->SaveAs(Form("%s/figs/fitting/%s/%s_onedsource_cent%s_%s_ifile%i_ievt%i_ikt%i_ich0_AVG%d%s.png", 
                                  path, frames[thisframe], isPathUrqmd, centleg[ICENT], energy, ifile, ievt, ikt, NEVT_AVG, 
                                  fitQualityTag));
              ikt_plotted[ikt] = true; // FIXME uncomment to only save first per ikt
              if(!is3Dfit)
              {
                canvas->Clear();
              }
            }
            if(is3Dfit)
            {
              leg->Clear();
              // Delete old pads before creating new ones
              for(int ipad = 0; ipad < NOSL; ipad++)
              {
                if(pad[ipad])
                {
                  pad[ipad]->Delete();
                  pad[ipad] = nullptr;
                }
              }
              // Create new pads
              pad[0] = new TPad("pad1","",0.,0.,0.37,1.);
              pad[1] = new TPad("pad2","",0.37,0.,0.68,1.);
              pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
              canvas->cd();
              for(int ipad = 0; ipad < NOSL; ipad++) pad[ipad]->Draw();
            }
          }

          // ---- end of plotting ----

          // DO NOT SAVE bad fits
          if(strcmp(fitQuality, "GOODFIT") != 0)
          {
            NDF=0; // prbably not needed, but just in case
            cout << "Bad fit, skipping saving. Reason: " << fitQuality << endl;
            
            continue;
          }
          // SAVING results
          // --------------
          //Ngoodfits++; // FIXME uncomment this if needed for custom histogram root-naming
          // Debug: check for under/overflow or non-finite values before filling
          if(!TMath::Finite(alpha) || !TMath::Finite(R) || !TMath::Finite(N)) {
            cerr << Form("DEBUG: Non-finite fit parameters alpha=%g R=%g N=%g; skipping fill (ifile=%d ievt=%d ikt=%d)\n", alpha, R, N, ifile, ievt, ikt);
          } else {
            double axmin = alpha_vs_R_all->GetXaxis()->GetXmin();
            double axmax = alpha_vs_R_all->GetXaxis()->GetXmax();
            double aymin = alpha_vs_R_all->GetYaxis()->GetXmin();
            double aymax = alpha_vs_R_all->GetYaxis()->GetXmax();
            if(alpha < axmin || alpha > axmax || R < aymin || R > aymax) {
              cerr << Form("DEBUG: alpha/R (%.4g, %.4g) outside alpha_vs_R_all range [%.4g,%.4g]x[%.4g,%.4g] (energy=%s ifile=%d ievt=%d ikt=%d)\n", alpha, R, axmin, axmax, aymin, aymax, energy, ifile, ievt, ikt);
            }
            alphahist[ikt]->Fill(alpha);
            Rhist[ikt]->Fill(R);
            if(is3Dfit)
            {
              R_out_hist[ikt]->Fill(R_out);
              R_side_hist[ikt]->Fill(R_side);
              R_long_hist[ikt]->Fill(R_long);
            }
            Nhist[ikt]->Fill(N);
            cerr << "Filling alpha vs R histogram ikt " << ikt << " with alpha,R: " << alpha << "," << R << endl;
            alpha_vs_R[ikt]->Fill(alpha,R);
            alpha_vs_R_all->Fill(alpha,R);
          }

          // Saving to arrays (then vectors):
          alpha_vec[ikt] = alpha;
          alpha_errdn_vec[ikt] = dalpha; // without asymmetric errors for now
          alpha_errup_vec[ikt] = dalpha;
          R_vec[ikt] = R;
          R_errdn_vec[ikt] = dR;
          R_errup_vec[ikt] = dR;
          if(is3Dfit)
          {
            R_out_vec[ikt] = R_out;
            R_out_errdn_vec[ikt] = dR_out;
            R_out_errup_vec[ikt] = dR_out;
            R_side_vec[ikt] = R_side;
            R_side_errdn_vec[ikt] = dR_side;
            R_side_errup_vec[ikt] = dR_side;
            R_long_vec[ikt] = R_long;
            R_long_errdn_vec[ikt] = dR_long;
            R_long_errup_vec[ikt] = dR_long;
          }
          N_vec[ikt] = N;
          N_errup_vec[ikt] = dN;
          N_errdn_vec[ikt] = dN;
          // Fill per-ikt and the global confidence histograms
          confidencehist[ikt]->Fill(conflev);
          confidencehist_all->Fill(conflev);
          //canvas->Clear();
        } // end of ikt loop
        
        if(averagingComplete) // MOVED OUTSIDE: only create graphs when averaging is complete
        {
          //cerr << "Arrived here as well." << endl; // debug averaging logic
          TGraphAsymmErrors* alpha_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, alpha_vec, xLow, xHigh, alpha_errdn_vec, alpha_errup_vec);
          alpha_vs_kt->SetTitle(Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#alpha",energy));
          alpha_vs_kt->SetName(Form("alpha_vs_kt_%d", Ngoodfits));
          TGraphAsymmErrors* R_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, R_vec, xLow, xHigh, R_errdn_vec, R_errup_vec);
          R_vs_kt->SetTitle(Form("#R(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R",energy));
          R_vs_kt->SetName(Form("R_vs_kt_%d", Ngoodfits));
          TGraphAsymmErrors* N_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, N_vec, xLow, xHigh, N_errdn_vec, N_errup_vec);
          N_vs_kt->SetTitle(Form("#N(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#N",energy));
          N_vs_kt->SetName(Form("N_vs_kt_%d", Ngoodfits));

          alpha_vs_KT_all.push_back(alpha_vs_kt);
          R_vs_KT_all.push_back(R_vs_kt);
          N_vs_KT_all.push_back(N_vs_kt);

          TGraphAsymmErrors* R_out_vs_kt = nullptr;
          TGraphAsymmErrors* R_side_vs_kt = nullptr;
          TGraphAsymmErrors* R_long_vs_kt = nullptr;
          if(is3Dfit)
          {
            R_out_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, R_out_vec, xLow, xHigh, R_out_errdn_vec, R_out_errup_vec);
            R_out_vs_kt->SetTitle(Form("#R_{out}(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R_{out}",energy));
            R_out_vs_kt->SetName(Form("R_out_vs_kt_%d", Ngoodfits));
            R_side_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, R_side_vec, xLow, xHigh, R_side_errdn_vec, R_side_errup_vec);
            R_side_vs_kt->SetTitle(Form("#R_{side}(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R_{side}",energy));
            R_side_vs_kt->SetName(Form("R_side_vs_kt_%d", Ngoodfits));
            R_long_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, R_long_vec, xLow, xHigh, R_long_errdn_vec, R_long_errup_vec);
            R_long_vs_kt->SetTitle(Form("#R_{long}(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R_{long}",energy));
            R_long_vs_kt->SetName(Form("R_long_vs_kt_%d", Ngoodfits));

            R_out_vs_KT_all.push_back(R_out_vs_kt);
            R_side_vs_KT_all.push_back(R_side_vs_kt);
            R_long_vs_KT_all.push_back(R_long_vs_kt);
          }
          
          Ngoodfits++;
        }
      } // end ievt
      cerr << "Finished all events in ifile " << ifile << endl;
      cout << "about to close input file." << endl;
      file->Close();
      delete file;
      cout << "input file closed, delete done." << endl;
    } // end ifile
    cerr << "Finished all input files, cleaning up." << endl;
    // Clean up histograms array to prevent ROOT from trying to clean them up twice
    for (int ikt = 0; ikt < NKT; ++ikt) {
      for (int iosl = 0; iosl < NOSL + 1; ++iosl) {
        if (histograms[ikt][iosl]) {
          delete histograms[ikt][iosl];
          histograms[ikt][iosl] = nullptr;
        }
      }
    }
    cerr << "histograms deleted." << endl;
  } // end of iframe loop

  
  // Clean up pads for 3D fit
  if(is3Dfit)
  {
    for(int ipad = 0; ipad < NOSL; ipad++)
    {
      if(pad[ipad])
      {
        pad[ipad]->Delete();
        pad[ipad] = nullptr;
      }
    }
    cerr << "pads deleted." << endl;
  }
  
  // Delete canvas directly without calling Clear/Close
  
  if(canvas && !is3Dfit) {
    cerr << "about to delete canvas." << endl;
    delete canvas;
    canvas = nullptr;
    cerr << "canvas deleted." << endl;
  }
  
  delete myLevy_reader;

  const char* rhofitmax_syst_label = (rho_fitmax_syst == 0) ? "" : (rho_fitmax_syst == 1) ? "_strictrhoFitMax" : "_looserhoFitMax";
  // Create output file
  TFile* file_output = new TFile(Form("./results/%s_onedfitresults_%s_cent%s_%s_AVG%d%s%s.root", 
                                      isPathUrqmd, frames[thisframe], centleg[ICENT], energy, NEVT_AVG,
                                      qlcms_syst_label, rhofitmax_syst_label), "RECREATE"); // Compare with Yan, add:  _Yan
  cout << "output file created." << endl;
  file_output->cd();
  cout << "output file cd() done." << endl;

  // Writing out to files
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    cout << "ikt: " << ikt << endl;
    //alphahist[ikt]->Write(); // these two will not be used as alpha_vs_R_all and the TGraphAsymmErrors vectors are written out
    //Rhist[ikt]->Write();
    if(is3Dfit)
    {
      R_out_hist[ikt]->Write();
      R_side_hist[ikt]->Write();
      R_long_hist[ikt]->Write();
    }
    alpha_vs_R[ikt]->Write();
    Nhist[ikt]->Write(); // I could put this with alpha and R together in TH3, but stick with this for now
    confidencehist[ikt]->Write();
  }
  alpha_vs_R_all->Write();
  confidencehist_all->Write();
  cout << "histograms written." << endl;

  for(size_t i=0; i<alpha_vs_KT_all.size(); i++)
  {
    //cout << "i: " << i << endl;
    // Maybe instead of writing out all, average them first? - not needed yet, the output are still way smaller files than input files
    alpha_vs_KT_all[i]->Write();
    R_vs_KT_all[i]->Write();
    if(is3Dfit)
    {
      R_out_vs_KT_all[i]->Write();
      R_side_vs_KT_all[i]->Write();
      R_long_vs_KT_all[i]->Write();
    }
    N_vs_KT_all[i]->Write();
  }
  cout << "vectors in form of TGraphAsymmErrors written." << endl;
  file_output->Close();
  return 0;
}

