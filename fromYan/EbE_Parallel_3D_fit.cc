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
#include "json.hpp"

using namespace std;

using json = nlohmann::json;
const int NKT = 8;
std::vector<double> ktbins;

// Global variables
const int NCUT = 1;//9986;
const int NEVT = 10000;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 4;
const int NPARS = 5;
const int NOSL = 3;
const double Mass2_pi = 0.019479835;

const int colors[10] = {632,416,600,432,616,400,8,9,42,46};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.0025,0.0025,0.04,0.,0.,0.04,0.,0.,0.04};
const float leftmargins[9]  = {0.2,0.0005,0.0005,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.03,0.03,0.03,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.16,0.16,0.16,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
//const float topmargins[4] = {0.05, 0.05, 0.05, 0.05};    
//const float bottommargins[4] = {0.15, 0.15, 0.15, 0.15}; 
//const float leftmargins[4] = {0.20, 0.10, 0.20, 0.10};   
//const float rightmargins[4] = {0.05, 0.15, 0.05, 0.15};
const char* frames[3] = {"LCMS","pcms","lab"};
const char* osl_labels[NOSL] = {"out","side","long"};//,"outlong"};


std::vector<TH1D*> histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

double rfitmax;
std::vector<std::vector< double> > rfitmax_vec = {
                                          {60.0, 40.0, 80.0},
                                          {60.0, 40.0, 80.0},
                                          {60.0, 40.0, 80.0},
                                          {60.0, 40.0, 80.0},
                                          {60.0, 40.0, 80.0},
                                          {50.0, 40.0, 80.0},
                                         };

std::vector<std::vector< double> > rfitmin_vec = {
                                          {1.0, 1.0, 1.0},
                                          {1.0, 1.0, 1.0},
                                          {1.0, 1.0, 1.0},
                                          {1.0, 1.0, 1.0},
                                          {1.0, 1.0, 1.0},
                                          {1.0, 1.0, 1.0}
                                         };
double rfitmin;

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

// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  NDF = 0;
  double logL = 0.0;
  double alpha  = params[0];
  double R_out  = params[1];
  double R_side = params[2];
  double R_long = params[3];
  //double R_outlong_2  = params[4];
  //double R_outlong_proj = std::sqrt(0.5 * R_out * R_out + 0.5 * R_long * R_long + R_outlong_2);
  double N  = params[4];
  for (int i = 0; i < NOSL; ++i)
  {
    double R = (i == 0) ? R_out : ((i == 1) ? R_side : R_long );
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
  NDF -= NPARS;
//  std::cout << "N: " << N << std::endl;
  return 2.0 * logL;
}

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
int main(int argc, char** argv)
{
  int ikt;
  if(argc < 3) {
    std::cerr << "Invalid input levy_fit <system> <kt num> <evt avg> <beammom>" << std::endl; 
    return 0;
  }
  std::string system = argv[1];
  ikt = atoi(argv[2]); 
  int evt_avg = atoi(argv[3]); 
	int beamMomentum = std::atoi(argv[4]);
  int firstfit = 1;
  int firstbadfit = 1;
  float snn_from_beam = TMath::Sqrt(2*0.9315*(beamMomentum+2*0.9315));
	std::ifstream infile;
  int ibeam = 0;
	switch (beamMomentum)
	{
		case 13:
			infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_13.json");
      ibeam = 5;
			break;
		case 19:
			infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_19.json");
      ibeam = 4;
			break;
		case 30:
      if(!system.compare("BeBe")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/BeBe_150.json");
      if(!system.compare("ArSc")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_30.json");
      if(!system.compare("XeLa")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/XeLa_150.json");
      ibeam = 3;
			break;
		case 40:
			infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_40.json");
      ibeam = 2;
			break;
		case 75:
			infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_75.json");
      ibeam = 1;
			break;
		case 150:
      if(!system.compare("BeBe")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/BeBe_150.json");
      if(!system.compare("ArSc")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/ArSc_150.json");
      if(!system.compare("XeLa")) infile.open("/afs/cern.ch/user/b/bporfy/private/ArSc_hbtmeasure/EventCutsArSc/HBTProject/hbtclass/inp/XeLa_150.json");
      ibeam = 0;
			break;
		default:
			break;
	}
	json jsondata;
	infile >> jsondata;
	int beamMom = jsondata["energy"].template get<int>();
	//if(beamMom != beamMomentum){
	//	std::cerr << "Beam momentum not matching, check JSON file" << std::endl;
	//	exit(1);
	//}
  if(!system.compare("BeBe")) ktbins = jsondata["data"]["final"]["kT_sides"].template get<std::vector<double>>();
  if(!system.compare("ArSc")) ktbins = jsondata["data"]["final"]["kT_sides"].template get<std::vector<double>>();
  if(!system.compare("XeLa")) ktbins = jsondata["data"]["test"]["kT_sides"].template get<std::vector<double>>();
	
  std::ostringstream path;
  path << "/eos/user/b/bporfy/www/plots/urqmd/ArSc" << beamMomentum << "/3D/plots/";

  cout << "about to create levy reader." << endl;
	myLevy_reader = new Levy_reader("/eos/user/b/bporfy/www/files/UrQMD/levy_values_20250322_ultimate_v2_last_final_cdwcb.dat");
  cout << "levy reader created." << endl;

  TH1D* alphahist[NCUT][NKT];
  TH1D* Rohist[NCUT][NKT];
  TH1D* Rshist[NCUT][NKT];
  TH1D* Rlhist[NCUT][NKT];
  TH1D* Ravghist[NCUT][NKT];
  TH1D* Rolhist[NCUT][NKT];
  TH1D* Lambdahist[NCUT][NKT];
  TH2D* alphaRavghist[NCUT][NKT];
  TH1D* conflevhist[NCUT][NKT];

  // Open the file and get the histograms
  TFile *file;
  if(!system.compare("XeLa"))
  {
    file = TFile::Open(Form("/eos/home-b/bporfy/Private/UrQMD_data/%s%i_qLCMScut_pTcut_mTlim_t10000_etadecay_test_mT_default.root", system.c_str(), beamMomentum));
  }
  else 
  {
    file = TFile::Open(Form("/eos/home-b/bporfy/Private/UrQMD_data/%s%i_PSDCent_pairCuts_qLCMScut_pTcut_mTlim_t10000_etadecay_final_mT_default.root", system.c_str(), beamMomentum));
  }
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }
	TString inputFile = file->GetName();
	std::cout << "Using " << inputFile << " input..." << std::endl;
	inputFile.Remove(0,63);
	inputFile.Remove(inputFile.Sizeof()-6, inputFile.Sizeof());

  std::cout << "output prefix: " << inputFile << std::endl;

//  TH3F* sourcehists[NKT][NFRAME];

  TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(35);
  Tl.SetNDC(kTRUE);
  TLegend *leg = new TLegend(0.5,0.3,0.7,0.5);
  for(int icut = 0; icut < NCUT; icut++)
  {
    //TCanvas* canvas = new TCanvas("c1", "", 900, 800);
    TCanvas* canvas = new TCanvas("c1", "", 1200, 600);
    TPad *pad[NOSL];
  	pad[0] = new TPad("pad1","",0.,0.,0.37,1.);
  	pad[1] = new TPad("pad2","",0.37,0.,0.68,1.);
  	pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
    canvas->cd();
    for(int ipad = 0; ipad < NOSL; ipad++)
      pad[ipad]->Draw();
    rfitmax = rfitmax_vec[ibeam][icut];
    rfitmin = rfitmin_vec[ibeam][icut];
    alphahist[icut][ikt] = new TH1D(Form("alphahist_icut%i_ikt%i",icut, ikt),"",200,0.,2.);
    alphahist[icut][ikt]->SetDirectory(nullptr);
    alphaRavghist[icut][ikt] = new TH2D(Form("alphaRavghist_icut%i_ikt%i",icut, ikt),"",200,0.,2.,200,2.,18.);
    alphaRavghist[icut][ikt]->SetDirectory(nullptr);
    Rohist[icut][ikt] = new TH1D(Form("Rohist_icut%i_ikt%i",icut, ikt),"",200,0.,20.);
    Rohist[icut][ikt]->SetDirectory(nullptr);
    Rshist[icut][ikt] = new TH1D(Form("Rshist_icut%i_ikt%i",icut, ikt),"",200,0.,20.);
    Rshist[icut][ikt]->SetDirectory(nullptr);
    Rlhist[icut][ikt] = new TH1D(Form("Rlhist_icut%i_ikt%i",icut, ikt),"",200,0.,20.);
    Rlhist[icut][ikt]->SetDirectory(nullptr);
    Ravghist[icut][ikt] = new TH1D(Form("Ravghist_icut%i_ikt%i",icut, ikt),"",200,0,15.);
    Ravghist[icut][ikt]->SetDirectory(nullptr);
    Rolhist[icut][ikt] = new TH1D(Form("Rolhist_icut%i_ikt%i",icut, ikt),"",600,-30.,30.);
    Rolhist[icut][ikt]->SetDirectory(nullptr);
    Lambdahist[icut][ikt] = new TH1D(Form("Lambdahist_icut%i_ikt%i",icut, ikt),"",200,0.0,2.0);
    Lambdahist[icut][ikt]->SetDirectory(nullptr);
    conflevhist[icut][ikt] = new TH1D(Form("clhist_icut%i_ikt%i",icut, ikt),"",1000,0.,1.);
    conflevhist[icut][ikt]->SetDirectory(nullptr);
    for(int ievt = 0; ievt < NEVT; ievt += evt_avg)
    {
      bool badfit = false;
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
          if (histograms[ikt][iframe].size() == NOSL) 
          {
//            std::cout << "histograms size is 4, deleting." << std::endl;
            for(int i = 0; i < NOSL; i++) delete histograms[ikt][iframe][i];
            histograms[ikt][iframe].clear();
          }
          for(int iosl = 0; iosl < NOSL; iosl++)
          {
            // Form the histogram name
            TH1D* tempHist_evt_avg = nullptr;
            TString histName;
            TString histName_pos;
            for(int ievTemp = ievt; ievTemp < ievt+evt_avg; ievTemp++)
            {
              histName = Form("D_%s_%s_ev%i_ch0_KT%i", osl_labels[iosl], frames[iframe], ievTemp, ikt);
              histName_pos = Form("D_%s_%s_ev%i_ch1_KT%i", osl_labels[iosl], frames[iframe], ievTemp, ikt);

              TH1D* tempHist_neg = (TH1D*)file->Get(histName_pos)->Clone();
              TH1D* tempHist_ch0 = (TH1D*)file->Get(histName)->Clone();
              tempHist_neg->Add(tempHist_ch0);
              tempHist_ch0->Delete();
              if(ievTemp == ievt) 
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
    
          if (histograms[ikt][iframe].size() != NOSL)
          {
            std::cerr << "Error: Expected 4 histograms" << std::endl;
            continue;
          }
          //too few hits
          //if(histograms[ikt][iframe][0]->GetEntries() < 2 * 200) continue;
    
          // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
		      //ROOT::Minuit2::Minuit2Minimizer minimizer ( ROOT::Minuit2::kCombined );
        
          // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
        
          // Create the function to be minimized
          ROOT::Math::Functor functor(&logLikelihood, NPARS);
//          ROOT::Math::Functor functor(&chiSquare, 5);
          minimizer->SetFunction(functor);
        
          // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.3, 0.01, 0.5, 2.0);
//          minimizer->SetFixedVariable(0, "alpha", 1.3);
          minimizer->SetLimitedVariable(1, "R_out", 4.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(2, "R_side", 5.0, 0.01, 0., 20.);
          minimizer->SetLimitedVariable(3, "R_long", 9.0, 0.01, 0., 20.);
//          minimizer->SetFixedVariable(3, "R_long", 20.0);
//          minimizer->SetLimitedVariable(4, "R_outlong_2", 5.0, 0.01, -900, 900.);
//          minimizer->SetFixedVariable(4, "N_out", 1.);
//          minimizer->SetFixedVariable(5, "N_side", 1.);
//          minimizer->SetFixedVariable(6, "N_long", 1.);
//          minimizer->FixVariable(...);
          minimizer->SetVariable(4, "N", 1., 0.01);
          //minimizer->FixVariable(5);
//          std::cout << "Numer of parameters: " << minimizer->NDim() << std::endl;
        
          // Minimize the function
          minimizer->Minimize();
        
          //cout << "minimizing done. " << endl;
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          double conflev = TMath::Prob(chi2val, NDF);
          //std::cout << "alpha: "  << results[0] << " +- " << errors[0] << std::endl;
          //std::cout << "R_out: "  << results[1] << " +- " << errors[1] <<std::endl;
          //std::cout << "R_side: " << results[2] << " +- " << errors[2] <<std::endl;
          //std::cout << "R_long: " << results[3] << " +- " << errors[3] <<std::endl;
          ////std::cout << "R_outlong_2: " << results[4] << " +- " << errors[4] <<std::endl;
          //std::cout << "N "       << results[4] << " +- " << errors[4] <<std::endl;
          //std::cout << "chi2: "   << chi2val << std::endl;
          conflevhist[icut][ikt]->Fill(conflev);
		      int fitstatus = minimizer->Status();
		      int fitcovstatus = minimizer->CovMatrixStatus();
		      if(fitstatus<0 || fitstatus>5) fitstatus=5;
		      if(fitcovstatus<-1 || fitcovstatus>3) fitcovstatus=-1;
          if((fitstatus == 0 || fitstatus == 3) && (fitcovstatus == 3 || fitcovstatus == 2 || fitcovstatus == 1) )
          {
            alphahist[icut][ikt]->Fill(results[0]);
            Rohist[icut][ikt]->Fill(results[1]);
            Rshist[icut][ikt]->Fill(results[2]);
            Rlhist[icut][ikt]->Fill(results[3]);
            double Ravg = sqrt( (results[1]*results[1] + results[2]*results[2] + results[3]*results[3]) / 3. );
            Ravghist[icut][ikt]->Fill(Ravg);
            //Rolhist[icut][ikt]->Fill(abs(results[4]));
            alphaRavghist[icut][ikt]->Fill(results[0],Ravg);
            Lambdahist[icut][ikt]->Fill(results[4]);
            badfit = false;
          }
          else badfit = true;

//          cout << "histograms filled." << endl;             
          if(donotdraw) {/*delete minimizer->*/ continue;}
          double alpha = results[0];
          double dalpha = errors[0];
          double N = results[4];
          double dN = errors[4];
          for(int iosl = 0; iosl < NOSL; iosl++)
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
						// 4 pad setting
            //double padht = pad[ipad]->GetAbsHNDC();
            //double xfactor = 2./(canvas->GetAbsWNDC()/padwdth);
            //double yfactor = 2./(canvas->GetAbsHNDC()/padht);
						// 3 pad setting
            double xfactor =  3./(canvas->GetAbsWNDC()/padwdth);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetRangeUser(0.2,200.);
            //histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleSize(0.07/xfactor);
            //histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleSize(0.07/yfactor);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleSize(0.08/xfactor);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleSize(0.1);
            histograms[ikt][iframe][iosl]->GetYaxis()->SetLabelSize(0.06);
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitle("D(#rho)");
            if(iframe == 0 && iosl == 0) histograms[ikt][iframe][iosl]->GetYaxis()->SetTitleOffset(0.9);
            if(iframe == 0) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitle(Form("#rho_{%s} [fm]",osl_labels[iosl]));
            if(iframe == 0) histograms[ikt][iframe][iosl]->GetXaxis()->SetTitleOffset(0.8*xfactor);
            histograms[ikt][iframe][iosl]->SetMinimum(5.e-7);
            histograms[ikt][iframe][iosl]->SetMaximum(3.);
            TF1* f_levyfunc = new TF1(Form("levyfunc%i%i",iframe,iosl), fitFunction, rfitmin, rfitmax, 3);
            TF1* f_levyfunc_full = new TF1(Form("levyfunc_full%i%i",iframe,iosl), fitFunction, 0.1, 5000., 3);
            f_levyfunc->SetParNames("alpha","R","norm");
            f_levyfunc->SetLineStyle(1);
            double Rosl = results[1+iosl];
            double dRosl = errors[1+iosl];
            //if(iosl == 3) Rosl = TMath::Sqrt(std::abs(results[1+iosl]));
            f_levyfunc_full->SetParameters(results[0], Rosl, results[4]); // alpha, (Rout or Rside or Rlong), N
            f_levyfunc->SetParameters(results[0],Rosl,N);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelSize(0.06/xfactor);
            histograms[ikt][iframe][iosl]->GetXaxis()->SetLabelOffset(-0.025/xfactor/xfactor/xfactor);
            histograms[ikt][iframe][iosl]->Draw("pe");


            f_levyfunc_full->SetLineStyle(2);
            f_levyfunc_full->Draw("same");
            f_levyfunc->SetLineStyle(1);
            f_levyfunc->Draw("same");

            if(iosl == 0) leg->AddEntry(histograms[ikt][iframe][iosl],"D(#rho)","PE");

            Tl.SetTextSize(30);
//            Tl.DrawLatex(0.58-0.05*iosl, 0.78+iframe*0.1, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.73+iframe*0.1, Form("conf.lev. = %1.5f", conflev));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.68+iframe*0.1, Form("R = (%1.2f #pm %1.2f) fm",Rosl, dRosl));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.63+iframe*0.1, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
//            Tl.DrawLatex(0.58-0.05*iosl, 0.58+iframe*0.1, Form("N = %1.2f #pm %1.2f",N, dN));
            if(iosl == 0) Tl.DrawLatex(0.24, 0.40, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
            //if(iosl == 0) Tl.DrawLatex(0.24, 0.56, Form("conf.lev. = %1.5f", conflev));
            if(iosl == 0) Tl.DrawLatex(0.25, 0.33, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
            if(iosl == 0) Tl.DrawLatex(0.24, 0.26, Form("^{*}#lambda = %1.2f #pm %1.2f",N, dN));
            //if(iosl == 0) Tl.DrawLatex(0.24, 0.38, Form("Fit status: %s", statuses[fitstatus]));
            //if(iosl == 0) Tl.DrawLatex(0.24, 0.32, Form("Cov. matrix: %s", covstatuses[fitcovstatus+1]));
            //if(iosl == 0) Tl.DrawLatex(0.24, 0.26, Form("Edm %1.3f", minimizer->Edm()));
            Tl.DrawLatex((iosl == 0 ? 0.24 : 0.18), 0.20, Form("R_{%s} = (%1.2f #pm %1.2f) fm^{2}", osl_labels[iosl], Rosl, dRosl));
//
          }
        // Clean up
        /*delete minimizer->*/
        }
        if(!donotdraw)
        {
          canvas->cd();
          Tl.SetTextSize(30);
          Tl.DrawLatex(0.08, 0.92, Form("UrQMD %i event",evt_avg) );
          Tl.DrawLatex(0.08, 0.86, "#pi#kern[-0.3]{{}^{#pm}}#pi#kern[-0.3]{{}^{#pm}} pair-source projections" );
          Tl.DrawLatex(0.4, 0.92, "0#minus10% ^{40}Ar#plus^{45}Sc" );
          Tl.DrawLatex(0.4, 0.86, Form("#sqrt{s_{NN}} #approx %1.2f GeV", snn_from_beam));
//          Tl.DrawLatex(0.71, 0.91, "p_{T} [GeV/c] = 0.15#minus1.0");
//          Tl.DrawLatex(0.71, 0.91, "#pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
          //Tl.DrawLatex(0.55, 0.94, "#pi#kern[-0.3]{{}^{#pm}}#pi#kern[-0.3]{{}^{#pm}}");//#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
          Tl.DrawLatex(0.70, 0.90, Form("K_{T} [GeV/c] = %1.2f #minus %1.2f",ktbins[ikt],ktbins[ikt+1]) );
          //Tl.DrawLatex(0.70, 0.87, Form("Fit range %1.0f - %1.0f",rfitmin, rfitmax) );
//          Tl.DrawLatex(0.71, 0.65, "|#eta| < 1.0"); 
//          Tl.DrawLatex(0.71, 0.59, "Q_{LCMS} < 0.5k_{T}"); 
//          Tl.DrawLatex(0.71, 0.59, "Q_{LCMS} < k_{T}"); 
//          Tl.DrawLatex(0.71, 0.85, Form("#LTm_{T}#GT = %1.3f",sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
          if(!badfit) {
            if(firstfit == 1) 
            {
//              canvas->Print(Form("%s%s/fits/%s_projections_icut%i_evtavg%i_ikt%i_ich0.pdf(", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
              canvas->Print(Form("%s%s/fits/%s_projections_icut%i_ievt%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, ievt, evt_avg, ikt));
            }
            else 
            {
//              canvas->Print(Form("%s%s/fits/%s_projections_icut%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
              canvas->Print(Form("%s%s/fits/%s_projections_icut%i_ievt%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, ievt, evt_avg, ikt));
            }
            firstfit++;
          }
          else {
            if(firstbadfit == 1) 
            {
//              canvas->Print(Form("%s%s/fits/%s_bad_fit_projections_icut%i_evtavg%i_ikt%i_ich0.pdf(", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
              canvas->Print(Form("%s%s/fits/%s_bad_fit_projections_icut%i_ievt%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, ievt, evt_avg, ikt));
            }
            else 
            {
//              canvas->Print(Form("%s%s/fits/%s_bad_fit_projections_icut%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
              canvas->Print(Form("%s%s/fits/%s_bad_fit_projections_icut%i_ievt%i_evtavg%i_ikt%i_ich0.pdf", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, ievt, evt_avg, ikt));
            }
            firstbadfit++;
          }
//          cout << "canvas saved." << endl;
//          for(int ipad = 0; ipad < 3; ipad++) pad[ipad]->Clear();
//          cout << "pads cleared." << endl;
          canvas->Clear();
          leg->Clear();
          pad[0] = new TPad("pad1","",0.,0.,0.37,1.);
          pad[1] = new TPad("pad2","",0.37,0.,0.68,1.);
          pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
					//4 pad setting
          //pad[0] = new TPad("pad1","",0.,0.5,0.5,1.);
          //pad[1] = new TPad("pad2","",0.5,0.5,1.0,1.);
          //pad[2] = new TPad("pad3","",0.0,0.,0.5,0.5);
          //pad[3] = new TPad("pad4","",0.5,0.,1.0,0.5);
          for(int ipad = 0; ipad < NOSL; ipad++) pad[ipad]->Draw();
//          cout << "canvas cleared." << endl;
        }
      /*
      }
      */
    }
//    canvas->Print(Form("%s%s/fits/%s_projections_icut%i_evtavg%i_ikt%i_ich0.pdf)", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
//    canvas->Print(Form("%s%s/fits/%s_bad_fit_projections_icut%i_evtavg%i_ikt%i_ich0.pdf)", path.str().c_str(), frames[thisframe],inputFile.Data(), icut, evt_avg, ikt), "pdf");
    canvas->Clear();
    if(canvas) delete canvas;
  }
  cout << "File loop ended." << endl;
  file->Close();
  if(file) delete file;
  if(writetofile)
  {
    TFile* file_output = new TFile(Form("out/fitresults_%s_qLCMS_cut_010cent_all.root",frames[thisframe]), "RECREATE");  
    file_output->cd();
    for(int icut = 0; icut < NCUT; icut++)
    {
      cout << ikt << endl;
      std::cout << alphahist[icut][ikt]->GetEntries() << std::endl;
      alphahist[icut][ikt]->Write();
      alphaRavghist[icut][ikt]->Write();
      Rohist[icut][ikt]->Write();
      Rshist[icut][ikt]->Write();
      Rlhist[icut][ikt]->Write();
      Ravghist[icut][ikt]->Write();
      Rolhist[icut][ikt]->Write();
      Lambdahist[icut][ikt]->Write();
      conflevhist[icut][ikt]->Write();
    }
    file_output->Close();
  }
  delete myLevy_reader;
  return 0;
}

