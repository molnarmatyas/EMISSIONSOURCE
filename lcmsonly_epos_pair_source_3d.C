// execute: root -b -q 'lcmsonly_epos_pair_source_3d.C(arg1, arg2...)'
// make sure there is an 'analysed' directory in the parent dir!

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include <iostream>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h"                         // Used in PID selection
#include "TRandom3.h"                       // Used in PID selection
#include "TComplex.h"
#include "TVector2.h"
#include <fstream>
#include "TNtuple.h"
#include "TTree.h"

int nmax= 100000000;
using namespace std;
const double Mass2_pi = 0.019479835;
//const double Mass2_pi = 0.13957*0.13957;
const double Mass2_ka = 0.24371698032;
bool file_ok = false;
TFile *EPOS_file;
TTree *tree;
int *id, *ist, *ity, *ior, *jor, np;
float bim;
float *zus, *px, *py, *pz, *m, *x, *y, *z, *t;
Int_t numberOfEvents, all_numberOfEvents;

const int NKTBIN = 10;
const int NCH = 2;
//const int NEVT = 10;
const int NCENT = 10; // number of centrality classes

const double kt_limits[NKTBIN+1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
const char * centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};
//icent = 0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70%, 8: 70-80%, 9: 80-90%

TFile *file_output;   //here you get output file field with histograms for different particles
TH1D *bimhist;
TH1D *centhist;
TH1D *etahist;
TH1D *nphist;
TH2D *npist_vs_np;
TH2D *mult_pions_vs_allcharged;
TH2D *bim_vs_np;
TH1F *pion_pair_source_avg_lab[NKTBIN][4];
TH1F *pion_pair_source_avg_lcms[NKTBIN][4];
//TH1F *pion_pair_source_avg_pcms[NKTBIN][NCH][3];

// To limit printout after a given count
int pair_count = 0;

const char* dirs[4] = {"out","side","long","sqrtrho2"};
int ifile;

inline bool exists_test0(TString name)
{
  ifstream f(name);
  return f.good();
}

void opentree(TString filename)
{
  zus = new float [nmax]; //private use
  px = new float [nmax];//  p_x of particle
  py = new float [nmax];//  p_y of particle
  pz = new float [nmax];//  p_z of particle
  m = new float [nmax]; //  mass of particle
  x = new float [nmax]; //  x component of formation point
  y = new float [nmax]; //  y component of formation point
  z = new float [nmax]; //  z component of formation point
  t = new float [nmax]; //  formation time
  id = new int [nmax];  //  particle id  (see array "idt" in function "idtrafo" in file "ids.f": first column "epos id" second one PDG id
  ist = new int [nmax]; //  particle status (hadron last generation (0) or not (1); other numbers refer to partons, Pomerons, etc)
  ity = new int [nmax]; //  type of particle origin (20-29 from soft strings, 30-39 from hard strings, 40-59 from remnants, 60 from fluid)
  ior = new int [nmax]; //  index of father  (resonance decay products have only a father)
  jor = new int [nmax]; //  index of mother  (mothers are needed for exemple for strings: the partons between ior and jor constitute the string)
//  np = new int ;        //  number of particles
//  bim = new float ;     //  impact parameter (usually; other choices are possible)

  EPOS_file = new TFile(filename);
  if(EPOS_file->IsZombie()) cout << "Zombie file" << endl;
  else
  {
    if(EPOS_file->Get("teposevent0"))
    {
      file_ok = true;
      tree = (TTree*)EPOS_file->Get("teposevent0");
      tree->SetBranchAddress("np",&np);
      tree->SetBranchAddress("bim",&bim);
      tree->SetBranchAddress("zus",zus);
      tree->SetBranchAddress("px",px);
      tree->SetBranchAddress("py",py);
      tree->SetBranchAddress("pz",pz);
      tree->SetBranchAddress("e",m);
      tree->SetBranchAddress("x",x);
      tree->SetBranchAddress("y",y);
      tree->SetBranchAddress("z",z);
      tree->SetBranchAddress("t",t);
      tree->SetBranchAddress("id", id);
      tree->SetBranchAddress("ist",ist);
      tree->SetBranchAddress("ity",ity);
      tree->SetBranchAddress("ior",ior);
      tree->SetBranchAddress("jor",jor);
    }
    else
    {
      file_ok = false;
      EPOS_file->Close();
    }
  }
}

//=========================centrality===========================
Int_t centrality( Float_t impactparameter )
{
  //Float_t ImpactParameterLimits[] = {0, 3.32, 4.73, 6.68, 8.17, 9.42, 10.55, 11.57, 12.48, 13.37};
  Float_t ImpactParameterLimits[] = {0, 3.38, 4.8, 6.68, 8.17, 9.42, 10.55, 11.57, 12.48, 13.37}; // 0-10% same as Yan Huang's
  Int_t MiddleBinID[]      = {0,    1,    2,    3,    4,    5,     6,     7,     8,     9};
  Int_t myCentrality ;
  if      ( impactparameter > ImpactParameterLimits[9] )  { myCentrality = MiddleBinID[0] ; }
  else if ( impactparameter > ImpactParameterLimits[8] )  { myCentrality = MiddleBinID[1] ; }
  else if ( impactparameter > ImpactParameterLimits[7] )  { myCentrality = MiddleBinID[2] ; }
  else if ( impactparameter > ImpactParameterLimits[6] )  { myCentrality = MiddleBinID[3] ; }
  else if ( impactparameter > ImpactParameterLimits[5] )  { myCentrality = MiddleBinID[4] ; }
  else if ( impactparameter > ImpactParameterLimits[4] )  { myCentrality = MiddleBinID[5] ; }
  else if ( impactparameter > ImpactParameterLimits[3] )  { myCentrality = MiddleBinID[6] ; }
  else if ( impactparameter > ImpactParameterLimits[2] )  { myCentrality = MiddleBinID[7] ; }
  else if ( impactparameter > ImpactParameterLimits[1] )  { myCentrality = MiddleBinID[8] ; }
  else                                                    { myCentrality = MiddleBinID[9] ; }
  return myCentrality;
}

Int_t centrality_mult( Float_t mult )
{
  Float_t MultiplicityLimits[] = {3806.67,3534.29,3210,2226.67,1640,1261.67,985,823.333,694,-50};
  Int_t MiddleBinID[]      = {0,    1,    2,    3,    4,    5,     6,     7,     8,     9};
  Int_t myCentrality ;
  if      ( mult > MultiplicityLimits[0] )  { myCentrality = MiddleBinID[0] ; }
  else if ( mult > MultiplicityLimits[1] )  { myCentrality = MiddleBinID[1] ; }
  else if ( mult > MultiplicityLimits[2] )  { myCentrality = MiddleBinID[2] ; }
  else if ( mult > MultiplicityLimits[3] )  { myCentrality = MiddleBinID[3] ; }
  else if ( mult > MultiplicityLimits[4] )  { myCentrality = MiddleBinID[4] ; }
  else if ( mult > MultiplicityLimits[5] )  { myCentrality = MiddleBinID[5] ; }
  else if ( mult > MultiplicityLimits[6] )  { myCentrality = MiddleBinID[6] ; }
  else if ( mult > MultiplicityLimits[7] )  { myCentrality = MiddleBinID[7] ; }
  else if ( mult > MultiplicityLimits[8] )  { myCentrality = MiddleBinID[8] ; }
  else                                                    { myCentrality = MiddleBinID[9] ; }
  return myCentrality;
}



void calculations(int ICENT)
{
//=====================================================
  numberOfEvents = tree->GetEntries();
// cout << "Number of events in file: " << numberOfEvents << endl;
  all_numberOfEvents += numberOfEvents;
//  cout << "Analized number of events: " << all_numberOfEvents << endl;
  int nevt_kept = 0;
  for(int ievt = 0; ievt < numberOfEvents; ievt++)
  {
//    if(nevt_kept > 0) break;
    tree->GetEntry(ievt);
    //=====================================================
    int icent = 9-centrality(bim);
    cout << "centrality bin of event: " << icent << endl;
    bimhist->Fill(bim);
    centhist->Fill(icent);
    int np_ist = 0;
    int np_pion_ist = 0;
    //nphist->Fill(np);
    
    // !!! centrality selection
    if(ICENT == NCENT)
    {
      cout << "Doing calculations for all centrality classes." << endl;
    }
    else if(ICENT == NCENT+1)
    {
      if(icent > 1) continue; // icent>1 ~ 0-10%
      cout << "Doing calculations for centrality class(es) 0-10 percent." << endl;
    }
    else
    {
      if(ICENT != icent) continue; // normal centrality selection
    }
    nevt_kept++;
    //=====================================================

    //==== TRACK LOOP ====//
    //=====================
//    cout << "event processing..." << endl;
    int npions = 0;
    for(int itrk = 0; itrk < np; itrk++)
    {
      if ( ist[itrk] != 0 ) continue; //core-corona urqmd
      /**/
      if(fabs(id[itrk]) == 120 || fabs(id[itrk]) == 130|| fabs(id[itrk]) == 1120){
        np_ist++; // counting multiplicity of ist==0 charged hadrons
      }
      /**/
      if(int(fabs(id[itrk])) != 120) continue; // This should be good as well, with or without type casting
      np_pion_ist++; // for counting multiplicity of ist==0 charged pions
      double pTi = TMath::Sqrt(px[itrk] * px[itrk] + py[itrk] * py[itrk]);
      double p1 = sqrt(pTi*pTi + pz[itrk]*pz[itrk]);
      double eta1 = 0.5 * TMath::Log(TMath::Abs((p1 + pz[itrk]) / (p1 - pz[itrk])));
      if(fabs(eta1) > 1.) continue;
      etahist->Fill(eta1);
      if(pTi < 0.15 || pTi > 1.) continue;
      //if(p1 < 0.15 || p1 > 1.) continue;
      npions++;
      for(int jtrk = itrk+1; jtrk < np; jtrk++)
      {
        if ( ist[jtrk] != 0 ) continue; //core-corona urqmd
        if(id[itrk] != id[jtrk]) continue;
        double pTj = TMath::Sqrt(px[jtrk] * px[jtrk] + py[jtrk] * py[jtrk]);
        if(pTj < 0.15 || pTj > 1.) continue;
        double p2 = sqrt(pTj*pTj + pz[jtrk]*pz[jtrk]);
        //if(p2 < 0.15 || p2 > 1.) continue;
        double eta2 = 0.5 * TMath::Log(TMath::Abs((p2 + pz[jtrk]) / (p2 - pz[jtrk])));
        if(fabs(eta2) > 1.) continue;

        double kTr = 0.5 * TMath::Sqrt( (px[itrk] + px[jtrk]) * (px[itrk] + px[jtrk]) + (py[itrk] + py[jtrk]) * (py[itrk] + py[jtrk]));
        //double kTr = 0.5 * TMath::Sqrt( (px[itrk] + px[jtrk]) * (px[itrk] + px[jtrk]) + (py[itrk] + py[jtrk]) * (py[itrk] + py[jtrk]) + (pz[itrk] + pz[jtrk]) * (pz[itrk] + pz[jtrk]));
        if(kTr < 0.175 || kTr > 0.675) continue;
        int ikt = (kTr - kt_limits[0]) / (kt_limits[1] - kt_limits[0]); //do not forget to cross-check!! only works for evenly spaced bins
        //int ikt = std::upper_bound(std::begin(kt_limits), std::end(kt_limits), kTr) - kt_limits - 1;
        //cout << "kTr, ikt, kt_limits[ikt]: " << kTr << ", " << ikt << ", " << kt_limits[ikt] << endl;

        double pxi = px[itrk];
        double pxj = px[jtrk];
        double pyi = py[itrk];
        double pyj = py[jtrk];
        double pzi = pz[itrk];
        double pzj = pz[jtrk];
        double Ei = sqrt(p1*p1 + Mass2_pi);
        double Ej = sqrt(p2*p2 + Mass2_pi);

        double qx = (pxi - pxj);
        double qy = (pyi - pyj);
        double qT2 = qx*qx+qy*qy;
        double qz2 = (4.*(pzi*Ej-pzj*Ei)*(pzi*Ej-pzj*Ei) / 
                     ((Ei+Ej)*(Ei+Ej)-(pzi+pzj)*(pzi+pzj)));
        double Qlcms = sqrt(qT2+qz2);
        //if(Qlcms > kt_limits[ikt+1]) continue;
        //if(Qlcms > kt_limits[ikt]) continue;
        if(Qlcms > kTr) continue;

        double ti = t[itrk];
        double tj = t[jtrk];
        double xi = x[itrk];
        double xj = x[jtrk];
        double yi = y[itrk];
        double yj = y[jtrk];
        double zi = z[itrk];
        double zj = z[jtrk];
        double rx = (xi - xj);
        double ry = (yi - yj);
        double rz = (zi - zj);
        double rt = (ti - tj);
        double K0 = Ei + Ej;
        double Kx = pxi + pxj;
        double Ky = pyi + pyj;
        double Kz = pzi + pzj;
        double Kt = sqrt(Kx*Kx+Ky*Ky);
        //double K2 = Kx * Kx + Ky * Ky + Kz * Kz;
        //double M = sqrt(K0*K0 - K2);
        double cosphi = Kx/Kt;
        double sinphi = Ky/Kt;
        //double betax = Kx/K0;
        //double betay = Ky/K0;
        double betaz = Kz/K0;

        double rhoout  = fabs( rx * cosphi + ry * sinphi - rt*Kt/K0);
        double rhoside = fabs(-rx * sinphi + ry * cosphi);
        double rholong = fabs(rz - betaz*rt);

        double rhoout_lcms  = fabs( rx * cosphi + ry * sinphi - (Kt/(K0*K0-Kz*Kz))*(K0*rt-Kz*rz) );
        //double rhoside_lcms = rhoside;
        double rholong_lcms = (K0*rz-Kz*rt)/sqrt(K0*K0-Kz*Kz);
        //double pcms_coeff = ((Kx * rx + Ky * ry + Kz * rz) * (K0 - M)/ K2 - rt) / M;
        //double rhoout_pcms  = fabs( rx * cosphi + ry * sinphi + pcms_coeff * Kt);
        //double rhoside_pcms = rhoside;
        //double rholong_pcms = fabs(rz + pcms_coeff * Kz);
        double rho = TMath::Sqrt(rhoout*rhoout + rhoside*rhoside + rholong*rholong);
        //double rho_lcms = TMath::Sqrt(rhoout_lcms*rhoout_lcms + rhoside_lcms*rhoside_lcms + rholong_lcms*rholong_lcms);
        double rho_lcms = TMath::Sqrt(rhoout_lcms*rhoout_lcms + rhoside*rhoside + rholong_lcms*rholong_lcms);

        pion_pair_source_avg_lab[ikt][0]->Fill(rhoout);
        pion_pair_source_avg_lcms[ikt][0]->Fill(rhoout_lcms);
        //pion_pair_source_avg_pcms[ikt][ich][0]->Fill(rhoout_pcms);
        pion_pair_source_avg_lab[ikt][1]->Fill(rhoside);
        pion_pair_source_avg_lcms[ikt][1]->Fill(rhoside);
        //pion_pair_source_avg_pcms[ikt][ich][1]->Fill(rhoside);
        pion_pair_source_avg_lab[ikt][2]->Fill(rholong);
        pion_pair_source_avg_lcms[ikt][2]->Fill(rholong_lcms);
        //pion_pair_source_avg_pcms[ikt][ich][2]->Fill(rholong_pcms);
        pion_pair_source_avg_lab[ikt][3]->Fill(rho);
        pion_pair_source_avg_lcms[ikt][3]->Fill(rho_lcms);

        // Print out rho of the 1st 100 particles
        if(pair_count < 100)
        {
          cout << "Particle i momentum coordinates: px = " << pxi << ", py = " << pyi << ", pz = " << pzi << endl;
          cout << "Particle j momentum coordinates: px = " << pxj << ", py = " << pyj << ", pz = " << pzj << endl;
          cout << "Particle rho_lcms = " << rho_lcms << endl;
          pair_count++;
        }
      }
    }//end of track loop
    //    cout << "end of track loop; number of pions analyzed: " << npions << endl;

    // Multiplicity histos
    //nphist->Fill(static_cast<double>(np_ist));
    nphist->Fill(static_cast<double>(np_pion_ist));
    npist_vs_np->Fill(static_cast<double>(np_ist), np);
    mult_pions_vs_allcharged->Fill(static_cast<double>(np_pion_ist), static_cast<double>(np_ist));
    bim_vs_np->Fill(bim, static_cast<double>(np_ist));

    file_output->cd();
    for(int ikt = 0; ikt < NKTBIN; ikt++)
      for(int iosl = 0; iosl < 4; iosl++)
      {
        pion_pair_source_avg_lcms[ikt][iosl]->SetName(Form("pion_pair_source_avg_lcms_ifile%i_ievt%i_ikt%i_%s", ifile, ievt, ikt, dirs[iosl]));
        pion_pair_source_avg_lcms[ikt][iosl]->Write();
        pion_pair_source_avg_lcms[ikt][iosl]->Reset("icesm");
        pion_pair_source_avg_lab[ikt][iosl]->SetName(Form("pion_pair_source_avg_lab_ifile%i_ievt%i_ikt%i_%s", ifile, ievt, ikt, dirs[iosl]));
        //pion_pair_source_avg_lab[ikt][iosl]->Write();
        pion_pair_source_avg_lab[ikt][iosl]->Reset("icesm");
      }
    //    cout << "event processed." << endl;
  }//////// ==========  END OF LOOP over EVENTS
  cout << "analysis done on " << nevt_kept << " events." << endl;
}

void delete_tab()
{
  delete [] zus;
  delete [] px;
  delete [] py;
  delete [] pz;
  delete [] m;
  delete [] x;
  delete [] y;
  delete [] z;
  delete [] t;
  delete [] id;
  delete [] ist;
  delete [] ity;
  delete [] ior;
  delete [] jor;
}

//========================= main =======================
void lcmsonly_epos_pair_source_3d(const int NFILEMAX=400, const int NEVT=25, const char* energy = "200GeV", int ICENT=11)//epos_correlations_eventbyevent(int orig_arg)
{
  // Not too efficient, as it needs full analysis runs for centrality classes instead of collecting them together, but this was the easier approach for development
  if(ICENT == NCENT || ICENT < 0)
  {
    ICENT = NCENT; // centleg element: "all"
    cout << "WARNING: all centrality classes together." << endl;
  }else if(ICENT > NCENT)
  {
    ICENT = NCENT + 1; // centleg element: "0-10%"
    cout << "WARNING: 0-10 percent centrality together." << endl;
  }else{
    cout << centleg[ICENT] << " percent centrality class to be analysed." << endl;
  }
  //file_output = new TFile(Form("./analysed/EPOS_3d_source_%scent_all_%s.root",centleg[ICENT], energy), "RECREATE");   //here you get output file field with histograms for different particles
  file_output = new TFile(Form("./analysed/EPOS_3d_source_%scent_all_%s.root",centleg[ICENT], energy), "RECREATE");   //here you get output file field with histograms for different particles
  bimhist = new TH1D("epos_impact_parameter_distribution","",2000,0.,20.);
  centhist = new TH1D("epos_centrality_class_distribution","",10,0.,10.);
  etahist = new TH1D("epos_pion_eta_distribution","",2000,-10.,10.);
  nphist = new TH1D("epos_pion_multiplicity_distribution","Pion npart, ist=0",100,0.,5000.); // 15000 60000 // !!! make sure title is correct
  npist_vs_np = new TH2D("npist_vs_np","",100,0.,5000.,800,0.,40000.);
  mult_pions_vs_allcharged = new TH2D("mult_pions_vs_allcharged","",100,0.,5000.,100,0.,5000.);
  bim_vs_np = new TH2D("bim_vs_np","",100,0.,20.,100,0.,40000.);

  //const int nbins = 100; //168; FIXME only for comparison with Yan
  //double xmin = 3.e-1;//1.e-1;
  //double xmax = 1.2e3;//1.e3;
  const int nbins = 168;
  double xmin = 1.e-1;
  double xmax = 1.e3;
  double binfactor = TMath::Power(xmax/xmin,1.0/nbins);
  double xbins[nbins+1];
  xbins[0] = xmin;
  for(int ibin = 1; ibin <= nbins; ibin++)
    xbins[ibin] = xmin * TMath::Power(binfactor,ibin);
  for(int ikt = 0; ikt < NKTBIN; ikt++)
    for(int iosl = 0; iosl < 4; iosl++)
    {
      pion_pair_source_avg_lab[ikt][iosl]  = new TH1F(Form("pion_pair_source_avg_lab_ikt%i_%s", ikt, dirs[iosl]), "", nbins, xbins);
//      pion_pair_source_avg_lcms[ikt][ich][iosl] = new TH1F(Form("pion_pair_source_avg_lcms_ikt%i_ich%i_%s", ikt, ich, dirs[iosl]), "", nbins, 0., 200.);
      pion_pair_source_avg_lcms[ikt][iosl] = new TH1F(Form("pion_pair_source_avg_lcms_ikt%i_%s", ikt, dirs[iosl]), "", nbins, xbins);
//      pion_pair_source_avg_pcms[ikt][ich][iosl] = new TH1F(Form("pion_pair_source_avg_pcms_ikt%i_ich%i_%s", ikt, ich, dirs[iosl]), "", nbins, xbins);
    }
//  for(int ikt = 0; ikt < NKTBIN; ikt++)
//    for(int ich = 0; ich < NCH; ich++)
//    {
//      pion_pair_source_avg_lab[ikt][ich] = new TH3F(Form("pion_pair_source_avg_lab_ikt%i_ich%i", ikt, ich), "", nbins, xbins, nbins, xbins, nbins, xbins);
//      pion_pair_source_avg_lcms[ikt][ich] = new TH3F(Form("pion_pair_source_avg_lcms_ikt%i_ich%i", ikt, ich), "", nbins, xbins, nbins, xbins, nbins, xbins);
//      pion_pair_source_avg_pcms[ikt][ich] = new TH3F(Form("pion_pair_source_avg_pcms_ikt%i_ich%i", ikt, ich), "", nbins, xbins, nbins, xbins, nbins, xbins);
//    }

////////////////////////////////////////////////////////////
//reading the files
  //==========================
  for(int i = 1; i < NFILEMAX+1; i++)
  {
    ifile = i;
    /*
    // Deprecated
    TString number;
    number.Form("%d",i);
    //TString path = "/phenix/u/kincsesd/plhf1/EPOS/files/z-gg2my59a-";
    TString path = "/home/starelte/OneDrive/KutatÁsás/epos4_simdata/200GeV_5nfull5nfreeze/output/z-prueba_EOS2_";
    //TString path = "/home/starelte/OneDrive/KutatÁsás/epos4_simdata/7.7GeV_100nfull10nfreeze/z-prueba_EOS2_";
    //TString path = "/mnt/c/Users/MolnarMatyas/OneDrive-elte.hu/KutatÁsás/epos4_simdata/7.7GeV_100nfull10nfreeze/z-prueba_EOS2_";
    path.Append(number);
    path.Append(".root");
    */
    //TString path = "/mnt/c/Users/MolnarMatyas/OneDrive-elte.hu/KutatÁsás/epos4_simdata/x3ff/"; // !!! change this to actual path
    //TString path = "/home/starelte/OneDrive/KutatÁsás/epos4_simdata/x3ff/";
    TString path = "../../epos4_simdata/x3ff/";
    //TString path = "/home/starelte/OneDrive/KutatÁsás/EPOS4/Stuff/z-200GeV_0-10cent_1evt.root"; // FIXME only for comparison with Yan, , first event enough
    path.Append(Form("%s/output/z-prueba_EOS2_%d.root",energy,i)); // FIXME
    if(exists_test0(path))
    {
      opentree(path);
      if(file_ok)
      {
        if(i % 10 == 1) cout << " ***** File " << path << "  opened ***" << endl;
        file_output->cd();
        calculations(ICENT);
        cout << "calculations done." << endl;
        delete_tab();
        cout << "delete done." << endl;
      }
      else
      file_ok = false;
      EPOS_file -> Close();
      cout << "input file closed." << endl;
    }
  }
  cout << "file loop ended." << endl;
  {
    cout << "Calculating centrality bin multiplicity thresholds..." << endl;
    // Define centrality binning
    const int nQuantiles = 10;
    double quantileValues[nQuantiles];  
    double probabilities[nQuantiles];  

    for (int i = 0; i < nQuantiles; i++) {
        probabilities[i] = 1.0 - (i + 1) * 0.1;  // 90%, 80%, ..., 10%, 0%
    }

    // Compute centrality thresholds
    nphist->GetQuantiles(nQuantiles, quantileValues, probabilities);
    // Print threshold values
    cout << "Centrality Thresholds:" << endl;
    for (int i = 0; i < nQuantiles; i++) {
        cout << "Centrality " << i * 10 << "-" << (i + 1) * 10 
              << "% threshold: " << quantileValues[i] << endl;
    }
    cout << "Threshold array:" << endl;
    cout << "{";
    for (int i=0; i<nQuantiles-1; i++)
    {
      cout << quantileValues[i] << ",";
    }
    cout << quantileValues[nQuantiles-1] << "}" << endl;

    /*
    pion-kaon-proton, ist==0
    {2338.89,1391.67,996.429,816.667,616.667,472.222,335.714,247.17,202.83,-50}
    pion:
    {1855,1100,708.333,550,387.5,268.085,149.242,82.1429,34.9265,-50}
    np:
    {3806.67,3534.29,3210,2226.67,1640,1261.67,985,823.333,694,-50}
    */
  }

  file_output->cd();
  bimhist->Write();
  centhist->Write();
  etahist->Write();
  nphist->Write();
  npist_vs_np->Write();
  bim_vs_np->Write();
  mult_pions_vs_allcharged->Write();
//  for(int ikt = 0; ikt < NKTBIN; ikt++)
//    for(int ich = 0; ich < NCH; ich++)
//      for(int iosl = 0; iosl < 3; iosl++)
//      {
////        pion_pair_source_avg_lab[ikt][ich][iosl]->Write();
//        pion_pair_source_avg_lcms[ikt][ich][iosl]->Write();
////        pion_pair_source_avg_pcms[ikt][ich][iosl]->Write();
//      }
  file_output -> Close();
  cout << "output file closed." << endl;
  return;
}//////////// End of Analysis
