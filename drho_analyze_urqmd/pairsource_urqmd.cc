// After the urqmd .f19 output files have been converted with convert_f19.cc (thus, creating the ...._tree.root files):
// Store compact LCMS pair separations; histograms are built by the fit/FFT consumers.
// TODO: add bim limits for centrality selection when minBias data will be available

#include <TH1.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TParameter.h>
#include <TNamed.h>

#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <fstream>

#include "../header_for_all_emissionsource.h"
// const int NKT = 10;
// const double ktbins[NKT+1]= {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
int barWidth = 70;

//const double Mass2_pi = 0.019479835;

//const double ktbins[NKT+1]= {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
//const double qmax[NKT] =       {0.05, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26}; // deprecated I guess
//const int NCENT = 10; // number of centrality classes
//const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};

int kTBinNum(double KT)
{
  // Keep both outer edges; discard pairs outside the configured kT range.
  if (!std::isfinite(KT) || KT < ktbins[0] || KT > ktbins[NKT])
    return -1;

  int ibin = 0;
  while (KT > ktbins[ibin + 1])
    ibin++;
  return ibin;
}

void progressbar(Int_t event, Int_t nEvents)
{
  float progress = event*1.0 / nEvents * 1.0;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << "% ( " << event << " / " << nEvents << " )" << "\r";
  std::cout.flush();
  
}

// Reimplemented in header_for_all_emissionsource.h:
// const char* qLCMS_cut[3] = {"default", "strict", "loose"};
// const double qLCMS_cut_values[3] = {0.15, 0.05, 0.25}; // in GeV/c
const char* qLCMS_outputtext[3] = {"", "_strictqLCMS", "_looseqLCMS"};


int main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "usage: pairsource_urqmd <energy_string> <qLCMS_cut> [charge_mode]" << std::endl;
    std::cerr << "  charge_mode: 0=both (default), 1=positive only, 2=negative only" << std::endl;
    return 0;
  }

  int qLCMS_cuttype = (int)atoi(argv[2]);
  if(qLCMS_cuttype < 0 || qLCMS_cuttype > 2)
  {
    std::cerr << "qLCMS_cuttype must be 0 (default), 1 (strict) or 2 (loose)" << std::endl;
    qLCMS_cuttype = 0;
    return 0;
  }
  std::cout << "Using qLCMS cut type: " << qLCMS_cut[qLCMS_cuttype] << std::endl;

  int charge_mode = 0;
  if(argc >= 4)
  {
    charge_mode = (int)atoi(argv[3]);
  }
  if(charge_mode < 0 || charge_mode > 2)
  {
    std::cerr << "charge_mode must be 0 (both), 1 (positive only) or 2 (negative only)" << std::endl;
    charge_mode = 0;
    return 0;
  }
  const char* charge_mode_label[3] = {"both", "positive", "negative"};
  const char* charge_mode_suffix[3] = {"", "_pos", "_neg"};
  std::cout << "Using charge selection: " << charge_mode_label[charge_mode] << std::endl;

  TFile *treefile = new TFile(Form("AuAu_%s_tree.root",argv[1]),"read");
  if(treefile->IsZombie()) {
    std::cerr << "Cannot open input event file" << std::endl;
    return 1;
  }
  TTree *event = dynamic_cast<TTree*>(treefile->Get("urqmd_tree"));
  if(!event) {
    std::cerr << "Input file has no urqmd_tree" << std::endl;
    return 1;
  }
  //TString inputFile = treefile->GetName(); 
  TFile* outFile = new TFile(Form("UrQMD_3d_source_%scent_all_%s%s%s%s.root",centleg[11],argv[1],
                                  qLCMS_outputtext[qLCMS_cuttype], charge_mode_suffix[charge_mode],
                                  useFFTmethod ? "_FFT" : ""),
                                  useFFTmethod ? "CREATE" : "RECREATE"); // FIXME make 11->ICENT for other centralities
  if(outFile->IsZombie()) {
    std::cerr << "Cannot create output file (FFT mode refuses to overwrite existing files)" << std::endl;
    return 1;
  }
  
  // Output TTree for pair source data (replaces per-event histograms for much smaller files)
  TTree *pairTree = new TTree("pairTree", "Pion pair source data in LCMS");
  Int_t out_ievent;
  Short_t out_ich, out_iKT;
  Float_t out_rho, out_rhoout, out_rhoside, out_rholong;
  UChar_t out_osl_signs;
  pairTree->Branch("ievent", &out_ievent, "ievent/I");
  pairTree->Branch("ich", &out_ich, "ich/S");
  pairTree->Branch("iKT", &out_iKT, "iKT/S");
  pairTree->Branch("rho", &out_rho, "rho/F");
  pairTree->Branch("rhoout", &out_rhoout, "rhoout/F");
  pairTree->Branch("rhoside", &out_rhoside, "rhoside/F");
  pairTree->Branch("rholong", &out_rholong, "rholong/F");
  if(useFFTmethod) {
    // Three sign bits recover the signed vector without duplicating its floats.
    pairTree->Branch("rho_osl_signs", &out_osl_signs, "rho_osl_signs/b");
  }

  TH1D *KTdist = new TH1D("ktdist","K_{T} distribution",100,0,10);
  TH1D *rapdist = new TH1D("rapdist","rapidity distribution",100,-10,10);
  TH1D *post_rapdist = new TH1D("post_rapdist","post-cut rapidity distribution",100,-10,10);
  TH1D *etadist = new TH1D("etadist","pseudorapidity distribution",100,-10,10);
  TH1D *post_etadist = new TH1D("post_etadist","post-cut pseudorapidity distribution",100,-10,10);

  std::cout << "File " << " loaded" << std::endl;
  std::vector<double> *PX = 0;
  std::vector<double> *PY = 0;
  std::vector<double> *PZ = 0;
  std::vector<double> *energy = 0;
  std::vector<double> *rx = 0; 
  std::vector<double> *ry = 0; 
  std::vector<double> *rz = 0; 
  std::vector<double> *time = 0;
  std::vector<int> *ID = 0;
  event->SetBranchAddress("pX",&PX);
  event->SetBranchAddress("pY",&PY);
  event->SetBranchAddress("pZ",&PZ);
  event->SetBranchAddress("E",&energy);
  event->SetBranchAddress("rX",&rx);
  event->SetBranchAddress("rY",&ry);
  event->SetBranchAddress("rZ",&rz);
  event->SetBranchAddress("t",&time);
  event->SetBranchAddress("PID",&ID);

  Int_t nEvent = event->GetEntries();
  for (Int_t ievent = 0; ievent < nEvent; ievent++) {
    event->GetEntry(ievent);
    progressbar(ievent+1.0, nEvent);
    
    int eventsize = PX->size();
    int pid = 0;
    for(int ich = 0; ich < NCH; ich++)
    {
      if(charge_mode == 1 && ich == 0) continue;
      if(charge_mode == 2 && ich == 1) continue;
      (ich == 0) ? pid = -211 : pid = 211;
      //std::cout << "Now pid: " << pid << std::endl;
      for (int j = 0; j < eventsize; j++)
      {
        if (ID->at(j) == pid /* || ID->at(j) == -211 */) {
          auto Ej = energy->at(j); //TMath::Sqrt(pj*pj + Mass2_pi);
          for (int k = j+1; k < eventsize; k++)
          {
            if (ID->at(k) == pid /* || ID->at(k) == -211 */) {
              //auto pTk = TMath::Sqrt(PX->at(k)*PX->at(k) + PY->at(k)*PY->at(k));
              //auto pk = TMath::Sqrt(pTk*pTk + PZ->at(k)*PZ->at(k));
              auto Ek = energy->at(k); //TMath::Sqrt(pk*pk + Mass2_pi);

              auto t = time->at(j) - time->at(k);
              auto Rx = rx->at(j) - rx->at(k);
              auto Ry = ry->at(j) - ry->at(k);
              auto Rz = rz->at(j) - rz->at(k);
                

              auto K0 = 0.5 * (Ej + Ek);
              auto Kx = 0.5 * (PX->at(j) + PX->at(k));
              auto Ky = 0.5 * (PY->at(j) + PY->at(k));
              auto Kz = 0.5 * (PZ->at(j) + PZ->at(k));
                
              auto Qx = PX->at(j) - PX->at(k);
              auto Qy = PY->at(j) - PY->at(k);

              auto qZLCMS2 = (PZ->at(j)*Ek-PZ->at(k)*Ej)*(PZ->at(j)*Ek-PZ->at(k)*Ej) / (K0*K0 - Kz*Kz);
              auto qLCMS = TMath::Sqrt(Qx*Qx + Qy*Qy + qZLCMS2);

              auto KT = TMath::Sqrt(Kx*Kx + Ky*Ky);
              KTdist->Fill(KT); // retain the full pre-cut diagnostic distribution
              int iKT = kTBinNum(KT);
              if(iKT < 0) continue;
              auto cosphi = Kx / KT;
              auto sinphi = Ky / KT;

              auto rhoout = Rx * cosphi + Ry * sinphi - (KT / (K0*K0 - Kz*Kz)) * (K0 * t - Kz * Rz);
              auto rhoside = - Rx * sinphi + Ry * cosphi;
              auto rholong = (K0 * Rz - Kz * t) / TMath::Sqrt(K0*K0 - Kz*Kz);

              auto rho2 = rhoout*rhoout + rhoside*rhoside + rholong*rholong;

              /*
              if((iKT == 0 || iKT == 15) && qLCMS < ktbins[iKT+1])
                std::cout << "qlcms: " << qLCMS << " kT: " << ktbins[iKT] << "-" << ktbins[iKT+1] << std::endl;
              */
              
              // CUTS
              double pTj = TMath::Sqrt(PX->at(j)*PX->at(j) + PY->at(j)*PY->at(j));
              double pTk = TMath::Sqrt(PX->at(k)*PX->at(k) + PY->at(k)*PY->at(k));
              if(pTj<0.15 || pTj>=1.0) continue;
              if(pTk<0.15 || pTk>=1.0) continue;
              double pzj = PZ->at(j);
              double pzk = PZ->at(k);
              double pj = TMath::Sqrt(pTj*pTj + pzj*pzj);
              double pk = TMath::Sqrt(pTk*pTk + pzk*pzk);
              double etaj = 0.5 * TMath::Log(TMath::Abs((pj + pzj) / (pj - pzj)));
              double etak = 0.5 * TMath::Log(TMath::Abs((pk + pzk) / (pk - pzk)));
              double yj = 0.5 * TMath::Log((Ej + pzj) / (Ej - pzj)) - 2.02172; // so far only manually, based on bporfy TODO boost automatically to CMS
              double yk = 0.5 * TMath::Log((Ek + pzk) / (Ek - pzk)) - 2.02172;
              if(k == j+1) // fill only once per pair
              {
                etadist->Fill(etaj);
                rapdist->Fill(yj);
              }
              if(fabs(etaj) > 1.) continue;
              if(fabs(etak) > 1.) continue;
              //if(fabs(yj) > 1.) continue;
              //if(fabs(yk) > 1.) continue;
              if(k == j+1) // fill only once per pair
              {
                post_etadist->Fill(etaj);
                post_rapdist->Fill(yj);
              }
              if(qLCMS < TMath::Sqrt(qLCMS_cut_values[qLCMS_cuttype] * TMath::Sqrt(KT*KT + Mass2_pi))) // 150 Mev * mT baseline; strict 50 MeV, loose 250 MeV
              //if(qLCMS < qmax[iKT]) 
              //if(qLCMS < ktbins[iKT+1]) 
              {
                // Fill TTree instead of histograms - much more compact storage
                out_ievent = ievent;
                out_ich = ich;
                out_iKT = iKT;
                out_rho = TMath::Sqrt(rho2);
                out_rhoout = TMath::Abs(rhoout);
                out_rhoside = TMath::Abs(rhoside);
                out_rholong = TMath::Abs(rholong);
                out_osl_signs = (rhoout < 0 ? 1 : 0) | (rhoside < 0 ? 2 : 0) | (rholong < 0 ? 4 : 0);
                pairTree->Fill();
              }
            } // only particles with matching PID
          } // Loop over pairs - loop for particle k
        } // only particles with matching PID
      } // end of track loop - loop for particle j
    } // end of charge loop
  } // end of event loop
  std::cout << std::endl;
  outFile->cd();
  // Counts include events with no accepted pairs, which cannot be inferred from pairTree.
  TParameter<Long64_t>("source_events", nEvent).Write();
  TParameter<int>("pair_source_schema", useFFTmethod ? 2 : 1).Write();
  TParameter<int>("source_charge_mode", charge_mode).Write();
  TParameter<int>("source_qlcms_cut", qLCMS_cuttype).Write();
  TParameter<double>("source_kt_min", ktbins[0]).Write();
  TParameter<double>("source_kt_max", ktbins[NKT]).Write();
  TNamed("source_energy", argv[1]).Write();
  TNamed("source_centrality", centleg[11]).Write();
  pairTree->Write();
  KTdist->Write();
  rapdist->Write();
  post_rapdist->Write();
  etadist->Write();
  post_etadist->Write();

  outFile->Close();

  return 0;
}
