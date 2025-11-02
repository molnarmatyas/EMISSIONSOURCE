// root.exe -b -q collect_alphavsR.cpp

#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

const int NENERGIES = 6;
const int NKT = 10;
const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
const char* energies[NENERGIES] = {"7.7GeV", "9.2GeV", "11.5GeV", "14.5GeV", "19.6GeV", "27GeV"};

void plotthem(const char* energy, int ikt)
{
  TFile *file = TFile::Open(Form("levyfit/results/onedfitresults_lcms_cent0-10_%s.root",energy));
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return;
  }

  TH2F* alpha_vs_R_ikt;
  if(ikt==NKT)
  {
    alpha_vs_R_ikt = (TH2F*)file->Get("alpha_vs_R_all");
    alpha_vs_R_ikt->SetTitle(Form("#alpha vs R, all K_T, #sqrt{s_{NN}}=%s", energy));
  }
  else
  {
    alpha_vs_R_ikt = (TH2F*)file->Get(Form("alpha_vs_R_ikt%i",ikt));
    alpha_vs_R_ikt->SetTitle(Form("#alpha vs R, %.3f<K_{T}<%.3f, #sqrt{s_{NN}}=%s",ktbins[ikt], ktbins[ikt + 1], energy));
  }

  if (!alpha_vs_R_ikt)
  {
    std::cerr << "Error getting histogram" << std::endl;
    return;
  }
 
  alpha_vs_R_ikt->GetXaxis()->SetTitle("#alpha");
  alpha_vs_R_ikt->GetYaxis()->SetTitle("R [fm]");
  //alpha_vs_R_ikt->GetXaxis()->SetRangeUser(0.6,1.8);
  //alpha_vs_R_ikt->GetYaxis()->SetRangeUser(2.5,11.);
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s",energy), 800, 600);
  alpha_vs_R_ikt->Draw("COLZ");
  if(ikt==NKT)
  {
    c->SaveAs(Form("figs/alpha_vs_R_cent0-10_%s_iktall.png",energy));
  }
  else
  {
    c->SaveAs(Form("figs/alpha_vs_R_cent0-10_%s_ikt%i.png",energy,ikt));
  }
  //c->SetGrid();
  delete c;
  file->Close();

}

void collect_alphavsR()
{
  for(int ienergy=0; ienergy<NENERGIES; ienergy++)
  {
    for(int ikt=0; ikt<NKT; ikt++)
    {
      plotthem(energies[ienergy], ikt);
    }
    // all kT bins
    plotthem(energies[ienergy], NKT);
  }
}