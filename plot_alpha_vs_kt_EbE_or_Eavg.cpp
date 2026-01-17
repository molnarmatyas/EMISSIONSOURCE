// root.exe -b -q plot_alpha_vs_kt_EbE_or_Eavg.cpp\(\"7.7GeV\"\)
// !!! note: not K_T but m_T

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>

#include "header_for_all_emissionsource.h"
//const int NCENT = 10; // number of centrality classes
//const int NKT = 10;
//double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
//const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};
//const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};

//const double Mass2_pi = 0.019479835;
//const double Mass2_ka = 0.24371698032;

// ---------  MAIN  -----------------------------------------
void plot_alpha_vs_kt_EbE_or_Eavg(const char* energy="9.2GeV", bool urqmd=false, int NEVT_AVG=1) 
{
  double mtbins[NKT+1] = {0};
  // Convert from k_T bin limits to m_T values
  for(int ii=0; ii<NKT+1; ii++)
  {
    mtbins[ii] = sqrt(ktbins[ii]*ktbins[ii] + Mass2_pi);
    cerr << "mtbin[" << ii << "] = " << mtbins[ii] << endl;
  }
  
  for(int icent=0; icent<NCENT+2; icent++)
  { 
    if(urqmd) {
      if(icent!=11) continue; // only 0-10% available for urqmd
    }
    const char* isurqmd = urqmd ? "UrQMD_" : "EPOS_"; // add "urqmd_" to the filename if urqmd is true
    //const char* uppercase = urqmd ? "LCMS" : "lcms"; // OK, now I really hate not building my own stuff from the scratch. FIXME with every other single thing before... HOORAY FIXED IT
    TFile *f = new TFile(Form("levyfit/results/%sonedfitresults_lcms_cent%s_%s_AVG%d.root", isurqmd, centleg[icent],energy,NEVT_AVG), "READ");

    if (!f || f->IsZombie()) {
      std::cerr << "Error opening file" << std::endl;
      return;
    }
    Double_t alpha_vs_kt[NKT] = {0};
    Double_t alpha_errdn_vs_kt[NKT] = {0};
    Double_t alpha_errup_vs_kt[NKT] = {0};
    Double_t R_vs_kt[NKT] = {0};
    Double_t R_errdn_vs_kt[NKT] = {0};
    Double_t R_errup_vs_kt[NKT] = {0};
    Double_t N_vs_kt[NKT] = {0}; 
    Double_t N_errdn_vs_kt[NKT] = {0};
    Double_t N_errup_vs_kt[NKT] = {0};
    for(int ikt=0; ikt<NKT; ikt++)
    {
      TH2F* alphaR = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", ikt));
      alpha_vs_kt[ikt] = alphaR->GetMean(1);
      alpha_errdn_vs_kt[ikt] = alphaR->GetStdDev(1);//GetMeanError(1);
      alpha_errup_vs_kt[ikt] = alphaR->GetStdDev(1);//GetMeanError(1);
      R_vs_kt[ikt] = alphaR->GetMean(2);
      R_errdn_vs_kt[ikt] = alphaR->GetStdDev(2);//GetMeanError(2);
      R_errup_vs_kt[ikt] = alphaR->GetStdDev(2);//GetMeanError(2);
      TH1F* Nhist = (TH1F*)f->Get(Form("Nhist_ikt%i", ikt));
      N_vs_kt[ikt] = Nhist->GetMean();
      N_errdn_vs_kt[ikt] = Nhist->GetStdDev();
      N_errup_vs_kt[ikt] = Nhist->GetStdDev();
    }
    
    double ktbin_centers[NKT];
    for(int ii=0; ii<NKT; ii++)
    {
      // K_T bin centers
      ktbin_centers[ii] = (mtbins[ii+1]+mtbins[ii])/2;
    }
    double binWidths[NKT];
    for(int ii=0; ii<NKT; ii++)
    {
      // K_T bin centers
      binWidths[ii] = mtbins[ii+1]-mtbins[ii];
    }
 
    // ....


    // Compute the asymmetric errors for TGraphAsymmErrors
    double xLow[NKT], xHigh[NKT];
    for (int i = 0; i < NKT; ++i) {
      xLow[i] = binWidths[i]/2;
      xHigh[i] = binWidths[i]/2;
    }

    // Create the TGraphAsymmErrors object
    TGraphAsymmErrors *graphalpha = new TGraphAsymmErrors(NKT, ktbin_centers, alpha_vs_kt, xLow, xHigh, alpha_errdn_vs_kt, alpha_errdn_vs_kt);
    TGraphAsymmErrors *graphR     = new TGraphAsymmErrors(NKT, ktbin_centers, R_vs_kt, xLow, xHigh, R_errdn_vs_kt, R_errup_vs_kt);
    TGraphAsymmErrors *graphN     = new TGraphAsymmErrors(NKT, ktbin_centers, N_vs_kt, xLow, xHigh, N_errdn_vs_kt, N_errup_vs_kt);

    // --- alpha ---
    // Style settings
    //graphalpha->SetTitle(Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s, %s%%;K_{T} (GeV/c);#alpha",energy, centleg[icent]));
    graphalpha->SetTitle(Form("#alpha(m_{T}), #sqrt{s_{NN}}=%s, %s%%;m_{T} (GeV/c);#alpha",energy, centleg[icent]));
    graphalpha->SetMarkerStyle(20);
    graphalpha->SetMarkerSize(1.2);
    graphalpha->SetLineWidth(2);

    // Draw the graph
    TCanvas *c = new TCanvas("c", Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s",energy), 800, 600);
    c->SetGrid();
    gStyle->SetOptStat(0);
    //graph->GetYaxis()->SetRangeUser(0,2);
    // Customize the axes
    graphalpha->GetXaxis()->SetLimits(0.0, 0.7); // Set X-axis range
    graphalpha->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y-axis range
    graphalpha->Draw("AP");

    // Save the plot to a file
    c->SaveAs(Form("figs/%salpha_vs_kt_%s_cent%s.png", isurqmd, energy, centleg[icent]));

    // --- R ---

    // Style settings
    //graphR->SetTitle(Form("R(K_{T}), #sqrt{s_{NN}}=%s, %s%%;K_{T} (GeV/c);R",energy, centleg[icent]));
    graphR->SetTitle(Form("R(m_{T}), #sqrt{s_{NN}}=%s, %s%%;m_{T} (GeV/c);R",energy, centleg[icent]));
    graphR->SetMarkerStyle(20);
    graphR->SetMarkerSize(1.2);
    graphR->SetLineWidth(2);

    // Draw the graph
    //TCanvas *c = new TCanvas("c", Form("#R(K_{T}), #sqrt{s_{NN}}=%s",energy), 800, 600);
    c->SetGrid();
    gStyle->SetOptStat(0);
    //graph->GetYaxis()->SetRangeUser(0,2);
    // Customize the axes
    graphR->GetXaxis()->SetLimits(0.0, 0.7); // Set X-axis range
    graphR->GetYaxis()->SetRangeUser(0.0, 10.0); // Set Y-axis range FIXME from 3.0 to 10.0?
    graphR->Draw("AP");

    // Save the plot to a file
    c->SaveAs(Form("figs/%sR_vs_kt_%s_cent%s.png", isurqmd, energy, centleg[icent]));

    
    // --- N ---

    // Style settings
    graphN->SetTitle(Form("N(K_{T}), #sqrt{s_{NN}}=%s, %s%%;m_{T} (GeV/c);N",energy, centleg[icent]));
    graphN->SetMarkerStyle(20);
    graphN->SetMarkerSize(1.2);
    graphN->SetLineWidth(2);

    // Draw the graph
    //TCanvas *c = new TCanvas("c", Form("N(K_{T}), #sqrt{s_{NN}}=%s",energy), 800, 600);
    c->SetGrid();
    gStyle->SetOptStat(0);
    //graph->GetYaxis()->SetNangeUser(0,2);
    // Customize the axes
    graphN->GetXaxis()->SetLimits(0.0, 0.7); // Set X-axis range
    //graphN->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y-axis range
    graphN->Draw("AP");

    // Save the plot to a file
    c->SaveAs(Form("figs/%sN_vs_kt_%s_cent%s.png", isurqmd, energy, centleg[icent]));
    

    // Save graphs to a .root file
    TFile* outFile = new TFile(Form("alphaNR_vs_kt/%salphaNR_graphs_%s_cent%s.root", isurqmd, energy, centleg[icent]), "RECREATE");
    graphalpha->Write(Form("graphalpha%s",centleg[icent]));
    graphN->Write(Form("graphN%s",centleg[icent]));
    graphR->Write(Form("graphR%s",centleg[icent]));
    outFile->Close();

    std::cout << Form("Graphs saved to %salphaNR_graphs_%s_cent%s.root", isurqmd, energy, centleg[icent]) << std::endl;

    // Clean up
    delete c;
    delete graphalpha;
    delete graphN;
    delete graphR;
  }
}
