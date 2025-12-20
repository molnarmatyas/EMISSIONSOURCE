// This ROOT macro collects the rho_out, rho_side, rho_long histograms from
// analysed/UrQMD*.root output file, performs the given averaging over events
// and prints the results to a three-panel canvas, saves it in figs/ folder.


#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "header_for_all_emissionsource.h"

const float rightmargins[9] = {0.0025,0.0025,0.04,0.,0.,0.04,0.,0.,0.04};
const float leftmargins[9]  = {0.2,0.0005,0.0005,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.03,0.03,0.03,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.16,0.16,0.16,0.,0.,0.,0.2,0.2,0.2};
const int NPARS = 5;
const int NOSL = 3;

// Simple text progress bar
static void PrintProgress(int current, int total)
{
  if (total <= 0) return;
  const int barWidth = 40;
  double fraction = static_cast<double>(current) / static_cast<double>(total);
  int pos = static_cast<int>(barWidth * fraction);
  std::cout << "[";
  for (int i = 0; i < barWidth; ++i)
  {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << std::fixed << std::setprecision(1) << (fraction * 100.0) << "% ("
            << current << "/" << total << ")\r";
  std::cout.flush();
}

// Convert a 1D rho histogram to a physical density D(rho)
// For directional projections (out/side/long):
//  - divide by bin width to get density per fm
//  - normalize to unit integral
// NOTE: Do NOT apply the 1/(4π r^2) factor here; that is only for radial D(rho).
static void MakeDensity1D(TH1* h, const char* xaxisTitle)
{
  if (!h) return;
  double integral_counts = h->Integral();
  if (integral_counts <= 0) return;

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    double c = h->GetBinContent(ib);
    double e = h->GetBinError(ib);
    double w = h->GetXaxis()->GetBinWidth(ib);
    if (w <= 0) continue;
    h->SetBinContent(ib, c / w);
    h->SetBinError(ib, e / w);
  }
  // Normalize to unit integral (over x)
  h->Scale(1.0 / integral_counts);

  h->SetStats(0);
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xaxisTitle);
  h->GetYaxis()->SetTitle("D(#rho)");
  h->GetXaxis()->SetRangeUser(0.2, 5000.);
  h->SetMinimum(5.e-13);
  h->SetMaximum(1.e-2);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetTitleOffset(1.2);
}

// Helper: load per-event, per-kT, per-OSL histogram and sum charges into one event histogram
static TH1D* LoadEventOSL(TFile* f, int ievt, int ikt, const char* osl_tag)
{
  // Hist name pattern from pairsource_urqmd.cc: D_<osl>_lcms_ev%i_ch%i_KT%i
  // where <osl> is out/side/long
  if (!f) return nullptr;
  TH1D* h_evt_sum = nullptr;
  int NCH = 2; // sum over charges
  for (int ich = 0; ich < NCH; ++ich)
  {
    TString hname = Form("D%s_lcms_ev%d_ch%d_KT%d", osl_tag, ievt, ich, ikt);
    TH1D* h = dynamic_cast<TH1D*>(f->Get(hname));
    if (!h) continue;
    if (!h_evt_sum)
    {
      h_evt_sum = dynamic_cast<TH1D*>(h->Clone(Form("evt_%s_ievt%d_ikt%d", osl_tag, ievt, ikt)));
      if (h_evt_sum) h_evt_sum->SetDirectory(nullptr);
      //if (h_evt_sum) h_evt_sum->Sumw2();
    }
    else
    {
      h_evt_sum->Add(h);
    }
  }
  return h_evt_sum;
}

int print_drho_osl3d(int nevt_avg=5000)
{
  TFile *f_in = TFile::Open("analysed/UrQMD_3d_source_0-10cent_all_30A.root", "read");
  if (!f_in || f_in->IsZombie()) {
    std::cerr << "Error: Could not open input file." << std::endl;
    return 1;
  }

  // Ensure output directory exists
  //gSystem->mkdir("figs", kTRUE);

  const char* osl_tags[3] = {"_out", "_side", "_long"};
  const char* xaxis_labels[3] = {"#rho_{out} [fm]", "#rho_{side} [fm]", "#rho_{long} [fm]"};

  for(int ikt=0; ikt<NKT; ikt++)
  {
    if(ikt!=6) continue; // Only kT bin 6 for now
    // Create a canvas with three pads
    TCanvas *c1 = new TCanvas(Form("c1_ikt%d", ikt), "Rho OSL 3D", 1200, 600);
    //c1->Divide(3, 1);
    TPad *pad[NOSL];
  	pad[0] = new TPad("pad1","",0.,0.,0.37,1.);
  	pad[1] = new TPad("pad2","",0.37,0.,0.68,1.);
  	pad[2] = new TPad("pad3","",0.68,0.,1.,1.);
    c1->cd();
    for(int ipad = 0; ipad < NOSL; ipad++)
      pad[ipad]->Draw();

    // Loop over the three OSL histograms
    for (int i = 0; i < 3; i++) {
      // c1->cd(i + 1);
      // gPad->SetLogx(1);
      // gPad->SetLogy(1);
      c1->cd();
      int ipad = i;
      gStyle->SetLabelSize(0.06,"Y");
      c1->SetLogx(1);
      c1->SetLogy(1);
      pad[ipad]->cd();
      gPad->SetRightMargin(rightmargins[ipad]);
      gPad->SetLeftMargin(leftmargins[ipad]);
      gPad->SetTopMargin(topmargins[ipad]);
      pad[ipad]->SetBottomMargin(bottommargins[ipad]);
      gPad->SetLogx(1);
      gPad->SetLogy(1);

      // Average over events: sum charges within each event, then average across events
      TH1D* h_avg = nullptr;
      int nEventsUsed = 0;
      for (int ievt = 0; ievt < nevt_avg; ++ievt)
      {
        TH1D* h_evt = LoadEventOSL(f_in, ievt, ikt, osl_tags[i]);
        if (!h_evt)
        {
          // If the very first event is missing, likely no such histograms for this kT
          if (ievt == 0) break;
          continue;
        }
        if (!h_avg)
        {
          h_avg = dynamic_cast<TH1D*>(h_evt->Clone(Form("avg%s_ikt%d", osl_tags[i], ikt)));
          if (h_avg) h_avg->SetDirectory(nullptr);
          //if (h_avg) h_avg->Sumw2();
        }
        else
        {
          h_avg->Add(h_evt);
        }
        nEventsUsed++;
        delete h_evt;

        // Progress display every 100 events and at the end
        if ( (ievt % 100 == 0) || (ievt == nevt_avg - 1) )
          PrintProgress(ievt + 1, nevt_avg);
      }

      // Newline after finishing this OSL component to keep output tidy
      if (nEventsUsed > 0)
        std::cout << std::endl;

      // Average the histogram
      if (h_avg && nEventsUsed > 0) {
        h_avg->Scale(1.0 / nEventsUsed);

        // Convert to physical D(rho) for directional projections
        MakeDensity1D(h_avg, xaxis_labels[i]);
        h_avg->GetXaxis()->SetRangeUser(0.2, 200.);
        h_avg->SetMinimum(5.e-7);
        h_avg->SetMaximum(3.0);

        // Set descriptive title (no axis text here)
        h_avg->SetTitle(Form("D(%s), k_{T} in [%.3g, %.3g] GeV/c; Averaged over %d events",
                             (i==0?"#rho_{out}":(i==1?"#rho_{side}":"#rho_{long}")),
                             ktbins[ikt], ktbins[ikt+1], nEventsUsed));
        h_avg->Draw("PE");
      } else {
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
        tl.DrawLatex(0.2, 0.5, "No histograms found for this selection.");
      }
    }

    // Save the canvas
    c1->SaveAs(Form("figs/rho_osl3d_ikt%d_avg%d.png", ikt, nevt_avg));
    delete c1;
  }




  return 0;
}
