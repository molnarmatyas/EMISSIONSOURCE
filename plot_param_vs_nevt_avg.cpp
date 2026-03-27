#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMarker.h>
#include <vector>
#include <string>
#include <cmath>
#include "header_for_all_emissionsource.h"

void circle_defaultnevtavg_point(TGraphErrors* param_vs_nevt_graph, int ienergy)
{
    if(!param_vs_nevt_graph || param_vs_nevt_graph->GetN() <= 0) {
        return;
    }

    double x = NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]];
    int ipoint = -1;
    
    // Find the default point in the graph
    for(int i=0; i<param_vs_nevt_graph->GetN(); i++)
    {
        if(std::abs(param_vs_nevt_graph->GetX()[i] - x) < 1e-9)
        {
            ipoint = i;
            break;
        }
    }

    if(ipoint < 0) {
        return;
    }

    double y = param_vs_nevt_graph->GetY()[ipoint];

    // Open-circle marker size is in screen space, so it stays circular on linear/log axes.
    TMarker* marker = new TMarker(x, y, 107); // 107 is an open circle marker, thicker
    marker->SetMarkerColor(kGreen-8);
    marker->SetMarkerSize(2.8);
    marker->Draw("same");
    
}

// ------------------ MAIN ------------------------
void plot_param_vs_nevt_avg(int ikt =-1) {
    // NEVT_AVG values defined in header_for_all_emissionsource.h as NEVT_AVGsyst array and NEVT_AVG_DEFAULT index
    std::vector<int> nevt_avg;
    for(int i=0; i<NEVTAVGS; i++)
    {
        nevt_avg.push_back(NEVT_AVGsyst[i]);
    }
        
    // Create vectors to store parameters for different energies
    std::vector<TGraphErrors*> alpha_graphs;
    std::vector<TGraphErrors*> R_graphs;
    std::vector<TGraphErrors*> N_graphs;
    std::vector<TGraphErrors*> conflev_graphs;
    
    // Colors for different energies
    int colors[NENERGIES] = {kRed, kBlue, kTeal+3, kMagenta, kOrange+7, kGreen+2, kCyan+2, kViolet+2, kPink+7, kGray+2, kAzure+7};
    
    for(size_t e = 0; e < NENERGIES; e++) {
        std::vector<double> alphas, alphaErrs;
        std::vector<double> Rs, RErrs;
        std::vector<double> Ns, NErrs;
        std::vector<double> conflevs, conflevErrs;
        std::vector<double> nevt_points, nevt_errs;
        
        for(int nevt : nevt_avg) {
            // Construct filename
            TString filename = Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root",
                                  energies[e], nevt);
            
            TFile* file = TFile::Open(filename);
            if(!file) {
                printf("Could not open file: %s\n", filename.Data());
                continue;
            }
            
            // Get the histogram (using alpha_vs_R_all for all kT bins combined)
            const char* histname = "alpha_vs_R_all";
            if(ikt>=0)
            {
                histname = Form("alpha_vs_R_ikt%i", ikt);
            }
            TH2F* h2 = (TH2F*)file->Get(histname);
            if(!h2) {
                printf("Could not find alpha_vs_R_all/alpha_vs_R_ikt# histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            // Check if histogram has entries
            if(h2->GetEntries() == 0) {
                printf("Histogram alpha_vs_R_all/alpha_vs_R_ikt# has no entries in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            
            // Get mean values and errors
            double alpha = h2->GetMean(1);  // X axis = alpha
            double alphaErr = h2->GetRMS(1) / sqrt(h2->GetEntries());
            double R = h2->GetMean(2);      // Y axis = R
            double RErr = h2->GetRMS(2) / sqrt(h2->GetEntries());
            
            // Same for lambda (in other convention, displayed here as N) and confidence level
            // N is for diff. ikt only
            int ikt_N = ikt;
            if(ikt_N < 0) ikt_N = 0; // just to get the first histogram, later will be added together
            TH1F* hN = (TH1F*)file->Get(Form("Nhist_ikt%i", ikt_N));
            if(!hN) {
                printf("Could not find Nhist_ikt# histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            // Check if combining all kT bins (original ikt was < 0)
            if(ikt < 0) {
                // Add all kT bins together
                for(int ikt_add = 1; ikt_add < NKT; ikt_add++) {
                    TH1F* hN_add = (TH1F*)file->Get(Form("Nhist_ikt%i", ikt_add));
                    if(hN_add) {
                        hN->Add(hN_add);
                    }
                }
                // Rename to avoid confusion
                hN->SetName("Nhist_allkt");
                hN->SetTitle("N, all kT bins combined");
            }
            
            // Check if hN has entries before computing mean/RMS
            if(hN->GetEntries() == 0) {
                printf("Nhist_ikt# is empty in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            
            double N = hN->GetMean();
            double NErr = hN->GetRMS() / sqrt(hN->GetEntries());
            // Confidence level is stored together and for all kT bins separately as well
            const char* conflev_histname = (ikt>=0) ? Form("confidencehist_ikt%i", ikt) : "confidencehist_all";
            TH1F* hConfLev = (TH1F*)file->Get(conflev_histname);
            if(!hConfLev) {
                printf("Could not find confidence histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            if(hConfLev->GetEntries() == 0) {
                printf("Confidence histogram is empty in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            double conflev = hConfLev->GetMean();
            double conflevErr = hConfLev->GetRMS() / sqrt(hConfLev->GetEntries());

            // Push all quantities only after all required histograms were found.
            alphas.push_back(alpha);
            alphaErrs.push_back(alphaErr);
            Rs.push_back(R);
            RErrs.push_back(RErr);
            Ns.push_back(N);
            NErrs.push_back(NErr);
            conflevs.push_back(conflev);
            conflevErrs.push_back(conflevErr);
            nevt_points.push_back(nevt);
            nevt_errs.push_back(0);
            
            file->Close();
        }
        
        // Create graphs
        TGraphErrors* g_alpha = nullptr;
        TGraphErrors* g_R = nullptr;
        TGraphErrors* g_N = nullptr;
        TGraphErrors* g_conflev = nullptr;
        if(!nevt_points.empty()) {
            g_alpha = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), alphas.data(), nevt_errs.data(), alphaErrs.data());
            g_R = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), Rs.data(), nevt_errs.data(), RErrs.data());
            g_N = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), Ns.data(), nevt_errs.data(), NErrs.data());
            g_conflev = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), conflevs.data(), nevt_errs.data(), conflevErrs.data());
        }
            
        // Use solid markers and slightly transparent, thicker lines to connect points
        int baseCol = colors[e];
        int transCol = TColor::GetColorTransparent(baseCol, 0.4); // 40% opacity

        if(g_alpha) {
            g_alpha->SetMarkerStyle(20);
            g_alpha->SetMarkerColor(baseCol);
            g_alpha->SetMarkerSize(1.0);
            g_alpha->SetLineColor(transCol);
            g_alpha->SetLineWidth(3);

            g_R->SetMarkerStyle(20);
            g_R->SetMarkerColor(baseCol);
            g_R->SetMarkerSize(1.0);
            g_R->SetLineColor(transCol);
            g_R->SetLineWidth(3);

            g_N->SetMarkerStyle(20);
            g_N->SetMarkerColor(baseCol);
            g_N->SetMarkerSize(1.0);
            g_N->SetLineColor(transCol);
            g_N->SetLineWidth(3);

            g_conflev->SetMarkerStyle(20);
            g_conflev->SetMarkerColor(baseCol);
            g_conflev->SetMarkerSize(1.0);
            g_conflev->SetLineColor(transCol);
            g_conflev->SetLineWidth(3);
        }
        
        alpha_graphs.push_back(g_alpha);
        R_graphs.push_back(g_R);
        N_graphs.push_back(g_N);
        conflev_graphs.push_back(g_conflev);
    }
    
    // Create and divide canvas
    TCanvas* c1 = new TCanvas("c1", "Parameter vs NEVT_AVG", 2400, 2400);
    c1->Divide(2,2);
    
    // Plot alpha
    c1->cd(1);
    gPad->SetLogx();
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    TLegend* leg1 = new TLegend(0.7, 0.1, 0.9, 0.4);
    bool first_alpha_drawn = false;
    
    for(size_t i = 0; i < alpha_graphs.size(); i++) {
        if(!alpha_graphs[i]) continue;
        if(!first_alpha_drawn) {
            alpha_graphs[i]->SetTitle("Levy #alpha vs NEVT_AVG");
            alpha_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            alpha_graphs[i]->GetYaxis()->SetTitle("#alpha");
            alpha_graphs[i]->GetYaxis()->SetRangeUser(1.3, 1.9);
            // assuming that lowest energies - plotted first - are highstat, this range setting should work
            //alpha_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg[NEVTAVGS-1]);
            alpha_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            alpha_graphs[i]->Draw("APE");
            // overlay a connecting line so points remain visible with error bars
            alpha_graphs[i]->Draw("L same");
            first_alpha_drawn = true;
        } else {
            alpha_graphs[i]->Draw("PE same");
            alpha_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(alpha_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg1->AddEntry(alpha_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg1->Draw();
    
    // Plot R
    c1->cd(2);
    gPad->SetLogx();
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    //TLegend* leg2 = new TLegend(0.1, 0.1, 0.3, 0.4);
    // Place R legend in the upper-right with same size as alpha legend (width=0.2, height=0.3)
    TLegend* leg2 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_R_drawn = false;
    
    for(size_t i = 0; i < R_graphs.size(); i++) {
        if(!R_graphs[i]) continue;
        if(!first_R_drawn) {
            R_graphs[i]->SetTitle("Levy R vs NEVT_AVG");
            R_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_graphs[i]->GetYaxis()->SetTitle("R [fm]");
            R_graphs[i]->GetYaxis()->SetRangeUser(3, 7); // R between 0 and 10 fm
            R_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            R_graphs[i]->Draw("APE");
            R_graphs[i]->Draw("L same");
            first_R_drawn = true;
        } else {
            R_graphs[i]->Draw("PE same");
            R_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(R_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg2->AddEntry(R_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg2->Draw();

    // Plot N
    c1->cd(3);
    gPad->SetLogx();
    // Place N legend in the upper-right with same size as alpha/R legend (width=0.2, height=0.3)
    TLegend* leg3 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_N_drawn = false;

    for(size_t i = 0; i < N_graphs.size(); i++) {
        if(!N_graphs[i]) continue;
        if(!first_N_drawn) {
            N_graphs[i]->SetTitle("Levy N vs NEVT_AVG");
            N_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            N_graphs[i]->GetYaxis()->SetTitle("N");
            N_graphs[i]->GetYaxis()->SetRangeUser(0.92, 1.12); // N between 0.98 and 1.12
            N_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            N_graphs[i]->Draw("APE");
            N_graphs[i]->Draw("L same");
            first_N_drawn = true;
        } else {
            N_graphs[i]->Draw("PE same");
            N_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(N_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg3->AddEntry(N_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg3->Draw();

    // Plot confidence level
    c1->cd(4);
    gPad->SetLogx();
    // Place confidence legend in the upper-right with same size as alpha/R legend
    TLegend* leg4 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_conflev_drawn = false;

    for(size_t i = 0; i < conflev_graphs.size(); i++) {
        if(!conflev_graphs[i]) continue;
        if(!first_conflev_drawn) {
            conflev_graphs[i]->SetTitle("Fit Confidence Level vs NEVT_AVG");
            conflev_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            conflev_graphs[i]->GetYaxis()->SetTitle("Confidence Level");
            conflev_graphs[i]->GetYaxis()->SetRangeUser(1e-15, 1.05); // Confidence level between 0 and 1
            conflev_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            // Set log scale for y-axis
            gPad->SetLogy();
            conflev_graphs[i]->Draw("APE");
            conflev_graphs[i]->Draw("L same");
            first_conflev_drawn = true;
        } else {
            conflev_graphs[i]->Draw("PE same");
            conflev_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(conflev_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg4->AddEntry(conflev_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg4->Draw();

    // Final, save
    const char* ikt_suffix = (ikt>=0) ? Form("_ikt%i", ikt) : "_allkt";
    c1->SaveAs(Form("figs/levy_params_vs_nevt_avg%s.png", ikt_suffix));
}
