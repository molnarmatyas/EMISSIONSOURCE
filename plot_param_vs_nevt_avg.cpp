#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <vector>
#include <string>
#include "header_for_all_emissionsource.h"

void plot_param_vs_nevt_avg(int ikt =-1) {
    // Define the NEVT_AVG values to process
    std::vector<int> nevt_avg = {10, 25, 50, 100, 200, 500, 1000, 5000, 10000}; // removed 1, not meaningful for low energies, runs too long
    std::vector<std::string> energies = {"3p0","3p2","3p5","3p9","4p5","7p7","9p2","11p5","14p5","19p6","27"};
    
    // Create vectors to store parameters for different energies
    std::vector<TGraphErrors*> alpha_graphs;
    std::vector<TGraphErrors*> R_graphs;
    std::vector<TGraphErrors*> N_graphs;
    std::vector<TGraphErrors*> conflev_graphs;
    
    // Colors for different energies
    int colors[] = {kRed, kBlue, kTeal+3, kMagenta, kOrange+7, kGreen+2, kCyan+2, kViolet+2, kPink+7, kGray+2, kAzure+7};
    
    for(size_t e = 0; e < energies.size(); e++) {
        std::vector<double> alphas, alphaErrs;
        std::vector<double> Rs, RErrs;
        std::vector<double> Ns, NErrs;
        std::vector<double> conflevs, conflevErrs;
        std::vector<double> nevt_points, nevt_errs;
        
        for(int nevt : nevt_avg) {
            // Construct filename
            TString filename = Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root",
                                  energies[e].c_str(), nevt);
            
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
            
            alphas.push_back(alpha);
            alphaErrs.push_back(alphaErr);
            Rs.push_back(R);
            RErrs.push_back(RErr);
            nevt_points.push_back(nevt);
            nevt_errs.push_back(0);

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
            if(ikt_N < 0) {
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
            
            double N = hN->GetMean();
            double NErr = hN->GetRMS() / sqrt(hN->GetEntries());
            Ns.push_back(N);
            NErrs.push_back(NErr);
            // Confidence level is stored together and for all kT bins separately as well
            const char* conflev_histname = (ikt>=0) ? Form("confidencehist_ikt%i", ikt) : "confidencehist_all";
            TH1F* hConfLev = (TH1F*)file->Get(conflev_histname);
            double conflev = hConfLev->GetMean();
            double conflevErr = hConfLev->GetRMS() / sqrt(hConfLev->GetEntries());
            conflevs.push_back(conflev);
            conflevErrs.push_back(conflevErr);
            
            file->Close();
        }
        
        // Create graphs
        TGraphErrors* g_alpha = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &alphas[0], &nevt_errs[0], &alphaErrs[0]);
        TGraphErrors* g_R = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &Rs[0], &nevt_errs[0], &RErrs[0]);
        TGraphErrors* g_N = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &Ns[0], &nevt_errs[0], &NErrs[0]);
        TGraphErrors* g_conflev = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &conflevs[0], &nevt_errs[0], &conflevErrs[0]);
            
        // Use solid markers and slightly transparent, thicker lines to connect points
        int baseCol = colors[e];
        int transCol = TColor::GetColorTransparent(baseCol, 0.4); // 40% opacity

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
    TLegend* leg1 = new TLegend(0.2, 0.2, 0.4, 0.5);
    
    for(size_t i = 0; i < alpha_graphs.size(); i++) {
        if(i == 0) {
            alpha_graphs[i]->SetTitle("Levy #alpha vs NEVT_AVG");
            alpha_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            alpha_graphs[i]->GetYaxis()->SetTitle("#alpha");
            alpha_graphs[i]->Draw("APE");
            // overlay a connecting line so points remain visible with error bars
            alpha_graphs[i]->Draw("L same");
        } else {
            alpha_graphs[i]->Draw("PE same");
            alpha_graphs[i]->Draw("L same");
        }
        leg1->AddEntry(alpha_graphs[i], Form("%s GeV", energies[i].c_str()), "LP");
    }
    leg1->Draw();
    
    // Plot R
    c1->cd(2);
    gPad->SetLogx();
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    TLegend* leg2 = new TLegend(0.2, 0.2, 0.4, 0.5);
    
    for(size_t i = 0; i < R_graphs.size(); i++) {
        if(i == 0) {
            R_graphs[i]->SetTitle("Levy R vs NEVT_AVG");
            R_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_graphs[i]->GetYaxis()->SetTitle("R [fm]");
            R_graphs[i]->Draw("APE");
            R_graphs[i]->Draw("L same");
        } else {
            R_graphs[i]->Draw("PE same");
            R_graphs[i]->Draw("L same");
        }
        leg2->AddEntry(R_graphs[i], Form("%s GeV", energies[i].c_str()), "LP");
    }
    leg2->Draw();

    // Plot N
    c1->cd(3);
    gPad->SetLogx();
    // Place N legend in the upper-right with same size as alpha/R legend (width=0.2, height=0.3)
    TLegend* leg3 = new TLegend(0.6, 0.6, 0.8, 0.9);

    for(size_t i = 0; i < N_graphs.size(); i++) {
        if(i == 0) {
            N_graphs[i]->SetTitle("Levy N vs NEVT_AVG");
            N_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            N_graphs[i]->GetYaxis()->SetTitle("N");
            N_graphs[i]->GetYaxis()->SetRangeUser(0.95, 1.12); // N between 0 and 1
            N_graphs[i]->Draw("APE");
            N_graphs[i]->Draw("L same");
        } else {
            N_graphs[i]->Draw("PE same");
            N_graphs[i]->Draw("L same");
        }
        leg3->AddEntry(N_graphs[i], Form("%s GeV", energies[i].c_str()), "LP");
    }
    leg3->Draw();

    // Plot confidence level
    c1->cd(4);
    gPad->SetLogx();
    // Place confidence legend in the upper-right with same size as alpha/R legend
    TLegend* leg4 = new TLegend(0.6, 0.6, 0.8, 0.9);

    for(size_t i = 0; i < conflev_graphs.size(); i++) {
        if(i == 0) {
            conflev_graphs[i]->SetTitle("Fit Confidence Level vs NEVT_AVG");
            conflev_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            conflev_graphs[i]->GetYaxis()->SetTitle("Confidence Level");
            conflev_graphs[i]->GetYaxis()->SetRangeUser(0, 1.05); // Confidence level between 0 and 1
            conflev_graphs[i]->Draw("APE");
            conflev_graphs[i]->Draw("L same");
        } else {
            conflev_graphs[i]->Draw("PE same");
            conflev_graphs[i]->Draw("L same");
        }
        leg4->AddEntry(conflev_graphs[i], Form("%s GeV", energies[i].c_str()), "LP");
    }
    leg4->Draw();

    // Final, save
    c1->SaveAs("figs/levy_params_vs_nevt_avg.png");
}