#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <string>

void plot_param_vs_nevt_avg() {
    // Define the NEVT_AVG values to process
    std::vector<int> nevt_avg = {1, 50, 100, 200, 500, 1000, 5000, 10000};
    std::vector<std::string> energies = {"27"}; //"3p0", 
    
    // Create vectors to store parameters for different energies
    std::vector<TGraphErrors*> alpha_graphs;
    std::vector<TGraphErrors*> R_graphs;
    
    // Colors for different energies
    int colors[] = {kRed, kBlue};
    
    for(size_t e = 0; e < energies.size(); e++) {
        std::vector<double> alphas, alphaErrs;
        std::vector<double> Rs, RErrs;
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
            TH2F* h2 = (TH2F*)file->Get("alpha_vs_R_all");
            if(!h2) {
                printf("Could not find histogram in file: %s\n", filename.Data());
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
            
            file->Close();
        }
        
        // Create graphs
        TGraphErrors* g_alpha = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &alphas[0], &nevt_errs[0], &alphaErrs[0]);
        TGraphErrors* g_R = new TGraphErrors(nevt_points.size(),
            &nevt_points[0], &Rs[0], &nevt_errs[0], &RErrs[0]);
            
        g_alpha->SetMarkerStyle(20);
        g_alpha->SetMarkerColor(colors[e]);
        g_alpha->SetLineColor(colors[e]);
        g_R->SetMarkerStyle(20);
        g_R->SetMarkerColor(colors[e]);
        g_R->SetLineColor(colors[e]);
        
        alpha_graphs.push_back(g_alpha);
        R_graphs.push_back(g_R);
    }
    
    // Create and divide canvas
    TCanvas* c1 = new TCanvas("c1", "Parameter vs NEVT_AVG", 1200, 600);
    c1->Divide(2,1);
    
    // Plot alpha
    c1->cd(1);
    gPad->SetLogx();
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    TLegend* leg1 = new TLegend(0.2, 0.2, 0.4, 0.3);
    
    for(size_t i = 0; i < alpha_graphs.size(); i++) {
        if(i == 0) {
            alpha_graphs[i]->SetTitle("Levy #alpha vs NEVT_AVG");
            alpha_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            alpha_graphs[i]->GetYaxis()->SetTitle("#alpha");
            alpha_graphs[i]->Draw("APE");
        } else {
            alpha_graphs[i]->Draw("PE same");
        }
        leg1->AddEntry(alpha_graphs[i], Form("%s GeV", energies[i].c_str()), "PE");
    }
    leg1->Draw();
    
    // Plot R
    c1->cd(2);
    gPad->SetLogx();
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    TLegend* leg2 = new TLegend(0.2, 0.2, 0.4, 0.3);
    
    for(size_t i = 0; i < R_graphs.size(); i++) {
        if(i == 0) {
            R_graphs[i]->SetTitle("Levy R vs NEVT_AVG");
            R_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_graphs[i]->GetYaxis()->SetTitle("R [fm]");
            R_graphs[i]->Draw("APE");
        } else {
            R_graphs[i]->Draw("PE same");
        }
        leg2->AddEntry(R_graphs[i], Form("%s GeV", energies[i].c_str()), "PE");
    }
    leg2->Draw();
    
    c1->SaveAs("figs/levy_params_vs_nevt_avg.png");
}