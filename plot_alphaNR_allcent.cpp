// root.exe -b -q plot_alphaNR_allcent.cpp
// Plots alpha(sqrt(sNN)) and R(sqrt(sNN)) for all centralities on one canvas each

#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>

const int NCENT = 10; // number of centrality classes
const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};
const int colors[NCENT+2] = {kRed, kGray, kGreen, kCyan, kOrange, kBlack, kViolet, kPink, kBlue, kAzure, kSpring, kMagenta}; // kYellow

const char* energies[] = {"7.7GeV", "9.2GeV", "11.5GeV", "14.5GeV", "19.6GeV", "27GeV", "39GeV", "62.4GeV", "130GeV", "200GeV"};
const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
const double energydouble[NENERGIES] = {7.7, 9.2, 11.5, 14.5, 19.6, 27., 39., 62.4, 130., 200.};
const int NKT = 10;

const char* energies_urqmd[] = {"3p0","3p2","3p5","3p9","4p5","7p7","9p2","11p5","14p5","19p6","27"};
const int NENERGIES_urqmd = sizeof(energies_urqmd) / sizeof(energies_urqmd[0]);
const double energydouble_urqmd[NENERGIES_urqmd] = {3.0, 3.2, 3.5, 3.9, 4.5, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0};

// !!! set these to correct values before running !!!
int NEVT_AVG = 1000;

int calculate_syst_avg(int iparam, int icent, TGraphAsymmErrors* graph) {
    for(int ienergy=0; ienergy<NENERGIES; ienergy++) {
        TFile* f = new TFile(Form("levyfit/results/onedfitresults_lcms_cent%s_%s.root", centleg[icent], energies[ienergy]), "READ");
        /*
        // Create a buffer histogram
        TH2F* alphaR = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 0));
        TH2F* alphaR_summed = (TH2F*)alphaR->Clone("alpha_vs_R_summed");
        for(int ikt=1; ikt<NKT; ikt++) {
            alphaR_summed->Add((TH2F*)f->Get(Form("alpha_vs_R_ikt%i", ikt)));
        }
        */
        // For only m_T = 300 GeV !!!
        TH2F* alphaR = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 1));
        //alphaR->Scale(0.75);
        //TH2F* alphaR_summed = (TH2F*)alphaR->Clone("alpha_vs_R_summed");
        //alphaR_summed->Add((TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 2)), 0.25);
        TH2F* alphaR_summed = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 2));
        //alphaR_summed->Scale(0.25);


        // Plotting for debugging etc.
        gStyle->SetOptStat(1);
        TCanvas* c1 = new TCanvas("c2", "c2", 800, 600);
        c1->cd();
        alphaR_summed->SetTitle(Form("#alpha vs R, %s%% cent, %s", centleg[icent], energies[ienergy]));
        alphaR_summed->Draw("COLZ");
        c1->SaveAs(Form("figs/alpha_vs_R_summed_%s_cent%s.png", energies[ienergy], centleg[icent]));
        c1->Close();
        gStyle->SetOptStat(0);
        // Average and stddev from the summed buffer histo
        if(iparam==0) {
            //graph->SetPoint(ienergy, energydouble[ienergy], alphaR_summed->GetMean(1));
            //graph->SetPointError(ienergy, 0, 0, alphaR_summed->GetStdDev(1), alphaR_summed->GetStdDev(1));
            graph->SetPoint(ienergy, energydouble[ienergy], 0.75*alphaR->GetMean(1) + 0.25*alphaR_summed->GetMean(1)); // only mT=300
            graph->SetPointError(ienergy, 0, 0, 0.25*alphaR->GetStdDev(1) + 0.75*alphaR_summed->GetStdDev(1), 0.25*alphaR->GetStdDev(1) + 0.75*alphaR_summed->GetStdDev(1)); // only mT=300
        } else if(iparam==1) {
            //graph->SetPoint(ienergy, energydouble[ienergy], alphaR_summed->GetMean(2));
            //graph->SetPointError(ienergy, 0, 0, alphaR_summed->GetStdDev(2), alphaR_summed->GetStdDev(2));
            graph->SetPoint(ienergy, energydouble[ienergy], 0.75*alphaR->GetMean(2) + 0.25*alphaR_summed->GetMean(2)); // only mT=300
            graph->SetPointError(ienergy, 0, 0, 0.25*alphaR->GetStdDev(2) + 0.75*alphaR_summed->GetStdDev(2), 0.25*alphaR->GetStdDev(2) + 0.75*alphaR_summed->GetStdDev(2)); // only mT=300
        }
        f->Close();
    }
    return 0;
}

int calculate_syst_avg_urqmd(int iparam, int icent, TGraphAsymmErrors* graph) {
    if(icent != 11) {
        std::cerr << "Error: Centrality class " << icent << " is not available for UrQMD." << std::endl;
        return -1;
    }
    for(int ienergy=0; ienergy<NENERGIES_urqmd; ienergy++) {
        TFile* f = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent%s_%s_AVG%d.root", centleg[icent], energies_urqmd[ienergy],NEVT_AVG), "READ");
        /*
        // Create a buffer histogram
        TH2F* alphaR = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 0));
        TH2F* alphaR_summed = (TH2F*)alphaR->Clone("alpha_vs_R_summed");
        for(int ikt=1; ikt<NKT; ikt++) {
            alphaR_summed->Add((TH2F*)f->Get(Form("alpha_vs_R_ikt%i", ikt)));
        }
        */
        // For only m_T = 300 GeV !!!
        TH2F* alphaR = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 1));
        //alphaR->Scale(0.75);
        //TH2F* alphaR_summed = (TH2F*)alphaR->Clone("alpha_vs_R_summed");
        //alphaR_summed->Add((TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 2)), 0.25);
        TH2F* alphaR_summed = (TH2F*)f->Get(Form("alpha_vs_R_ikt%i", 2));
        //alphaR_summed->Scale(0.25);


        // Plotting for debugging etc.
        gStyle->SetOptStat(1);
        TCanvas* c1 = new TCanvas("c2", "c2", 800, 600);
        c1->cd();
        alphaR_summed->SetTitle(Form("#alpha vs R, %s%% cent, %s", centleg[icent], energies_urqmd[ienergy]));
        alphaR_summed->Draw("COLZ");
        c1->SaveAs(Form("figs/urqmd_alpha_vs_R_summed_%s_cent%s.png", energies_urqmd[ienergy], centleg[icent]));
        c1->Close();
        gStyle->SetOptStat(0);
        // Average and stddev from the summed buffer histo
        if(iparam==0) {
            //graph->SetPoint(ienergy, energydouble_urqmd[ienergy], alphaR_summed->GetMean(1));
            //graph->SetPointError(ienergy, 0, 0, alphaR_summed->GetStdDev(1), alphaR_summed->GetStdDev(1));
            graph->SetPoint(ienergy, energydouble_urqmd[ienergy], 0.75*alphaR->GetMean(1) + 0.25*alphaR_summed->GetMean(1)); // only mT=300
            graph->SetPointError(ienergy, 0, 0, 0.25*alphaR->GetStdDev(1) + 0.75*alphaR_summed->GetStdDev(1), 0.25*alphaR->GetStdDev(1) + 0.75*alphaR_summed->GetStdDev(1)); // only mT=300
        } else if(iparam==1) {
            //graph->SetPoint(ienergy, energydouble_urqmd[ienergy], alphaR_summed->GetMean(2));
            //graph->SetPointError(ienergy, 0, 0, alphaR_summed->GetStdDev(2), alphaR_summed->GetStdDev(2));
            graph->SetPoint(ienergy, energydouble_urqmd[ienergy], 0.75*alphaR->GetMean(2) + 0.25*alphaR_summed->GetMean(2)); // only mT=300
            graph->SetPointError(ienergy, 0, 0, 0.25*alphaR->GetStdDev(2) + 0.75*alphaR_summed->GetStdDev(2), 0.25*alphaR->GetStdDev(2) + 0.75*alphaR_summed->GetStdDev(2)); // only mT=300
        }
        f->Close();
    }
    return 0;
}

// ----------------------  MAIN  ----------------------------------------
void plot_alphaNR_allcent() {
    gStyle->SetOptStat(0);
    //gStyle->SetTitleAlign(13);
    //gStyle->SetTitleX(0.2);
    //gStyle->SetTitleY(0.95);

    for (int iparam = 0; iparam < 2; ++iparam) {
        const char* param = (iparam == 0) ? "alpha" : "R";
        TCanvas* c = new TCanvas(Form("c_%s_final", param), Form("Final %s vs sqrt(sNN)", param), 800, 600);
        TLegend* leg = new TLegend(0.6, 0.18, 0.88, 0.38);//(0.6, 0.6, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        bool first = true;

        for (int icent = 11; icent < NCENT+2; ++icent) { // !!! icent=10 : only "all" and "0-10%"
            c->cd();
            // For "statistical error"
            // TODO? rename to another variable, plot with error bars (+small horiz. lines at the end) instead of filled uncertainty area
            /*
            TString fname = Form("alphaNR_vs_kt/kt_averaged_%s_vs_energy_cent%s.root", param, centleg[icent]);
            TFile* f = TFile::Open(fname);
            if (!f || f->IsZombie()) {
                std::cerr << "Could not open file: " << fname << std::endl;
                continue;
            }

            TGraphAsymmErrors* graph = (TGraphAsymmErrors*)f->Get("Graph");
            if (!graph) {
                std::cerr << "Graph not found in file: " << fname << std::endl;
                f->Close();
                continue;
            }
            // If uncommenting, do not forget to uncomment f->Close below
            */
            /**/
            // For "systematic error"
            TGraphAsymmErrors* graph = new TGraphAsymmErrors(NENERGIES);
            int sysavg_status = calculate_syst_avg(iparam, icent, graph);
            if(sysavg_status != 0) {
                std::cerr << "Error calculating systematic average for parameter " << param << " and centrality " << centleg[icent] << std::endl;
                continue;
            }
            /**/

            /*
            // Not needed since 200 GeV simulation is already available with EPOS4
            if(iparam==0 && icent==11)
            {
                int n=graph->GetN();
                graph->SetPoint(n, 200, 1.595);
                graph->SetPointError(n, 0, 0, 0.085, 0.085);
            }
            */


            graph->SetMarkerStyle(20);// + icent); // Different marker style for each centrality
            graph->SetMarkerColor(colors[icent]);
            graph->SetLineColor(colors[icent]);
            graph->SetFillColorAlpha(colors[icent], 0.3);

            // Make it similar to STAR Preliminary alpha(sNN) plots
            graph->GetXaxis()->SetLimits(1,400);//RangeUser 9999
            //graph->GetYaxis()->SetLimits(0.1,2.1);
            c->SetLogx(1);
            graph->SetTitle("EPOS4 & UrQMD 0-10% #pi^{#pm}#pi^{#pm}");
            if(iparam==0) graph->GetYaxis()->SetRangeUser(0.1, 2.1);//(1.3,1.8);
            else if(iparam==1) graph->GetYaxis()->SetRangeUser(1.5, 14.5);//(4.5, 6.5);


            if (first) {
                //graph->Draw("AP");
                graph->Draw("A3");
                graph->Draw("PX SAME");
                graph->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
                graph->GetYaxis()->SetTitle((std::string(param) == "alpha") ? "#alpha" : "R [fm]");
                graph->GetYaxis()->SetTitleOffset(1.2);
                first = false;
            } else {
                //graph->Draw("P SAME");
                graph->Draw("3 SAME");
                graph->Draw("PX SAME");
            }

            leg->AddEntry(graph, Form("%s%% cent, EPOS4", centleg[icent]), "lp");
            
            /*
            f->Close();
            //delete graph;
            */


            // UrQMD
            // -----
            TGraphAsymmErrors* graph_urqmd = new TGraphAsymmErrors(NENERGIES_urqmd);
            sysavg_status = calculate_syst_avg_urqmd(iparam, icent, graph_urqmd);
            if(sysavg_status != 0) {
                std::cerr << "Error calculating systematic average for parameter " << param << " and centrality " << centleg[icent] << std::endl;
                continue;
            }
            
            graph_urqmd->SetMarkerStyle(21);
            graph_urqmd->SetMarkerColor(kTeal);
            graph_urqmd->SetLineColor(kTeal);
            graph_urqmd->SetFillColorAlpha(kTeal, 0.3);
            graph_urqmd->SetTitle("EPOS4 & UrQMD 0-10% #pi^{#pm}#pi^{#pm}");
            graph_urqmd->GetXaxis()->SetLimits(1,400);//RangeUser 9999
            if(iparam==0) graph_urqmd->GetYaxis()->SetRangeUser(0.1,2.1);//(1.3,1.8);
            else if(iparam==1) graph_urqmd->GetYaxis()->SetRangeUser(1.5, 14.5);//(4.5, 6.5)


            if (first) {
                graph_urqmd->Draw("3 SAME");
                graph_urqmd->Draw("PX SAME");
                //graph_urqmd->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
                //graph_urqmd->GetYaxis()->SetTitle((std::string(param) == "alpha") ? "#alpha" : "R [fm]");
                //graph_urqmd->GetYaxis()->SetTitleOffset(1.2);
                //first = false;
            } else {
                //graph->Draw("P SAME");
                graph_urqmd->Draw("3 SAME");
                graph_urqmd->Draw("PX SAME");
            }

            leg->AddEntry(graph_urqmd, Form("%s%% cent, UrQMD", centleg[icent]), "lp");
        }

        leg->Draw();

        /*
        TLatex title;
        title.SetNDC();
        title.SetTextSize(0.04);
        title.DrawLatex(0.15, 0.93, Form("Averaged %s vs #sqrt{s_{NN}}", param));
        */

        c->SaveAs(Form("figs/final_%s_vs_energy.pdf", param));
        //c->SaveAs(Form("figs/final_%s_vs_energy.pdf", param));
        delete c;
    }
}
