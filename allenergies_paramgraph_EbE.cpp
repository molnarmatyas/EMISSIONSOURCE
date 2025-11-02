// root.exe -b -q allenergies_paramgraph_EbE.cpp

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <vector>
#include <iostream>

const int NCENT = 10; // number of centrality classes
const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};

// Function to overlay graphs for a specific parameter
void overlay_for_parameter(const std::vector<TString>& energies, const std::vector<int>& colors, const char* param, int icent) {
    // Create a canvas
    TCanvas* c = new TCanvas(Form("c_%s", param), Form("%s vs K_{T} for all energies, %s%%", param, centleg[icent]), 800, 600);
    c->SetGrid();

    // Create a legend
    //TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position: top-right corner
    TLegend* legend = new TLegend(0.15,0.2); // Position: auto
    legend->SetTextSize(0.03);
    legend->SetBorderSize(1);

    // Load and draw graphs
    bool firstGraph = true;
    for (size_t i = 0; i < energies.size(); ++i) {
        TString fileName = Form("alphaNR_vs_kt/alphaNR_graphs_%s_cent%s.root", energies[i].Data(), centleg[icent]);
        TFile* file = TFile::Open(fileName, "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << fileName << std::endl;
            continue;
        }

        TString graphName = Form("graph%s%s", param, centleg[icent]);
        TGraphAsymmErrors* graph = (TGraphAsymmErrors*)file->Get(graphName);
        if (!graph) {
            std::cerr << "Error: Could not find " << graphName << " in " << fileName << std::endl;
            file->Close();
            continue;
        }

        // Style the graph
        graph->SetLineColor(colors[i]);
        graph->SetMarkerColor(colors[i]);
        graph->SetMarkerStyle(20 + i); // Different marker style for each energy

        // Add to legend
        legend->AddEntry(graph, energies[i], "lp");

        // Plotting range
        double yMin = 0.0, yMax = 2.0; // Default range
        if (strcmp(param, "N") == 0) {
            yMin = 0.2;
            yMax = 1.5;
        } else if (strcmp(param, "R") == 0) {
            yMin = 3.0;
            yMax = 12.0;
        } // Keep default for alpha (0.0 to 2.0)

        // Draw the graph
        if (firstGraph) {
            graph->SetTitle(Form("%s vs K_{T} for all energies, %s%%;K_{T} (GeV/c);%s", param, centleg[icent], param));
            graph->GetXaxis()->SetLimits(0.0, 0.7); // Set X-axis range
            graph->GetYaxis()->SetRangeUser(yMin, yMax); // Adjust Y-axis range as needed
            graph->Draw("AP");
            firstGraph = false;
        } else {
            graph->Draw("P SAME");
        }

        file->Close();
    }

    // Draw the legend
    legend->Draw();

    // Save the plot
    c->SaveAs(Form("figs/overlay_%s_vs_kt_cent%s.png", param, centleg[icent]));
    c->SaveAs(Form("alphaNR_vs_kt/overlay_%s_vs_kt_cent%s.root", param, centleg[icent]));

    std::cout << "Overlay plot for " << param << " saved as overlay_" << param << "_vs_kt.png and .root" << std::endl;

    // Clean up
    delete legend;
    delete c;
}

// Main function to plot all parameters
void allenergies_paramgraph_EbE() {
    for(int icent=0; icent<NCENT+2; icent++)
    {
        // List of energies
        std::vector<TString> energies = {"7.7GeV", "9.2GeV", "11.5GeV", "14.5GeV", "19.6GeV", "27GeV", "39GeV", "62.4GeV", "130GeV", "200GeV"};
        std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kGray, kCyan, kOrange, kBlack, kViolet, kPink, kTeal, kAzure, kSpring}; // Colors for each graph

        // Parameters to process
        std::vector<TString> parameters = {"#alpha", "R"}; //"N", 

        double ymin = 0.;
        double ymax = 2.0;
        for (const auto& param : parameters) {
            overlay_for_parameter(energies, colors, param.Data(), icent);
        }
    }
}

