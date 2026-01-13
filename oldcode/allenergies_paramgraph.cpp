// root.exe -b -q allenergies_paramgraph.cpp

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <vector>
#include <iostream>

// Function to overlay graphs for a specific parameter
void overlay_for_parameter(const std::vector<TString>& energies, const std::vector<int>& colors, const char* param) {
    // Create a canvas
    TCanvas* c = new TCanvas(Form("c_%s", param), Form("#%s vs K_{T} for all energies", param), 800, 600);
    c->SetGrid();

    // Create a legend
    //TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position: top-right corner
    TLegend* legend = new TLegend(0.15,0.2); // Position: auto
    legend->SetTextSize(0.03);
    legend->SetBorderSize(1);

    // Load and draw graphs
    bool firstGraph = true;
    for (size_t i = 0; i < energies.size(); ++i) {
        TString fileName = Form("alphaNR_graphs_%s.root", energies[i].Data());
        TFile* file = TFile::Open(fileName, "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << fileName << std::endl;
            continue;
        }

        TString graphName = Form("graph%s", param);
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
            graph->SetTitle(Form("#%s vs K_{T} for all energies;K_{T} (GeV/c);#%s", param, param));
            graph->GetXaxis()->SetLimits(0.0, 0.6); // Set X-axis range
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
    c->SaveAs(Form("figs/overlay_%s_vs_kt.png", param));
    c->SaveAs(Form("overlay_%s_vs_kt.root", param));

    std::cout << "Overlay plot for " << param << " saved as overlay_" << param << "_vs_kt.png and .root" << std::endl;

    // Clean up
    delete legend;
    delete c;
}

// Main function to plot all parameters
void allenergies_paramgraph() {
    // List of energies
    std::vector<TString> energies = {"7.7GeV", "9.2GeV", "19.6GeV", "200GeV"}; // 27 GeV to be applied
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan}; // Colors for each graph

    // Parameters to process
    std::vector<TString> parameters = {"alpha", "N", "R"};

    double ymin = 0.;
    double ymax = 2.0;
    for (const auto& param : parameters) {
        overlay_for_parameter(energies, colors, param.Data());
    }
}

