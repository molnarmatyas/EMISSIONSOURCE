// root.exe -b -q averaged_paramgraph.cpp

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

// Function to calculate K_T-averaged value and errors
void calculate_kt_averaged(const TGraphAsymmErrors* graph, double& avg, double& errLow, double& errHigh) {
    double sum = 0.0, weightSum = 0.0;
    double errSumLow = 0.0, errSumHigh = 0.0;

    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y, errYLow, errYHigh;
        graph->GetPoint(i, x, y);
        errYLow = graph->GetErrorYlow(i);
        errYHigh = graph->GetErrorYhigh(i);

        double weight = 1.0 / ((errYLow + errYHigh) / 2.0); // Use average error as weight
        sum += y * weight;
        weightSum += weight;

        errSumLow += pow(errYLow, 2); // Propagate errors
        errSumHigh += pow(errYHigh, 2);
    }

    avg = sum / weightSum;
    errLow = sqrt(errSumLow) / graph->GetN();  // Divide by sqrt(N) for average
    errHigh = sqrt(errSumHigh) / graph->GetN();
}

void plot_kt_averaged_vs_energy(const std::vector<TString>& energies, const std::vector<int>& colors, const char* param) {
    std::vector<double> sqrtSNN, avgValues, errLowValues, errHighValues;

    for (size_t i = 0; i < energies.size(); ++i) {
        TString fileName = Form("alphaNR_graphs_%sGeV.root", energies[i].Data());
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

        // Calculate K_T-averaged value and errors
        double avg, errLow, errHigh;
        calculate_kt_averaged(graph, avg, errLow, errHigh);

        // Add results to vectors
        sqrtSNN.push_back(energies[i].Atof()); // Convert energy string to numeric
        avgValues.push_back(avg);
        errLowValues.push_back(errLow);
        errHighValues.push_back(errHigh);

        file->Close();
    }

    const int nPoints = sqrtSNN.size();
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(nPoints);

    for (int i = 0; i < nPoints; ++i) {
        graph->SetPoint(i, sqrtSNN[i], avgValues[i]);
        graph->SetPointError(i, 0.0, 0.0, errLowValues[i], errHighValues[i]); // No X-axis errors
    }

    // Create canvas and draw graph
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetTitleSize(0.04, "XYZ");
    TCanvas* c = new TCanvas(Form("c_avg_%s", param), Form("K_{T}-Averaged %s vs Energy", param), 1080, 720);
    c->SetGrid();
    c->SetLogx();
    

    graph->SetTitle(Form("K_{T}-averaged #%s vs energy;#sqrt{s_{NN}} (GeV);<#%s>_{K_{T}}", param, param));
    graph->SetMarkerStyle(21); // Larger solid circles
    graph->SetMarkerSize(1.8);
    graph->SetLineWidth(2);
    graph->SetLineColor(kBlack);  // Solid black line for better visibility
    graph->SetMarkerColor(kBlue); // Blue solid markers
    graph->SetFillColor(kBlue);   // Fill the marker with solid color
    
    // Set the fill color with alpha
    int fillColor = TColor::GetColorTransparent(kRed, 0.5);
    graph->SetFillColor(fillColor);

    // Draw the graph with filled error bars
    graph->Draw("A3");

    // Draw the graph points and lines
    graph->Draw("P");
    //graph->Draw("AP");

    // Save the plot
    c->SaveAs(Form("figs/kt_averaged_%s_vs_energy.png", param));
    c->SaveAs(Form("kt_averaged_%s_vs_energy.root", param));

    std::cout << "K_{T}-averaged plot for " << param << " saved as kt_averaged_" << param << "_vs_energy.png and .root" << std::endl;

    delete c;
}

void plot_all_kt_averaged(const std::vector<TString>& energies, const std::vector<int>& colors) {
    plot_kt_averaged_vs_energy(energies, colors, "alpha");
    plot_kt_averaged_vs_energy(energies, colors, "N");
    plot_kt_averaged_vs_energy(energies, colors, "R");
}

void averaged_paramgraph() {
    std::vector<TString> energies = {"7.7", "9.2", "19.6", "200"}; // "27"
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kOrange};
    plot_all_kt_averaged(energies, colors);
}