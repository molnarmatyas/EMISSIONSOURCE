// root.exe -b -q averaged_paramgraph_EbE.cpp

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

const int NCENT = 10; // number of centrality classes
const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};

TGraphAsymmErrors* RemovePoint(const TGraphAsymmErrors* original, int removeIndex) {
    int n = original->GetN();
    TGraphAsymmErrors* newGraph = new TGraphAsymmErrors(n - 1);

    int j = 0;
    for (int i = 0; i < n; ++i) {
        if (i == removeIndex) continue;

        double x, y;
        original->GetPoint(i, x, y);
        newGraph->SetPoint(j, x, y);
        newGraph->SetPointError(j,
            original->GetErrorXlow(i),
            original->GetErrorXhigh(i),
            original->GetErrorYlow(i),
            original->GetErrorYhigh(i));
        ++j;
    }

    return newGraph;
}


// Function to calculate K_T-averaged value and errors
void calculate_kt_averaged(const TGraphAsymmErrors* graph, double& avg, double& errLow, double& errHigh) {
    double sum = 0.0, weightSum = 0.0;
    double errSumLow = 0.0, errSumHigh = 0.0;

    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y, errYLow, errYHigh;
        graph->GetPoint(i, x, y);
        if (y <= 0 || std::isnan(y) || std::isinf(y)) //continue;
        {
            graph = RemovePoint(graph, i);
            continue;
        }
        errYLow = graph->GetErrorYlow(i);
        errYHigh = graph->GetErrorYhigh(i);

        double weight = 1.0 / ((errYLow + errYHigh) / 2.0); // Use average error as weight
        sum += y * weight;
        weightSum += weight;

        errSumLow += pow(errYLow, 2); // Propagate errors
        errSumHigh += pow(errYHigh, 2);
    }

    avg = sum / weightSum;
    cerr << "avg: " << avg << endl;
    errLow = sqrt(errSumLow) / graph->GetN();  // Divide by sqrt(N) for average
    errHigh = sqrt(errSumHigh) / graph->GetN();
}

void plot_kt_averaged_vs_energy(const std::vector<TString>& energies, const std::vector<int>& colors, const char* param, int icent) {
    std::vector<double> sqrtSNN, avgValues, errLowValues, errHighValues;

    for (size_t i = 0; i < energies.size(); ++i) {
        TString fileName = Form("alphaNR_vs_kt/alphaNR_graphs_%sGeV_cent%s.root", energies[i].Data(),centleg[icent]);
        //cerr << energies[i].Data() << endl;
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
    

    graph->SetTitle(Form("K_{T}-averaged %s vs energy, %s%%;#sqrt{s_{NN}} (GeV);<%s>_{K_{T}}", param, centleg[icent], param));
    graph->SetMarkerStyle(21); // Larger solid circles
    graph->SetMarkerSize(1.8);
    graph->SetLineWidth(2);
    graph->SetLineColor(kBlack);  // Solid black line for better visibility
    graph->SetMarkerColor(kBlue); // Blue solid markers
    graph->SetFillColor(kBlue);   // Fill the marker with solid color
    
    // Set the fill color with alpha
    int fillColor = TColor::GetColorTransparent(kRed, 0.5);
    graph->SetFillColor(fillColor);

    //if(strcmp(param, "alpha") == 0) graph->GetYaxis()->SetRangeUser(1.3, 1.7);
    //else if(strcmp(param, "R") == 0) graph->GetYaxis()->SetRangeUser(4.5, 6.5);

    // Draw the graph with filled error bars
    //graph->Draw("A3");

    // Draw the graph points and lines
    //graph->Draw("P");
    graph->Draw("AP");

    // Save the plot
    c->SaveAs(Form("figs/kt_averaged_%s_vs_energy_cent%s.png", param, centleg[icent])); 
    TFile* outfile = new TFile(Form("alphaNR_vs_kt/kt_averaged_%s_vs_energy_cent%s.root", param, centleg[icent]), "RECREATE");
    outfile->cd();
    graph->Write();
    outfile->Close();

    std::cout << "K_{T}-averaged plot for " << param << " saved as kt_averaged_" << param << "_vs_energy.png and .root" << std::endl;

    delete c;
}

void plot_all_kt_averaged(const std::vector<TString>& energies, const std::vector<int>& colors) {
    for(int icent=0; icent<NCENT+2; icent++)
    {
        plot_kt_averaged_vs_energy(energies, colors, "alpha", icent);
        //plot_kt_averaged_vs_energy(energies, colors, "N");
        plot_kt_averaged_vs_energy(energies, colors, "R", icent);
    }
}

void averaged_paramgraph_EbE() {
    std::vector<TString> energies = {"7.7", "9.2", "11.5", "14.5", "19.6", "27", "39", "62.4", "130", "200"};
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kGray, kCyan, kOrange, kBlack, kViolet, kPink, kTeal, kAzure, kSpring}; // Colors for each graph
    plot_all_kt_averaged(energies, colors);
}
