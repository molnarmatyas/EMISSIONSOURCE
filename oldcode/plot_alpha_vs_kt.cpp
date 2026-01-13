// ./extract_params.sh 9.2GeV > extracted_params.h &&  root.exe -b -q 'plot_alpha_vs_kt.cpp("9.2GeV")'

#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include "extracted_params.h"

// ---------  MAIN  -----------------------------------------
void plot_alpha_vs_kt(const char* energy="9.2GeV") {
    const int nPoints = 10;//5;
    double ktBinEdges[nPoints+1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};//{0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    double ktBins[nPoints];
    for(int ii=0; ii<nPoints; ii++)
    {
        // K_T bin centers
        ktBins[ii] = (ktBinEdges[ii+1]+ktBinEdges[ii])/2;
    }
    double binWidths[nPoints];
    for(int ii=0; ii<nPoints; ii++)
    {
        // K_T bin centers
        binWidths[ii] = ktBinEdges[ii+1]-ktBinEdges[ii];
    }

    // Compute the asymmetric errors for TGraphAsymmErrors
    double xLow[nPoints], xHigh[nPoints];
    for (int i = 0; i < nPoints; ++i) {
        xLow[i] = binWidths[i]/2;
        xHigh[i] = binWidths[i]/2;
    }

    // Create the TGraphAsymmErrors object
    TGraphAsymmErrors *graphalpha = new TGraphAsymmErrors(nPoints, ktBins, alphaValues, xLow, xHigh, alphaErrLow, alphaErrHigh);
    TGraphAsymmErrors *graphN     = new TGraphAsymmErrors(nPoints, ktBins, NValues, xLow, xHigh, NErrLow, NErrHigh);
    TGraphAsymmErrors *graphR     = new TGraphAsymmErrors(nPoints, ktBins, RValues, xLow, xHigh, RErrLow, RErrHigh);

    // --- alpha ---
    // Style settings
    graphalpha->SetTitle(Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#alpha",energy));
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
    c->SaveAs(Form("figs/alpha_vs_kt_%s.png",energy));

    // --- R ---

    // Style settings
    graphR->SetTitle(Form("#R(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R",energy));
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
    //graphR->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y-axis range
    graphR->Draw("AP");

    // Save the plot to a file
    c->SaveAs(Form("figs/R_vs_kt_%s.png",energy));

    // --- N ---

    // Style settings
    graphN->SetTitle(Form("#N(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#N",energy));
    graphN->SetMarkerStyle(20);
    graphN->SetMarkerSize(1.2);
    graphN->SetLineWidth(2);

    // Draw the graph
    //TCanvas *c = new TCanvas("c", Form("#N(K_{T}), #sqrt{s_{NN}}=%s",energy), 800, 600);
    c->SetGrid();
    gStyle->SetOptStat(0);
    //graph->GetYaxis()->SetNangeUser(0,2);
    // Customize the axes
    graphN->GetXaxis()->SetLimits(0.0, 0.7); // Set X-axis range
    //graphN->GetYaxis()->SetRangeUser(0.0, 2.0); // Set Y-axis range
    graphN->Draw("AP");

    // Save the plot to a file
    c->SaveAs(Form("figs/N_vs_kt_%s.png",energy));

    // Save graphs to a .root file
    TFile* outFile = new TFile(Form("alphaNR_graphs_%s.root",energy), "RECREATE");
    graphalpha->Write("graphalpha");
    graphN->Write("graphN");
    graphR->Write("graphR");
    outFile->Close();

    std::cout << Form("Graphs saved to alphaNR_graphs_%s.root",energy) << std::endl;

    // Clean up
    delete c;
    delete graphalpha;
    delete graphN;
    delete graphR;
}