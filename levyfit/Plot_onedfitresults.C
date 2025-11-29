#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <cmath>
#include <iostream>

#define SQR(x) ((x) * (x))

const double Mass2_pi = 0.019479835;
const int NKT = 10;

const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
const char* path = "~/work/figs/EPOS_project/3d_proj/";
const char* frames[3] = {"lcms","pcms","lab"};
const int thisframe = 0;

void Plot_onedfitresults() 
{
    // Open the fitresults.root file
    TFile *file = TFile::Open(Form("onedfitresults_%s.root",frames[thisframe]));
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file fitresults.root" << std::endl;
        return;
    }

    // Create arrays to store the histograms and TGraphErrors
    TH1D* alphahist[NKT];
    TH1D* Rhist[NKT];

    TGraphErrors *alphaGraph = new TGraphErrors(NKT);
    TGraphErrors *RGraph = new TGraphErrors(NKT);

    // Loop over the kT indices
    for (int ikt = 0; ikt < NKT; ++ikt) {
        // Form the histogram names and retrieve them from the file
        alphahist[ikt] = (TH1D*)file->Get(Form("alphahist_ikt%d", ikt));
        Rhist[ikt] = (TH1D*)file->Get(Form("Rhist_ikt%d", ikt));

        if (!alphahist[ikt] || !Rhist[ikt] ) {
            std::cerr << "Error retrieving histograms for kT index " << ikt << std::endl;
            continue;
        }

        // Calculate the x value for the graph
        double kt_mid = 0.5 * (ktbins[ikt] + ktbins[ikt + 1]);
        double x_value = std::sqrt(Mass2_pi + SQR(kt_mid));

        // Extract mean and stddev values and set points in TGraphErrors
        double mean_alpha = alphahist[ikt]->GetMean();
        double stddev_alpha = alphahist[ikt]->GetStdDev();
        alphaGraph->SetPoint(ikt, x_value, mean_alpha);
        alphaGraph->SetPointError(ikt, 0, stddev_alpha);

        double mean_R = Rhist[ikt]->GetMean();
        double stddev_R = Rhist[ikt]->GetStdDev();
        RGraph->SetPoint(ikt, x_value, mean_R);
        RGraph->SetPointError(ikt, 0, stddev_R);
    }

    // Create a canvas for the alpha kT dependency plot
    TCanvas *c3 = new TCanvas("alpha_ktdep", "Alpha kT Dependency", 800, 600);
    alphaGraph->SetTitle("Alpha kT Dependency; kT Value; Mean Value");
    alphaGraph->SetMarkerStyle(20);
    alphaGraph->SetMinimum(1.);
    alphaGraph->SetMaximum(2.);
    alphaGraph->Draw("AP");
    c3->SaveAs(Form("%s%s/alpha1d_ktdep.pdf",path,frames[thisframe]));

    // Create a canvas for the R kT dependency plot
    TCanvas *c4 = new TCanvas("R_ktdep", "R kT Dependency", 800, 600);
    RGraph->SetLineColor(kRed);
    RGraph->SetMarkerColor(kRed);

    RGraph->SetTitle("R kT Dependency; kT Value; Mean Value");
    RGraph->SetMarkerStyle(20);
    RGraph->SetMinimum(3.8);
    RGraph->SetMaximum(12.2);
    RGraph->Draw("AP");

    c4->SaveAs(Form("%s%s/R1d_ktdep.pdf",path,frames[thisframe]));

    // Close the file
    file->Close();
    return;
}
