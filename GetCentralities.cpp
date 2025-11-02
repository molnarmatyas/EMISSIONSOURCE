#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>

int GetCentralities()
{
    // Open ROOT file and retrieve histogram
    TFile *file = TFile::Open("EPOS_3d_source_010cent_all.root");
    TH1D *mult = (TH1D*)file->Get("epos_pion_multiplicity_distribution");

    if (!mult) {
        cout << "Error: Histogram not found!" << endl;
        return -1;
    }

    // Define centrality binning
    const int nQuantiles = 10;
    double quantileValues[nQuantiles];  
    double probabilities[nQuantiles] = {0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.0};

    /*
    for (int i = 0; i < nQuantiles; i++) {
        probabilities[i] = 1.0 - (i + 1) * 0.1;  // 90%, 80%, ..., 10%, 0%
    }
    */

    // Compute centrality thresholds
    mult->GetQuantiles(nQuantiles, quantileValues, probabilities);

    // Print threshold values
    cout << "Centrality Thresholds:" << endl;
    for (int i = 0; i < nQuantiles; i++) {
        cout << "Centrality " << probabilities[nQuantiles-1-i] *100. << "-" << probabilities[nQuantiles-2-i] *100. 
             << "% threshold: " << quantileValues[i] << endl;
    }

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Multiplicity Distribution", 800, 600);
    gPad->SetLogy(); // Set logarithmic y-axis

    // Draw histogram
    mult->SetLineColor(kBlue);
    mult->SetLineWidth(2);
    mult->Draw();

    // Overlay vertical lines and labels
    for (int i = 0; i < nQuantiles; i++) {
        // Draw vertical line
        TLine *line = new TLine(quantileValues[i], 1, quantileValues[i], mult->GetMaximum());
        line->SetLineColor(kRed);
        line->SetLineStyle(2);  // Dashed line
        line->SetLineWidth(2);
        line->Draw();

        // Add text label for centrality class
        TLatex *label = new TLatex(quantileValues[i] * 1.05, mult->GetMaximum() / 3, 
                                   Form("%.0f-%.0f%%", probabilities[nQuantiles-1-i] *100., probabilities[nQuantiles-2-i] *100.));
        label->SetTextSize(0.03);
        label->SetTextColor(kRed);
        label->SetTextAngle(90);  // Vertical orientation
        label->Draw();
    }

    // Save the canvas
    c1->SaveAs("figs/Multiplicity_Centrality.png");

    return 0;
}

