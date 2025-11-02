// root.exe -b -q Print_bim_allenergies.cpp

const int colors[] = {kRed, kGray, kGreen, kCyan, kOrange, kBlack, kViolet, kPink, kBlue, kAzure, kSpring, kMagenta}; // kYellow

int Print_bim_allenergies()
{
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  double toright = 0.0;//0.48;
  TLegend* leg = new TLegend(0.2+toright, 0.88, 0.4+toright, 0.6);

  const char* energies[] = {"7.7GeV", "9.2GeV", "11.5GeV", "14.5GeV", "19.6GeV", "27GeV", "62.4GeV", "130GeV", "200GeV"};
  const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
  const double energydouble[NENERGIES] = {7.7, 9.2, 11.5, 14.5, 19.6, 27., 62.4, 130., 200.};

  // Define centrality binning
  const int nQuantiles = 10;
  double quantileValues[nQuantiles];  
  double probabilities[nQuantiles] = {0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05,0.0};
  //double probabilities[nQuantiles] = {0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
  
  for(int ienergy=0; ienergy<NENERGIES; ienergy++) //+NENERGIES-1
  {
    TFile* f = new TFile(Form("analysed/EPOS_3d_source_allcent_all_%s.root", energies[ienergy]),"READ");
    TH1F* bimhist = (TH1F*)f->Get("epos_impact_parameter_distribution");
    bimhist->SetTitle("EPOS4 impact parameter distribution, lower cent. thresholds");
    bimhist->Rebin(50);
    bimhist->SetMarkerStyle(20+ienergy);
    bimhist->SetMarkerSize(0.8);
    bimhist->SetMarkerColor(colors[ienergy]);
    bimhist->SetLineColor(colors[ienergy]);
    bimhist->GetYaxis()->SetRangeUser(0, 130);
    //bimhist->GetXaxis()->SetRangeUser(12,20);
    bimhist->GetXaxis()->SetTitle("b [fm]");
    bimhist->Draw("PL SAME");
    leg->AddEntry(bimhist, Form("#sqrt{s_{NN}}=%s",energies[ienergy]));

    // Compute centrality thresholds
    bimhist->GetQuantiles(nQuantiles, quantileValues, probabilities);
    // Print threshold values
    std::cout << " ** Centrality Thresholds for energy " << energies[ienergy] << endl;
    for (int i = 0; i < nQuantiles; i++) {
        std::cout << "Centrality " << probabilities[nQuantiles-1-i] *100. << "-" << probabilities[nQuantiles-2-i] *100. 
             << "% lower threshold: " << quantileValues[nQuantiles-1-i] << endl; // !!! nQuantiles-1-i for descending order
    }
    // Overlay vertical lines and labels
    for (int i = 0; i < nQuantiles; i++) {
        // Draw vertical line
        TLine *line = new TLine(quantileValues[i], 1, quantileValues[i], bimhist->GetMaximum());
        line->SetLineColor(kRed);
        line->SetLineStyle(2);  // Dashed line
        line->SetLineWidth(2);
        line->Draw();

        // Add text label for centrality class
        TLatex *label = new TLatex(quantileValues[i] * 1.05, bimhist->GetMaximum() / 3, 
                                   Form("%.0f-%.0f%%", probabilities[i-1] *100., probabilities[i] *100.));
        label->SetTextSize(0.03);
        label->SetTextColor(kRed);
        label->SetTextAngle(90);  // Vertical orientation
        label->Draw();
    }
  }

  leg->Draw();
  c1->Print("figs/centrality/bimdist_allenergies.png");

  // --------------------------------------
  // Now, do a merged one with all energies
  c1->Clear();
  TFile* f0 = new TFile(Form("analysed/EPOS_3d_source_allcent_all_%s.root", energies[0]),"READ");
  TH1F* bimhist_merged = (TH1F*)f0->Get("epos_impact_parameter_distribution");
  bimhist_merged->SetTitle("EPOS4 impact parameter distribution, lower cent. thresholds");
  bimhist_merged->SetMarkerStyle(20+0);
  bimhist_merged->SetMarkerSize(0.8);
  bimhist_merged->SetMarkerColor(colors[0]);
  bimhist_merged->SetLineColor(colors[0]);
  bimhist_merged->GetXaxis()->SetTitle("b [fm]");

  // Loop over all energies and merge the histograms
  for(int ienergy=1; ienergy<NENERGIES; ienergy++)
  {
    TFile* f = new TFile(Form("analysed/EPOS_3d_source_allcent_all_%s.root", energies[ienergy]),"READ");
    TH1F* bimhist = (TH1F*)f->Get("epos_impact_parameter_distribution");
    bimhist_merged->Add(bimhist);
  }

  // Set properties for the merged histogram
  bimhist_merged->Rebin(50);
  bimhist_merged->Draw("PL");

  // Compute centrality thresholds
  bimhist_merged->GetQuantiles(nQuantiles, quantileValues, probabilities);
  std::cout << " --- Centrality Thresholds from merged bim histogram:" << endl;
  for (int i = 0; i < nQuantiles; i++) {
      std::cout << "Centrality " << probabilities[nQuantiles-1-i] *100. << "-" << probabilities[nQuantiles-2-i] *100. 
            << "% lower threshold: " << quantileValues[nQuantiles-1-i] << endl; // !!! nQuantiles-1-i for descending order
  }
  // Overlay vertical lines and labels
  for (int i = 0; i < nQuantiles; i++) {
      // Draw vertical line
      TLine *line = new TLine(quantileValues[i], 1, quantileValues[i], bimhist_merged->GetMaximum());
      line->SetLineColor(kRed);
      line->SetLineStyle(2);  // Dashed line
      line->SetLineWidth(2);
      line->Draw();

      // Add text label for centrality class
      TLatex *label = new TLatex(quantileValues[i] * 1.05, bimhist_merged->GetMaximum() / 3, 
                                  Form("%.0f-%.0f%%", probabilities[i-1] *100., probabilities[i] *100.));
      label->SetTextSize(0.03);
      label->SetTextColor(kRed);
      label->SetTextAngle(90);  // Vertical orientation
      label->Draw();
  }

  c1->Print("figs/centrality/bimdist_allenergies_merged.png");


  return 0;
}
