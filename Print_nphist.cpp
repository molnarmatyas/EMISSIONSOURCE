// root.exe -b -q  Print_nphist.cpp\(\"9.2GeV\"\)

int Print_nphist(const char* energy = "7.7GeV")
{
    TFile *f = new TFile(Form("analysed/EPOS_3d_source_allcent_all_%s.root",energy), "READ");
    gStyle->SetOptStat(1111); // Enable statistics box (1111 = show all stats)
    gStyle->SetStatX(0.4);    // (default is 0.78)
    gStyle->SetStatY(0.9);    // (default is 0.85)
    TCanvas *c1 = new TCanvas();
    TH1F* np = (TH1F*)f->Get("epos_pion_multiplicity_distribution");

    np->SetTitle(Form("EPOS4 pion multiplicity, #sqrt{s_{NN}} = %s", energy));
    np->GetXaxis()->SetRangeUser(0,2500);
    
    c1->SetLogy(1);
    np->Draw();
    c1->Print(Form("figs/np_%s.png",energy));
    
    return 0;
}