#ifndef PI
#define PI 3.14159265358979323846
#endif
int Real_bim()
{
  TFile *f = new TFile("7.7GeV_EPOS_3d_source_010cent_all.root", "READ");
  TH1D *bim = (TH1D*)f->Get("epos_impact_parameter_distribution");
  TH1D *cent = (TH1D*)f->Get("epos_centrality_class_distribution");
  bim->Sumw2();
  cent->Sumw2();
  bim->Rebin(100); // so it will be smoother
  int nbins = bim->GetNbinsX();
  cerr << "bim nbins " << nbins << endl;
  /*
  For all histogram types: nbins, xlow, xup

  bin = 0; underflow bin
  bin = 1; first bin with low-edge xlow INCLUDED
  bin = nbins; last bin with upper-edge xup EXCLUDED
  bin = nbins+1; overflow bin
  */
  for(int i=0+1; i<nbins+1; i++)
  {
    double bimcenter = bim->GetBinCenter(i);
    double mult = (double)(bim->GetBinContent(i));
    double norm = 2 * PI * bimcenter;
    bim->SetBinContent(i,mult/norm);
  }

  gStyle->SetOptStat(1111); // Enable statistics box (1111 = show all stats)
  gStyle->SetStatX(0.4);    // (default is 0.78)
  gStyle->SetStatY(0.9);    // (default is 0.85)
  //gStyle->SetStatW(0.2);    // Width of the stat box
  //gStyle->SetStatH(0.15);   // Height of the stat box

  TCanvas *c1 = new TCanvas();
  bim->Draw();
  c1->Print("figs/binwidthnorm_bimdist.png");
  //f->Close();

  //cerr << cent->GetBinContent(1) << endl;
  //cerr << cent->GetBinContent(10) << endl;
  cent->SetBinContent(1, 2*cent->GetBinContent(1));
  cent->SetBinContent(10, 0.5*cent->GetBinContent(10));
  //cent->Sumw2();
  cent->Draw();
  c1->Print("figs/binwidthnorm_centdist.png");
  
  
  return 0;
}

