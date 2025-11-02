TFile *file = TFile::Open("yourfile.root");
TTree *tree = (TTree*)file->Get("teposevent0");

// Create a histogram from the "bim" branch
tree->Draw("bim>>hBim");

// Retrieve the created histogram
TH1F *hBim = (TH1F*)gDirectory->Get("hBim");

// Apply rebinning
hBim->Rebin(10);

// Draw the rebinned histogram
hBim->Draw();

