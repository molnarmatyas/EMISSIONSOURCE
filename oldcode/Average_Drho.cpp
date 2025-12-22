// run: root.exe -b -q 'Average_Drho.cpp(true)'
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <iostream>
#include <TString.h>

const char* urqmd_energies[] = {"3p0","3p2","3p5","3p9","4p5","7p7"};
const int NKT=10;
const int NCH=2;

int Average_Drho(int avgN=100, bool urqmd=true, int NEVT=1000)
{
  int engyarr_size;
  const char** energyarray;
  if(urqmd)
  {
    engyarr_size = sizeof(urqmd_energies)/sizeof(urqmd_energies[0]);
    energyarray = urqmd_energies;
  }
  else
  {
    cerr << "EPOS averaging not implemented yet." << endl;
    return 0; // FIXME
    //energyarray = energies;
  }

  cout << "Averaging for " << avgN << " events." << endl;

  for(int ienergy=0; ienergy<engyarr_size; ienergy++)
  {
    const char* energy = energyarray[ienergy];
    cout << "Processing energy: " << energy << endl;
    // Read the file with given energy
    TFile *file = TFile::Open(Form("drho_analyze/AuAu_%s_drho.root",energy));
    if (!file)
    {
      std::cerr << "Error opening file" << std::endl;
      return 1;
    }
    TFile *file_output = new TFile(Form("drho_analyze/AuAu_%s_drho_avg.root",energy), "RECREATE");
    if (!file_output)
    {
      std::cerr << "Error creating output file" << std::endl;
      return 1;
    }

    for(int ich=0; ich<NCH; ich++)
    {
      cout << "Processing charge (-- = 0 or ++ = 1): " << ich << endl;
      for(int ikt=0; ikt<NKT; ikt++)
      {
        cout << "Processing kt bin: " << ikt << endl;
        // Loop over events and average histograms
        int newevt = 0;
        for(int ievt=0;ievt<NEVT;ievt+=avgN)
        {
          TH1D* hist = (TH1D*)file->Get(Form("D_LCMS_ev%d_ch%d_KT%d", ievt,ich,ikt));
          if (!hist)
          {
            std::cerr << "Error: Histogram D_LCMS_ev" << ievt << "_ch" << ich << "_KT" << ikt <<" not found in the file." << std::endl;
            continue;
          }
          hist->Reset();
          hist->SetName(Form("D_LCMS_ev%d_ch%d_KT%d", newevt,ich,ikt));
          newevt++;
          for(int jevt=ievt+1; jevt<ievt+avgN; jevt++)
          {
            if( (TH1D*)file->Get(Form("D_LCMS_ev%d_ch%d_KT%d", jevt,ich,ikt)) )
            {
              hist->Add((TH1D*)file->Get(Form("D_LCMS_ev%d_ch%d_KT%d", jevt,ich,ikt)));
            }
            else
            {
              std::cerr << "Error: Histogram D_LCMS_ev" << jevt << "_ch" << ich << "_KT" << ikt << " not found in the file." << std::endl;
              break;
            }
          }
          file_output->cd();
          hist->Write();
        } // end of event loop
      } // end of kt loop
    } // end of ch loop

    file_output->Close();
    file->Close();

  } // end of energy loop
  cout << "All done!" << endl;

  return 0;
}
