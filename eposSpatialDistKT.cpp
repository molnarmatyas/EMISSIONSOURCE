/*
Original from:
Santiago Bernal Langarica
October, 2024

Modified by:
Molnar, Matyas
January 2025
*/
// to be run like
// root.exe -b -q eposSpatialDistKT.cpp\(0,100,1,0,\"200GeV\"\)

// TODO: centrality filtering to be implemented

#include <TH1.h>
#include <TFile.h>
#include <TMath.h>

#include <stdio.h>
#include <stdlib.h>

#include <fstream>

#include <iostream>
#include <cstring>

// -----
#ifndef __CINT__
#include <zlib.h>  // from distPoints.cpp, what is this for???
#endif
char fbuffer[256];
//const Double_t pi = TMath::Pi();
// -----

using namespace std;

const Int_t NbinsGtr = 46, Nbinsless = 6;
const Double_t pi = TMath::Pi();
const double Mass2_pi = 0.019479835;

// !!! cuts
const auto pTmin = 0.15; //0.1
const auto pTmax = 1.0;  //2.0
const auto etamax = 1.0;
double qLCMSmax = 0.1; // will change acc. to given KT bin

const int nkT = 10;
const double kT_bounds[nkT+1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};

void WRONGNORM_rdivision(TH1D *h, Int_t N)
{
    long int Integral = h->GetEntries();
    cerr << "Integral = " << Integral << endl;
    Int_t nbins = N;
    for (Int_t i = 0; i <= nbins; i++)
    {
        auto BinValue = h->GetBinContent(i);
        auto BinCenter = h->GetBinCenter(i);
        auto BinWidth = h->GetBinWidth(i);
        auto BinError = h->GetBinError(i);
        BinError = sqrt(BinValue); //h->GetBinError(i);
        
        auto Scale = 4. * pi * BinCenter * BinCenter * BinWidth * Integral;
        h->SetBinContent(i, BinValue / Scale);
        h->SetBinError(i, BinError / Scale);
    }
}

void rdivision(TH1D *h, Int_t N)
{
  //
}
void rdivisionSTOP(TH1D *h, Int_t N)
{
    Int_t nbins = N;
    h->Sumw2();
    //h->ComputeIntegral();
    //cerr << "Integral = " << h->Integral() << endl;
    double Integral = h->Integral();
    //cerr << "GetEntries pre: = " << h->GetEntries() << endl;
    for (Int_t i = 0; i <= nbins; i++)
    {
        auto BinValue = h->GetBinContent(i);
        auto BinCenter = h->GetBinCenter(i);
        auto BinWidth = h->GetBinWidth(i);
        auto BinError = h->GetBinError(i);
        BinError = sqrt(BinValue); //h->GetBinError(i);
        
        auto Scale = 4. * pi * BinCenter * BinCenter * BinWidth;
        h->SetBinContent(i, BinValue / Scale);
        h->SetBinError(i, BinError / Scale);
    }
    //h->ComputeIntegral();
    //double Integral = h->Integral();

    //h->Print("all");
    h->Scale(1.0/Integral);
    cerr << "Integral = " << Integral << endl;
}

void distPoints_(Int_t nkT)
{
    /* Implemented from distPoints.cpp */
    cout << "----------------------------" << endl;
    cout << "Creating individual kT files" << endl;
    TFile *inputfile = TFile::Open("distance_dists_EPOS_KT.root");
    TH1D *spatdistLCMS[nkT];
    for(int iFile=0; iFile<nkT; iFile++)
    {
      spatdistLCMS[iFile] = (TH1D*) inputfile->Get(Form("spatdistLCMS%d",iFile));
      Int_t Nbins = spatdistLCMS[iFile]->GetXaxis()->GetNbins();
    
      ofstream myfile;
      myfile.open(Form("distance_dists_EPOS_KT%d.txt",iFile));
      for (int i = 1; i <= Nbins; i++)
      {
          auto bincenter = spatdistLCMS[iFile]->GetBinCenter(i);
          auto binvalue = spatdistLCMS[iFile]->GetBinContent(i);
          auto binerr = spatdistLCMS[iFile]->GetBinError(i);
          myfile << bincenter << "  " << binvalue << "  " << binerr << "\n";
      }
      myfile.close();
      cout << "File " << iFile << " ready!" << endl;
    }
}


// -------------------- MAIN ----------------------
void eposSpatialDistKT(Int_t ini = 0, Int_t nevt = 100, Int_t nfiles = 1, Int_t evtfrom = 0, const char* energy = "9.2GeV") //nfiles tot = 1334
{
    TStopwatch timer;
    timer.Start();

    // TFile out("distance_dists_EPOS_KT-10_files1.root","RECREATE");
    TFile out("distance_dists_EPOS_KT.root","RECREATE");
    const int Nbins = 168;
    double xmin = 1.e-1;
    double xmax = 1.e3;
    double binfactor = TMath::Power(xmax/xmin,1.0/Nbins);
    double logbins[Nbins+1];
    logbins[0] = xmin;
    for(int ibin = 1; ibin <= Nbins; ibin++)
        logbins[ibin] = xmin * TMath::Power(binfactor,ibin);
    TH1D* spatdistLCMSKT[nkT];
    for(int iKT=0; iKT<nkT; iKT++)
    {
        spatdistLCMSKT[iKT] = new TH1D(Form("spatdistLCMS%d",iKT)
                                      ,Form("Spatial size of charged pion source in the LCMS for kT #in [%.3f,%.3f] GeV",kT_bounds[iKT],kT_bounds[iKT+1])
                                      ,Nbins,logbins);
    }
    TH1D *KTdist = new TH1D("ktdist","K_{T} distribution",100,0,10);
    // do not Sumw2() the dists yet!

    Int_t NoFilesNotWorking = 0, NoFilesNotExisting = 0, NoZombieFiles = 0, NoEmptyFiles = 0;
    std::vector< int > FilesNotWorking, FilesNotExisting, ZombieFiles, EmptyFiles;

    Int_t EventCounter =  0;

    for (Int_t nf = ini; nf < (ini + nfiles); nf++)
    {
        char inputfilename[100];
        snprintf(inputfilename,100,"z-prueba_EOS2_%s.root",energy);
        //snprintf(inputfilename,100,"/phenix/u/kincsesd/plhf1/EPOS/files/z-gg2my59a-1201.root",energy);
        //snprintf(inputfilename,100,"/storage/mexnica/Mate/EPOS4/AuAu_7.7GeV_100evts_100files/z-prueba_EOS2.%d.root",nf);

        if (gSystem->AccessPathName(inputfilename))
        {
            cout << "File " << nf << " does not exist in that folder" << endl;
            NoFilesNotWorking++;
            FilesNotWorking.push_back(nf);
            NoFilesNotExisting++;
            FilesNotExisting.push_back(nf);
            continue;
        }

        TFile *inputfile = TFile::Open(inputfilename);

        if (!inputfile || inputfile->IsZombie())
        {
            cout << "File " << nf << " is empty or zombie" << endl;
            NoFilesNotWorking++;
            FilesNotWorking.push_back(nf);
            NoZombieFiles++;
            ZombieFiles.push_back(nf);
            continue;
        }

        if (gDirectory->FindObject("teposevent"))
        {
            cout << "File " << nf << " is empty or does not contain the event tree" << endl;
            NoFilesNotWorking++;
            FilesNotWorking.push_back(nf);
            NoEmptyFiles++;
            EmptyFiles.push_back(nf);
            continue;
        }

        if(inputfile->FindKeyAny("teposevent") == nullptr)
        {
            cout << "File " << nf << " is empty or does not contain the event tree" << endl;
            NoFilesNotWorking++;
            FilesNotWorking.push_back(nf);
            NoEmptyFiles++;
            EmptyFiles.push_back(nf);
            continue;
        }
        
        
        // TTree *header = (TTree*) inputfile->Get("teposhead");
        // Int_t entriesHeader = header->GetEntries();
        TTree *event = (TTree*) inputfile->Get("teposevent");
        Int_t entriesEvent = event->GetEntries();
        if(entriesEvent>nevt) entriesEvent=nevt;

        cout << "File " << nf << " loaded" << endl;

        // Int_t iversn, laproj, latarg, maproj, matarg, nfull;
        // Float_t engy;

        Int_t NEvTot = 150000;
        Float_t Px[NEvTot], Py[NEvTot], Pz[NEvTot], ene[NEvTot], x[NEvTot], y[NEvTot], z[NEvTot], Ti[NEvTot], bim;
        Int_t id[NEvTot], ist[NEvTot], np;

        // header->SetBranchAddress("iversn",&iversn);
        // header->SetBranchAddress("laproj",&laproj);
        // header->SetBranchAddress("latarg",&latarg);
        // header->SetBranchAddress("maproj",&maproj);
        // header->SetBranchAddress("matarg",&matarg);
        // header->SetBranchAddress("engy",&engy);
        // header->SetBranchAddress("nfull",&nfull);

        event->SetBranchAddress("px",Px);
        event->SetBranchAddress("py",Py);
        event->SetBranchAddress("pz",Pz);
        event->SetBranchAddress("e",ene);
        event->SetBranchAddress("x",x);
        event->SetBranchAddress("y",y);
        event->SetBranchAddress("z",z);
        event->SetBranchAddress("t",Ti);
        event->SetBranchAddress("id",id);
        event->SetBranchAddress("np",&np);
        event->SetBranchAddress("ist",&ist);
        event->SetBranchAddress("bim",&bim);

        // header->GetEntry(0);

        // auto vsn = iversn;
        // auto zproj = laproj;
        // auto ztarg = latarg;
        // auto Mproj = maproj;
        // auto Mtarg = matarg;
        // auto sqrtsNN = engy;
        // auto TotEv = nfull;

        for (Int_t ievent = evtfrom; ievent < entriesEvent; ievent++)
        {
            event->GetEntry(ievent);
            cerr << "Processing event index " << ievent << endl;

            auto Npart = np;
            auto ENE = ene;
            auto PX = Px;
            auto PY = Py;
            auto PZ = Pz;
            auto time = Ti;
            auto rx = x;
            auto ry = y;
            auto rz = z;
            auto ID = id;
            auto IST = ist;
            // auto b = bim;
            
            EventCounter++;
            

            for (int itrack = 0; itrack < Npart; itrack++)
            {
                if (fabs(ID[itrack]) == 120)
                {
                    // !!! ist cut for particle => final stage
                    if(IST[itrack] != 0) continue;
                    auto pTi = sqrt(PX[itrack]*PX[itrack] + PY[itrack]*PY[itrack]);
                    auto pi = sqrt(pTi*pTi + PZ[itrack]*PZ[itrack]);
                    auto Ei = sqrt(pi*pi + Mass2_pi); //ENE[itrack];
                    auto etai = 0.5 * TMath::Log(TMath::Abs((pi+PZ[itrack])/(pi-PZ[itrack])));
                    // !!! pT cut & eta cut
                    if(pTi>pTmax || pTi<pTmin || TMath::Abs(etai)>etamax) continue;
                    for (int jtrack = itrack+1; jtrack < Npart; jtrack++)
                    {
                        if (ID[jtrack] == ID[itrack])
                        {
                            // !!! ist cut for pair => final stage
                            if(IST[jtrack] != 0) continue;
                            auto pTj = sqrt(PX[jtrack]*PX[jtrack] + PY[jtrack]*PY[jtrack]);
                            auto pj = sqrt(pTj*pTj + PZ[jtrack]*PZ[jtrack]);
                            auto Ej = sqrt(pj*pj + Mass2_pi);//ENE[jtrack]; 
                            auto etaj = 0.5 * TMath::Log(TMath::Abs((pj+PZ[jtrack])/(pj-PZ[jtrack])));
                            if(pTj>pTmax || pTj<pTmin || TMath::Abs(etaj)>etamax) continue;

                            auto t = time[itrack] - time[jtrack];
                            auto Rx = rx[itrack] - rx[jtrack];
                            auto Ry = ry[itrack] - ry[jtrack];
                            auto Rz = rz[itrack] - rz[jtrack];
                            
                            auto s2 = t*t - (Rx*Rx + Ry*Ry + Rz*Rz);
                            auto r2 = Rx*Rx + Ry*Ry + Rz*Rz;

                            auto K0 = 0.5 * (Ei + Ej);//(ENE[itrack] + ENE[jtrack]);
                            auto Kx = 0.5 * (PX[itrack] + PX[jtrack]);
                            auto Ky = 0.5 * (PY[itrack] + PY[jtrack]);
                            auto Kz = 0.5 * (PZ[itrack] + PZ[jtrack]);
                            
                            auto Q0 = (Ei - Ej);//ENE[itrack] - ENE[jtrack];
                            auto Qx = PX[itrack] - PX[jtrack];
                            auto Qy = PY[itrack] - PY[jtrack];
                            auto Qz = PZ[itrack] - PZ[jtrack];

                            auto qZLCMS2 = (PZ[itrack]*Ej-PZ[jtrack]*Ei)*(PZ[itrack]*Ej-PZ[jtrack]*Ei) / (K0*K0 - Kz*Kz); // 4* not needed, defined in K#
                            auto qLCMS = TMath::Sqrt(Qx*Qx + Qy*Qy + qZLCMS2);

                            auto KT = TMath::Sqrt(Kx*Kx + Ky*Ky);
                            auto cosphi = Kx / KT;
                            auto sinphi = Ky / KT;

                            auto rhoout = Rx * cosphi + Ry * sinphi - (KT / (K0*K0 - Kz*Kz)) * (K0 * t - Kz * Rz);
                            auto rhoside = - Rx * sinphi + Ry * cosphi;
                            auto rholong = (K0 * Rz - Kz * t) / TMath::Sqrt(K0*K0 - Kz*Kz);

                            auto rho2 = rhoout*rhoout + rhoside*rhoside + rholong*rholong;

                            KTdist->Fill(KT);
                          
                            double sqrtrho2 = TMath::Sqrt(rho2);
                            for(int iKT=0; iKT<nkT; iKT++)
                            {
                                // !!! qLCMS cut included
                                double cut_multiplier = 1.; // e.g. 1., 0.8, 2.0 etc. or could even be iKT-dependent
                                //qLCMSmax = kT_bounds[iKT+1];
                                 qLCMSmax = KT;
                                qLCMSmax *= cut_multiplier;
                                if(KT>kT_bounds[iKT] && KT<kT_bounds[iKT+1] && qLCMS<=qLCMSmax)
                                {
                                    spatdistLCMSKT[iKT]->Fill(sqrtrho2);
                                }
                            }

                        } // pions, EPOS id=120 with...
                    } // Loop over pairs
                } // ...pions, EPOS id=120
            } // Loop over particles
        } // Loop over events
    } // Loop over files

    // do not Sumw2() even here!
    // wrong rdivision: spatdistLCMSKT1->Scale(1. / (4. * pi * (spatdistLCMSKT1->GetEntries())));
    //                  spatdistLCMSKT1->Scale(1. / (4. * pi * EventCounter * (spatdistLCMSKT1->GetEntries())));

    for(int iKT=0; iKT<nkT; iKT++)
    {
        rdivisionSTOP(spatdistLCMSKT[iKT],Nbins);
    }

    KTdist->Sumw2();

    for(int iKT=0; iKT<nkT; iKT++)
    {
        out.WriteObject(spatdistLCMSKT[iKT], Form("spatdistLCMS%d",iKT));
    }
    out.WriteObject(KTdist,"KTdist");

    out.Close();

    cout << "Non working files: ";
    copy(FilesNotWorking.begin(), FilesNotWorking.end(),ostream_iterator<int>(cout, " "));
    cout << endl;
    
    cout << "Non existing files: ";
    copy(FilesNotExisting.begin(), FilesNotExisting.end(),ostream_iterator<int>(cout, " "));
    cout << endl;

    cout << "Zombie files: ";
    copy(ZombieFiles.begin(), ZombieFiles.end(),ostream_iterator<int>(cout, " "));
    cout << endl;

    cout << "Empty files: ";
    copy(EmptyFiles.begin(), EmptyFiles.end(),ostream_iterator<int>(cout, " "));
    cout << endl;

    cout << "Number of files not working: " << NoFilesNotWorking << " / " << nfiles << endl;
    cout << "Number of non-existing files: " << NoFilesNotExisting << " / " << NoFilesNotWorking << endl;
    cout << "Number of zombie files: " << NoZombieFiles << " / " << NoFilesNotWorking << endl;
    cout << "Number of empty files: " << NoEmptyFiles << " / " << NoFilesNotWorking << endl;


    // ------------------------------------------
    //     IMPLEMENTED FROM distPoints.cpp
    // ------------------------------------------
    distPoints_(nkT);

    timer.Print();
    exit(0);
}
