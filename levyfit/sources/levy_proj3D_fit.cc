//20250516: draw three dirs in one canvas
//20250516: draw three dirs of different kt bins for a single event
//20250519: draw (alpha, Rout) TH2D
//20250520: read kT infomations from input txt, draw the alpha and Rout... vs mT
//20250523: read bim and centbin from the output info file from RhoDis.C 
//20250527: run for all files in the list
//20250529: save (alpha, R) histograms for one job. to run with batch jobs
//20250530: run for different energies
//20250619: C.L cut modified for Nevts check
//20250619: put EvtNum into the inputfile dir
//20250623: Include the overflow bins when integral 
//20250623: N param[4] analysis
//20250924: N lambda save to Tgrapha
//20251008: Nch used for centrality


// Standard Library
    #include <fstream>
    #include <iostream>
    #include <sstream>
    #include <string>
    #include <vector>
// ROOT Core
    #include <TCanvas.h>
    #include <TFile.h>
    #include <TF1.h>
    #include <TGraph2D.h>
    #include <TGraphErrors.h>
    #include <TH1.h>
    #include <TH2.h>
    #include <TLine.h>
    #include <TMath.h>
    #include <TVector3.h>
// ROOT Graphics and UI
    #include <TLatex.h>
    #include <TLegend.h>
    #include <TEllipse.h>
// ROOT Math
    #include <Math/Factory.h>
    #include <Math/Functor.h>
// Levy reader
    #include "Levy_proj_reader.h"
    #include "infoDataRead.h"

    using namespace std;

    Levy_reader* myLevy_reader;

    double fit_min = 1.0;//1 default
    double fit_max = 30.0;//20 default for Nevts Merged
    double Xaxis_min = 0.2;
    double Xaxis_max = 200;
   // double Xaxis_max = 1.0e17;
    int linewid = 2;

  //  const int NDIR = 3;
    const int NPAR = NDIR + 2; //alpha, Rout, Rside, Rlong, N
    int NDF; // For tracking the number of degrees of freedom in the fit


//centtrality bins define
    //const int kCentBin = 6;
    //const char* centName[6] = {"0-10%","10-20%","20-30%","30-40%","40-50%","50-100%"};


    const int kCentBin = 10;
    const char* centName[10] = {"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};
    //const char* centName[10] = {"0-10%","10-20%","20-30%","30-40%","40-50%","50-100%","0-10%","35-40%","40-100%"};//from the centrality set, 1: 0-5%; 2: 5-10%;
         //centrality bin
      int centbin = 0;
      int centbinNameId=centbin;

//defines for kT bins
      //const int kKtBin = 10;
      Float_t  ktbd[kKtBin+1]={0.0, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 1.0};//17


     // Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.775, 1.0};//merge high kTs, 11
     // Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 1.0};//merge high kTs, 6


   // Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675};

//files directory
    std::string Dir = "/users/huang/EPOS4_Analysis/AuAu_x3ff/RhoDis";//Yan input
    TString outDir("3DFits");//Yan output

    //const std::string Dir = "/home/yan/Levy/3D_Fits/Compare_EPOS3/EPOS3_Evt180/";//EPOS3 compare input
    //TString outDir("/home/yan/Levy/3D_Fits/EPOS4_Sim/LevyFit/Fits/EPOS3_data_180");

    // TString outDir("/home/yan/Levy/3D_Fits/Compare_EPOS3/Fits");//EPOS3_Compare output

    //infile to read the event information 


    //int EvtNum = 10;

    std::vector<TH1*> histograms(3); // 预分配3个位置
    double LevyProj1DFunc(const double *x, const double *par)
    {
      double alpha = par[0];
      double R = par[1];
      double N = par[2];
      double Rcc = (R*pow(2.,1./alpha));
      return (2.*N/Rcc)*(myLevy_reader->getValue_1d(alpha, x[0]/Rcc));
    }

    double logLikelihood(const double *params)
    {
      NDF = 0;
      double logL = 0.0;
      double alpha  = params[0];
      double N  = params[4];
      for (int idir = 0; idir < NDIR; idir++)
      {
        double *thispar = new double[NDIR];
        thispar[0] = alpha;
        thispar[1] = params[idir+1]; // idir = 0,1,2 -> Rout,Rside,Rlong; but Rout=param[1], Rside=param[2], Rlong=param[3];
        thispar[2] = N;
        double integral = histograms[idir]->Integral(1,histograms[idir]->GetNbinsX()+1);
        //cout<<"Nbins: "<<histograms[idir]->GetNbinsX()<<endl;
        for (int ibin = 1; ibin <= histograms[idir]->GetNbinsX(); ++ibin)
        {
          double x = histograms[idir]->GetXaxis()->GetBinCenter(ibin);
          if(x > fit_max || x < fit_min) continue;
          double binwidth = histograms[idir]->GetXaxis()->GetBinWidth(ibin);
          double observed = histograms[idir]->GetBinContent(ibin);
          double expected = LevyProj1DFunc(&x, thispar)*binwidth*integral;
          if(expected <= 0.) continue; // Avoid log(0) or negative expected values
          if(observed != 0.)
          logL += expected + observed*log(observed/expected) - observed;
          NDF++;
          //cout<<"idir: "<<idir<<" NDF: "<<NDF<<endl;
        }
        //std::cout << "this params: " << thispar[0] << ", " << thispar[1] <<", "<<thispar[2] << std::endl;
        delete thispar;
      }
      NDF -= 5;
      return 2.0 * logL;
    }


//Main function here
    int main(int argc, char** argv)
    {

  // Initialize Levy reader object
      myLevy_reader = new Levy_reader("levy_proj3D_values.dat"); // Download this file from https://csanad.web.elte.hu/phys/levy_proj3D_values.dat
      int Nevts = std::atof(argv[1]);//Nevts, how many single events are combined, for the input directory
      int EvtNum = std::atof(argv[2]);//Number of final merged events generated.
      double Ecm = std::atof(argv[3]);
      
      if(Ecm == 62.4) {fit_max = 30;}
 
  //files directory 
      //input
      std::ostringstream oee;
      oee << Form("/%.1fGeV/Nevts%d", Ecm, Nevts);// 自动转换为字符串，比如 "/Ecm27.0"
      std::string Dirh = Dir + "/out_histos" + oee.str();
      std::string DirkT = Dir + "/kTMean" + oee.str();
      //if the out directory not exits, create it
      //output
      outDir += Form("/%.1fGeV/Nevts%d", Ecm, Nevts);
      TString outFits = outDir + Form("/Fits");



      if (gSystem->AccessPathName(outFits)) { // 检查路径是否存在:ml-citation{ref="3" data="citationList"}
        if (gSystem->mkdir(outFits, kTRUE) == -1) {
          Error("SavePDF", "目录创建失败: %s (错误码: %d)", outFits.Data(), TSystem::GetErrno());
          return -1;
        }
      }



//for all events in all files
    char ktName[kKtBin][25];
    for(int i=0;i<kKtBin;i++) sprintf(ktName[i],"%.3f-%.3f",ktbd[i],ktbd[i+1]);
    const char* rhoDirString[NDIR] = {"out", "side", "long"};
//(alpha, R) 2D histograms for all events
    TH2D* hAlpha_R[NDIR][kKtBin];
    TString hNameA_R[NDIR][kKtBin];
    TString hTitleA_R[NDIR][kKtBin];
    for (int idR = 0; idR < NDIR; idR++){
      for(Int_t jIM=0; jIM<kKtBin; jIM++){
        hNameA_R[idR][jIM] = Form("h#alpha_R_%s_%d", rhoDirString[idR], jIM);
        hTitleA_R[idR][jIM] = Form("%s, k_{T} %s GeV ;  #alpha; R_{%s} [fm]",centName[centbinNameId], ktName[jIM], rhoDirString[idR]);//R vs alpha
        hAlpha_R[idR][jIM] = new TH2D(hNameA_R[idR][jIM], hTitleA_R[idR][jIM], 200, 0.0, 2.0, 200, 1.0, 15);//R out
      }//kT bins loop
    }//three directions, out, side, long
    TH1D* hAlpha[kKtBin];
    TH1D* hRout[kKtBin];
    TH1D* hRside[kKtBin];
    TH1D* hRlong[kKtBin];
    TH1D* hN_lambda[kKtBin];
    //C.L cheeck
    TH1D* hFitStatus[kKtBin];
    TH1D* hCovStatus[kKtBin];
    TH2D* hMinStatus[kKtBin];
    for(Int_t jIM=0; jIM<kKtBin; jIM++){
      TString hNameA = Form("halpha_%d", jIM);
      TString hTitleA = Form("%s, k_{T} %s GeV ;  #alpha",centName[centbinNameId], ktName[jIM]);
      TString hNameRo = Form("R_out_%d", jIM);
      TString hTitleRo = Form("%s, k_{T} %s GeV ;  R_{out} [fm]",centName[centbinNameId], ktName[jIM]);
      TString hNameRs = Form("R_side_%d", jIM);
      TString hTitleRs = Form("%s, k_{T} %s GeV ;  R_{side} [fm]",centName[centbinNameId], ktName[jIM]);
      TString hNameRl = Form("R_long_%d", jIM);
      TString hTitleRl = Form("%s, k_{T} %s GeV ;  R_{long} [fm]",centName[centbinNameId], ktName[jIM]);
      TString hNameN = Form("hN_lambda_%d", jIM);
      TString hTitleN = Form("%s, k_{T} %s GeV ;  N(#lambda)",centName[centbinNameId], ktName[jIM]);

      hAlpha[jIM] = new TH1D(hNameA, hTitleA, 200, 0.0, 2.0);
      hRout[jIM] = new TH1D(hNameRo, hTitleRo, 200, 1.0, 15);
      hRside[jIM] = new TH1D(hNameRs, hTitleRs, 200, 1.0, 15);
      hRlong[jIM] = new TH1D(hNameRl, hTitleRl, 200, 1.0, 15);
      hN_lambda[jIM] = new TH1D(hNameN, hTitleN, 200, 0.0, 2.0);
      //C.L cheeck
      hFitStatus[jIM] = new TH1D(Form("hFitStatus_%d",jIM), Form("FitStatus for kTbin %d; Fit Status",jIM), 5, 0, 5);
      hCovStatus[jIM] = new TH1D(Form("hCovStatus_%d",jIM), Form("CovMatrixStatus for kTbin %d; CovMatrix Status",jIM), 5, 0, 5);
      hMinStatus[jIM] = new TH2D(Form("hMinStatus_%d",jIM), Form("CovMatrixStatus vs FitStatus for kTbin %d; Fit Status; CovMatrixStatus",jIM), 5, 0, 5, 5, 0, 5);
    }


  //histograms for all events
    TH1D *hbim = new TH1D("hbim", "impact parameter", 100, 0, 20);
    double bim = -1;

    int nEvtsTotal = 0;
    int nEvtsCent = 0;
  //conflev for all events
    TH1D *hCL = new TH1D("hCL", "confidence level", 100, -0.2, 1.2);
    TH1D *hCLBad = new TH1D("hCLBad", "CL cut", 100, -0.00001, 0.00011);
    TH2D *hCL_kT = new TH2D("hCL_kT", "CL vs kT", 100, 0.0, 1.0, 100, -0.2, 1.2);
    TH2D *hCL_kTBad = new TH2D("hCL_kTBad", "CL vs kT cut", 100, 0.0, 1.0, 100, -0.00001, 0.00011);

  //Read kT, mT mean and error, for all files
    kTMTData data = ReadkTandMTValues(DirkT,kKtBin);



  //find the mini and max fileID in a batch job
    int minFileID = INT_MAX;
    int maxFileID = INT_MIN;

   for(int fi=4; fi<argc;++fi){
     int fileID = std::atoi(argv[fi]);
     if (fileID < minFileID) minFileID = fileID;
     if (fileID > maxFileID) maxFileID = fileID;



     //root files initialization
     std::vector<std::string> root_files;
     for(int i = 0; i <EvtNum; ++i) {
       root_files.push_back(Dirh + "/D_rho_Cent" + std::to_string(centbin) + "_File" + std::to_string(fileID) + "_Evt" + std::to_string(i) + ".root");
     }



     // TString infilename = Form("%s/kTBin_results_%d.txt", DirkT.c_str(), fileID);
      //std::ifstream infileInfo(infilename);


  // Iterating through the histograms and perform the fitting on all of them
      for (int iEvt = 0; iEvt < EvtNum; iEvt++) {
        cout<<"============================================================================================================="<<endl;
        cout<<"======================================"<<" Proceeding Evt: "<<iEvt<<"====================================================="<<endl;
        cout<<"============================================================================================================="<<endl;
        cout<<"-------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<" Proceeding: "<<root_files[iEvt]<<endl;
        cout<<"-------------------------------------------------------------------------------------------------------------"<<endl;
    // Obtain histograms from data file
        //TFile* infile = new TFile("D_rho_0.root"); // An example file
        TFile *infile = new TFile(root_files[iEvt].c_str(), "READ");
        if (!infile || infile->IsZombie()) {
          std::cerr << "Error: Failed to open the file: " << root_files[iEvt] << std::endl;
          delete infile;
          continue;
        }



             //read from outfile ktBin_results from RhoDis.C
        //-------------------read from infile Start---------------------
      /*  if (!infileInfo.is_open()) {
          std::cerr << "文件打开失败: " << strerror(errno) <<",   "<<infilename<< std::endl; // 输出系统错误信息
        }
        std::vector<EventInfo> EvtInfo;

        if (GetEventInfo(infileInfo, iEvt, EvtInfo)) {
          cout<<"-------info read from file: -------"<<endl;
         // for(const auto& evt : EvtInfo){
          //for(Int_t jIM=0; jIM<kKtBin; jIM++){
            const auto& evt = EvtInfo[0];
            cout<<"eventIndex: "<< evt.evtInd<<", bim = "<< evt.bim<<", centbin = "<< evt.centBin<<", ktbin: "<<
            evt.ktbin<<", ktmean = "<< evt.ktmean<<", kterror = "<< evt.kterror<<", npairs = "<< evt.nPairs<<endl;
            //cout<<"eventIndex: "<< EvtInfo[jIM].evtInd<<", bim = "<< EvtInfo[jIM].bim<<", centbin = "<< EvtInfo[jIM].centBin<<", ktbin: "<<
            //EvtInfo[jIM].ktbin<<", ktmean = "<< EvtInfo[jIM].ktmean<<", kterror = "<< EvtInfo[jIM].kterror<<", npairs = "<< EvtInfo[jIM].nPairs<<endl;
           // }
           centbin = evt.centBin;
           bim = evt.bim;
           EvtInfo.clear();
        }*/


      
        cout<<"centbin = "<<centbin<<endl;
        /*if(centbin != 0 && centbin != 1) continue;//0-10%
        nEvtsCent++;
        hbim->Fill(bim);*/


    //output pdf to store all the fits
        //TString pdfFileName = TString(Form("%s/3D_Fits_%d_%d.pdf",outFits.Data(),fileID, iEvt));
        TString pdfFileName = TString(Form("%s/%.1fGeV_3D_Fits_%d_%d.pdf",outFits.Data(),Ecm,fileID, iEvt));

    // Create canvas for plotting
        TCanvas *c1 = new TCanvas("fits", "", 700, 500);
        //TCanvas *c2 = new TCanvas("c2", "", 1000, 1000);
        //c1->SetLogy();
        //c1->SetLogx();
        //c1->SetLeftMargin(0.13);
        c1->Divide(3,1,0,0);
        c1->Print(pdfFileName + "[");  
        



    //loop for different kt bins
        for(Int_t jIM=0; jIM<kKtBin; jIM++){

          //  cout<<"==========================Processing kT bin："<<jIM<<"===================="<<endl;
      //the histo names for different kt bins
          TString hName("hDrho_");
          TString hNameOut = hName + "Out_" + "Evt" + Form("%d",iEvt)  + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
          TString hNameSide = hName + "Side_" + "Evt" + Form("%d",iEvt)  + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
          TString hNameLong = hName + "Long_" + "Evt" + Form("%d",iEvt) + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
          cout<<"hNameout: "<<hNameOut<<" hNmaeSide: "<<hNameSide<<" hNameLong: "<<hNameLong<<endl;
          histograms[0] = (TH1*)infile->Get(hNameOut);  // 索引0对应out
          histograms[1] = (TH1*)infile->Get(hNameSide); // 索引1对应side
          histograms[2] = (TH1*)infile->Get(hNameLong); // 索引2对应long
      // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
      // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
      // Create the function to be minimized
          ROOT::Math::Functor functor(&logLikelihood, NPAR);
          minimizer->SetFunction(functor);
      // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.5, 0.01, 0.5, 2.0);
          minimizer->SetLimitedVariable(1, "Rout",  5.0, 0.01, 1.0, 14.0);
          minimizer->SetLimitedVariable(2, "Rside", 5.0, 0.01, 1.0, 14.0);
          minimizer->SetLimitedVariable(3, "Rlong", 5.0, 0.01, 1.0, 14.0);
          minimizer->SetLimitedVariable(4, "N", 1.0, 0.01, 0.0, 1.0);
          //minimizer->SetVariable(4, "N", 1.0, 0.01);
      // Minimize the function
          minimizer->Minimize();
          minimizer->PrintResults();

          int status_code = minimizer->Status();
          int cov_status = minimizer->CovMatrixStatus();
          cout<<"minimizer status_code = "<<status_code<<",  CovMatrixStatus = "<<cov_status<<endl;
          hFitStatus[jIM]->Fill(status_code);
          hCovStatus[jIM]->Fill(cov_status);
          hMinStatus[jIM]->Fill(status_code, cov_status);



        // Save the results
          const double *params = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          double CL = TMath::Prob(chi2val, NDF);//confidence level
          cerr << "chi^2/NDF = " << chi2val << "/" << NDF << " -> C.L. = " << CL << endl;
          hCL->Fill(CL);
          hCL_kT->Fill(data.kTMean[jIM], CL);
          if(CL <= 0.0001) {
            hCLBad->Fill(CL);
            hCL_kTBad->Fill(data.kTMean[jIM], CL);
          //  continue;
          }
          



        // ROOT functions for drawing the fit
          TF1* f_levy_func_fitted[NDIR];
          TF1* f_levy_func_fullrange[NDIR];
        // Actual plotting in each direction
          int idir;
          for(idir = 0; idir<NDIR; idir++)
          {
            c1->cd(idir+1);
            gPad->SetLogx();
            gPad->SetLogy();
            gPad->SetGridx();
            gStyle->SetGridColor(kGray);
            gStyle->SetGridStyle(2);
            //gPad->SetMargins(0.15, 0.1, 0.1, 0.1); // 调整边距
            if(idir==0) gPad->SetLeftMargin(0.13);
            //gPad->SetLeftMargin(0.13);
        // Rescale histogram for plotting
            double integral = histograms[idir]->Integral(1,histograms[idir]->GetNbinsX()+1);//with overflow bins
            cout<<"integral: "<<integral<<endl;
            for (int ibin = 1; ibin <= histograms[idir]->GetNbinsX(); ibin++)
            {
              double content = histograms[idir]->GetBinContent(ibin);
              double error = histograms[idir]->GetBinError(ibin);
              double binwidth = histograms[idir]->GetXaxis()->GetBinWidth(ibin);
              histograms[idir]->SetBinContent(ibin, content / binwidth);
              histograms[idir]->SetBinError(ibin, error / binwidth);
            }
            histograms[idir]->Scale(1.0 / integral);

        // Set histogram plotting properties and draw it
            histograms[idir]->SetTitle("");
            histograms[idir]->GetXaxis()->SetTitle(Form("|#rho_{LCMS}^{%s}|",rhoDirString[idir]));
            histograms[idir]->GetYaxis()->SetTitle(Form("D(#rho_{LCMS}^{%s})",rhoDirString[idir]));
            histograms[idir]->GetYaxis()->SetTitleOffset(1.5);
            histograms[idir]->SetStats(0);
            //histograms[idir]->GetYaxis()->SetRangeUser(1e-6, 1e-0);
            //histograms[idir]->GetXaxis()->SetRangeUser(1, 200.);
            histograms[idir]->GetYaxis()->SetRangeUser(3e-7, 3);
            //histograms[idir]->GetYaxis()->SetRangeUser(3e-25, 1);
            histograms[idir]->GetXaxis()->SetRangeUser(Xaxis_min, Xaxis_max);
            histograms[idir]->SetMarkerColor(kBlue);
            // histograms[idir]->SetMarkerSize(0.3);
            //histograms[idir]->SetMarkerStyle(8);
            histograms[idir]->SetLineColor(kBlue);
            histograms[idir]->GetYaxis()->SetTitleSize(0.04);
            histograms[idir]->GetXaxis()->SetTitleSize(0.04);
            histograms[idir]->Draw("pe");
            // if(idir>0) histograms[idir]->Draw("SAME NOY");

        // Create a new function with the fitted parameters
            f_levy_func_fitted[idir] = new TF1(Form("f_levy_func_fitted_dir%d",idir), LevyProj1DFunc, fit_min, fit_max, NPAR-2);
            f_levy_func_fitted[idir]->SetParameters(params[0], params[idir+1], params[4]); // alpha, (Rout or Rside or Rlong), N
            f_levy_func_fitted[idir]->SetLineStyle(1); 
            f_levy_func_fitted[idir]->SetLineWidth(linewid);
            f_levy_func_fitted[idir]->SetLineColor(kRed);
            f_levy_func_fitted[idir]->Draw("same");
        // Another function to be drawin in the full range
            f_levy_func_fullrange[idir] = new TF1(Form("f_levy_func_fullrange_dir%d",idir), LevyProj1DFunc, Xaxis_min, Xaxis_max, NPAR-2);
            f_levy_func_fullrange[idir]->SetParameters(params[0], params[idir+1], params[4]); // alpha, (Rout or Rside or Rlong), N
            f_levy_func_fullrange[idir]->SetLineStyle(2);
            f_levy_func_fullrange[idir]->SetLineWidth(linewid);
            f_levy_func_fullrange[idir]->SetLineColor(kRed);
            f_levy_func_fullrange[idir]->Draw("same");

        // Draw parameters on the plot
            TLatex l;
            l.SetNDC();
            l.SetTextSize(0.055);
            if(idir == 0){//alpha
              l.DrawLatex(0.2,0.60,  Form("%s,k_{T} %s GeV", centName[centbin], ktName[jIM]));
              l.DrawLatex(0.2,0.55,  Form("N(#lambda) = %.2f #pm %.2f",params[4],errors[4]));
              l.DrawLatex(0.2,0.5,  Form("#alpha = %.2f #pm %.2f",params[0],errors[0]));
              l.DrawLatex(0.2,0.45, Form("#chi^{2} / NDF = %.2f / %d",chi2val,NDF));
              l.DrawLatex(0.2,0.40, Form("C.L.: %.2f%%",CL * 100));
            } 
            l.DrawLatex(0.2,0.35, Form("R_{%s} = (%.2f #pm %.2f) fm",rhoDirString[idir], params[idir+1], errors[idir+1]));//R

      //2D (alpha, Rout) Fill
            //if(CL > 0.0001)
            if(status_code == 0 && cov_status == 3) {
              hAlpha_R[idir][jIM]->Fill(params[0], params[idir+1]);
            }
            f_levy_func_fitted[idir]->SetBit(kCanDelete);
            f_levy_func_fullrange[idir]->SetBit(kCanDelete);   
          }//loop for each direction

             hAlpha[jIM]->Fill(params[0]);//alpha
             hRout[jIM]->Fill(params[1]);//Rout
             hRside[jIM]->Fill(params[2]);//R
             hRlong[jIM]->Fill(params[3]);//R
             hN_lambda[jIM]->Fill(params[4]);//N_lambda
          c1->Print(pdfFileName);
        }//loop for different kT bins
        c1->Print(pdfFileName + "]");
        delete c1;
	      infile->Close();
        delete infile;

      }//for all events
      nEvtsTotal += EvtNum;
     // infileInfo.close();

    }//for all files in the list
 
    cout<<"=================Finishing files from "<<minFileID<<" to "<<maxFileID<<", "<<nEvtsCent<<Form(" %s Events", centName[centbinNameId])<<" among "<<nEvtsTotal<<" in total ==============="<<endl;


    //========================2D R_Alpha[kKtBin] histograms for all events in all files==================================
      TCanvas *c2 = new TCanvas("c2", "R vs #alpha", 1200, 1600);
      //output pdf to store the 2D (alpha,Rout...) 
      TString pdfNameAlpha_R = TString(Form("%s/Alpha_R_fID_%d_%d.pdf",outDir.Data(), minFileID, maxFileID));
      c2->Print(pdfNameAlpha_R + "["); 
      //c2->cd();
     // hbim->Draw();
      //c2->Print(pdfNameAlpha_R);
      for(Int_t jkt=0; jkt<kKtBin; jkt++){
        c2->Clear();
        c2->Divide(2,2);
        for(int idir = 0; idir<NDIR; idir++){
          c2->cd(idir+1);
          hAlpha_R[idir][jkt]->Draw("colz");
        }//dir loop
        c2->Print(pdfNameAlpha_R);
      }//kT loop
      c2->Print(pdfNameAlpha_R + "]");  
      delete c2;

  //====================================alpha, Rout... vs mT for all events in all files================================
  //save output of alpha, Rout... vs mT to rootfile and pdf
      TFile *outputRoot = new TFile(Form("%s/Alpha_R_mT_fID_%d_%d.root",outDir.Data(), minFileID, maxFileID), "RECREATE");
      TGraphErrors *graph_R_mT[NDIR];
      TGraphErrors *graph_Alpha_mT[NDIR];
      TCanvas *Canvas_Alpha_R = new TCanvas("alpha_R", "#alpha and R vs m_{T}", 1200, 1600);
      TString pdf_R_alpha_mT = TString(Form("%s/Alpha_R_mT_fID_%d_%d.pdf",outDir.Data(), minFileID, maxFileID));

  //Read mT mean and error
      //kTMTData data = ReadkTandMTValues(Dir,kKtBin);
  //Get the alpha and Rout... mean and error from 2D histograms projections
      AlphaRStats alpha_R_stats = ComputeAlphaRStats(hAlpha_R, NDIR, kKtBin);

      //Lambda N mean and stdDev
      double lambdaMean[kKtBin];
      double lambdaStdDev[kKtBin];
      for (int j = 0; j < kKtBin; ++j){
        lambdaMean[j] = hN_lambda[j]->GetMean();
        lambdaStdDev[j] = hN_lambda[j]->GetStdDev();
      }
      TGraphErrors *graph_Lambda_mT;
      graph_Lambda_mT = new TGraphErrors(kKtBin, data.mTMean.data(), lambdaMean , data.mTError.data(), lambdaStdDev);//mT
      graph_Lambda_mT->SetName("graph_mT_lambda");




  //TGraphaErrors for alpha vs mT, Rout... vs mT
      Canvas_Alpha_R->Divide(3,2,0,0); 
      TPad* pad;
      double padHeight, fontScale;
      double mTMax= 1.2;
      for (int dj=0; dj<NDIR; ++dj){
    //alpha
        if(dj==0) {
          Canvas_Alpha_R->cd(dj+1);
          pad = (TPad*)gPad;
          padHeight = pad->GetWh();
          fontScale = padHeight / 1000.0; // 600是你设定的参考pad大小
          gPad->SetLeftMargin(0.15);
          //gPad->SetRightMargin(0.01);   // 收紧右边距
          gPad->SetFrameLineWidth(1);
          graph_Alpha_mT[dj] = new TGraphErrors(kKtBin, data.mTMean.data(), alpha_R_stats.mean_alpha[dj].data(), data.mTError.data(), alpha_R_stats.stdev_alpha[dj].data());//mT
          graph_Alpha_mT[dj]->SetName("graph_mT_alpha");
          graph_Alpha_mT[dj]->SetTitle("");
          // graph_Alpha_mT[dj]->SetTitle(Form("%s, #alpha vs <m_{T}>;<m_{T}> [GeV/c];#alpha",centName[centbinNameId])); // 标题和坐标轴标签
          graph_Alpha_mT[dj]->SetMarkerStyle(20);   
          graph_Alpha_mT[dj]->SetMarkerColor(kRed); 
          graph_Alpha_mT[dj]->SetFillColor(kRed-10);
          graph_Alpha_mT[dj]->SetFillStyle(3001);
          graph_Alpha_mT[dj]->SetLineColor(kBlue); 
          graph_Alpha_mT[dj]->Draw("APE3");  
          graph_Alpha_mT[dj]->GetXaxis()->SetLimits(0.1, mTMax);   // X compared with EPOS3 drawing (0.1,0.8)
          graph_Alpha_mT[dj]->GetYaxis()->SetRangeUser(0.9, 2.2);
          graph_Alpha_mT[dj]->GetXaxis()->SetTitleSize(0.05 * fontScale);
          graph_Alpha_mT[dj]->GetYaxis()->SetTitleSize(0.05 * fontScale);
          graph_Alpha_mT[dj]->GetXaxis()->SetTitleOffset(1.0 / fontScale);
          graph_Alpha_mT[dj]->GetYaxis()->SetTitleOffset(1.0 / fontScale);
          graph_Alpha_mT[dj]->GetYaxis()->SetTitle("#alpha");
          graph_Alpha_mT[dj]->GetXaxis()->SetNdivisions(505);
          graph_Alpha_mT[dj]->Write();
          graph_Lambda_mT->Write();
          TLatex latexA;
          latexA.SetNDC(); // 使用归一化坐标 (0~1)
          latexA.SetTextSize(0.09);
          latexA.DrawLatex(0.20, 0.18, Form("%s, #alpha", centName[centbinNameId]));
        }//only draw alpha for the first 2D (alpha,Rout)

        if(dj==1){//draw the board 
          Canvas_Alpha_R->cd(2);
          //gPad->SetLeftMargin(0.05);
          gPad->SetFillStyle(4000); gPad->SetFrameLineColor(0);
        }

    //Rout, Rside, Rlong
        Canvas_Alpha_R->cd(dj+4);
        pad = (TPad*)gPad;
        padHeight = pad->GetWh();
        fontScale = padHeight / 1000.0; 
        if(dj == 0) gPad->SetLeftMargin(0.15);
        graph_R_mT[dj] = new TGraphErrors(kKtBin, data.mTMean.data(), alpha_R_stats.mean_R[dj].data(), data.mTError.data(), alpha_R_stats.stdev_R[dj].data());//mT
        graph_R_mT[dj]->SetName("graph_mT_R");
        graph_R_mT[dj]->SetTitle(""); 
        //graph_kT_R = new TGraphErrors(kKtBin, kTMean, mean_R,kTError, RErr);//running only one event
        //graph_R_mT[dj]->SetTitle(Form("%s, R_{%s} vs <m_{T}>;<m_{T}> [GeV/c];R [fm]", centName[centbinNameId], rhoDirString[dj])); 
        //graph_kT_R->SetTitle(Form("%s, R vs <m_{T}>;<m_{T}> [GeV/c];R [fm]", centNameMerged[centbin])); 
        graph_R_mT[dj]->SetMarkerStyle(20);   
        graph_R_mT[dj]->SetMarkerColor(kBlue); 
        graph_R_mT[dj]->SetFillColor(kBlue-10);
        graph_R_mT[dj]->SetFillStyle(3001);   
        graph_R_mT[dj]->SetLineColor(kRed); 
        graph_R_mT[dj]->Draw("APE3");  
        graph_R_mT[dj]->GetXaxis()->SetLimits(0.1, mTMax); 
        cout<<"mtMax: "<<mTMax<<endl; 
        graph_R_mT[dj]->GetYaxis()->SetRangeUser(3.5, 11);//Y Range compared to EPOS3 drawing
        graph_R_mT[dj]->GetXaxis()->SetTitleSize(0.05 * fontScale);
        graph_R_mT[dj]->GetYaxis()->SetTitleSize(0.05 * fontScale);
        graph_R_mT[dj]->GetXaxis()->SetTitleOffset(1.0 / fontScale);
        graph_R_mT[dj]->GetYaxis()->SetTitleOffset(1.0 / fontScale);
        if(dj == 0) graph_R_mT[dj]->GetYaxis()->SetTitle("R [fm]");
        if(dj == 2) graph_R_mT[dj]->GetXaxis()->SetTitle("<m_{T}> [GeV/c]");
        graph_R_mT[dj]->GetXaxis()->SetNdivisions(505);
        graph_R_mT[dj]->Write();
        TLatex latex;
        latex.SetNDC(); // 使用归一化坐标 (0~1)
        latex.SetTextSize(0.09);
        latex.DrawLatex(0.20, 0.18, Form("%s, R_{%s}", centName[centbinNameId], rhoDirString[dj]));
      }//different direction loop

      //CL plots
      TCanvas* cCL = new TCanvas("cCL","conf", 1600, 1200);
      cCL->Divide(2,2);
      cCL->cd(1);
      hCL->Draw();
      cCL->cd(3);
      hCLBad->Draw();
      cCL->cd(2);
      hCL_kT->Draw("colz");
      cCL->cd(4);
      hCL_kTBad->Draw("colz");
      cCL->SaveAs(Form("%s/CL.png", outDir.Data()));
      //FitStatus plots
      TCanvas* cstatus = new TCanvas("cstatus","cstatus", 2500, 1000);
      cstatus->Divide(5,2);
      for(int js=0; js<kKtBin; ++js){
        cstatus->cd(js+1);
        hMinStatus[js]->Draw("colz");
      }
      cstatus->SaveAs(Form("%s/FitStatus.png", outDir.Data()));


  //save pdf
      Canvas_Alpha_R->SaveAs(pdf_R_alpha_mT);
      delete Canvas_Alpha_R;
      //save and delete histograms
      for (int i = 0; i < NDIR; ++i) {
        for (int j = 0; j < kKtBin; ++j) {
          hAlpha_R[i][j]->Write();
          delete hAlpha_R[i][j];
        }
      }
      for (int jk = 0; jk < kKtBin; ++jk) {
          hAlpha[jk]->Write();
          hRout[jk]->Write();
          hRside[jk]->Write();
          hRlong[jk]->Write();
          hN_lambda[jk]->Write();
          hFitStatus[jk]->Write();
          hCovStatus[jk]->Write();
          hMinStatus[jk]->Write();
          delete hAlpha[jk];
          delete hRout[jk];
          delete hRside[jk];
          delete hRlong[jk];
        }
      //CL check
      hCL->Write();
      hCLBad->Write();
      hCL_kT->Write();
      hCL_kTBad->Write();
      hbim->Write();
      outputRoot->Close();
      delete myLevy_reader;
      return 0;

    }
