//20250401:space and momentum distributions draw
//20250406:draw np and bim, centrality check.
//20250406:R and alpha, updated for levy_fit
//20250420: npartproj+nparttarg, ngl check,and 2D plots with bim
//20250420:kT bins
//20250422:centbin, kTbinz--> only ktBin,for each event
//20250430: kt bin mean print out,info file saved, including event, centbin, ktbin, ktMean,
//to do:Histos/D_rho_eventi.root generated for each event 
//20250502 kt Rebinning, to check the p=0.0 issue
//20250505 kt histogram for all events
//20250508 QLCMS cut applied
//20250509 QLCMS distribution
//20250513 kT rebins, remove the first and last bin
//20250520: improved calculation from (t2-t1)/2 to t2-t1
//20250521: compare D_rho for a single event, binning, weight check, t2-t1
//20250521: weight modified, weight when rebinning
//20250522: eta! mass! double rho! debug
//20250523: 1D hist no weight, no rebinning and integral, do this when fit
//20250530: analysis for different energies
//20250616: Nevents merged for low energies 
//20250617: Nevents merged for the same centbin
//20250619: run Nevts batch, add Nevts to output path
//20250701: kt bins extented to larger for 200GeV lambda check 
//20250911: Qlcms cut (Longitudinally Co-Moving System)
//20251005: merge high kTs, draw Nch vs b
//20251008: hNch used to define the centrality


#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <iostream>
#include <cmath>
#include <map>

// Define constants for cuts
const float PtMin = 0.15; // GeV/c
const float PtMax = 1.0;  // GeV/c
const float EtaMax = 1.0; // Pseudorapidity cut
const int npMax = 150000;//Size of the arrays, Maximum number of particles in any event npM is 83535//84000
Int_t npM = 0; // Maximum number of particles in any event, initialize npM to 0
double MassPi = 0.13957039;//mass of charged pions in GeV/c2 

const int kCentBin = 10;
//Float_t centbd[kCentBin+1] = {0, 3.34, 4.75, 6.71, 8.21, 9.46, 14.89};

//std::vector<float> centbd = {0, 3.34, 4.75, 6.71, 8.21, 9.46, 14.89};
//std::vector<float> centbd = {0, 4.75, 6.71, 8.21, 9.46, 10.57, 30};//backup 
//std::vector<float> centbd = {0, 4.60, 6.71, 8.21, 9.46, 10.57, 30};//ToBeFIX   {4.60, 6.60, 8.00, 9.40, 10.40, 11.40, 12.40, 13.40, 14.20, 15.00}

std::vector<float> centbd = {0, 4.70, 6.70, 8.10, 9.40, 10.40, 11.40, 12.40, 13.40, 14.30, 18.40};//ToBeFIX   {4.60, 6.60, 8.00, 9.40, 10.40, 11.40, 12.40, 13.40, 14.20, 15.00}
std::vector<float> nchbd = {4000.00, 2203.00, 1570.00, 1128.00, 817.00, 599.00, 424.00, 285.00, 214.00, 182.00, 0.0};//ToBeFIX   {3535.00, 2203.00, 1570.00, 1128.00, 817.00, 599.00, 424.00, 285.00, 214.00, 182.00, 175.00}
std::vector<float> nchbd_hist = {0.0, 182.0, 214.0, 285.0, 424.0, 599.0, 817.0, 1128.0, 1570.0, 2203.0, 4000.0};
//const char* centName[kCentBin] = {"0-10%",10-20%","20-30%","30-40%","40-50%","50-100%"};

//Float_t centbd[kCentBin+1] = {0,3.32,4.73,6.68,8.17,9.42,15.0};
//const char* centName[kCentBin] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-100%"};
Int_t centbin = -1;//from the centrality set, 1: 0-5%; 2: 5-10%;
// 每个 centbin 一个事件缓存队列
std::map<int, std::vector<int>> centbin_buffer;
// 每个 centbin 一个 merge 组编号
std::map<int, int> merge_index;

//const int kKtBin = 9;
//Float_t  ktbd[kKtBin+1]={0.0, 0.2, 0.25, 0.3, 0.35, 0.40, 0.45, 0.5, 0.6,1.0};//10
//Float_t  ktbd[kKtBin+1]={0.0, 0.2, 0.3, 0.35, 0.40, 0.45, 0.5, 1.0};//7


const int kKtBin = 17;
//Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675};//11
Float_t  ktbd[kKtBin+1]={0.0, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 1.0};//17

//==============kT bins Merge Check=====================
//Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.775, 1.0};//merge high kTs, 11
//Float_t  ktbd[kKtBin+1]={0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 1.0};//merge high kTs, 6



int NumRun = 10;
bool test = false;//true, false

bool Matyas = false;



bool isCharged(int id) {
	//main particles for each species
    //static const int charged_ids[] = {120, -120, 130, -130, 1120, 1130, 2230, 2330, 121, -121, 131, -131, 12, -12, 11, -11, 14, -14, 13, -13, 16, -16, 15, -15, -240, 240, 2231, 2331, 3331};
    //only:charged kaons pions protons
    static const int charged_ids[] = {120, -120, 130, -130, 1120, 1130, 2230, 2330, 121, -121, 131, -131, 12, -12, 11, -11, 14, -14, 13, -13, 16, -16, 15, -15, -240, 240, 2231, 2331, 3331};
    return std::find(std::begin(charged_ids), std::end(charged_ids), id) != std::end(charged_ids);
}




// Struct to save the info of a single particle
struct Particle {
	Int_t id, ist;
	Float_t t, x, y, z;  
	Float_t energy, px, py, pz, pt, eta;  
}; 


//Function to calculate the phase space factor
double PhaseSpaceFactor(double rho, int NPairs) {
	// return  4 * TMath::Pi() * rho * rho / (NPairs * (NPairs-1)/2);
	return  4 * TMath::Pi() * rho * rho / NPairs ;
	//double weight = 1.0 / (4 * TMath::Pi() * rho * rho) / (EvtNum*particlesStr.size()*(particlesStr.size()-1)/2);
}



void RhoDis(const std::vector<std::string>& inputFiles, const std::vector<int>& fileIDs, double Ecm, int Nevents)//Main function to draw the rho distribution
{
    TChain tree("teposevent");

    for (const auto& fname : inputFiles) {
        tree.Add(fname.c_str());
        std::cout << " Added file: " << fname << std::endl;
    }




	//TFile* in_File = new TFile(inputFile);
	// TFile* in_File = new TFile("z-200GeV_0-10cent_1evt.root");

	//if(Matyas) {in_File = new TFile("Merged_Matyas_Cent0_Evt1000.root");}
	//if(Matyas) {in_File = new TFile("z-200GeV_0-10cent_1evt.root");}

	//TTree* tree = (TTree*)in_File->Get("teposevent");

//centrality bins
	//if(Ecm == 27) {centbd = {0, 3.38, 4.65, 6.69, 8.15, 9.45, 15.13};}
	//if(Ecm == 130) {centbd = {0, 3.34, 4.75, 6.73, 8.24, 9.52, 15.17};}
	//if(Ecm == 200) {centbd = {0, 3.37, 4.78, 6.75, 8.27, 9.54, 15.23};}
	if(Ecm == 11.5) {
		centbd = {0.0, 4.80, 6.80, 8.20, 9.40, 10.60, 11.60, 12.50, 13.40, 14.30, 19.40};
		nchbd_hist = {0.0, 195.00, 241.00, 329.00, 486.00, 722.00, 1033.00, 1486.00, 2065.00, 2991.00, 5000.00};}

	if(Ecm == 14.5) {
		centbd = {0.0, 4.70, 6.60, 8.10, 9.30, 10.50, 11.50, 12.50, 13.30, 14.30, 19.80};
		nchbd_hist = {0.0, 202.00, 261.00, 380.00, 560.00, 832.00, 1212.00, 1715.00, 2425.00, 3486.00, 5700.00};}

	if(Ecm == 19.6) {
		centbd = {0.0, 4.80, 6.60, 8.10, 9.50, 10.50, 11.60, 12.50, 13.40, 14.40, 19.70};
		nchbd_hist = {0.0, 213.00, 291.00, 444.00, 661.00, 1007.00, 1452.00, 2120.00, 3014.00, 4334.00, 7000.00};}

	if(Ecm == 27) {
	  centbd = {0.0, 4.80, 6.60, 8.20, 9.50, 10.50, 11.50, 12.40, 13.30, 14.20, 19.10};
	  nchbd_hist = {0.0, 225.00, 316.00, 493.00, 780.00, 1189.00, 1700.00, 2489.00, 3572.00, 5076.00, 8200.00};}//for 0-10% merged

	if(Ecm == 39) {
		centbd = {0.0, 4.70, 6.80, 8.30, 9.50, 10.60, 11.70, 12.50, 13.40, 14.30, 20.30};
		nchbd_hist = {0.0, 229.00, 327.00, 530.00, 836.00, 1286.00, 1923.00, 2850.00, 4143.00, 6038.00, 9800.00};}

	if(Ecm == 62.4) {
		centbd = {0.0, 4.60, 6.50, 8.10, 9.40, 10.60, 11.60, 12.50, 13.30, 14.40, 18.40};
		nchbd_hist = {0.0, 271.00, 444.00, 705.00, 1096.00, 1676.00, 2606.00, 3827.00, 5360.00, 7568.00, 13000.00};}

	if(Ecm == 130) {
		centbd = {0.0, 4.70, 6.80, 8.30, 9.50, 10.70, 11.70, 12.60, 13.50, 14.50, 19.20};
		nchbd_hist = {0.0, 313.00, 503.00, 857.00, 1432.00, 2261.00, 3453.00, 5018.00, 7115.00, 10616.00, 18000.00};}

	if(Ecm == 200) {
		//centbd = {0, 4.78, 4.78, 6.75, 8.27, 9.54, 15.23};
		centbd = {0, 4.60, 6.60, 8.20, 9.50, 10.70, 11.60, 12.60, 13.40, 14.40, 20.00};
		nchbd_hist = {0.0, 342.00, 592.00, 1062.00, 1674.00, 2693.00, 4217.00, 6048.00, 8898.00, 12607.00, 20000.00};}


     //Qlcms cut
	double Qcut_c = 10;
	if(Ecm == 7.7 || Ecm == 9.2 || Ecm == 11.5 || Ecm == 14.5) {Qcut_c = 12;}
	if(Ecm == 19.6 || Ecm == 27 || Ecm == 39 || Ecm == 62.4 || Ecm == 130 || Ecm == 200) {Qcut_c = 10;}

	cout<<"Ecm = "<<Ecm<<" Qcut_c = "<<Qcut_c<<endl;



	Int_t np = -1;  // Number of particles in the event, initialized to 0
	Int_t npartproj = -1; 
	Int_t nparttarg = -1;
	Int_t ngl = -1;
	Float_t bim=-100.0;//impact parameter


//define important variables for calculation
	Int_t  ist[npMax],id[npMax];
	Float_t t[npMax], x[npMax], y[npMax], z[npMax], mass[npMax], energy[npMax], px[npMax], py[npMax], pz[npMax], pt[npMax], eta[npMax];
	//initialize raw arrays
	for (int j = 0; j < npMax; j++) {
		ist[j] = 0.0;
		id[j] = 0.0;
		t[j] = 0.0;
		x[j] = 0.0;
		y[j] = 0.0;
		z[j] = 0.0;
		mass[j] = 0.0;
		energy[j] = 0.0;
		px[j] = 0.0;
		py[j] = 0.0;
		pz[j] = 0.0;
		pt[j] = 0.0;
		eta[j] = 0.0;
	}
	//structure to save selected particles
	std::vector<Particle> particlesStr;

// Set up branch addresses
	/*tree->SetBranchAddress("np", &np);  // Link the "np" branch to the np variable
	tree->SetBranchAddress("npartproj", &npartproj); 
	tree->SetBranchAddress("nparttarg", &nparttarg); 
	tree->SetBranchAddress("ngl", &ngl); 
	tree->SetBranchAddress("bim", &bim);
	tree->SetBranchAddress("ist", &ist);
	tree->SetBranchAddress("id", &id);
	tree->SetBranchAddress("t", &t);
	tree->SetBranchAddress("x", &x);
	tree->SetBranchAddress("y", &y);
	tree->SetBranchAddress("z", &z);
	tree->SetBranchAddress("e", &mass);
	tree->SetBranchAddress("px", &px);
	tree->SetBranchAddress("py", &py);
	tree->SetBranchAddress("pz", &pz);*/

    tree.SetBranchAddress("np", &np);  // Link the "np" branch to the np variable
	tree.SetBranchAddress("npartproj", &npartproj); 
	tree.SetBranchAddress("nparttarg", &nparttarg); 
	tree.SetBranchAddress("ngl", &ngl); 
	tree.SetBranchAddress("bim", &bim);
	tree.SetBranchAddress("ist", &ist);
	tree.SetBranchAddress("id", &id);
	tree.SetBranchAddress("t", &t);
	tree.SetBranchAddress("x", &x);
	tree.SetBranchAddress("y", &y);
	tree.SetBranchAddress("z", &z);
	tree.SetBranchAddress("e", &mass);
	tree.SetBranchAddress("px", &px);
	tree.SetBranchAddress("py", &py);
	tree.SetBranchAddress("pz", &pz);

	gStyle->SetLineWidth(2);
	gStyle->SetHistLineWidth(2);
	gROOT->ForceStyle();


	TString path(Form("out_histos/%.1fGeV/Nevts%d", Ecm, Nevents));
	TString pathMean(Form("kTMean/%.1fGeV/Nevts%d", Ecm, Nevents));
	if (gSystem->AccessPathName(path)) {
		if (gSystem->mkdir(path, kTRUE) != 0) {
		std::cerr << "failed to creat the new fold [" << path << "]" << std::endl;
		}
	}
//creat pathMean
	if (gSystem->AccessPathName(pathMean)) {
		if (gSystem->mkdir(pathMean, kTRUE) != 0) {
		std::cerr << "failed to creat the new fold [" << pathMean << "]" << std::endl;
		}
	}



	// Creating histogram for D(rho) with log-binning
	const int nbinsRho = 168;//60 default
	double rho_min = 0.1;
	//double rho_max = 1000;//300.0;//1200
	double rho_max = 1.0e17;//300.0;//1200
	//const int nbinsRho = 16800;
	double binfactor = TMath::Power(rho_max/rho_min,1.0/nbinsRho);
	double rho_bins[nbinsRho+1];
	rho_bins[0] = rho_min;
	for(int ibin = 1; ibin <= nbinsRho; ibin++){
		rho_bins[ibin] = rho_min * TMath::Power(binfactor,ibin);
		// cout<<ibin<<", rho_bin: "<<rho_bins[ibin]<<endl;
	}


	TH1F * hmKt = new TH1F("hmKt","Kt bin finder",kKtBin,ktbd);
	//const int kCentBin = 6;
	TH1F *hCent = new TH1F("hCent", "Centrality Total",kCentBin,centbd.data());
	hCent->Sumw2();


	TH1F *hNchCent = new TH1F("hNchCent", "Nch and Centrality Total",kCentBin,nchbd_hist.data());
	hNchCent->Sumw2();
	// 设置横坐标只显示 centbd 的边界值
	/*for (int ibd = 1; ibd <= kCentBin; ibd++) {
		TString label = Form("%.2f-%.2f", centbd[ibd-1], centbd[ibd]);
		hCent->GetXaxis()->SetBinLabel(ibd, label.Data());
	}*/
	char ktName[kKtBin][25];
	for(int i=0;i<kKtBin;i++) sprintf(ktName[i],"%.3f-%.3f",ktbd[i],ktbd[i+1]);


//histograms for distributions
	int nBins = 200;
	TH1F *hx = new TH1F("hx", "x distribution", nBins, -100, 100);
	TH1F *hy = new TH1F("hy", "y distribution", nBins, -100, 100);
	TH1F *hz = new TH1F("hz", "z distribution", nBins, -100, 100);
	TH1F *ht = new TH1F("ht", "time distribution", nBins, 0.0, 100);
	TH1F *hpx = new TH1F("hpx", "px distribution", 100, -2, 2);
	TH1F *hpy = new TH1F("hpy", "py distribution", 100, -2, 2);
	TH1F *hpz = new TH1F("hpz", "pz distribution", 100, -2, 2);
	TH1F *hpt = new TH1F("hpt", "pt distribution", 100, 0, 2);
	TH1F *hpMom = new TH1F("hpMom", "momentum distribution", 100, 0, 2);
	TH1F *henergy = new TH1F("henergy", "energy distribution", 100, 0.0, 2);


	//distributions for particle pairs all events
	TH1F *hrx = new TH1F("hrx", "rx distribution", nBins, -100, 100);
	TH1F *hry = new TH1F("hry", "ry distribution", nBins, -100, 100);
	TH1F *hrz = new TH1F("hrz", "rz distribution", nBins, -100, 100);
	TH1F *hkT = new TH1F("hkT", "kT distribution", 100, 0.0, 2);
	//TH1F *hkTAll = new TH1F("hkTAll", "kT distribution all", 100, 0.0, 2);
	TH1F *hQlcms = new TH1F("hQlcms", "QLCMS", 100, 0.0, 2);
	TH1F *hQlcmsCut = new TH1F("hQlcmsCut", "QLCMS after cut", 100, 0.0, 2);
	TH2D *hQlcms_mT = new TH2D("hQlcms_mT","QLCMS vs sqrtmT", 100, 0.0, 2, 100, 0.0, 2);
	TH2D *hQlcms_mTCut = new TH2D("hQlcms_mT_Cut","QLCMS vs sqrtmT after Cut", 100, 0.0, 2, 100, 0.0, 2);
	TH2D *hQlcms_kT = new TH2D("hQlcms_kT","QLCMS vs kT", 100, 0.0, 2, 100, 0.0, 2);
	//histograms for multiplicity and centrality
	TH1F *hnp = new TH1F("hnp", "multiplicity", 100, 0, 100000);
	TH1F *hnpAll = new TH1F("hnpAll", "npartproj+nparttarg", 500, 0, 5000);//npartproj+nparttarg
	TH1F *hngl = new TH1F("hngl", "ngl", 2000, 0, 15000);
	TH1F *hbim = new TH1F("hbim", "impact parameter", 100, 0, 20);
	TH1F *hbim0 = new TH1F("hbim0", "impact parameter for 0-10%", 100, 0, 20);
	//TH2F *hNchBim = new TH2F("hNchBim", "Nch vs impact parameter", 5000, 0, 5000, 100, 0, 20);//Nch vs b default
	TH2F *hNchBim = new TH2F("hNchBim", "Nch vs impact parameter", 20000, 0, 20000, 300, 0, 30);//Nch vs b

	TH2F *hMulCent = new TH2F("hMulCent", "impact parameter vs multiplicity", 100, 0, 20, 2000, 0, 100000);
	TH2F *hBimNAll = new TH2F("hBimNAll", "bim vs npartAll",100, 0, 20, 500, 0, 5000);
	TH2F *hBimNgl = new TH2F("hBimNgl", "bim vs ngl",100, 0, 20, 2000, 0, 15000);


	//TH1F *correctionHist[kKtBin]= nullptr;
	//kt histograms
	TH1D * kt_histAll[kKtBin] = {};// for all events to get the kTMean, new TH1D("kt_histAll","k_{T} distributions for all events", 100, 0.0, 2);

	for(Int_t jM=0; jM<kKtBin; jM++){
		TString hNamekTAll = Form("hkTAll_Cent%d_kT%d", centbin, jM);
		TString hTitlekTAll = Form("kT distribution for all events Cent%d kT%d", centbin, jM);//rho
		kt_histAll[jM] = new TH1D(hNamekTAll, hTitlekTAll, 100, 0.0, 2);
		kt_histAll[jM]->SetDirectory(0);
	}


	//out file to save the general informations, eventIndex, bim,centbin, ktbin, ktMean(from kt histograms)
	/*std::ofstream outfile(Form("%s/kTBin_results_Cent%d_%d.txt", pathMean.Data(), centbin, fileID));
	outfile << "iMergeIndex  " <<"  bim  "<<"  CentBin  "<<"  kTBin  "<<"  kTMean  "<<"  kTError  "<<"  NPairsMerged"<<endl;
	std::ofstream outfilekT(Form("%s/kTMean_Cent%d_%d.txt",pathMean.Data(), centbin, fileID));
	outfilekT <<"kTBin  "<<"  kTMean  "<<"  kTError  "<<endl;
	std::ofstream outfilekTRec(Form("%s/kTAll_Cent%d_%d.txt",pathMean.Data(), centbin, fileID));*/
	//outfilekT <<"kTBin  "<<"  kTMean  "<<"  kTError  "<<endl;
   int fileID = 1; 

   cout<<"treeIndex default: "<<fileID<<endl;
		//out file to save the general informations, eventIndex, bim,centbin, ktbin, ktMean(from kt histograms)
	std::ofstream outfile(Form("%s/kTBin_results_%d.txt", pathMean.Data(), fileID));
	outfile << "EventIndex  " <<"  bim  "<<"  CentBin  "<<"  kTBin  "<<"  kTMean  "<<"  kTError  "<<"  NPairs"<<endl;
	std::ofstream outfilekT(Form("%s/kTMean_%d.txt",pathMean.Data(), fileID));
	outfilekT <<"kTBin  "<<"  kTMean  "<<"  kTError  "<<endl;
	std::ofstream outfilekTRec(Form("%s/kTAll_%d.txt",pathMean.Data(), fileID));

	int EvtNum = tree.GetEntries();
	cout<<"EvtNum:  "<<EvtNum<<endl;
	if(test) EvtNum = NumRun;
//Loop over events in ttree
	
	int Nevts = Nevents;//check Nevents effect
	//int Nch = 0;

	
	//int MergedEvtNum = (EvtNum + Nevts - 1) / Nevts;
	//cout<<"Running events number in this file = "<<EvtNum<<",  MergedEvtNum = "<<MergedEvtNum<<endl;
for (int iEvt = 0; iEvt < EvtNum; ++iEvt){

	    tree.GetEntry(iEvt);
		//Get the Maximun of the particles number, npMax
		if (np > npM) {
		npM = np; // Update npMax if the current np is larger
		}
		//plot of np and bim
		hnp->Fill(np);
		hbim->Fill(bim);
		hnpAll->Fill(npartproj+nparttarg);
		hngl->Fill(ngl);
		hCent->Fill(bim);

		hMulCent->Fill(bim,np);
		hBimNAll->Fill(bim,npartproj+nparttarg);
		hBimNgl->Fill(bim,ngl);

		//centrality from centbd
		//if(Matyas) {
		//centbin = hCent->FindBin(bim)-1;
		//if(centbin==kCentBin) centbin = -1;
		//cout<<"centrality bin = "<<centbin<<" , bim is: "<<bim<<" , np = "<<np<<endl;

		int Nch = 0;
		for (Int_t j = 0; j < np; j++) {//single particle
			if (isCharged(id[j])) Nch++;
		}
		hNchBim->Fill(Nch, bim);

		hNchCent->Fill(Nch);
		int bin = hNchCent->FindBin(Nch)-1;
		int centbin = kCentBin - 1 - bin;
		//if(centbin==kCentBin) centbin = -1;
		//if(centbin<0) centbin = -1;

		cout<<"centrality bin from Nch = "<<centbin<<" , bim is: "<<bim<<" , Nch = "<<Nch<<endl;
 
  //int fileID = tree.GetTreeNumber(); 

   cout<<"treeIndex: "<<tree.GetTreeNumber()<<endl;


 // 加入该 centbin 的缓冲池
    centbin_buffer[centbin].push_back(iEvt);
    int NPairsMerged[kKtBin] = {0};//for Nevts
	int NPairsMergedAll;//for Nevts all kt bins
	int iMerge;
	TH1D * kt_hist[kKtBin] = {};//kt for different kt for a single or Nevts events in a group
	for(Int_t jIM=0; jIM<kKtBin; jIM++){
		TString hNamekT = Form("hkTEvt%dCent%dkT%d",iMerge,centbin,jIM);
		TString hTitlekT = Form("kT distribution for event%d centrality bin%dkT%d",iMerge, centbin,jIM);//rho
		kt_hist[jIM] = new TH1D(hNamekT, hTitlekT, 100, 0.0, 2);
		kt_hist[jIM]->SetDirectory(0);
    }

	TH1D * rho_hist[kKtBin] = {};
	TH1D * hist_out[kKtBin] = {};
	TH1D * hist_side[kKtBin] = {};
	TH1D * hist_long[kKtBin] = {};

    

    // 检查是否达到了 Nevts 个
    if (centbin_buffer[centbin].size() >= Nevts) {//Merge Nevts, check if this group alreay has Nevts
        iMerge = merge_index[centbin]++;
        std::vector<int> group = centbin_buffer[centbin];
        


//for (int iMerge = 0; iMerge < MergedEvtNum; ++iMerge) {
		//histograms for each merged event
		for(Int_t jIM=0; jIM<kKtBin; jIM++){
			TString hName("hDrho_");
			TString hTitle("Rho distribution for ");
			TString hNameRho = hName + "Evt" + Form("%d",iMerge) + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
			TString hTitleRho = hTitle + "Event" + Form("%d",iMerge)+ "centrality bin " + Form("%d",centbin) + " kt " + Form("%d",jIM);//rho
			TString hNameOut = hName + "Out_" + "Evt" + Form("%d",iMerge)  + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
			TString hTitleOut = hTitle + "Out direction" + "Evt" + Form("%d",iMerge)  + "centrality bin " + Form("%d",centbin) + "kt " + Form("%d",jIM);//rho out
			TString hNameSide = hName + "Side_" + "Evt" + Form("%d",iMerge)  + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
			TString hTitleSide = hTitle + "Side direction" + "Evt" + Form("%d",iMerge)  + "centrality bin "  + Form("%d",centbin) + "kt " + Form("%d",jIM);//rho side
			TString hNameLong = hName + "Long_" + "Evt" + Form("%d",iMerge) + "Cent" + Form("%d",centbin) + "kt" + Form("%d",jIM);
			TString hTitleLong = hTitle + "Long direction " + "Evt" + Form("%d",iMerge) + "centrality bin " + Form("%d",centbin) + "kt " + Form("%d",jIM);//rho long

			rho_hist[jIM] = new TH1D(hNameRho, hTitleRho, nbinsRho, rho_bins);
			hist_out[jIM] = new TH1D(hNameOut, hTitleOut,  nbinsRho, rho_bins);
			hist_side[jIM] = new TH1D(hNameSide, hTitleSide,  nbinsRho, rho_bins);
			hist_long[jIM] = new TH1D(hNameLong, hTitleLong,  nbinsRho, rho_bins);
			rho_hist[jIM]->Sumw2();
			rho_hist[jIM]->SetDirectory(0);
			hist_out[jIM]->SetDirectory(0);
			hist_side[jIM]->SetDirectory(0);
			hist_long[jIM]->SetDirectory(0);


			//kt distributions
			//TString hNamekT = Form("hkTEvt%dCent%dkT%d",iMerge,centbin,jIM);
			//TString hTitlekT = Form("kT distribution for event%d centrality bin%dkT%d",iMerge, centbin,jIM);//rho
			//kt_hist[jIM] = new TH1D(hNamekT, hTitlekT, 100, 0.0, 2);
			//kt_hist[jIM]->SetDirectory(0);
		}




	  for (size_t local_idx = 0; local_idx < group.size(); ++local_idx) {
	    int EvtId = group[local_idx];
	    tree.GetEntry(EvtId);
	    cout<<"===============================Running event "<<local_idx<<" in the current "<<iMerge<<"th  Merged events for centbin "<<centbin<<". The "<<EvtId<<"th event in this file"<<"==============================="<<endl;

		//  if(bim>4.73) continue;//centrality 0-10%

	//Loop over particles in each event
		particlesStr.clear(); //initialize structure for each event
		//int Nch = 0;
		for (Int_t j = 0; j < np; j++) {//single particle
			// Apply the single-particle cuts
			// if (id[j] != 120) continue;          // Only pi+, PID = 120

			//if (isCharged(id[j])) Nch++;
			if (fabs(id[j]) != 120) continue;          //  pi+ and pi-, PID = +-120
			if (ist[j] != 0) continue;           // Only final detected particles

			//calculate energy from mass and p
			// energy[j] = sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j] + mass[j]*mass[j]);
			energy[j] = sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j] + MassPi*MassPi);

			//cout<<"mass: "<<mass[j]<<endl;

			// Calculate pT and eta
			Float_t pt, eta, pMom;
			pt = sqrt(px[j] * px[j] + py[j] * py[j]);
			pMom = sqrt(px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]);
			eta = 0.5 * log((pMom + pz[j]) / (pMom - pz[j])); 

			// pt and eta cuts
			if (pt < PtMin || pt > PtMax) continue;
			if (fabs(eta) > EtaMax) continue;

			// Store valid particles
			Particle p;
			p.id = id[j]; 
			p.ist = ist[j]; 
			p.t = t[j];
			p.x = x[j];
			p.y = y[j];
			p.z = z[j];
			p.energy = energy[j];
			p.px = px[j];
			p.py = py[j];
			p.pz = pz[j];
			p.pt = pt;
			p.eta = eta;
			particlesStr.push_back(p);

			//histograms for space and momentum distributions
			hx->Fill(x[j]);
			hy->Fill(y[j]);
			hz->Fill(z[j]);
			henergy->Fill(energy[j]);
			ht->Fill(t[j]);
			hpx->Fill(px[j]);
			hpy->Fill(py[j]);
			hpz->Fill(pz[j]);
			hpt->Fill(pt);
			hpMom->Fill(pMom);
			//cout<<"px: "<<px[j]<<" py: "<<py[j]<<" pz: "<<pz[j]<<" enrgy: "<<energy[j]<<endl;    
		}//loop single particle

		//hNchBim->Fill(Nch, bim);
		//hNchCent->Fill(Nch);
		//centbin = hNchCent->FindBin(Nch)-1;
		//if(centbin==kCentBin) centbin = -1;
		//cout<<"centrality bin from Nch = "<<centbin<<" , bim is: "<<bim<<" , Nch = "<<Nch<<endl;




        cout<<"EvtInd: "<<iEvt<<" number of saved particles: "<<particlesStr.size()<<endl;
		//outfile<<"--------------EvtInd: "<<iEvt<<" number of saved particles: "<<particlesStr.size()<<"----------------"<<endl;


		//number of pairs for each centbin and kTbin.
		int NPairs[kKtBin] = {0};
		Int_t ktbin;

	// Loop over all particle pairs in the current event
		for (size_t iPar = 0; iPar < particlesStr.size(); iPar++) {//p1
			Particle &p1 = particlesStr[iPar]; 
			// cout<<"==========particle1 ID: "<<p1.id<<", ist = "<<p1.ist<<",  px = "<<p1.px<<",  py = "<<p1.py<<",  pz = "<<p1.pz<<",  pt ="<<p1.pt<<",  eta = "<<p1.eta<<"=========="<<endl;
			for (size_t jPar = iPar + 1; jPar < particlesStr.size(); jPar++) 
			{
				Particle &p2 = particlesStr[jPar];
				if(p1.id != p2.id) continue;//pi+pi+, and pi-pi-
				//calcalute rho for pair
				double rx = p2.x - p1.x;
				double ry = p2.y - p1.y;
				double rz = p2.z - p1.z;
				//K components
				double Kx = (p1.px + p2.px) / 2;
				double Ky = (p1.py + p2.py) / 2;
				double Kz = (p1.pz + p2.pz) / 2;
				double K0 = (p1.energy + p2.energy) / 2;
				double kT2 = Kx * Kx + Ky * Ky;
				double kT = sqrt(kT2);
				double mT = sqrt(kT2 + MassPi * MassPi);
				double sqrt_mT = sqrt(mT);
				//QLCMS<kT cut
				double q_Num = 4*(p1.pz*p2.energy-p2.pz*p1.energy) * (p1.pz*p2.energy-p2.pz*p1.energy);
				double q_den = (p1.energy+p2.energy)*(p1.energy+p2.energy) - (p1.pz+p2.pz)*(p1.pz+p2.pz);
				double q_z_lcms2 = q_Num/q_den;
				double Q_LCMS = sqrt((p1.px-p2.px)*(p1.px-p2.px) + (p1.py-p2.py)*(p1.py-p2.py) + q_z_lcms2);
				hQlcms->Fill(Q_LCMS);
				hQlcms_kT->Fill(kT, Q_LCMS);
				hQlcms_mT->Fill(sqrt_mT,Q_LCMS);
				//Qlcms Cut
				if (Q_LCMS > Qcut_c/sqrt(1000) * sqrt_mT) continue; //c=12 default
				//if (Q_LCMS > kT) continue;
				hQlcmsCut->Fill(Q_LCMS);
				hQlcms_mTCut->Fill(sqrt_mT,Q_LCMS);
				ktbin = hmKt->FindBin(kT)-1;
				if(ktbin==kKtBin) ktbin = -1;
				if(ktbin<0) continue;
				//cout<<"kt = "<<kT<<", ktbin is:  "<<ktbin<<endl;
				kt_hist[ktbin]->Fill(kT);
				kt_histAll[ktbin]->Fill(kT);

				//t, use the deltaT/2
				double t = p2.t - p1.t;
				// Calculate cos(phi) and sin(phi) 
				double cos_phi = Kx / kT;
				double sin_phi = Ky / kT;
				//Compute rhoOut (equation13)
				double rhoOut = rx * cos_phi + ry * sin_phi - kT / (K0 * K0 - Kz * Kz) * (K0 * t - Kz * rz);
				// Compute rhoSide (equation 14)
				double rhoSide = -rx * sin_phi + ry * cos_phi;
				// Compute rhoLong (equation 15)
				double rhoLong = (K0 * rz - Kz * t) / sqrt(K0 * K0 - Kz * Kz);
				// Compute the total rho
				double rho = sqrt(rhoOut * rhoOut + rhoSide * rhoSide + rhoLong * rhoLong);
				// To get a normalized histogram dividing with 4*pi*rho^2 and (number of particle pairs) is needed
				//double weight = 1.0 / (4 * TMath::Pi() * rho * rho) / (particlesStr.size()*(particlesStr.size()-1)/2);
				//double weight = 1.0 / (4 * TMath::Pi() * rho * rho);

				//rho_hist[ktbin]->Fill(rho,weight);
				rho_hist[ktbin]->Fill(rho);
				hist_out[ktbin]->Fill(fabs(rhoOut));//,1.0/EvtNum);
				hist_side[ktbin]->Fill(fabs(rhoSide));//,1.0/EvtNum);
				hist_long[ktbin]->Fill(fabs(rhoLong));//,1.0/EvtNum);
				//hist->Fill(rho,weight);
				NPairs[ktbin]++;
				NPairsMerged[ktbin]++;
				NPairsMergedAll++;
				//histograms for space and momentum distributions
				hrx->Fill(rx);
				hry->Fill(ry);
				hrz->Fill(rz);
				hkT->Fill(kT);
				outfilekTRec<<kT<<endl;

				//cout<<"------particle2 ID: "<<p2.id<<", ist = "<<p2.ist<<",  px = "<<p2.px<<",  py = "<<p2.py<<",  pz = "<<p2.pz<<",  pt ="<<p2.pt<<",  eta = "<<p2.eta<<endl;
				//cout<<"rho =  "<<rho<<",   kt = "<<kT<<",   Q_LCMS = "<<Q_LCMS<<endl;
				//std::cout << "\n";
				//cout<<"rho pair: "<<rho<<"  factor: "<<phaseSpaceFactor<<endl;
			}//p2
		}//p1
}//loop over the merged Nevts=10 events, save as one single event
      centbin_buffer[centbin].clear();
//}//check if the group has Nevts already

//===========================================save the merged Nevts to iMerge th event output================


		//cout<<"iMerge EvtInd: "<<merge_index<<" number of saved particles: "<<particlesStr.size()<<endl;
		//outfile<<"--------------iMergeEvt: "<<iMerge<<" number of saved particles: "<<particlesStr.size()<<"----------------"<<endl;
        cout<<"iMerge EvtInd: "<<iMerge<<" number of saved particles: "<<NPairsMergedAll<<endl;
		outfile<<"--------------iMergeEvt: "<<iMerge<<" number of saved particles: "<<NPairsMergedAll<<"----------------"<<endl;


	//kt histograms for each event analysis, kt mean and error
	//for each kt bin, summarize the info
		for(int iBin=0;iBin<kKtBin;iBin++){//for each kT bin
			cout<<"kT bin: "<<iBin<<"  number of pairs: "<<NPairsMerged[iBin]<<endl;
			//for each kt bin histogram, find the mean of kt
			double weighted_sum = 0.0;
			double total_content = 0.0;
			// double error_sum = 0.0;
			for(int iBinkT=1; iBinkT<=kt_hist[iBin]->GetNbinsX();++iBinkT){
				const double bin_center = kt_hist[iBin]->GetBinCenter(iBinkT);
				const double bin_content = kt_hist[iBin]->GetBinContent(iBinkT);
				weighted_sum += bin_center * bin_content;
				total_content += bin_content; 
			}//get the mean <kt> for each kt bin
			double mean = (total_content > 0) ? weighted_sum / total_content : 0.0;

			//error bar calculated
			double error_sum = 0.0;
			for(int iBinkT=1; iBinkT<=kt_hist[iBin]->GetNbinsX();++iBinkT){
				const double bin_center = kt_hist[iBin]->GetBinCenter(iBinkT);
				const double bin_content = kt_hist[iBin]->GetBinContent(iBinkT);
				// error_sum += pow(bin_content * bin_center - mean, 2);//error bar sum
				error_sum += bin_content * pow(bin_center - mean, 2);  // 正确累加方式:ml-citation{ref="5" data="citationList"}
			}//get the error bar of <kt>
			double std_dev = sqrt(error_sum / total_content);
			//print out into file
			outfile <<iMerge<<"    " <<bim<<"    "<<centbin<<"    "<<iBin<<"    "<<mean<<"    "<<std_dev<<"    "<<NPairsMerged[iBin]<<endl;
		}//for each ktbin 

    //rho_hist Normalization
	/*	for(Int_t jN=0; jN<kKtBin; jN++){
			//Scale by integral
			double integral = rho_hist[jN]->Integral();
			cout<<"kTbin: "<<jN<<" Integral: "<<integral<<endl;
			for (int ih = 1; ih <= rho_hist[jN]->GetNbinsX(); ih++) {
				double bin_content = rho_hist[jN]->GetBinContent(ih);
				double bin_width = rho_hist[jN]->GetBinWidth(ih);
				double bin_error = rho_hist[jN]->GetBinError(ih);
				double bin_center = rho_hist[jN]->GetBinCenter(ih);
				double phasespace = PhaseSpaceFactor(bin_center, 1);
				rho_hist[jN]->SetBinContent(ih, bin_content / bin_width / phasespace);///bin_width
				rho_hist[jN]->SetBinError(ih, bin_error / bin_width / phasespace); // Error normalization / bin_width
			}
			rho_hist[jN]->Scale(1.0 / integral);
		}*/

		// Save the histogram to a file
		//output file for each event
		TFile *outputFile = new TFile(Form("%s/D_rho_Cent%d_File%d_Evt%d.root",path.Data(),centbin, fileID, iMerge), "RECREATE");//fileID_EntryEvtId
		// outputFile->cd();
		//hCent->Write();
		//hbim->Write();
		for(Int_t jds=0; jds<kKtBin; jds++){
		rho_hist[jds]->Write();
		hist_out[jds]->Write();
		hist_side[jds]->Write();
		hist_long[jds]->Write();
		kt_hist[jds]->Write();
		//rho_hist[jds]->Delete();
		}
		// outputFile->Write();
		outputFile->Close();
		delete outputFile;    // 显式释放文件对象:ml-citation{ref="6" data="citationList"}
 //}//loop over event, Nevts=10 merged together
//}//loop over iMerge events
}//check if the group has Nevts already, loop over all Merged group
}//loop over all Evts


		//print out npM; for all events
		std::cout << "Maximum number of particles in any event (npM): " << npM << endl;

//kt histogram analysis for all events, to get the mean and error of kT
	for(int iBinA=0;iBinA<kKtBin;iBinA++){
		double weighted_sumAll = 0.0;
		double total_contentAll = 0.0;
		// double error_sum = 0.0;
		for(int iBinkT=1; iBinkT<=kt_histAll[iBinA]->GetNbinsX();++iBinkT){
			const double bin_center = kt_histAll[iBinA]->GetBinCenter(iBinkT);
			const double bin_content = kt_histAll[iBinA]->GetBinContent(iBinkT);
			weighted_sumAll += bin_center * bin_content;
			total_contentAll += bin_content; 
		}//get the mean <kt> for each kt bin
		double kt_mean = (total_contentAll > 0) ? weighted_sumAll / total_contentAll : 0.0;

		//error bar calculated
		double error_sumAll = 0.0;
		for(int iBinkT=1; iBinkT<=kt_histAll[iBinA]->GetNbinsX();++iBinkT){
			const double bin_center = kt_histAll[iBinA]->GetBinCenter(iBinkT);
			const double bin_content = kt_histAll[iBinA]->GetBinContent(iBinkT);
			// error_sum += pow(bin_content * bin_center - mean, 2);//error bar sum
			error_sumAll += bin_content * pow(bin_center - kt_mean, 2);  // 正确累加方式:ml-citation{ref="5" data="citationList"}
		}//get the error bar of <kt>
		double std_devAll = sqrt(error_sumAll / total_contentAll); 
		//print out into file
		outfilekT <<iBinA<<"    "<<kt_mean<<"    "<<std_devAll<<"    "<<endl;
	}

	outfile.close();
	outfilekT.close();
	outfilekTRec.close();
	std::cout << "Result saved to kTBin_result.txt" << std::endl;


//output root file for all events
    TFile *outputFileAll = new TFile(Form("%s/DisHistos_%d.root",path.Data(),fileID), "RECREATE");//fileID all events histograms
    hCent->Write();
    hNchCent->Write();
    hbim->Write();
    hpt->Write();
    hpMom->Write();
    hkT->Write();
    hQlcms->Write();
    hQlcmsCut->Write();
    hQlcms_mT->Write();
    hQlcms_mTCut->Write();
    hQlcms_kT->Write();
    hx->Write();
    hy->Write();
    hz->Write();
    hrx->Write();
    hry->Write();
    hrz->Write();
    hnp->Write();
    hnpAll->Write();
    hngl->Write();
    hBimNgl->Write();
    hBimNAll->Write();
    hMulCent->Write();
    hNchBim->Write();
    outputFileAll->Close();
}
