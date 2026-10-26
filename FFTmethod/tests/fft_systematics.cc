// Deterministic ROOT fixtures for the final calculator, separate from real analysis data.
#include "../fft_common.h"
#include <TGraphAsymmErrors.h>
#include <TAxis.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
namespace fs=std::filesystem;
fs::path directory, program;
const char* names[]={"lambda","R_out","R_side","R_long","alpha"};
const double truth[]={0.7,4.,5.,6.,1.3}, shifts[]={0.,.1,-.2,.3,.5,-.4,-.6};
void require(bool value, const std::string& message) { if(!value) throw std::runtime_error(message); }
void near(double actual,double expected) { require(std::isfinite(actual) && std::abs(actual-expected)<1.e-10,"Unexpected systematic value"); }
std::string quote(std::string value) {
  std::string result="'";
  for(char c:value) result+=c=='\'' ? "'\\''" : std::string(1,c);
  return result+"'";
}
void run(std::vector<std::string> args,int expected=0) {
  std::string command="cd "+quote(directory.string())+" && "+quote(program.string());
  for(auto arg:args) command+=" "+quote(arg);
  command+=" > "+quote((directory/"run.log").string())+" 2>&1";
  int status=std::system(command.c_str());
  if(status<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=expected) {
    std::ifstream log(directory/"run.log"); std::cerr<<log.rdbuf();
    throw std::runtime_error("Unexpected calculator exit code");
  }
}
std::string filename(const char* energy,int variation) {
  std::string base=std::string("UrQMD_3d_source_0-10cent_all_")+energy;
  const char* suffix[]={"_Cq_fitresults","_strictqLCMS_Cq_fitresults","_looseqLCMS_Cq_fitresults",
    "_Cq_fitresults_narrowqFit","_Cq_fitresults_wideqFit","_Cq_coarseFFT_fitresults","_Cq_fineFFT_fitresults"};
  return base+suffix[variation]+"_FFT.root";
}
void makeFit(const char* energy,int variation,double energyShift) {
  TFile file((directory/"analysed"/filename(energy,variation)).c_str(),"RECREATE");
  const int cuts[]={0,1,2,0,0,0,0}, ranges[]={0,0,0,1,2,0,0}, grids[]={0,0,0,0,0,1,2};
  TNamed("source_energy",energy).Write(); TNamed("source_centrality","0-10").Write();
  TParameter<Long64_t>("source_events",10000).Write();
  TParameter<int>("source_qlcms_cut",cuts[variation]).Write();
  TParameter<int>("cqfit_range_variation",ranges[variation]).Write();
  TParameter<int>("fft_binning_variation",grids[variation]).Write();
  TParameter<double>("fft_source_half_range",fftSourceHalfRange).Write();
  TParameter<int>("fft_fit_complete",1).Write();
  int kt,bins=fftSourceBinsSyst[grids[variation]]; bool valid=true; double qmax;
  Long64_t firstEvent=0,events=-1;
  TTree fits("fit_results_FFT","fixture"), summary("fft_sources","fixture");
  fits.Branch("kt",&kt); fits.Branch("valid",&valid); fits.Branch("qmax",&qmax);
  summary.Branch("kt",&kt); summary.Branch("bins",&bins);
  summary.Branch("first_event",&firstEvent); summary.Branch("events",&events);
  TGraphAsymmErrors graphs[5];
  for(kt=NKT-1; kt>=0; --kt) {
    qmax=fftQmaxFit[kt]*fftQmaxFitScale[ranges[variation]];
    fits.Fill(); summary.Fill();
    for(int p=0; p<5; ++p) {
      int point=graphs[p].GetN();
      graphs[p].SetPoint(point,kT_center[kt],truth[p]+energyShift+.01*kt+shifts[variation]);
      graphs[p].SetPointError(point,kT_center[kt]-ktbins[kt],ktbins[kt+1]-kT_center[kt],.02,.03);
    }
  }
  fits.Write(); summary.Write();
  for(int p=0; p<5; ++p) graphs[p].Write(Form("%s_vs_KT_FFT",names[p]));
}
int main(int argc,char** argv) {
  if(argc!=2) return 1;
  program=fs::absolute(argv[1]);
  char temporary[]="/tmp/emissionsource_fft_systematics_XXXXXX";
  directory=mkdtemp(temporary);
  try {
    fs::create_directories(directory/"analysed");
    for(int v=0; v<7; ++v) { makeFit("3p0",v,0.); makeFit("27",v,.1); }
    {
      TFile old((directory/"drho.root").c_str(),"RECREATE");
      TParameter<int>("untouched",42).Write();
      for(int p=1; p<5; ++p) {
        TGraphAsymmErrors graph;
        for(int kt=0; kt<NKT; ++kt) {
          graph.SetPoint(kt,std::sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi),kt==2 ? 0. : truth[p]);
          graph.SetPointError(kt,.01,.01,.02,.02);
        }
        graph.Write(Form("%s_syserr_3p0",names[p]));
      }
    }
    run({"--energies","27,3p0","--output","results/summary.root","--plot-dir","figures/final","--compare-drho","drho.root"});
    {
      TFile result((directory/"results/summary_FFT.root").c_str());
      require(result.Get<TParameter<int>>("fft_systematics_complete")->GetVal()==1,"Complete systematics marked incomplete");
      require(result.Get<TParameter<int>>("fft_systematics_points")->GetVal()==2*NKT,"Energy/kT points merged or lost");
      for(int p=0; p<5; ++p) {
        auto* stat=result.Get<TGraphAsymmErrors>(Form("%s_stat_3p0_FFT",names[p]));
        auto* syst=result.Get<TGraphAsymmErrors>(Form("%s_syserr_3p0_FFT",names[p]));
        auto* total=result.Get<TGraphAsymmErrors>(Form("%s_total_3p0_FFT",names[p]));
        require(stat && syst && total && stat->GetN()==NKT,"Missing final graph");
        require(std::string(stat->GetName())==std::string(names[p])+"_stat_3p0_FFT","Graph key and object name differ");
        require(std::string(stat->GetXaxis()->GetTitle())=="m_{T} [GeV/c^{2}]","Missing mT axis units");
        for(int kt=0; kt<NKT; ++kt) {
          near(stat->GetPointX(kt),std::sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi));
          near(stat->GetPointY(kt),truth[p]+.01*kt); near(stat->GetErrorYlow(kt),.02); near(stat->GetErrorYhigh(kt),.03);
          near(syst->GetErrorYlow(kt),std::sqrt(.29)); near(syst->GetErrorYhigh(kt),std::sqrt(.17));
          near(total->GetErrorYlow(kt),std::sqrt(.29+.02*.02)); near(total->GetErrorYhigh(kt),std::sqrt(.17+.03*.03));
        }
        auto* qcut=result.Get<TGraphAsymmErrors>(Form("%s_qLCMS_3p0_FFT",names[p]));
        auto* range=result.Get<TGraphAsymmErrors>(Form("%s_qFit_3p0_FFT",names[p]));
        auto* grid=result.Get<TGraphAsymmErrors>(Form("%s_binning_3p0_FFT",names[p]));
        near(qcut->GetErrorYlow(0),.2); near(qcut->GetErrorYhigh(0),.1);
        near(range->GetErrorYlow(0),0.); near(range->GetErrorYhigh(0),.4);
        near(grid->GetErrorYlow(0),.5); near(grid->GetErrorYhigh(0),0.);
        auto* energy=result.Get<TGraphAsymmErrors>(Form("%s_sNN_stat_kt0_FFT",names[p]));
        require(energy && energy->GetN()==2,"Missing energy plot graph");
        require(std::string(energy->GetName())==std::string(names[p])+"_sNN_stat_kt0_FFT","Energy graph object name differs from key");
        near(energy->GetPointX(0),3.); near(energy->GetPointX(1),27.); near(energy->GetPointY(1),truth[p]+.1);
      }
      require(fs::file_size(directory/"figures/final/3p0/mT_alpha_FFT.png")>1000,"Missing final mT plot");
      require(fs::file_size(directory/"figures/final/sNN_R_out_kt0_FFT.pdf")>1000,"Missing energy plot");
      require(fs::file_size(directory/"figures/final/sNN_R_out_kt4_FFT.pdf")>1000,"Missing second selected energy plot");
      require(!fs::exists(directory/"figures/final/sNN_R_out_kt1_FFT.pdf") &&
              !fs::exists(directory/"figures/final/sNN_alpha_kt9_FFT.png"),"Wrote energy figures for unselected kT bins");
      require(result.Get("alpha_sNN_stat_kt9_FFT"),"Plot selection must retain numerical results for every kT bin");
      require(fs::file_size(directory/"figures/final/3p0/contributions_alpha_FFT.png")>1000,"Missing source-contribution plot");
      require(fs::exists(directory/"figures/final/3p0/comparison_Drho_alpha_FFT.pdf"),"Missing optional D(rho) comparison");
    }
    {
      std::ifstream csv(directory/"results/summary_breakdown_FFT.csv");
      std::string line; int rows=0;
      while(std::getline(csv,line)) { require(std::count(line.begin(),line.end(),',')==19,"Malformed breakdown row"); ++rows; }
      require(rows==1+2*NKT*5,"Breakdown lost rows");
    }
    {
      TFile result((directory/"results/summary_FFT.root").c_str());
      std::ifstream text(directory/"results/summary_graphs_FFT.txt");
      std::string name,header,line; int blocks=0;
      while(std::getline(text,name)) {
        if(name.empty()) continue;
        std::getline(text,header);
        auto* graph=result.Get<TGraphAsymmErrors>(name.c_str());
        require(graph && header.find(" y err_up err_dn")!=std::string::npos,"Malformed graph text block");
        for(int i=0; i<graph->GetN(); ++i) {
          std::getline(text,line); std::istringstream values(line); double x,y,up,down;
          require(bool(values>>x>>y>>up>>down),"Malformed graph values");
          near(x,graph->GetPointX(i)); near(y,graph->GetPointY(i));
          near(up,graph->GetErrorYhigh(i)); near(down,graph->GetErrorYlow(i));
        }
        std::getline(text,line); require(line.empty(),"Missing graph block separator"); ++blocks;
      }
      require(blocks==2*6*5+2*NKT*5,"Missing graph text tables");
    }
    {
      std::ifstream csv(directory/"results/summary_percent_breakdown_FFT.csv");
      std::string line; std::getline(csv,line); int rows=0;
      const std::vector<std::string> parameters={"lambda","R_out","R_side","R_long","alpha"};
      const std::vector<std::string> sources={"qLCMS","qFit","binning","total"};
      const double down[]={.2,0.,.5,std::sqrt(.29)},up[]={.1,.4,0.,std::sqrt(.17)};
      while(std::getline(csv,line)) {
        std::istringstream record(line); std::vector<std::string> fields; std::string field;
        while(std::getline(record,field,',')) fields.push_back(field);
        require(fields.size()==7 && std::stoi(fields[4])==NKT,"Malformed percent summary");
        int p=std::find(parameters.begin(),parameters.end(),fields[2])-parameters.begin();
        int s=std::find(sources.begin(),sources.end(),fields[3])-sources.begin();
        require(p<5 && s<4,"Unknown summary parameter/source");
        double expectedDown=0.,expectedUp=0.;
        for(int kt=0; kt<NKT; ++kt) {
          double value=truth[p]+(fields[0]=="27" ? .1 : 0.)+.01*kt;
          expectedDown+=100.*down[s]/value/NKT; expectedUp+=100.*up[s]/value/NKT;
        }
        near(std::stod(fields[5]),expectedDown); near(std::stod(fields[6]),expectedUp); ++rows;
      }
      require(rows==2*5*4,"Missing percent summary rows");
      std::ifstream latex(directory/"results/summary_percent_table_FFT.tex");
      std::string table((std::istreambuf_iterator<char>(latex)),std::istreambuf_iterator<char>());
      require(table.find("\\begin{table}")!=std::string::npos && table.find("FFT binning")!=std::string::npos && table.find("nan")==std::string::npos,"Malformed LaTeX summary");
    }
    {
      TFile varied((directory/"analysed"/filename("3p0",4)).c_str(),"UPDATE");
      auto* graph=varied.Get<TGraphAsymmErrors>("alpha_vs_KT_FFT");
      graph->RemovePoint(3); graph->Write("alpha_vs_KT_FFT",TObject::kOverwrite); // kT 6 missing, reverse order
    }
    run({"--energies","3p0,27","--output","results/sparse_FFT.root","--no-plots"},2);
    {
      TFile result((directory/"results/sparse_FFT.root").c_str());
      require(result.Get<TParameter<int>>("fft_systematics_complete")->GetVal()==0,"Incomplete systematics hidden");
      auto* graph=result.Get<TGraphAsymmErrors>("alpha_stat_3p0_FFT");
      require(graph->GetN()==NKT-1,"Missing variation became a zero error");
      for(int i=0; i<graph->GetN(); ++i) require(std::abs(graph->GetPointX(i)-std::sqrt(kT_center[6]*kT_center[6]+Mass2_pi))>1.e-8,"Sparse point mismatched by index");
    }
    {
      TFile varied((directory/"analysed"/filename("27",1)).c_str(),"UPDATE");
      TParameter<Long64_t>("source_events",9999).Write("source_events",TObject::kOverwrite);
    }
    run({"--energies","3p0,27","--output","results/events_FFT.root","--no-plots"},2);
    {
      TFile result((directory/"results/events_FFT.root").c_str());
      require(result.Get<TGraphAsymmErrors>("alpha_stat_27_FFT")->GetN()==0,"Different event samples were combined");
      TFile old((directory/"drho.root").c_str());
      require(old.Get<TParameter<int>>("untouched")->GetVal()==42,"D(rho) comparison file changed");
    }
    run({"--energies","invalid","--no-plots"},1);
    std::cout << "PASS FFT directional systematics, quadrature, sparse kT matching, metadata, tables and final figures\nFixtures: " << directory << '\n';
    return 0;
  } catch(const std::exception& error) { std::cerr << "FAIL: " << error.what() << "\nFixtures: " << directory << '\n'; return 1; }
}
