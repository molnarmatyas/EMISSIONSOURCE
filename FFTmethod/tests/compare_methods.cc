// Small final-plotter fixtures: no production fits or source files are changed.
#include "../fft_common.h"
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TParameter.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
using namespace std;
namespace fs=std::filesystem;
fs::path directory,program;

void require(bool condition,string message) { if(!condition) throw runtime_error(message); }
void near(double actual,double expected) { require(isfinite(actual) && abs(actual-expected)<1.e-9,"Wrong comparison value/error"); }
string quote(string value) {
  string result="'";
  for(char c:value) result+=c=='\'' ? "'\\''" : string(1,c);
  return result+"'";
}
void run(vector<string> args,int expected=0) {
  string command=quote(program.string());
  for(auto arg:args) command+=" "+quote(arg);
  command+=" > "+quote((directory/"run.log").string())+" 2>&1";
  int status=system(command.c_str());
  if(status<0 || !WIFEXITED(status) || WEXITSTATUS(status)!=expected) {
    ifstream log(directory/"run.log"); cerr << log.rdbuf();
    throw runtime_error("Unexpected comparison exit code");
  }
}
void makeFinal(int method,bool baseline=true) {
  TFile file((directory/(method==0 ? "drho.root" : "fft.root")).c_str(),"RECREATE");
  const char* parameters[]={"alpha","R_out","R_side","R_long","normalization"};
  const double values[]={1.3,4.,5.,6.,.9};
  const double sysDown[]={.02,.2,.4,.6,.01},sysUp[]={.03,.3,.5,.7,.02};
  const double statDown[]={.01,.1,.3,.5,.005},statUp[]={.02,.2,.4,.6,.006};
  for(const char* energy : {"3p0","27"}) for(int q=0; q<5; ++q) {
    TGraphAsymmErrors sys,stat;
    for(int kt=NKT-1; kt>=0; --kt) {
      if(method==1 && q==1 && kt==2) continue; // An omitted FFT point.
      if(q==3 && kt==4) continue; // Ratio/difference do not need Rlong.
      double value=values[q]+(method==1 ? .1 : 0.);
      if(q==1 && kt==1) value=5.+(method==1 ? .1 : 0.); // Keep zero squared difference.
      if(method==0 && q==1 && kt==3) value=0.; // A legacy zero placeholder.
      int i=sys.GetN(); double mt=sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi);
      sys.SetPoint(i,mt,value); sys.SetPointError(i,.01,.01,sysDown[q],sysUp[q]);
      stat.SetPoint(i,mt,value); stat.SetPointError(i,.01,.01,statDown[q],statUp[q]);
    }
    TString parameter=q==4 ? (method==0 ? "N" : "lambda") : parameters[q];
    TString suffix=method==0 ? "" : "_FFT";
    sys.Write(parameter+"_syserr_"+energy+suffix);
    if(baseline) stat.Write(parameter+"_stat_"+energy+suffix);
  }
  if(method==1) TParameter<int>("fft_systematics_complete",1).Write();
  TParameter<int>("input_untouched",42).Write();
}
TGraphAsymmErrors* graph(TFile& file,const char* name) {
  auto* result=file.Get<TGraphAsymmErrors>(name);
  require(result,"Missing graph: "+string(name)); return result;
}
void checkPanels(TFile& file,const char* name,int n,int methods=2) {
  auto* canvas=file.Get<TCanvas>(name); require(canvas,"Missing canvas: "+string(name));
  int count=0; double ymin=0.,ymax=0.;
  for(auto* object : *canvas->GetListOfPrimitives()) if(auto* pad=dynamic_cast<TPad*>(object)) {
    auto* frame=dynamic_cast<TH1*>(pad->GetPrimitive(Form("frame_%d",count)));
    require(frame,"Missing panel frame");
    if(count==0) { ymin=frame->GetMinimum(); ymax=frame->GetMaximum(); }
    else { near(frame->GetMinimum(),ymin); near(frame->GetMaximum(),ymax); }
    int bands=0,curves=0;
    for(auto* link=pad->GetListOfPrimitives()->FirstLink(); link; link=link->Next()) {
      auto* graph=dynamic_cast<TGraphAsymmErrors*>(link->GetObject());
      if(!graph) continue;
      TString option=link->GetOption(),graphName=graph->GetName();
      option.ToUpper();
      if(graphName.EndsWith("_syserr")) { require(option.Contains("3"),"Uncertainties must be connected bands"); ++bands; }
      if(graphName.EndsWith("_points")) { require(option.Contains("L"),"Parameter points must form a continuous curve"); ++curves; }
    }
    require(bands==methods && curves==methods,"Missing method overlay in panel");
    ++count;
  }
  require(count==n,"Wrong number of panels");
}
int main(int argc,char** argv) {
  if(argc!=2) return 1;
  program=fs::absolute(argv[1]);
  char temporary[]="/tmp/emissionsource_compare_XXXXXX"; directory=mkdtemp(temporary);
  try {
    makeFinal(0); makeFinal(1);
    vector<string> args={"--drho",(directory/"drho.root").string(),"--fft",(directory/"fft.root").string(),
      "--energies","27,3p0","--output",(directory/"results/comparison.root").string(),
      "--plot-dir",(directory/"figures").string()};
    auto oldSize=fs::file_size(directory/"drho.root"),fftSize=fs::file_size(directory/"fft.root");
    run(args);
    {
      TFile file((directory/"results/comparison_FFT.root").c_str());
      auto* out=graph(file,"R_out_FFT_kT_3p0_points"); require(out->GetN()==NKT-1,"Lost sparse FFT points");
      near(out->GetPointX(2),kT_center[3]); near(out->GetPointY(2),4.1);
      out=graph(file,"R_out_Drho_kT_3p0_points"); require(out->GetN()==NKT-1,"Kept legacy zero placeholder");
      auto* ratio=graph(file,"Rout_over_Rside_Drho_kT_3p0_syserr");
      near(ratio->GetPointY(0),.8);
      near(ratio->GetErrorYlow(0),sqrt(.008)); near(ratio->GetErrorYhigh(0),sqrt(.007696));
      auto* difference=graph(file,"Rout2_minus_Rside2_Drho_kT_3p0_syserr");
      near(difference->GetPointY(0),-9.); near(difference->GetPointY(1),0.);
      near(difference->GetErrorYlow(0),hypot(8.*.2,10.*.5));
      near(difference->GetErrorYhigh(0),hypot(8.*.3,10.*.4));
      auto* rms=graph(file,"R_RMS_Drho_kT_3p0_points");
      double r=sqrt(77./3.); near(rms->GetPointY(0),r);
      near(rms->GetErrorYlow(0),sqrt(16.*.01+25.*.09+36.*.25)/(3.*r));
      require(rms->GetN()==NKT-2 && ratio->GetN()==NKT-1,"Wrong derived radius requirements");
      auto* energy=graph(file,"R_out_FFT_sNN_kt0_points");
      near(energy->GetPointX(0),3.); near(energy->GetPointX(1),27.);
      require(out->GetLineStyle()==1 && energy->GetLineStyle()==2,"Methods need different line styles");
      require(ratio->GetFillStyle()==3004 && graph(file,"Rout_over_Rside_FFT_kT_3p0_syserr")->GetFillStyle()==3005,"Methods need different fill patterns");
      checkPanels(file,"kT_alpha_3p0_comparison_FFT",1);
      checkPanels(file,"kT_radii_3p0_comparison_FFT",3);
      require(!file.Get("kT_R_out_3p0_comparison_FFT"),"Kept the separate two-panel radius plots");
      checkPanels(file,"sNN_radii_kt0_comparison_FFT",3);
      checkPanels(file,"sNN_R_RMS_kt0_comparison_FFT",1);
      checkPanels(file,"sNN_normalization_kt0_comparison_FFT",2,1);
      checkPanels(file,"sNN_alpha_kt4_comparison_FFT",1);
      require(!file.Get("sNN_radii_kt1_comparison_FFT") && !file.Get("sNN_radii_kt9_comparison_FFT"),"Energy plots must default to kT bins 0 and 4 only");
      require(fs::file_size(directory/"figures/3p0/kT_alpha_3p0_comparison_FFT.pdf")>1000,"Missing kT PDF");
      require(fs::file_size(directory/"figures/3p0/kT_radii_3p0_comparison_FFT.pdf")>1000,"Missing three-panel kT radii PDF");
      require(!fs::exists(directory/"figures/3p0/kT_R_out_3p0_comparison_FFT.pdf"),"Wrote an obsolete radius figure");
      require(fs::file_size(directory/"figures/sNN_radii_kt0_comparison_FFT.png")>1000,"Missing energy PNG");
    }
    require(fs::file_size(directory/"drho.root")==oldSize && fs::file_size(directory/"fft.root")==fftSize,"Changed comparison inputs");
    {
      ifstream csv(directory/"results/comparison_FFT.csv"); string line; int rows=0;
      while(getline(csv,line)) { require(count(line.begin(),line.end(),',')==12,"Malformed CSV"); ++rows; }
      require(rows>200,"Missing CSV data");
    }
    auto custom=args;
    custom.insert(custom.end(),{"--kt","1","--no-plots","--output",(directory/"results/custom.root").string()});
    run(custom);
    {
      TFile file((directory/"results/custom_FFT.root").c_str());
      require(file.Get("alpha_FFT_sNN_kt1_points") && !file.Get("alpha_FFT_sNN_kt0_points"),"Ignored explicit energy-plot kT override");
    }
    makeFinal(0,false); // Older final-plotter output with systematic errors only.
    args.push_back("--no-plots"); run(args);
    {
      TFile file((directory/"results/comparison_FFT.root").c_str());
      require(file.Get<TParameter<int>>("comparison_missing_drho_baseline")->GetVal()==1,"Missing baseline errors not marked");
      auto* g=graph(file,"R_out_Drho_kT_3p0_points"); near(g->GetErrorYhigh(0),0.);
      ifstream csv(directory/"results/comparison_FFT.csv"); string content((istreambuf_iterator<char>(csv)),{});
      require(content.find(",0,nan,nan,")!=string::npos,"Missing baseline errors presented as measured zeros");
    }
    run({"--drho",(directory/"missing.root").string(),"--fft",(directory/"fft.root").string(),"--no-plots"},1);
    run({"--drho",(directory/"drho.root").string(),"--fft",(directory/"fft.root").string(),"--energies","5p2","--no-plots"},1);
    run({"--kt","10"},1);
    cout << "Method comparison checks passed in " << directory << endl;
  } catch(const exception& error) { cerr << error.what() << "\nFixtures: " << directory << endl; return 1; }
}
