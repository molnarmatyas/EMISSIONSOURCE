// Seven fits per energy: baseline, qLCMS x2, fit range x2, source binning x2.
#include "fft_common.h"
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>
using namespace std;

const int NPAR = 5, NVAR = 7, NSYST = 3;
const char* parameters[NPAR] = {"lambda", "R_out", "R_side", "R_long", "alpha"};
const char* titles[NPAR] = {"#lambda", "R_{out} [fm]", "R_{side} [fm]", "R_{long} [fm]", "#alpha"};
const char* sources[NSYST] = {"qLCMS", "qFit", "binning"};
// Same default convention as calc_and_plot_syserr.cpp for D(rho).
bool use_larger_same_side_shift = false;
bool makePlots = true;

struct FitData {
  bool valid[NKT] = {};
  double value[NKT][NPAR] = {}, low[NKT][NPAR] = {}, high[NKT][NPAR] = {};
  Long64_t events = -1;
  TString centrality;
};

double energyValue(TString energy) { energy.ReplaceAll("p", "."); return energy.Atof(); }

bool readFit(TString filename, TString energy, int cut, int range, int grid, FitData& data)
{
  TFile input(filename, "READ");
  if(input.IsZombie()) { cerr << "Missing fit file: " << filename << endl; return false; }
  auto* complete = input.Get<TParameter<int>>("fft_fit_complete");
  auto* cutInfo = input.Get<TParameter<int>>("source_qlcms_cut");
  auto* rangeInfo = input.Get<TParameter<int>>("cqfit_range_variation");
  auto* gridInfo = input.Get<TParameter<int>>("fft_binning_variation");
  auto* halfRange = input.Get<TParameter<double>>("fft_source_half_range");
  auto* energyInfo = input.Get<TNamed>("source_energy");
  auto* eventInfo = input.Get<TParameter<Long64_t>>("source_events");
  auto* centrality = input.Get<TNamed>("source_centrality");
  auto* fits = input.Get<TTree>("fit_results_FFT");
  auto* summary = input.Get<TTree>("fft_sources");
  if(!complete || complete->GetVal()!=1 || !cutInfo || cutInfo->GetVal()!=cut ||
     !rangeInfo || rangeInfo->GetVal()!=range || !gridInfo || gridInfo->GetVal()!=grid ||
     !halfRange || abs(halfRange->GetVal()-fftSourceHalfRange)>1.e-8 ||
     !energyInfo || energy!=energyInfo->GetTitle() || !eventInfo || !centrality || !fits || !summary) {
    cerr << "Incomplete or inconsistent settings in " << filename << "; regenerate with the current FFT programs.\n";
    return false;
  }
  data.events=eventInfo->GetVal(); data.centrality=centrality->GetTitle();
  int kt, bins, binsByKT[NKT] = {};
  Long64_t firstEvent, events;
  summary->SetBranchAddress("kt", &kt); summary->SetBranchAddress("bins", &bins);
  summary->SetBranchAddress("first_event", &firstEvent); summary->SetBranchAddress("events", &events);
  for(Long64_t row=0; row<summary->GetEntries(); ++row) {
    summary->GetEntry(row);
    if(firstEvent!=0 || events!=-1) { cerr << "Systematics expect all events: " << filename << endl; return false; }
    if(kt>=0 && kt<NKT) binsByKT[kt]=bins;
  }
  summary->ResetBranchAddresses();
  bool valid; double qmax;
  fits->SetBranchAddress("kt", &kt); fits->SetBranchAddress("valid", &valid); fits->SetBranchAddress("qmax", &qmax);
  for(Long64_t row=0; row<fits->GetEntries(); ++row) {
    fits->GetEntry(row);
    if(kt>=0 && kt<NKT) data.valid[kt]=valid && binsByKT[kt]==fftSourceBinsSyst[grid] &&
      abs(qmax-fftQmaxFit[kt]*fftQmaxFitScale[range])<1.e-8;
  }
  fits->ResetBranchAddresses();
  for(int p=0; p<NPAR; ++p) {
    auto* graph=input.Get<TGraphAsymmErrors>(Form("%s_vs_KT_FFT",parameters[p]));
    for(kt=0; kt<NKT; ++kt) {
      int point=-1;
      if(graph) for(int i=0; i<graph->GetN(); ++i)
        if(abs(graph->GetPointX(i)-kT_center[kt])<1.e-8) point=i;
      if(point<0) { data.valid[kt]=false; continue; }
      data.value[kt][p]=graph->GetPointY(point);
      data.low[kt][p]=graph->GetErrorYlow(point); data.high[kt][p]=graph->GetErrorYhigh(point);
      if(!isfinite(data.value[kt][p]) || !isfinite(data.low[kt][p]) || !isfinite(data.high[kt][p]) ||
         data.low[kt][p]<0 || data.high[kt][p]<0) data.valid[kt]=false;
    }
  }
  return true;
}

void systematicShift(double center, double first, double second, double& down, double& up)
{
  double a=first-center, b=second-center;
  down=up=0.;
  if((a>=0 && b>=0) || (a<0 && b<0)) {
    double shift=use_larger_same_side_shift ? (abs(a)>=abs(b) ? a : b) : 0.5*(a+b);
    if(shift>=0) up=shift; else down=-shift;
  } else {
    up=max(a,b); down=-min(a,b);
  }
}

// Statistical bars and systematic boxes. Used for per-energy and overlay plots.
void drawPlot(const vector<TGraphAsymmErrors*>& stat, const vector<TGraphAsymmErrors*>& syst,
              const vector<TString>& labels, TString heading, const char* xTitle, const char* yTitle,
              TString filename, bool logx=false, bool drawPoints=true)
{
  if(!makePlots) return;
  double xmin=1.e30, xmax=-1.e30, ymin=1.e30, ymax=-1.e30;
  for(size_t i=0; i<stat.size(); ++i) for(int j=0; j<stat[i]->GetN(); ++j) {
    double x=stat[i]->GetPointX(j), y=stat[i]->GetPointY(j);
    xmin=min(xmin,x-syst[i]->GetErrorXlow(j)); xmax=max(xmax,x+syst[i]->GetErrorXhigh(j));
    ymin=min(ymin,y-max(stat[i]->GetErrorYlow(j),syst[i]->GetErrorYlow(j)));
    ymax=max(ymax,y+max(stat[i]->GetErrorYhigh(j),syst[i]->GetErrorYhigh(j)));
  }
  if(xmin>xmax) return;
  double dx=max(0.01,0.08*(xmax-xmin)), dy=max(0.02,0.18*(ymax-ymin));
  TCanvas canvas("final_plot_FFT", "", 1000,750);
  canvas.SetLeftMargin(0.13); canvas.SetRightMargin(0.035); canvas.SetBottomMargin(0.13); canvas.SetTopMargin(0.08);
  canvas.SetTicks(1,1); canvas.SetLogx(logx);
  auto* frame=canvas.DrawFrame(logx ? xmin/1.15 : xmin-dx,ymin-dy,logx ? xmax*1.15 : xmax+dx,ymax+2.*dy);
  frame->SetTitle(Form(";%s;%s",xTitle,yTitle));
  for(TAxis* axis : {frame->GetXaxis(),frame->GetYaxis()}) {
    axis->SetTitleSize(0.045); axis->SetLabelSize(0.035); axis->SetTitleOffset(1.25);
  }
  const int colors[]={kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kViolet+1,kTeal+3,kPink+7};
  int columns=labels.size()>5 ? 2 : 1;
  double height=0.048*((labels.size()+columns-1)/columns);
  TLegend legend(0.58,max(0.40,0.88-height),0.95,0.88);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.028);
  legend.SetNColumns(columns);
  for(size_t i=0; i<stat.size(); ++i) {
    if(stat[i]->GetN()==0) continue;
    int color=colors[i%10];
    syst[i]->SetFillColorAlpha(color,0.18); syst[i]->SetLineColor(color);
    syst[i]->Draw("2 SAME");
    stat[i]->SetMarkerStyle(20+i%10); stat[i]->SetMarkerSize(0.8); stat[i]->SetMarkerColor(color); stat[i]->SetLineColor(color);
    if(drawPoints) stat[i]->Draw("PZ SAME");
    legend.AddEntry(drawPoints ? stat[i] : syst[i],labels[i],drawPoints ? "p" : "f");
  }
  legend.Draw();
  TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.032);
  text.DrawLatex(0.14,0.95,heading);
  text.SetTextSize(0.024);
  text.DrawLatex(0.15,0.85,drawPoints ? "Bars: statistical; boxes: systematic" : "Signed relative systematic shifts");
  canvas.SaveAs(filename+".png"); canvas.SaveAs(filename+".pdf");
}

int main(int argc, char** argv)
{
  TString inputDir="analysed", outputName="FFTmethod/results/syserr_results_FFT.root", plotDir="figs/syserr/FFT", drhoName="";
  vector<TString> selected;
  for(int i=1; i<argc; ++i) {
    string option=argv[i];
    if(option=="--no-plots") { makePlots=false; continue; }
    if(option=="--help") {
      cout << "Usage: calc_and_plot_syserr_FFT.exe [--input-dir DIR] [--energies 3p0,27]\n"
           << "  [--output file_FFT.root] [--plot-dir DIR] [--no-plots] [--compare-drho syserr_results.root]\n";
      return 0;
    }
    if(i+1==argc) { cerr << "Missing value for " << option << endl; return 1; }
    TString value=argv[++i];
    if(option=="--input-dir") inputDir=value;
    else if(option=="--output") outputName=fftOutputName(value,"");
    else if(option=="--plot-dir") plotDir=value;
    else if(option=="--compare-drho") drhoName=value;
    else if(option=="--energies") {
      istringstream list(value.Data()); string item;
      while(getline(list,item,',')) if(!item.empty() && find(selected.begin(),selected.end(),TString(item))==selected.end()) selected.push_back(item);
    } else { cerr << "Unknown option: " << option << endl; return 1; }
  }
  if(selected.empty()) for(const char* energy : energies)
    if(!gSystem->AccessPathName(Form("%s/UrQMD_3d_source_0-10cent_all_%s_Cq_fitresults_FFT.root",inputDir.Data(),energy))) selected.push_back(energy);
  if(selected.empty()) { cerr << "No baseline FFT fits found.\n"; return 1; }
  for(TString energy : selected) if(energyValue(energy)<=0) { cerr << "Invalid energy: " << energy << endl; return 1; }
  sort(selected.begin(),selected.end(),[](TString a,TString b){return energyValue(a)<energyValue(b);});
  for(TString folder : {TString(gSystem->DirName(outputName)),plotDir})
    if(gSystem->AccessPathName(folder) && gSystem->mkdir(folder,true)!=0) { cerr << "Cannot create " << folder << endl; return 1; }
  gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TFile output(outputName,"RECREATE");
  if(output.IsZombie()) return 1;
  TString tableName=outputName; tableName.Resize(tableName.Length()-9); tableName+="_breakdown_FFT.csv";
  ofstream csv(tableName.Data());
  if(!csv) { cerr << "Cannot create " << tableName << endl; return 1; }
  csv << setprecision(12) << "energy,kt,kt_low,kt_high,mT,parameter,value,stat_down,stat_up,qLCMS_down,qLCMS_up,qFit_down,qFit_up,binning_down,binning_up,sys_down,sys_up,total_down,total_up,status\n";
  TFile* drho=nullptr;
  if(drhoName!="") {
    drho=TFile::Open(drhoName,"READ");
    if(!drho || drho->IsZombie()) { cerr << "Cannot read D(rho) comparison file.\n"; return 1; }
  }
  vector<TGraphAsymmErrors*> allStat[NPAR], allSyst[NPAR];
  vector<TString> energyLabels;
  bool incomplete=false;
  int accepted=0;
  for(TString energy : selected) {
    TString base=inputDir+"/UrQMD_3d_source_0-10cent_all_"+energy;
    TString nominal=base+"_Cq_FFT.root";
    TString files[NVAR]={fftOutputName(nominal,"_fitresults"),
      fftOutputName(base+"_strictqLCMS_Cq_FFT.root","_fitresults"),
      fftOutputName(base+"_looseqLCMS_Cq_FFT.root","_fitresults"),
      fftOutputName(nominal,(TString("_fitresults")+fftFitRangeSuffix[1]).Data()),
      fftOutputName(nominal,(TString("_fitresults")+fftFitRangeSuffix[2]).Data()),
      fftOutputName(base+"_Cq"+fftBinningSuffix[1]+"_FFT.root","_fitresults"),
      fftOutputName(base+"_Cq"+fftBinningSuffix[2]+"_FFT.root","_fitresults")};
    const int cuts[]={0,1,2,0,0,0,0}, ranges[]={0,0,0,1,2,0,0}, grids[]={0,0,0,0,0,1,2};
    FitData data[NVAR]; bool loaded[NVAR];
    for(int v=0; v<NVAR; ++v) loaded[v]=readFit(files[v],energy,cuts[v],ranges[v],grids[v],data[v]);
    for(int v=1; v<NVAR; ++v) if(loaded[v] && loaded[0] &&
        (data[v].events!=data[0].events || data[v].centrality!=data[0].centrality)) {
      cerr << "Event count/centrality differs in " << files[v] << endl; loaded[v]=false;
    }
    TGraphAsymmErrors stat[NPAR], syst[NPAR], total[NPAR], contribution[NSYST][NPAR], relative[NSYST][NPAR];
    for(int kt=0; kt<NKT; ++kt) {
      bool ready=true;
      for(int v=0; v<NVAR; ++v) if(!loaded[v] || !data[v].valid[kt]) ready=false;
      double mt=sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi);
      double xlow=mt-sqrt(ktbins[kt]*ktbins[kt]+Mass2_pi), xhigh=sqrt(ktbins[kt+1]*ktbins[kt+1]+Mass2_pi)-mt;
      if(!ready) { incomplete=true; cerr << energy << ", kT bin " << kt << ": missing/invalid systematic fit; omitted.\n"; }
      else ++accepted;
      for(int p=0; p<NPAR; ++p) {
        double numbers[13]; fill(numbers,numbers+13,numeric_limits<double>::quiet_NaN());
        if(ready) {
          double value=data[0].value[kt][p], down2=0.,up2=0.;
          int point=stat[p].GetN();
          for(int s=0; s<NSYST; ++s) {
            double down,up;
            systematicShift(value,data[1+2*s].value[kt][p],data[2+2*s].value[kt][p],down,up);
            down2+=down*down; up2+=up*up;
            contribution[s][p].SetPoint(point,mt,value); contribution[s][p].SetPointError(point,xlow,xhigh,down,up);
            if(abs(value)>1.e-12) {
              int relativePoint=relative[s][p].GetN();
              relative[s][p].SetPoint(relativePoint,mt,0.);
              relative[s][p].SetPointError(relativePoint,xlow,xhigh,100.*down/abs(value),100.*up/abs(value));
            }
            numbers[3+2*s]=down; numbers[4+2*s]=up;
          }
          double down=sqrt(down2),up=sqrt(up2), low=data[0].low[kt][p],high=data[0].high[kt][p];
          stat[p].SetPoint(point,mt,value); stat[p].SetPointError(point,0.,0.,low,high);
          syst[p].SetPoint(point,mt,value); syst[p].SetPointError(point,xlow,xhigh,down,up);
          total[p].SetPoint(point,mt,value); total[p].SetPointError(point,xlow,xhigh,hypot(low,down),hypot(high,up));
          numbers[0]=value; numbers[1]=low; numbers[2]=high; numbers[9]=down; numbers[10]=up;
          numbers[11]=hypot(low,down); numbers[12]=hypot(high,up);
        }
        csv << energy << ',' << kt << ',' << ktbins[kt] << ',' << ktbins[kt+1] << ',' << mt << ',' << parameters[p];
        for(double number : numbers) csv << ',' << number;
        csv << ',' << (ready ? "complete" : "missing_or_invalid_fit") << '\n';
      }
    }
    TString folder=plotDir+"/"+energy; gSystem->mkdir(folder,true);
    TString label=Form("#sqrt{s_{NN}} = %g GeV",energyValue(energy));
    TString heading="UrQMD Au+Au, "+data[0].centrality+"%, FFT; "+label;
    output.cd();
    for(int p=0; p<NPAR; ++p) {
      stat[p].Write(Form("%s_stat_%s_FFT",parameters[p],energy.Data()));
      syst[p].Write(Form("%s_syserr_%s_FFT",parameters[p],energy.Data()));
      total[p].Write(Form("%s_total_%s_FFT",parameters[p],energy.Data()));
      for(int s=0; s<NSYST; ++s) contribution[s][p].Write(Form("%s_%s_%s_FFT",parameters[p],sources[s],energy.Data()));
      drawPlot({&stat[p]},{&syst[p]},{"FFT C(q)"},heading,"m_{T} [GeV/c^{2}]",titles[p],folder+"/mT_"+parameters[p]+"_FFT");
      drawPlot({&relative[0][p],&relative[1][p],&relative[2][p]},
               {&relative[0][p],&relative[1][p],&relative[2][p]},
               {"Analysis q_{LCMS}","Fit range","FFT binning"},heading,"m_{T} [GeV/c^{2}]","Relative systematic shift [%]",
               folder+"/contributions_"+parameters[p]+"_FFT",false,false);
      // The old D(rho) output stores systematic bands, but no statistical bars.
      if(drho && p!=0) {
        auto* old=drho->Get<TGraphAsymmErrors>(Form("%s_syserr_%s",parameters[p],energy.Data()));
        if(old) {
          TGraphAsymmErrors band(*old);
          for(int i=band.GetN()-1; i>=0; --i)
            if(!isfinite(band.GetPointY(i)) || band.GetPointY(i)<=0.) band.RemovePoint(i);
          TGraphAsymmErrors points(band);
          for(int i=0; i<points.GetN(); ++i) points.SetPointError(i,0.,0.,0.,0.);
          drawPlot({&stat[p],&points},{&syst[p],&band},{"FFT C(q)","D(#rho), systematic only"},heading,
                   "m_{T} [GeV/c^{2}]",titles[p],folder+"/comparison_Drho_"+parameters[p]+"_FFT");
        }
      }
      allStat[p].push_back((TGraphAsymmErrors*)stat[p].Clone());
      allSyst[p].push_back((TGraphAsymmErrors*)syst[p].Clone());
    }
    energyLabels.push_back(label);
  }
  for(int p=0; p<NPAR; ++p) {
    drawPlot(allStat[p],allSyst[p],energyLabels,"UrQMD Au+Au, FFT","m_{T} [GeV/c^{2}]",titles[p],plotDir+"/mT_overlay_"+parameters[p]+"_FFT");
    for(int kt=0; kt<NKT; ++kt) {
      TGraphAsymmErrors stat,syst;
      double mt=sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi);
      for(size_t e=0; e<selected.size(); ++e) for(int i=0; i<allStat[p][e]->GetN(); ++i) {
        if(abs(allStat[p][e]->GetPointX(i)-mt)>1.e-8) continue;
        int point=stat.GetN(); double value=allStat[p][e]->GetPointY(i);
        stat.SetPoint(point,energyValue(selected[e]),value);
        stat.SetPointError(point,0.,0.,allStat[p][e]->GetErrorYlow(i),allStat[p][e]->GetErrorYhigh(i));
        syst.SetPoint(point,energyValue(selected[e]),value);
        double width=0.025*energyValue(selected[e]);
        syst.SetPointError(point,width,width,allSyst[p][e]->GetErrorYlow(i),allSyst[p][e]->GetErrorYhigh(i));
      }
      output.cd(); stat.Write(Form("%s_sNN_stat_kt%d_FFT",parameters[p],kt)); syst.Write(Form("%s_sNN_syserr_kt%d_FFT",parameters[p],kt));
      drawPlot({&stat},{&syst},{"FFT C(q)"},Form("UrQMD Au+Au, %.3f < k_{T} < %.3f GeV/c",ktbins[kt],ktbins[kt+1]),
               "#sqrt{s_{NN}} [GeV]",titles[p],plotDir+Form("/sNN_%s_kt%d_FFT",parameters[p],kt),true);
    }
    for(auto* graph : allStat[p]) delete graph;
    for(auto* graph : allSyst[p]) delete graph;
  }
  output.cd();
  TNamed("fft_systematics_settings",Form("all events; qLCMS indices 0/1/2; qmax scales %.6g/%.6g/%.6g; source bins %d/%d/%d; half-range %.6g fm",
    fftQmaxFitScale[0],fftQmaxFitScale[1],fftQmaxFitScale[2],fftSourceBinsSyst[0],fftSourceBinsSyst[1],fftSourceBinsSyst[2],fftSourceHalfRange)).Write();
  TNamed("fft_systematics_convention",use_larger_same_side_shift ? "Same-side: larger shift; opposite sides: separate; sources combined in quadrature" :
         "Same-side: mean shift; opposite sides: separate; sources combined in quadrature").Write();
  TParameter<int>("fft_systematics_complete",!incomplete).Write();
  TParameter<int>("fft_systematics_points",accepted).Write();
  output.Close(); if(drho) { drho->Close(); delete drho; }
  cout << "Wrote " << accepted << " complete kT/energy points to " << outputName << "; table: " << tableName << endl;
  return incomplete ? 2 : 0;
}
