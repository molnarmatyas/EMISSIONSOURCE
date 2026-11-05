// Compare the graphs saved by the D(rho) and FFT final plotters.
// No fits or systematic variations are repeated here.
#include "fft_common.h"
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
using namespace std;

const int NQUANT=8;
const char* quantities[NQUANT]={"alpha","R_out","R_side","R_long","R_RMS",
                               "Rout_over_Rside","Rout2_minus_Rside2","normalization"};
const char* titles[NQUANT]={"#alpha","R_{out} [fm]","R_{side} [fm]","R_{long} [fm]",
                          "R_{RMS} [fm]","R_{out}/R_{side}","R_{out}^{2}-R_{side}^{2} [fm^{2}]","#lambda^{*}"};
const char* methods[]={"Drho","FFT"};
const char* methodLabels[]={"D(#rho)","FFT C(q)"};
// Display only: baseline errors remain in the ROOT graphs and CSV table.
bool showStatErrors=false;

struct Point {
  bool valid=false, hasStat=false;
  double value=0., statDown=0., statUp=0., sysDown=0., sysUp=0.;
};
struct EnergyData {
  TString name;
  double energy;
  Point points[2][NQUANT][NKT];
};
struct Series {
  int method;
  int kt=-1, group=0;
  TGraphAsymmErrors points, band;
};
struct Panel {
  TString label, yTitle;
  vector<Series> series;
};

double energyValue(TString name) { name.ReplaceAll("p","."); return name.Atof(); }

// Both final plotters save mT coordinates. Match by coordinate, never by index:
// invalid kT bins can be absent, and older files can contain zero placeholders.
int findPoint(TGraphAsymmErrors* graph, double mt)
{
  if(graph) for(int i=0; i<graph->GetN(); ++i)
    if(abs(graph->GetPointX(i)-mt)<1.e-6) return i;
  return -1;
}

bool validErrors(TGraphAsymmErrors* graph, int i)
{
  return i>=0 && isfinite(graph->GetPointY(i)) &&
    isfinite(graph->GetErrorYlow(i)) && isfinite(graph->GetErrorYhigh(i)) &&
    graph->GetErrorYlow(i)>=0. && graph->GetErrorYhigh(i)>=0.;
}

Point readPoint(TGraphAsymmErrors* syst, TGraphAsymmErrors* stat, int kt)
{
  Point p;
  double mt=sqrt(kT_center[kt]*kT_center[kt]+Mass2_pi);
  int i=findPoint(syst,mt), j=findPoint(stat,mt);
  if(!validErrors(syst,i) || syst->GetPointY(i)<=0.) return p;
  p.valid=true; p.value=syst->GetPointY(i);
  p.sysDown=syst->GetErrorYlow(i); p.sysUp=syst->GetErrorYhigh(i);
  if(validErrors(stat,j) && abs(stat->GetPointY(j)-p.value)<1.e-6*max(1.,abs(p.value))) {
    p.hasStat=true; p.statDown=stat->GetErrorYlow(j); p.statUp=stat->GetErrorYhigh(j);
  }
  return p;
}

// First-order propagation of the saved marginal errors, without covariances.
// A negative derivative swaps up/down: increasing Rside decreases Rout/Rside.
Point derive(const Point* radii, int quantity)
{
  Point result;
  int n=quantity==4 ? 3 : 2;
  for(int i=0; i<n; ++i) if(!radii[i].valid) return result;
  double out=radii[0].value, side=radii[1].value, lon=radii[2].value;
  double gradient[3]={};
  if(quantity==4) {
    result.value=sqrt((out*out+side*side+lon*lon)/3.);
    for(int i=0; i<3; ++i) gradient[i]=radii[i].value/(3.*result.value);
  } else if(quantity==5) {
    result.value=out/side; gradient[0]=1./side; gradient[1]=-out/(side*side);
  } else {
    result.value=out*out-side*side; gradient[0]=2.*out; gradient[1]=-2.*side;
  }
  result.valid=true; result.hasStat=true;
  double statDown2=0.,statUp2=0.,sysDown2=0.,sysUp2=0.;
  for(int i=0; i<n; ++i) {
    const auto& r=radii[i]; double g=gradient[i];
    result.hasStat=result.hasStat && r.hasStat;
    statDown2+=pow(g*(g>=0. ? r.statDown : r.statUp),2);
    statUp2+=pow(g*(g>=0. ? r.statUp : r.statDown),2);
    sysDown2+=pow(g*(g>=0. ? r.sysDown : r.sysUp),2);
    sysUp2+=pow(g*(g>=0. ? r.sysUp : r.sysDown),2);
  }
  if(result.hasStat) { result.statDown=sqrt(statDown2); result.statUp=sqrt(statUp2); }
  result.sysDown=sqrt(sysDown2); result.sysUp=sqrt(sysUp2);
  return result;
}

void append(Series& series, const Point& p, double x, double xlow, double xhigh)
{
  if(!p.valid) return;
  int i=series.points.GetN();
  series.points.SetPoint(i,x,p.value);
  series.points.SetPointError(i,0.,0.,p.statDown,p.statUp);
  series.band.SetPoint(i,x,p.value);
  series.band.SetPointError(i,xlow,xhigh,p.sysDown,p.sysUp);
}

void style(Series& series)
{
  const int ktColors[]={kRed+1,kBlue+1,kGreen+2,kMagenta+1};
  int color=series.kt<0 ? (series.method==0 ? kRed+1 : kBlue+1) : ktColors[series.group%4];
  int marker=series.kt<0 ? 20 : (series.group%2==0 ? 20 : 21);
  for(auto* graph : {&series.points,&series.band}) {
    graph->SetLineColor(color); graph->SetMarkerColor(color); graph->SetFillColor(color);
    graph->SetLineWidth(2); graph->SetLineStyle(series.method==0 ? 1 : 2);
    graph->SetFillStyle(series.method==0 ? 3004 : 3005);
    graph->SetMarkerStyle(marker+(series.method==0 ? 0 : 4)); graph->SetMarkerSize(1.15);
  }
}

void writeSeries(Series& series, TString name, const char* xTitle, const char* yTitle)
{
  style(series);
  series.points.SetNameTitle(name+"_points",Form(";%s;%s",xTitle,yTitle));
  series.band.SetNameTitle(name+"_syserr",Form(";%s;%s",xTitle,yTitle));
  series.points.Write(); series.band.Write();
}

// One common y range is calculated from every panel, including error extents.
void drawPanels(vector<Panel>& panels, TString name, TString heading, const char* xTitle,
                double xmin, double xmax, bool energyAxis, TString plotDir)
{
  double ymin=1.e30,ymax=-1.e30;
  for(auto& panel : panels) for(auto& s : panel.series) for(int i=0; i<s.points.GetN(); ++i) {
    double y=s.points.GetPointY(i);
    ymin=min(ymin,y-max(showStatErrors ? s.points.GetErrorYlow(i) : 0.,s.band.GetErrorYlow(i)));
    ymax=max(ymax,y+max(showStatErrors ? s.points.GetErrorYhigh(i) : 0.,s.band.GetErrorYhigh(i)));
  }
  if(ymin>ymax) return;
  double dy=max(.03,ymax-ymin); ymin-=.08*dy; ymax+=.40*dy;
  int n=panels.size();
  TCanvas canvas(name,"",n==1 ? 900 : 1550,680);
  for(int p=0; p<n; ++p) {
    canvas.cd();
    auto* pad=new TPad(Form("panel_%d",p),"",double(p)/n,0.,double(p+1)/n,1.);
    pad->Draw(); pad->cd();
    pad->SetLeftMargin(p==0 ? .17 : 0.); pad->SetRightMargin(p==n-1 ? .045 : 0.);
    pad->SetBottomMargin(.16); pad->SetTopMargin(.045); pad->SetTicks(1,1); pad->SetLogx(energyAxis);
    auto* frame=pad->DrawFrame(xmin,ymin,xmax,ymax);
    frame->SetName(Form("frame_%d",p));
    frame->SetTitle(Form(";%s;%s",xTitle,p==0 ? panels[p].yTitle.Data() : ""));
    for(auto* axis : {frame->GetXaxis(),frame->GetYaxis()}) {
      axis->SetTitleFont(43); axis->SetLabelFont(43);
      axis->SetTitleSize(27); axis->SetLabelSize(n==3 ? 20 : 23); axis->SetTitleOffset(1.25);
    }
    frame->GetYaxis()->SetTitleOffset(n==3 ? 1.45 : 1.25);
    if(p>0) frame->GetYaxis()->SetLabelSize(0.);
    if(!energyAxis) frame->GetXaxis()->SetNdivisions(505);
    if(energyAxis) { frame->GetXaxis()->SetMoreLogLabels(); frame->GetXaxis()->SetNoExponent(); }
    if(energyAxis && p<n-1) {
      double decade=pow(10.,floor(log10(xmax)));
      frame->GetXaxis()->ChangeLabelByValue(floor(xmax/decade)*decade,-1,0.);
    }
    for(auto& series : panels[p].series) if(series.band.GetN()>0)
      series.band.Draw(series.band.GetN()>1 ? "3 SAME" : "2 SAME");
    for(auto& series : panels[p].series) if(series.points.GetN()>0)
      series.points.Draw(showStatErrors ? "LPZ SAME" : "LPX SAME");
    double left=pad->GetLeftMargin(),width=1.-left-pad->GetRightMargin();
    auto x=[&](double value){return left+value*width;};
    auto y=[&](double value){return pad->GetBottomMargin()+value*(1.-pad->GetBottomMargin()-pad->GetTopMargin());};
    TLatex text; text.SetNDC(); text.SetTextFont(62); text.SetTextSize(.035);
    if(p==0) {
      TString first=heading,second;
      if(n>1 && heading.Contains(",")) { int split=heading.First(','); first=heading(0,split); second=heading(split+1,heading.Length()); }
      text.DrawLatex(x(.035),y(.945),first);
      if(!second.IsNull()) text.DrawLatex(x(.035),y(.885),second);
    }
    // The middle radius panel has room for one shared legend.
    if((n==1 || p==1) && panels[p].series.size()>1) {
      bool grouped=energyAxis && panels[p].series.front().kt>=0;
      auto* legend=new TLegend(x(.04),y(grouped ? .75 : .80),x(.97),y(.90));
      legend->SetFillColor(kWhite); legend->SetFillStyle(1001); legend->SetBorderSize(1);
      legend->SetLineColor(kGray+1); legend->SetLineWidth(1); legend->SetTextFont(43); legend->SetTextSize(n==1 ? 20 : 17);
      legend->SetNColumns(2);
      if(grouped) legend->SetHeader("m_{T} [MeV/c^{2}]");
      for(auto& series : panels[p].series) if(series.points.GetN()>0) {
        TString entry=methodLabels[series.method];
        if(grouped) {
          double low=1000.*sqrt(ktbins[series.kt]*ktbins[series.kt]+Mass2_pi);
          double high=1000.*sqrt(ktbins[series.kt+1]*ktbins[series.kt+1]+Mass2_pi);
          entry=Form("%s, %.0f#minus%.0f",series.method==0 ? "D(#rho)" : "FFT",low,high);
        }
        legend->AddEntry(&series.band,entry,"lfp");
      }
      legend->Draw();
    }
    TString label=panels[p].label; label.ReplaceAll(" [fm]","");
    label.ReplaceAll("R_{",""); label.ReplaceAll("}","");
    if(!label.IsNull()) { text.SetTextSize(.040); text.DrawLatex(x(.75),y(.945),label); }
    bool empty=true;
    for(auto& series : panels[p].series) if(series.points.GetN()>0) empty=false;
    if(empty) text.DrawLatex(x(.25),y(.50),"No valid points");
  }
  canvas.cd();
  canvas.Write();
  canvas.SaveAs(plotDir+"/"+name+".png"); canvas.SaveAs(plotDir+"/"+name+".pdf");
}

int main(int argc,char** argv)
{
  TString drhoName="syserr_results.root", fftName="FFTmethod/results/syserr_results_FFT.root";
  TString outputName="FFTmethod/results/comparison_Drho_FFT.root", plotDir="figs/comparison/Drho_FFT";
  vector<TString> selected;
  vector<int> ktBins;
  bool makePlots=true;
  for(int i=1; i<argc; ++i) {
    string option=argv[i];
    if(option=="--no-plots") { makePlots=false; continue; }
    if(option=="--help") {
      cout << "Usage: compare_methods.exe [--drho syserr_results.root] [--fft syserr_results_FFT.root]\n"
           << "  [--energies 3p0,5p2,27] [--kt 0,4] [--output comparison_FFT.root]\n"
           << "  [--plot-dir DIR] [--no-plots]\n"
           << "Default: common energies; energy plots use ikt_indices_to_plot in the common header.\n"
           << "kT plots always include all available bins.\n";
      return 0;
    }
    if(i+1==argc) { cerr << "Missing value for " << option << endl; return 1; }
    TString value=argv[++i];
    if(option=="--drho") drhoName=value;
    else if(option=="--fft") fftName=value;
    else if(option=="--output") outputName=fftOutputName(value,"");
    else if(option=="--plot-dir") plotDir=value;
    else if(option=="--energies" || option=="--kt") {
      istringstream list(value.Data()); string item;
      while(getline(list,item,',')) {
        if(option=="--energies") {
          if(find(begin(energies),end(energies),TString(item))==end(energies)) {
            cerr << "Unknown energy: " << item << endl; return 1;
          }
          if(find(selected.begin(),selected.end(),TString(item))==selected.end()) selected.push_back(item);
        } else {
          int kt=-1; istringstream number(item);
          if(!(number>>kt) || !number.eof() || kt<0 || kt>=NKT) { cerr << "Invalid kT index: " << item << endl; return 1; }
          if(find(ktBins.begin(),ktBins.end(),kt)==ktBins.end()) ktBins.push_back(kt);
        }
      }
    } else { cerr << "Unknown option: " << option << endl; return 1; }
  }
  if(ktBins.empty()) ktBins.assign(begin(ikt_indices_to_plot),end(ikt_indices_to_plot));
  sort(ktBins.begin(),ktBins.end());
  TFile drho(drhoName,"READ"), fft(fftName,"READ");
  if(drho.IsZombie() || fft.IsZombie()) { cerr << "Run both final plotters first, or supply their ROOT files with --drho and --fft.\n"; return 1; }
  auto outputPath=filesystem::weakly_canonical(outputName.Data());
  if(outputPath==filesystem::weakly_canonical(drhoName.Data()) || outputPath==filesystem::weakly_canonical(fftName.Data())) {
    cerr << "Comparison output must differ from the input files.\n"; return 1;
  }
  auto* complete=fft.Get<TParameter<int>>("fft_systematics_complete");
  if(complete && complete->GetVal()!=1) cerr << "FFT systematics are incomplete; plotting available points only.\n";
  if(selected.empty()) for(const char* e : energies)
    if(drho.Get(Form("alpha_syserr_%s",e)) && fft.Get(Form("alpha_syserr_%s_FFT",e))) selected.push_back(e);
  sort(selected.begin(),selected.end(),[](TString a,TString b){return energyValue(a)<energyValue(b);});
  vector<EnergyData> data;
  bool missingStat=false;
  int common=0;
  for(TString energy : selected) {
    EnergyData d; d.name=energy; d.energy=energyValue(energy);
    for(int m=0; m<2; ++m) for(int q : {0,1,2,3,7}) {
      TString parameter=q==7 ? (m==0 ? "N" : "lambda") : quantities[q];
      TString suffix=m==0 ? "" : "_FFT";
      TFile& file=m==0 ? drho : fft;
      auto* syst=file.Get<TGraphAsymmErrors>(parameter+"_syserr_"+energy+suffix);
      auto* stat=file.Get<TGraphAsymmErrors>(parameter+"_stat_"+energy+suffix);
      for(int kt=0; kt<NKT; ++kt) {
        auto& point=d.points[m][q][kt]; point=readPoint(syst,stat,kt);
        if(m==0 && point.valid && !point.hasStat) missingStat=true;
        if(m==1 && point.valid && !point.hasStat) {
          cerr << "Missing/inconsistent FFT baseline errors for " << parameter << ", " << energy << ", kT bin " << kt << endl;
          return 1;
        }
      }
    }
    for(int m=0; m<2; ++m) for(int kt=0; kt<NKT; ++kt) {
      Point radii[]={d.points[m][1][kt],d.points[m][2][kt],d.points[m][3][kt]};
      for(int q=4; q<=6; ++q) d.points[m][q][kt]=derive(radii,q);
    }
    for(int q=0; q<4; ++q) for(int kt=0; kt<NKT; ++kt)
      if(d.points[0][q][kt].valid && d.points[1][q][kt].valid) ++common;
    data.push_back(d);
  }
  if(common==0) { cerr << "No common valid alpha/radius points found. Check the energy selection and input graphs.\n"; return 1; }
  if(missingStat) cerr << "Some D(rho) baseline errors are unavailable and recorded as nan in the CSV.\n";
  cout << "Plot display: connected systematic bands; baseline bars " << (showStatErrors ? "shown" : "hidden")
       << "; stored graph errors and CSV values are retained.\n"
       << "Derived uncertainties use independent-radius propagation. Normalization overlays use the original source/correlation fit amplitudes.\n";
  for(TString folder : {TString(gSystem->DirName(outputName)),plotDir})
    if(gSystem->AccessPathName(folder) && gSystem->mkdir(folder,true)!=0) { cerr << "Cannot create " << folder << endl; return 1; }
  gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TFile output(outputName,"RECREATE");
  TString csvName=outputName; csvName.ReplaceAll(".root",".csv");
  ofstream csv(csvName.Data());
  if(output.IsZombie() || !csv) { cerr << "Cannot create comparison outputs.\n"; return 1; }
  csv << setprecision(12) << "method,energy,energy_GeV,kt,kt_low,kt_high,quantity,value,baseline_error_available,baseline_down,baseline_up,sys_down,sys_up\n";
  for(const auto& d : data) {
    TString folder=plotDir+"/"+d.name;
    if(makePlots) gSystem->mkdir(folder,true);
    vector<Panel> radii;
    for(int q=1; q<=3; ++q) radii.push_back({titles[q],"R [fm]",{}});
    for(int q=0; q<NQUANT; ++q) {
      vector<Panel> panels={{"",titles[q],{}}};
      for(int m=0; m<2; ++m) {
        Series series; series.method=m;
        for(int kt=0; kt<NKT; ++kt) {
          const auto& point=d.points[m][q][kt];
          append(series,point,kT_center[kt],kT_center[kt]-ktbins[kt],ktbins[kt+1]-kT_center[kt]);
          if(!point.valid) continue;
          csv << methods[m] << ',' << d.name << ',' << d.energy << ',' << kt << ',' << ktbins[kt] << ',' << ktbins[kt+1] << ','
              << (q==7 ? (m==0 ? "N" : "lambda") : quantities[q]) << ',' << point.value << ',' << point.hasStat << ',';
          if(point.hasStat) csv << point.statDown << ',' << point.statUp;
          else csv << "nan,nan";
          csv << ',' << point.sysDown << ',' << point.sysUp << '\n';
        }
        writeSeries(series,Form("%s_%s_kT_%s",quantities[q],methods[m],d.name.Data()),"k_{T} [GeV/c]",titles[q]);
        auto& panel=(q>=1 && q<=3) ? radii[q-1] : panels[0];
        panel.series.push_back(series);
      }
      if(makePlots && (q==0 || q>=4)) drawPanels(panels,Form("kT_%s_%s_comparison_FFT",quantities[q],d.name.Data()),
        Form("UrQMD Au+Au 0#minus10%%, #sqrt{s_{NN}} = %g GeV",d.energy),"k_{T} [GeV/c]",ktbins[0]-.015,ktbins[NKT]+.015,
        false,folder);
    }
    if(makePlots) drawPanels(radii,Form("kT_radii_%s_comparison_FFT",d.name.Data()),
      Form("UrQMD Au+Au 0#minus10%%, #sqrt{s_{NN}} = %g GeV",d.energy),"k_{T} [GeV/c]",ktbins[0]-.015,ktbins[NKT]+.015,
      false,folder);
  }
  vector<Panel> radii;
  for(int q=1; q<=3; ++q) radii.push_back({titles[q],"R [fm]",{}});
  for(int q=0; q<NQUANT; ++q) {
    vector<Panel> single={{"",titles[q],{}}};
    for(size_t group=0; group<ktBins.size(); ++group) {
      int kt=ktBins[group];
      for(int m=0; m<2; ++m) {
        auto& panel=(q>=1 && q<=3) ? radii[q-1] : single[0];
        Series series; series.method=m; series.kt=kt; series.group=group;
        for(const auto& d : data) append(series,d.points[m][q][kt],d.energy,.018*d.energy,.018*d.energy);
        writeSeries(series,Form("%s_%s_sNN_kt%d",quantities[q],methods[m],kt),"#sqrt{s_{NN}} [GeV]",titles[q]);
        panel.series.push_back(series);
      }
    }
    if(makePlots && (q==0 || q>=4)) drawPanels(single,Form("sNN_%s_comparison_FFT",quantities[q]),
      "UrQMD Au+Au 0#minus10%","#sqrt{s_{NN}} [GeV]",data.front().energy/1.15,data.back().energy*1.15,true,plotDir);
  }
  if(makePlots) drawPanels(radii,"sNN_radii_comparison_FFT","UrQMD Au+Au 0#minus10%","#sqrt{s_{NN}} [GeV]",
    data.front().energy/1.15,data.back().energy*1.15,true,plotDir);
  TNamed("comparison_drho_input",drhoName.Data()).Write(); TNamed("comparison_fft_input",fftName.Data()).Write();
  TNamed("comparison_derived_convention","R_RMS=sqrt((Rout^2+Rside^2+Rlong^2)/3), computed from each method's reported radii; ratio=Rout/Rside; difference=Rout^2-Rside^2; no kT averaging").Write();
  TNamed("comparison_error_convention","Baseline bars and systematic bands kept separate; derived errors use first-order propagation without radius covariances, swapping asymmetric up/down errors for negative derivatives; missing baseline errors are nan in CSV").Write();
  TNamed("comparison_normalizations","Source-fit N and correlation-fit lambda are overlaid as lambda*; input amplitudes and their errors are used directly").Write();
  TNamed("comparison_baseline_convention","D(rho): existing final-plotter baseline spreads; FFT: residual-rescaled fit errors; these error estimates have different meanings").Write();
  TParameter<int>("comparison_missing_drho_baseline",missingStat).Write();
  TParameter<int>("comparison_show_baseline_errors",showStatErrors).Write();
  output.Close();
  cout << "Compared " << data.size() << " energies; wrote " << outputName << " and " << csvName << endl;
  if(makePlots) cout << "Figures: " << plotDir << endl;
  return 0;
}
