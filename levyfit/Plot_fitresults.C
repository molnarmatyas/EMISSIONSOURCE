#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <cmath>
#include <iostream>

#define SQR(x) ((x) * (x))

const double Mass2_pi = 0.019479835;
const int NKT = 10;
//const double ktbins[NKT + 1] = {0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575};
const double ktbins[NKT + 1] = {0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675};

const char* path = "~/work/figs/EPOS_project/3d_proj/";
const char* frames[3] = {"lcms","pcms","lab"};
const int thisframe = 0;


void Plot_fitresults() 
{
    // Open the fitresults.root file
    TFile *file = TFile::Open(Form("fitresults_%s_010cent_all.root",frames[thisframe]));
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file fitresults.root" << std::endl;
        return;
    }

    TGraphErrors* Rosl_prediction[3];
    for(int iosl = 0; iosl < 3; iosl++)
    {
      Rosl_prediction[iosl] = new TGraphErrors(1);
      Rosl_prediction[iosl]->SetMarkerStyle(20+iosl);
      Rosl_prediction[iosl]->SetMarkerSize(1.5);
    }

    Rosl_prediction[0]->SetPoint(0,sqrt(Mass2_pi+0.22*0.22),8.434523);
    Rosl_prediction[0]->SetPointError(0,0.,0.081998);
    Rosl_prediction[1]->SetPoint(0,sqrt(Mass2_pi+0.22*0.22),7.461424);
    Rosl_prediction[1]->SetPointError(0,0.,0.061037);
    Rosl_prediction[2]->SetPoint(0,sqrt(Mass2_pi+0.22*0.22),9.520307);
    Rosl_prediction[2]->SetPointError(0,0.,0.098702);

    // Create arrays to store the histograms and TGraphErrors
    TH1D* alphahist[NKT];
    TH1D* Rohist[NKT];
    TH1D* Rshist[NKT];
    TH1D* Rlhist[NKT];
    TH1D* Ravghist[NKT];
    TH1D* clhist[NKT];

    TGraphErrors *alphaGraph = new TGraphErrors(NKT);
    TGraphErrors *RoGraph = new TGraphErrors(NKT);
    TGraphErrors *RsGraph = new TGraphErrors(NKT);
    TGraphErrors *RlGraph = new TGraphErrors(NKT);
    TGraphErrors *RavgGraph = new TGraphErrors(NKT);

    // Define custom colors
    TColor* color0 = new TColor(1500, 0.121568, 0.466666, 0.705882); // #1f77b4
    TColor* color1 = new TColor(1501, 1.000000, 0.498039, 0.054901); // #ff7f0e
    TColor* color2 = new TColor(1502, 0.172549, 0.627451, 0.172549); // #2ca02c
    TColor* color3 = new TColor(1503, 0.839215, 0.152941, 0.156862); // #d62728
    TColor* color4 = new TColor(1504, 0.580392, 0.403921, 0.741176); // #9467bd
    TColor* color5 = new TColor(1505, 0.549019, 0.337254, 0.294117); // #8c564b
    TColor* color6 = new TColor(1506, 0.890196, 0.466666, 0.760784); // #e377c2
    TColor* color7 = new TColor(1507, 0.498039, 0.498039, 0.498039); // #7f7f7f
    TColor* color8 = new TColor(1508, 0.737254, 0.741176, 0.133333); // #bcbd22
    TColor* color9 = new TColor(1509, 0.090196, 0.745098, 0.811764); // #17becf

    // Loop over the kT indices
    for (int ikt = 0; ikt < NKT; ++ikt) {
        // Form the histogram names and retrieve them from the file
        alphahist[ikt] = (TH1D*)file->Get(Form("alphahist_ikt%d", ikt));
        Rohist[ikt] = (TH1D*)file->Get(Form("Rohist_ikt%d", ikt));
        Rshist[ikt] = (TH1D*)file->Get(Form("Rshist_ikt%d", ikt));
        Rlhist[ikt] = (TH1D*)file->Get(Form("Rlhist_ikt%d", ikt));
        Ravghist[ikt] = (TH1D*)file->Get(Form("Ravghist_ikt%d", ikt));
        clhist[ikt] = (TH1D*)file->Get(Form("clhist_ikt%d", ikt));

        if (!alphahist[ikt] || !Rohist[ikt] || !Rshist[ikt] || !Rlhist[ikt]) {
            std::cerr << "Error retrieving histograms for kT index " << ikt << std::endl;
            continue;
        }

        // Create a canvas and plot the R histograms on top of each other
        TCanvas *c1 = new TCanvas(Form("c1_ikt%d", ikt), Form("c1_ikt%d", ikt), 800, 600);
        Rshist[ikt]->SetLineColor(kRed);
        Rohist[ikt]->SetLineColor(kBlue);
        Rlhist[ikt]->SetLineColor(kGreen);

        Rshist[ikt]->Draw();
        Rohist[ikt]->Draw("SAME");
        Rlhist[ikt]->Draw("SAME");

        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(Rohist[ikt], "R_{o}", "l");
        legend->AddEntry(Rshist[ikt], "R_{s}", "l");
        legend->AddEntry(Rlhist[ikt], "R_{l}", "l");
        legend->Draw();

        c1->SaveAs(Form("%s%s/R_ikt%d.pdf", path, frames[thisframe], ikt));
        delete c1;

        // Plot alpha histogram separately
        TCanvas *c2 = new TCanvas(Form("c2_ikt%d", ikt), Form("c2_ikt%d", ikt), 800, 600);
        alphahist[ikt]->Draw();
        c2->SaveAs(Form("%s%s/alpha_ikt%d.pdf", path, frames[thisframe], ikt));
        c2->Clear();
        c2->SetLogx(1);
        clhist[ikt]->Draw();
        c2->SaveAs(Form("%s%s/cl_ikt%d.pdf", path, frames[thisframe], ikt));
        c2->SetLogx(0);
        delete c2;
        // Calculate the x value for the graph
        double kt_mid = 0.5 * (ktbins[ikt] + ktbins[ikt + 1]);
        double x_value = std::sqrt(Mass2_pi + SQR(kt_mid));

        // Extract mean and stddev values and set points in TGraphErrors
        double mean_alpha = alphahist[ikt]->GetMean();
        double stddev_alpha = alphahist[ikt]->GetStdDev();
        alphaGraph->SetPoint(ikt, x_value, mean_alpha);
        alphaGraph->SetPointError(ikt, 0, stddev_alpha);

        double mean_Ro = Rohist[ikt]->GetMean();
        double stddev_Ro = Rohist[ikt]->GetStdDev();
        RoGraph->SetPoint(ikt, x_value, mean_Ro);
        RoGraph->SetPointError(ikt, 0, stddev_Ro);

        double mean_Rs = Rshist[ikt]->GetMean();
        double stddev_Rs = Rshist[ikt]->GetStdDev();
        RsGraph->SetPoint(ikt, x_value, mean_Rs);
        RsGraph->SetPointError(ikt, 0, stddev_Rs);

        double mean_Rl = Rlhist[ikt]->GetMean();
        double stddev_Rl = Rlhist[ikt]->GetStdDev();
        RlGraph->SetPoint(ikt, x_value, mean_Rl);
        RlGraph->SetPointError(ikt, 0, stddev_Rl);

        double mean_Ravg = Ravghist[ikt]->GetMean();
        double stddev_Ravg = Ravghist[ikt]->GetStdDev();
        RavgGraph->SetPoint(ikt, x_value, mean_Ravg);
        RavgGraph->SetPointError(ikt, 0, stddev_Ravg);
    }
    // Close the file
    file->Close();
    delete file;
    // Open the fitresults.root file
    TFile *file2 = TFile::Open("phenixdata.root");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error opening file2 phenixdata.root" << std::endl;
        return;
    }
    TGraphAsymmErrors* alphaphenix = (TGraphAsymmErrors*)file2->Get("fit_graph_avgdet_cut0_par2_start0_cent0");
    alphaphenix->SetMarkerStyle(20);
    alphaphenix->SetMarkerColor(1);
    alphaphenix->RemovePoint(0);
    TGraphAsymmErrors* Rphenix = (TGraphAsymmErrors*)file2->Get("fit_graph_avgdet_cut0_par1_start0_cent0");
    Rphenix->SetMarkerStyle(20);
    Rphenix->SetMarkerColor(1);
    Rphenix->RemovePoint(0);

    file2->Close();
    delete file2;

    TFile *file3 = TFile::Open("phenixsyst.root");
    if (!file3 || file3->IsZombie()) {
        std::cerr << "Error opening file3 phenixsyst.root" << std::endl;
        return;
    }
    TGraph* asymmsystalphaup = (TGraph*)file3->Get("alpha_up");
    TGraph* asymmsystalphadn = (TGraph*)file3->Get("alpha_down");
    TGraph* asymmsystRup = (TGraph*)file3->Get("R_up");
    TGraph* asymmsystRdn = (TGraph*)file3->Get("R_down");

    TGraphAsymmErrors* alphasyst = new TGraphAsymmErrors(alphaphenix->GetN());
    TGraphAsymmErrors* Rsyst = new TGraphAsymmErrors(Rphenix->GetN());
   
    for(int ibin = 0; ibin < alphaphenix->GetN(); ibin++)
    {
      double mt = alphaphenix->GetX()[ibin];
      double val = alphaphenix->GetY()[ibin];
      double errup = val * asymmsystalphaup->GetY()[ibin+1];
      double errdn = val * asymmsystalphadn->GetY()[ibin+1];
      alphasyst->SetPoint(ibin,mt,val);
      alphasyst->SetPointError(ibin,0.01,0.01,errdn,errup);
    }
    alphasyst->SetLineWidth(1);
    alphasyst->SetFillStyle(0);
    alphasyst->SetLineStyle(1);
    alphasyst->SetLineColor(1);

    for(int ibin = 0; ibin < Rphenix->GetN(); ibin++)
    {
      double mt = Rphenix->GetX()[ibin];
      double val = Rphenix->GetY()[ibin];
      double errup = val * asymmsystRup->GetY()[ibin+1];
      double errdn = val * asymmsystRdn->GetY()[ibin+1];
      Rsyst->SetPoint(ibin,mt,val);
      Rsyst->SetPointError(ibin,0.01,0.01,errdn,errup);
    }
    Rsyst->SetLineWidth(1);
    Rsyst->SetFillStyle(0);
    Rsyst->SetLineStyle(1);
    Rsyst->SetLineColor(1);

    file3->Close();
    delete file3;
    // Create a canvas for the alpha kT dependency plot
    TCanvas *c3 = new TCanvas("alpha_ktdep", "", 800, 600);
    c3->SetTopMargin(0.015);
    c3->SetRightMargin(0.01);
    c3->SetBottomMargin(0.14);
    c3->SetLeftMargin(0.12);
    gStyle->SetLabelSize(0.06,"x");
    gStyle->SetLabelSize(0.06,"y");
    gStyle->SetTitleSize(0.06,"x");
    gStyle->SetTitleSize(0.06,"y");
    gStyle->SetTitleOffset(0.95,"y");
    double mtmin = 0.18;
    double mtmax = 0.72;
    TLatex Tl;
    Tl.SetNDC(true);
//    Tl.SetTextSize(20);
    TLegend *legend1 = new TLegend(0.15, 0.80, 0.95, 0.90);
    legend1->SetNColumns(2);
    legend1->SetTextSize(0.05);
    legend1->SetBorderSize(0.);
    legend1->SetFillStyle(0);

    alphaGraph->SetTitle(";m_{T} [GeV/c];#alpha");
    alphaGraph->SetMarkerStyle(20);
    alphaGraph->SetFillColorAlpha(1503,0.5);
    alphaGraph->SetFillStyle(1001);
    alphaGraph->SetLineStyle(2);
    alphaGraph->SetLineColor(1503);
    alphaGraph->SetLineWidth(3);
    alphaGraph->SetMinimum(0.95);
    alphaGraph->SetMaximum(2.05);
    alphaGraph->GetXaxis()->SetLimits(mtmin,mtmax);
    alphaGraph->Draw("A3");
    alphaGraph->Draw("LX");
    alphaphenix->Draw("pe same");
    alphasyst->Draw("e2");
    for(int ibin = 0; ibin < alphaGraph->GetN(); ibin++)
      cout << alphaGraph->GetX()[ibin] << "\t" << alphaGraph->GetY()[ibin] << "\t" << alphaGraph->GetErrorY(ibin) << endl;
    legend1->AddEntry(alphaphenix, "#alpha PHENIX [2407.08586]", "pe");
    legend1->AddEntry(alphaGraph, "#LT#alpha#GT EPOS3", "l");
    legend1->Draw("same");
    Tl.DrawLatex(0.18,0.92,"0#minus10% Au+Au @ #sqrt{s_{NN}} = 200 GeV, #pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
//    c3->SaveAs(Form("%s%s/alpha_ktdep.pdf",path,frames[thisframe]));
    c3->Clear();

    // Create a canvas for the alpha kT dependency plot
    TLegend *legend2 = new TLegend(0.5, 0.70, 0.85, 0.88);
    legend2->SetTextSize(0.05);
    legend2->SetBorderSize(0.);
    legend2->SetFillStyle(0);

    RavgGraph->SetTitle(";m_{T} [GeV/c];R [fm]");
    RavgGraph->SetMarkerStyle(20);
    RavgGraph->SetMinimum(3.8);
    RavgGraph->SetMaximum(12.8);
    RavgGraph->GetXaxis()->SetLimits(mtmin,mtmax);
    RavgGraph->SetFillColorAlpha(1504,0.5);
    RavgGraph->SetFillStyle(1001);
    RavgGraph->SetLineStyle(2);
    RavgGraph->SetLineColor(1504);
    RavgGraph->SetLineWidth(3);
    RavgGraph->Draw("A3");
    RavgGraph->Draw("LX");
    Rphenix->Draw("pe same");
    Rsyst->Draw("e2");
    for(int ibin = 0; ibin < RavgGraph->GetN(); ibin++)
      cout << RavgGraph->GetX()[ibin] << "\t" << RavgGraph->GetY()[ibin] << "\t" << RavgGraph->GetErrorY(ibin) << endl;
    Tl.DrawLatex(0.52,0.62,"#bar{R} = #sqrt{(R_{o}^{2}+R_{s}^{2}+R_{l}^{2})/3}");
    legend2->AddEntry(Rphenix, "R_{1d} PHENIX [2407.08586]", "pe");
    legend2->AddEntry(RavgGraph, "#LT#bar{R}#GT EPOS3", "l");
    legend2->Draw("same");
    Tl.DrawLatex(0.18,0.92,"0#minus10% Au+Au @ #sqrt{s_{NN}} = 200 GeV, #pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");
//    c3->SaveAs(Form("%s%s/Ravg_ktdep.pdf",path,frames[thisframe]));
    c3->Clear();

    TLegend *legend3 = new TLegend(0.5, 0.60, 0.85, 0.88);
    legend3->SetTextSize(0.05);
    legend3->SetBorderSize(0.);
    legend3->SetFillStyle(0);
    // Create a canvas for the R kT dependency plot
    RoGraph->SetLineColor(1500);
    RsGraph->SetLineColor(1501);
    RlGraph->SetLineColor(1502);
    RoGraph->SetMarkerColor(1500);
    RsGraph->SetMarkerColor(1501);
    RlGraph->SetMarkerColor(1502);

    RsGraph->SetTitle(";m_{T} [GeV/c];R [fm]");
    RoGraph->SetMarkerStyle(20);
    RsGraph->SetMarkerStyle(20);
    RlGraph->SetMarkerStyle(20);
    RoGraph->SetFillColorAlpha(1500,0.5);
    RoGraph->SetFillStyle(1001);
    RsGraph->SetFillColorAlpha(1501,0.5);
    RsGraph->SetFillStyle(1001);
    RlGraph->SetFillColorAlpha(1502,0.5);
    RlGraph->SetFillStyle(1001);
    RsGraph->SetMinimum(3.8);
    RsGraph->SetMaximum(12.8);
    RsGraph->GetXaxis()->SetLimits(mtmin,mtmax);
    RsGraph->Draw("A3");
    RoGraph->Draw("3 SAME");
    RlGraph->Draw("3 SAME");
    RoGraph->SetLineStyle(2);
    RsGraph->SetLineStyle(3);
    RlGraph->SetLineStyle(4);
    RoGraph->SetLineWidth(3);
    RsGraph->SetLineWidth(3);
    RlGraph->SetLineWidth(3);
    RsGraph->Draw("lx");
    RoGraph->Draw("lx");
    RlGraph->Draw("lx");
    Rphenix->Draw("pe same");
    Rsyst->Draw("e2");
    for(int ibin = 0; ibin < RoGraph->GetN(); ibin++)
      cout << RoGraph->GetX()[ibin] << "\t" << RoGraph->GetY()[ibin] << "\t" << RoGraph->GetErrorY(ibin) << endl;
    for(int ibin = 0; ibin < RsGraph->GetN(); ibin++)
      cout << RsGraph->GetX()[ibin] << "\t" << RsGraph->GetY()[ibin] << "\t" << RsGraph->GetErrorY(ibin) << endl;
    for(int ibin = 0; ibin < RlGraph->GetN(); ibin++)
      cout << RlGraph->GetX()[ibin] << "\t" << RlGraph->GetY()[ibin] << "\t" << RlGraph->GetErrorY(ibin) << endl;
//    for(int iosl = 0; iosl < 3; iosl++) Rosl_prediction[iosl]->Draw("pe same");
    Tl.DrawLatex(0.18,0.92,"0#minus10% Au+Au @ #sqrt{s_{NN}} = 200 GeV, #pi#kern[-0.3]{{}^{#plus}}#pi#kern[-0.3]{{}^{#plus}}#plus #pi#kern[-0.3]{{}^{#minus}}#pi#kern[-0.3]{{}^{#minus}}");

    legend3->AddEntry(Rphenix, "R_{1d} PHENIX [2407.08586]", "pe");
    legend3->AddEntry(RoGraph, "#LTR_{o}#GT EPOS3", "l");
    legend3->AddEntry(RsGraph, "#LTR_{s}#GT EPOS3", "l");
    legend3->AddEntry(RlGraph, "#LTR_{l}#GT EPOS3", "l");
    legend3->Draw("same");

//    c3->SaveAs(Form("%s%s/R_ktdep.pdf",path,frames[thisframe]));

    return;
}
