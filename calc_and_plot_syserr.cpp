// Collecting, dealing with, calculating systematic uncertainties:
//   qlcms
//	 rhofitmax
//	 nevt_avg
//	 mT choice
// Then, finally plotting them with result points on m_T vs Lévy parameter and sqrt(sNN) plots
// Do not forget that the point markers here should be almost the same size as the lines connecting them; the uncertainty bars should be filled areas instead of normal error bars

#include "header_for_all_emissionsource.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TSystem.h>
#include <vector>
#include <algorithm>
#include <sstream>

// To be set
int NEVT_AVG_DEFAULT = 5; // index
int NEVT_AVGsyst[] = {10, 25, 50, 100, 200, 500, 1000, 5000, 10000}; // removed 1, not meaningful for low energies, runs too long


int calc_and_plot_syserr()
{
    // Get m_T bin centers
    double mtbin_centers[NKT];
    for(int ii=0; ii<NKT; ii++)
    {
        mtbin_centers[ii] = 0.5 * (ktbins[ii] + ktbins[ii+1]);
        mtbin_centers[ii] = sqrt(mtbin_centers[ii]*mtbin_centers[ii] + Mass2_pi);
    }
    // x-error arrays (half-bin widths) - set to zero for now
    double xerr_low[NKT] = {0};
    double xerr_high[NKT] = {0};

    // CONTAINERS
    // Default systematics -> Lévy params
    TGraphAsymmErrors* alpha_default[NENERGIES];
    TGraphAsymmErrors* R_default[NENERGIES];
    TGraphAsymmErrors* N_default[NENERGIES];
    // qLCMS systematics
    TGraphAsymmErrors* alpha_qlcms[NENERGIES][3]; // 0: default, 1: strict, 2: loose
    TGraphAsymmErrors* R_qlcms[NENERGIES][3];
    TGraphAsymmErrors* N_qlcms[NENERGIES][3];
    // rhofitmax systematics
    TGraphAsymmErrors* alpha_rhofitmax[NENERGIES][3]; // 0: default, 1: strict, 2: loose
    TGraphAsymmErrors* R_rhofitmax[NENERGIES][3];
    TGraphAsymmErrors* N_rhofitmax[NENERGIES][3];
    // NEVT_AVG systematics
    TGraphAsymmErrors* alpha_nevtavg[NENERGIES][3]; // 0: default, 1: reasonably smaller, 2: reasonably larger
    TGraphAsymmErrors* R_nevtavg[NENERGIES][3];
    TGraphAsymmErrors* N_nevtavg[NENERGIES][3];

    // Main mT vs param graph
    TGraphAsymmErrors* alpha_syserr[NENERGIES];
    TGraphAsymmErrors* R_syserr[NENERGIES];
    TGraphAsymmErrors* N_syserr[NENERGIES];

    for(int ienergy = 0; ienergy < NENERGIES; ienergy++)
    {
        // Default
        TFile* f_default = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root", 
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT]), "READ");
        double alpha_vals[NKT] = {0};
        double alpha_errdn[NKT] = {0};
        double alpha_errup[NKT] = {0};
        double R_vals[NKT] = {0};
        double R_errup[NKT] = {0};
        double R_errdn[NKT] = {0};
        double N_vals[NKT] = {0};
        double N_errup[NKT] = {0};
        double N_errdn[NKT] = {0};
        for(int ikt=0; ikt<NKT; ikt++)
        {
            TH2F* alphaR = (TH2F*)f_default->Get(Form("alpha_vs_R_ikt%i", ikt));
            alpha_vals[ikt] = alphaR->GetMean(1);
            alpha_errdn[ikt] = alphaR->GetStdDev(1);
            alpha_errup[ikt] = alphaR->GetStdDev(1);
            R_vals[ikt] = alphaR->GetMean(2);
            R_errdn[ikt] = alphaR->GetStdDev(2);
            R_errup[ikt] = alphaR->GetStdDev(2);
            N_vals[ikt] = ((TH1F*)f_default->Get(Form("Nhist_ikt%i", ikt)))->GetMean();
            N_errdn[ikt] = ((TH1F*)f_default->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
            N_errup[ikt] = ((TH1F*)f_default->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
        }        
        alpha_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
        R_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
        N_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
        f_default->Close();
        // Reset arrays for next use
        for(int ikt=0; ikt<NKT; ikt++)
        {
            alpha_vals[ikt] = 0;
            alpha_errdn[ikt] = 0;
            alpha_errup[ikt] = 0;
            R_vals[ikt] = 0;
            R_errup[ikt] = 0;
            R_errdn[ikt] = 0;
            N_vals[ikt] = 0;
            N_errup[ikt] = 0;
            N_errdn[ikt] = 0;
        }

        // qLCMS
        const char* qlcms_labels[3] = {"","_strictqLCMS","_looseqLCMS"};
        for(int iq=0; iq<3; iq++)
        {
            TFile* f_qlcms = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d%s.root", 
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT], qlcms_labels[iq]), "READ");
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_qlcms->Get(Form("alpha_vs_R_ikt%i", ikt));
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                N_vals[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetMean();
                N_errdn[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
                N_errup[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
            }
            alpha_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            f_qlcms->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
            }
        }

        // rhofitmax
        const char* rhofitmax_labels[3] = {"","_strictrhoFitMax","_looserhoFitMax"};
        for(int ir=0; ir<3; ir++)
        {
            TFile* f_qlcms = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d%s.root", 
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT], rhofitmax_labels[ir]), "READ");
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_qlcms->Get(Form("alpha_vs_R_ikt%i", ikt));
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                N_vals[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetMean();
                N_errdn[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
                N_errup[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
            }
            alpha_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            f_qlcms->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
            }
        }

        // NEVT_AVG
        for(int in=0; in<3; in++)
        {
            int index = 0;
            if(in==0) index = NEVT_AVG_DEFAULT; // default
            else if(in==1) index = NEVT_AVG_DEFAULT - 1; // reasonably smaller
            else if(in==2) index = NEVT_AVG_DEFAULT + 1; // reasonably larger
            TFile* f_qlcms = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root", 
                                            energies[ienergy], NEVT_AVGsyst[index]), "READ");
            
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_qlcms->Get(Form("alpha_vs_R_ikt%i", ikt));
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                N_vals[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetMean();
                N_errdn[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
                N_errup[ikt] = ((TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt)))->GetStdDev();
            }
            alpha_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            f_qlcms->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
            }
        }

        // Now, having all systematics collected for this energy, calculate the total systematic uncertainties for each mT bin and Lévy parameter
        double alpha_syserr_up[NKT] = {0};
        double alpha_syserr_dn[NKT] = {0};
        double R_syserr_up[NKT] = {0};
        double R_syserr_dn[NKT] = {0};
        double N_syserr_up[NKT] = {0};
        double N_syserr_dn[NKT] = {0};
        for(int ikt=0; ikt<NKT; ikt++)
        {
            // Collect deviations from default for each systematic source
            // FIXME probably check the direction of deviations (up/down) later
            alpha_syserr_up[ikt] += pow(alpha_qlcms[ienergy][1]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // strict qLCMS
            alpha_syserr_dn[ikt] += pow(alpha_qlcms[ienergy][2]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // loose qLCMS
            alpha_syserr_up[ikt] += pow(alpha_rhofitmax[ienergy][1]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // strict rhofitmax
            alpha_syserr_dn[ikt] += pow(alpha_rhofitmax[ienergy][2]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // loose rhofitmax
            alpha_syserr_up[ikt] += pow(alpha_nevtavg[ienergy][1]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // smaller nevt_avg
            alpha_syserr_dn[ikt] += pow(alpha_nevtavg[ienergy][2]->GetY()[ikt] - alpha_default[ienergy]->GetY()[ikt], 2); // larger nevt_avg
            alpha_syserr_up[ikt] = sqrt(alpha_syserr_up[ikt]);
            alpha_syserr_dn[ikt] = sqrt(alpha_syserr_dn[ikt]);

            R_syserr_up[ikt] += pow(R_qlcms[ienergy][1]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // strict qLCMS
            R_syserr_dn[ikt] += pow(R_qlcms[ienergy][2]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // loose qLCMS
            R_syserr_up[ikt] += pow(R_rhofitmax[ienergy][1]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // strict rhofitmax
            R_syserr_dn[ikt] += pow(R_rhofitmax[ienergy][2]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // loose rhofitmax
            R_syserr_up[ikt] += pow(R_nevtavg[ienergy][1]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // smaller nevt_avg
            R_syserr_dn[ikt] += pow(R_nevtavg[ienergy][2]->GetY()[ikt] - R_default[ienergy]->GetY()[ikt], 2); // larger nevt_avg
            R_syserr_up[ikt] = sqrt(R_syserr_up[ikt]);
            R_syserr_dn[ikt] = sqrt(R_syserr_dn[ikt]);

            N_syserr_up[ikt] += pow(N_qlcms[ienergy][1]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // strict qLCMS
            N_syserr_dn[ikt] += pow(N_qlcms[ienergy][2]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // loose qLCMS
            N_syserr_up[ikt] += pow(N_rhofitmax[ienergy][1]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // strict rhofitmax
            N_syserr_dn[ikt] += pow(N_rhofitmax[ienergy][2]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // loose rhofitmax
            N_syserr_up[ikt] += pow(N_nevtavg[ienergy][1]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // smaller nevt_avg
            N_syserr_dn[ikt] += pow(N_nevtavg[ienergy][2]->GetY()[ikt] - N_default[ienergy]->GetY()[ikt], 2); // larger nevt_avg
            N_syserr_up[ikt] = sqrt(N_syserr_up[ikt]);
            N_syserr_dn[ikt] = sqrt(N_syserr_dn[ikt]);
        }
        // Create TGraphAsymmErrors for total systematic uncertainties
        alpha_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_default[ienergy]->GetY(), xerr_low, xerr_high, alpha_syserr_dn, alpha_syserr_up);
        R_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, R_default[ienergy]->GetY(), xerr_low, xerr_high, R_syserr_dn, R_syserr_up);
        N_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, N_default[ienergy]->GetY(), xerr_low, xerr_high, N_syserr_dn, N_syserr_up);
         
    } // end energy loop

    // NOTE: From now on, even though collected, I won't deal with "statistical" uncertainties here (StdDev)
    
    ////////////////////////////////////////////////////
    /////////////  PLOTTING   /////////////////////////
    ////////////////////////////////////////////////////
    
    // SYSTEMATICS TO mT vs PARAM

    // SYSTEMATICS TO sqrt(sNN) vs PARAM
    
    ////////////////////////////////////////////////////
    ////////  PLOTTING mT vs parameter with sys bands  
    ////////////////////////////////////////////////////
    auto parseEnergy = [](const char* s)->double{
        std::string str(s);
        std::size_t pos = str.find_first_not_of("0123456789.");
        std::string num = (pos==std::string::npos)?str:str.substr(0,pos);
        try { return std::stod(num); } catch(...) { return 0.0; }
    };

    double cutoff = 7.7; // GeV dividing low/high energies; treat 7.7 as LOW side

    // Parameter loop: 0=alpha,1=R,2=N
    // ensure output dir exists
    gSystem->mkdir("figs/syserr", kTRUE);

    for(int iparam=0; iparam<3; iparam++)
    {
        TCanvas* can = new TCanvas(Form("can_param_%d", iparam), "", 1400, 700);
        can->Divide(2,1);

        // Compute a global y-range for this parameter so left/right panels are directly comparable
        double global_ymin = 1e9, global_ymax = -1e9;
        for(int ie=0; ie<NENERGIES; ie++) {
            TGraphAsymmErrors* gdef_tmp = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
            TGraphAsymmErrors* gsys_tmp = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
            if(!gdef_tmp || !gsys_tmp) continue;
            for(int i=0;i<gdef_tmp->GetN();i++){
                double y = gdef_tmp->GetY()[i];
                double errup = gsys_tmp->GetEYhigh()[i];
                double errdn = gsys_tmp->GetEYlow()[i];
                global_ymin = std::min(global_ymin, y - errdn);
                global_ymax = std::max(global_ymax, y + errup);
            }
        }
        if(global_ymin > global_ymax){ global_ymin = 0.; global_ymax = 1.; }

        for(int side=0; side<2; side++) // 0: low, 1: high
        {
            can->cd(side+1);
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            // Use global y-range computed above
            double ymin = global_ymin;
            double ymax = global_ymax;
            double ypad = 0.12*(ymax - ymin);
            ymin -= ypad; ymax += ypad;

            // Draw an empty frame with axis titles
            double xmin = mtbin_centers[0];
            double xmax = mtbin_centers[NKT-1];
            TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
            const char* ytitle = (iparam==0)?"#alpha":(iparam==1)?"R [fm]":"N";
            frame->GetXaxis()->SetTitle("m_{T} (GeV/c^{2})");
            frame->GetYaxis()->SetTitle(ytitle);

            // Legend
            // Legend: smaller and moved to upper-left
            TLegend* leg = new TLegend(0.12,0.72,0.42,0.92);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.035);

            // Colors
            int colors[] = {kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+2,kOrange+7,kViolet+1,kCyan+1};

            int colorIdx=0;
            for(int ie=0; ie<NENERGIES; ie++)
            {
                double e = parseEnergy(energies[ie]);
                // treat 7.7 as LOW side: low = e <= cutoff, high = e > cutoff
                if((side==0 && e>cutoff) || (side==1 && e<=cutoff)) continue;
                // graphs
                TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
                TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
                if(!gdef || !gsys) continue;
                int col = colors[colorIdx % (sizeof(colors)/sizeof(int))];
                // filled band: use sys graph draw option 3
                int fillcol = TColor::GetColorTransparent(col, 0.35);
                gsys->SetFillColor(fillcol);
                gsys->SetFillStyle(1001);
                gsys->SetLineColor(col);
                gsys->SetLineWidth(1);
                gsys->Draw("3 same");
                // thin connecting line for central values
                gdef->SetLineColor(col);
                gdef->SetLineWidth(1);
                gdef->SetMarkerStyle(20);
                gdef->SetMarkerSize(0.7);
                gdef->Draw("LP same");

                std::stringstream label;
                label<<energies[ie];
                leg->AddEntry(gdef, label.str().c_str(), "l");
                colorIdx++;
            }
            leg->Draw();
        }

        // Save canvas
        can->SaveAs(Form("figs/syserr/mT_vs_param_%d.png", iparam));
        delete can;
    }

    // OUTPUT FILE to store results (keeps previous behavior)
    TFile* outfile = new TFile("syserr_results.root", "RECREATE");
    outfile->cd();
    for(int ie=0; ie<NENERGIES; ie++)
    {
        alpha_syserr[ie]->SetName(Form("alpha_syserr_%s", energies[ie]));
        R_syserr[ie]->SetName(Form("R_syserr_%s", energies[ie]));
        N_syserr[ie]->SetName(Form("N_syserr_%s", energies[ie]));
        alpha_syserr[ie]->Write();
        R_syserr[ie]->Write();
        N_syserr[ie]->Write();
    }
    outfile->Close();

    return 0;
}
