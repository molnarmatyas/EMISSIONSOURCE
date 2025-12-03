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

const char* levy_params[3] = {"alpha","R","N"};

// To be set
int NEVT_AVG_DEFAULT = 5; // index
int NEVT_AVGsyst[] = {10, 25, 50, 100, 200, 500, 1000, 5000, 10000}; // removed 1, not meaningful for low energies, runs too long

void correct_syserr_direction(double* param_uncert_up, double* param_uncert_dn, 
                              TGraphAsymmErrors* systcheckgraph, TGraphAsymmErrors* defaultgraph,
                              int ikt)
{
    // func(*param*_syserr_up[ikt], *param*_syserr_dn[ikt], *param*_*systcheck*[ienergy][1/2], *param*_default[ienergy])
    double diff = systcheckgraph->GetY()[ikt] - defaultgraph->GetY()[ikt];
    if(diff>=0.)
    {
        // systcheck datapoint above default
        param_uncert_up[ikt] += pow(diff, 2);
    }
    else
    {
        // systcheck datapoint below default
        param_uncert_dn[ikt] += pow(diff, 2);
    }
}


// -------- MAIN FUNCTION --------
int calc_and_plot_syserr(int energy_to_plot=-1)
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
            // Collect deviations from default for each systematic source - check the general direction of deviations (up/down)
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_qlcms[ienergy][1], alpha_default[ienergy], ikt); // strict qLCMS
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_qlcms[ienergy][2], alpha_default[ienergy], ikt); // loose qLCMS
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_rhofitmax[ienergy][1], alpha_default[ienergy], ikt); // strict rhofitmax
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_rhofitmax[ienergy][2], alpha_default[ienergy], ikt); // loose rhofitmax
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_nevtavg[ienergy][1], alpha_default[ienergy], ikt); // smaller nevt_avg
            correct_syserr_direction(alpha_syserr_up, alpha_syserr_dn, alpha_nevtavg[ienergy][2], alpha_default[ienergy], ikt); // larger nevt_avg
            alpha_syserr_up[ikt] = sqrt(alpha_syserr_up[ikt]);
            alpha_syserr_dn[ikt] = sqrt(alpha_syserr_dn[ikt]);
            // cerr << "alpha_syserr_up: " << alpha_syserr_up[ikt] << ", alpha_syserr_dn: " << alpha_syserr_dn[ikt] << endl;

            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_qlcms[ienergy][1], R_default[ienergy], ikt); // strict qLCMS
            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_qlcms[ienergy][2], R_default[ienergy], ikt); // loose qLCMS
            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_rhofitmax[ienergy][1], R_default[ienergy], ikt); // strict rhofitmax
            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_rhofitmax[ienergy][2], R_default[ienergy], ikt); // loose rhofitmax
            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_nevtavg[ienergy][1], R_default[ienergy], ikt); // smaller nevt_avg
            correct_syserr_direction(R_syserr_up, R_syserr_dn, R_nevtavg[ienergy][2], R_default[ienergy], ikt); // larger nevt_avg
            R_syserr_up[ikt] = sqrt(R_syserr_up[ikt]);
            R_syserr_dn[ikt] = sqrt(R_syserr_dn[ikt]);
            // cerr << "R_syserr_up: " << R_syserr_up[ikt] << ", R_syserr_dn: " << R_syserr_dn[ikt] << endl;

            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_qlcms[ienergy][1], N_default[ienergy], ikt); // strict qLCMS
            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_qlcms[ienergy][2], N_default[ienergy], ikt); // loose qLCMS
            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_rhofitmax[ienergy][1], N_default[ienergy], ikt); // strict rhofitmax
            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_rhofitmax[ienergy][2], N_default[ienergy], ikt); // loose rhofitmax
            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_nevtavg[ienergy][1], N_default[ienergy], ikt); // smaller nevt_avg
            correct_syserr_direction(N_syserr_up, N_syserr_dn, N_nevtavg[ienergy][2], N_default[ienergy], ikt); // larger nevt_avg
            N_syserr_up[ikt] = sqrt(N_syserr_up[ikt]);
            N_syserr_dn[ikt] = sqrt(N_syserr_dn[ikt]);
            // cerr << "N_syserr_up: " << N_syserr_up[ikt] << ", N_syserr_dn: " << N_syserr_dn[ikt] << endl;
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
    
    
    ////////////////////////////////////////////////////
    ////////  PLOTTING mT vs parameter with sys bands  
    ////////////////////////////////////////////////////
    auto parseEnergy = [](const char* s)->double{
        std::string str(s);
        std::size_t pos = str.find_first_not_of("0123456789.");
        std::string num = (pos==std::string::npos)?str:str.substr(0,pos);
        try { return std::stod(num); } catch(...) { return 0.0; }
    };

    // Helper: Find the least-crowded quadrant for legend placement (NDC coords)
    auto findBestLegendPos = [&parseEnergy](const std::vector<TGraphAsymmErrors*>& graphs, double xmin, double xmax, double ymin, double ymax,
                                             int iparam, double cutoff, int NENERGIES, int side) 
        -> std::tuple<double,double,double,double> {
        double xmid = 0.5*(xmin + xmax);
        double ymid = 0.5*(ymin + ymax);
        int count[4] = {0,0,0,0}; // TL,TR,BL,BR
        
        for(int ie=0; ie<NENERGIES; ie++) {
            double e = parseEnergy(energies[ie]);
            if((side==0 && e>cutoff) || (side==1 && e<=cutoff)) continue;
            TGraphAsymmErrors* g = graphs[ie];
            if(!g) continue;
            for(int i=0; i<g->GetN(); i++) {
                double x = g->GetX()[i];
                double y = g->GetY()[i];
                int quad = (x < xmid ? 0 : 1) + (y > ymid ? 0 : 2); // 0=TL,1=TR,2=BL,3=BR
                count[quad]++;
            }
        }
        // Find least-crowded quadrant
        int minQuad = 0, minCount = count[0];
        for(int q=1; q<4; q++) {
            if(count[q] < minCount) { minQuad = q; minCount = count[q]; }
        }
        // Return legend coords: width=0.20, height=0.20
        double legW = 0.20, legH = 0.20;
        double x1, y2;
        if(minQuad == 0) { x1 = 0.15; y2 = 0.72; } // TL
        else if(minQuad == 1) { x1 = 0.65; y2 = 0.72; } // TR
        else if(minQuad == 2) { x1 = 0.15; y2 = 0.25; } // BL
        else { x1 = 0.65; y2 = 0.25; } // BR
        return std::make_tuple(x1, y2, x1+legW, y2+legH);
    };

    double cutoff = 4.5; // GeV dividing low/high energies; treat 7.7 as LOW side

    // Parameter loop: 0=alpha,1=R,2=N
    // ensure output dir exists
    gSystem->mkdir("figs/syserr", kTRUE);

    for(int iparam=0; iparam<3; iparam++)
    {
        TCanvas* can = new TCanvas(Form("can_param_%d", iparam), "", 2800, 1400);
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
            double ymins[] = {1.3, 3, 0.97};
            double ymaxs[] = {2.1, 8., 1.13};
            double ymin = ymins[iparam];//global_ymin;
            double ymax = ymaxs[iparam];//global_ymax;
            double ypad = 0.12*(ymax - ymin);
            //ymin -= ypad; ymax += ypad;

            // Draw an empty frame with axis titles
            double xmin = mtbin_centers[0]-0.05;
            double xmax = mtbin_centers[NKT-1]+0.05;
            TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
            const char* ytitle = (iparam==0)?"#alpha":(iparam==1)?"R [fm]":"N";
            frame->GetXaxis()->SetTitle("m_{T} (GeV/c^{2})");
            frame->GetYaxis()->SetTitle(ytitle);

            // Legend: auto-position in least-crowded quadrant
            auto [leg_x1, leg_y2, leg_x2, leg_y1] = findBestLegendPos(
                (iparam==0)? std::vector<TGraphAsymmErrors*>(alpha_default, alpha_default+NENERGIES) :
                (iparam==1)? std::vector<TGraphAsymmErrors*>(R_default, R_default+NENERGIES) :
                             std::vector<TGraphAsymmErrors*>(N_default, N_default+NENERGIES),
                xmin, xmax, ymin, ymax, iparam, cutoff, NENERGIES, side);
            TLegend* leg = new TLegend(leg_x1, leg_y2, leg_x2, leg_y1);
            leg->SetFillColor(0);
            leg->SetBorderSize(1); // set to 0 to remove border
            leg->SetTextSize(0.023);

            // Colors
            int colors[] = {kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+2,kOrange+7,kViolet+1,kCyan+1};

            int colorIdx=0;
            for(int ie=0; ie<NENERGIES; ie++)
            {
                if(energy_to_plot>-1 && ie!=energy_to_plot) continue; // plot only selected energy if specified
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
                // gdef->SetLineColor(col);
                // gdef->SetLineWidth(3);
                // gdef->SetMarkerStyle(20);
                // gdef->SetMarkerSize(0.0);
                // gdef->Draw("LP same");
                gsys->SetLineColor(col);
                gsys->SetLineWidth(3);
                gsys->SetMarkerStyle(20);
                gsys->SetMarkerSize(0.0);
                gsys->Draw("LPX same");
                //gStyle->SetLegendTextSize(0.015);

                std::stringstream label;
                label<<energies[ie];
                // leg->AddEntry(gdef, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "l");
                leg->AddEntry(gsys, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "l");
                colorIdx++;
            }
            leg->Draw();
        }

        // Save canvas
        const char* energyplotted = (energy_to_plot>-1)? Form("_energy_%s", energies[energy_to_plot]) : "_allenergies";
        can->SetTitle(Form("m_{T} vs %s systematic uncertainties%s",
                            (iparam==0)?"#alpha":(iparam==1)?"R":"N", energyplotted));
        //can->SaveAs(Form("figs/syserr/mT_vs_param_%d%s.png", iparam, energyplotted));
        can->SaveAs(Form("figs/syserr/mT_vs_param_%s%s.png", levy_params[iparam], energyplotted));
        delete can;
    }

    
    // SYSTEMATICS TO sqrt(sNN) vs PARAM
    

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
