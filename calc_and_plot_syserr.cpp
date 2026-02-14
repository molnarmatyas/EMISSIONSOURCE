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
#include <fstream>
#include <iomanip>
#include <TF1.h>

const char* levy_params[3] = {"alpha","R","N"};

// To be set
//int NEVT_AVG_DEFAULT = 5; // single index

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

// I do not know if this makes any sense, this will return the average syst. difference in percent
double calc_syst_contribution(TGraphAsymmErrors* defaultgraph, double* errups, double* errdns, bool up=true)
{
    double percent = 0.;
    int points = 0;
    for(int ikt=0; ikt<NKT; ikt++)
    {
        double defy = defaultgraph->GetY()[ikt];
        if(defy < 1e-12) continue; // skip tiny/invalid points
        if(up)
        {
            // errups is an absolute upward uncertainty (magnitude)
            double mag = errups ? errups[ikt] : 0.;
            percent += mag / defy;
            points++;
        }
        else
        {
            // errdns is an absolute downward uncertainty (magnitude)
            double mag = errdns ? errdns[ikt] : 0.;
            percent += mag / defy;
            points++;
        }
    }
    percent /= points;
    percent *= 100.; // to percent
    return percent;
}

// Return median percent contribution (same inputs as calc_syst_contribution)
double calc_syst_contrib_median(TGraphAsymmErrors* defaultgraph, double* errups, double* errdns, bool up=true)
{
    if(!defaultgraph) return 0.;
    std::vector<double> vals;
    vals.reserve(NKT);
    for(int ikt=0; ikt<NKT; ikt++){
        double defy = defaultgraph->GetY()[ikt];
        if(defy < 1e-12) continue;
        double mag = 0.;
        if(up){ mag = errups ? errups[ikt] : 0.; }
        else { mag = errdns ? errdns[ikt] : 0.; }
        // store percent for this bin
        vals.push_back( (mag / defy) * 100. );
    }
    if(vals.empty()) return 0.;
    std::sort(vals.begin(), vals.end());
    size_t n = vals.size();
    if(n % 2 == 1) return vals[n/2];
    // even -> average two middle
    return 0.5 * (vals[n/2 - 1] + vals[n/2]);
}

Double_t myfunc(double x){return 0.85+pow(x,-0.14);}


// -------- MAIN FUNCTION --------
int calc_and_plot_syserr(int energy_to_plot=-1)
{
    // Get m_T bin centers
    double mtbin_centers[NKT];
    for(int ii=0; ii<NKT; ii++)
    {
        //mtbin_centers[ii] = 0.5 * (ktbins[ii] + ktbins[ii+1]);
        mtbin_centers[ii] = kT_center[ii];
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

    // CSV lines to collect percent breakdown per energy
    std::vector<std::string> csv_lines;

    // Arrays to store per-energy percent contributions (median) for reuse in sqrt(sNN) plots
    double alpha_qlcms_pct_up_arr[NENERGIES] = {0}, alpha_qlcms_pct_dn_arr[NENERGIES] = {0};
    double alpha_rhofit_pct_up_arr[NENERGIES] = {0}, alpha_rhofit_pct_dn_arr[NENERGIES] = {0};
    double alpha_nevt_pct_up_arr[NENERGIES] = {0}, alpha_nevt_pct_dn_arr[NENERGIES] = {0};
    double R_qlcms_pct_up_arr[NENERGIES] = {0}, R_qlcms_pct_dn_arr[NENERGIES] = {0};
    double R_rhofit_pct_up_arr[NENERGIES] = {0}, R_rhofit_pct_dn_arr[NENERGIES] = {0};
    double R_nevt_pct_up_arr[NENERGIES] = {0}, R_nevt_pct_dn_arr[NENERGIES] = {0};
    double N_qlcms_pct_up_arr[NENERGIES] = {0}, N_qlcms_pct_dn_arr[NENERGIES] = {0};
    double N_rhofit_pct_up_arr[NENERGIES] = {0}, N_rhofit_pct_dn_arr[NENERGIES] = {0};
    double N_nevt_pct_up_arr[NENERGIES] = {0}, N_nevt_pct_dn_arr[NENERGIES] = {0};

    // Also store mT-averaged default values and the ikt=2 default points per energy
    double avg_alpha[NENERGIES] = {0}, avg_R[NENERGIES] = {0}, avg_N[NENERGIES] = {0};
    double val_alpha_ikt2[NENERGIES] = {0}, val_R_ikt2[NENERGIES] = {0}, val_N_ikt2[NENERGIES] = {0};
    // And the m_T choice syst uncert
    double alpha_mTchoice_syst_up[NENERGIES] = {0}, alpha_mTchoice_syst_dn[NENERGIES] = {0};
    double R_mTchoice_syst_up[NENERGIES] = {0}, R_mTchoice_syst_dn[NENERGIES] = {0};
    double N_mTchoice_syst_up[NENERGIES] = {0}, N_mTchoice_syst_dn[NENERGIES] = {0};

    for(int ienergy = 0; ienergy < NENERGIES; ienergy++)
    {
        // Default
        TFile* f_default = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root", 
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]]), "READ");
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
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]], qlcms_labels[iq]), "READ");
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
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]], rhofitmax_labels[ir]), "READ");
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
            if(in==0) index = NEVT_AVG_DEFAULT[ienergy]; // default
            else if(in==1) index = NEVT_AVG_DEFAULT[ienergy] - 1; // reasonably smaller
            else if(in==2) index = NEVT_AVG_DEFAULT[ienergy] + 1; // reasonably larger
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

        // Give one-one number about the systematic uncertainty contributions from each source
        /*
        // Average
        double alpha_qlcms_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double alpha_rhofitmax_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double alpha_nevtavg_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double R_qlcms_percent_up = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double R_rhofitmax_percent_up = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double R_nevtavg_percent_up = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double N_qlcms_percent_up = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn);
        double N_rhofitmax_percent_up = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn);
        double N_nevtavg_percent_up = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn);

        double alpha_qlcms_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double alpha_rhofitmax_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double alpha_nevtavg_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double R_qlcms_percent_dn = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double R_rhofitmax_percent_dn = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double R_nevtavg_percent_dn = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double N_qlcms_percent_dn = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        double N_rhofitmax_percent_dn = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        double N_nevtavg_percent_dn = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        */
        // Median, might make more sense
        double alpha_qlcms_percent_up = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double alpha_rhofitmax_percent_up = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double alpha_nevtavg_percent_up = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double R_qlcms_percent_up = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double R_rhofitmax_percent_up = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double R_nevtavg_percent_up = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double N_qlcms_percent_up = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn);
        double N_rhofitmax_percent_up = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn);
        double N_nevtavg_percent_up = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn);

        double alpha_qlcms_percent_dn = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double alpha_rhofitmax_percent_dn = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double alpha_nevtavg_percent_dn = calc_syst_contrib_median(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double R_qlcms_percent_dn = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double R_rhofitmax_percent_dn = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double R_nevtavg_percent_dn = calc_syst_contrib_median(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double N_qlcms_percent_dn = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        double N_rhofitmax_percent_dn = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        double N_nevtavg_percent_dn = calc_syst_contrib_median(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);

        // Append CSV block lines: energy header, parameter sub-header, then lines for alpha,R,N
        std::ostringstream h1;
        h1 << "[" << energies[ienergy] << " GeV]";
        csv_lines.push_back(h1.str());

        csv_lines.push_back("Levy parameter,qlcms cut [%],rhofitmax [%],nevtavg [%]");
        std::ostringstream arow;
        arow << "alpha," << std::fixed << std::setprecision(3)
            << " +" << alpha_qlcms_percent_up << " / -" << alpha_qlcms_percent_dn << ","
            << " +" << alpha_rhofitmax_percent_up << " / -" << alpha_rhofitmax_percent_dn << ","
            << " +" << alpha_nevtavg_percent_up << " / -" << alpha_nevtavg_percent_dn;
        csv_lines.push_back(arow.str());
        std::ostringstream rrow;
        rrow << "R," << std::fixed << std::setprecision(3)
            << " +" << R_qlcms_percent_up << " / -" << R_qlcms_percent_dn << ","
            << " +" << R_rhofitmax_percent_up << " / -" << R_rhofitmax_percent_dn << ","
            << " +" << R_nevtavg_percent_up << " / -" << R_nevtavg_percent_dn;
        csv_lines.push_back(rrow.str());
        std::ostringstream nrow;
        nrow << "N," << std::fixed << std::setprecision(3)
            << " +" << N_qlcms_percent_up << " / -" << N_qlcms_percent_dn << ","
            << " +" << N_rhofitmax_percent_up << " / -" << N_rhofitmax_percent_dn << ","
            << " +" << N_nevtavg_percent_up << " / -" << N_nevtavg_percent_dn;
        csv_lines.push_back(nrow.str());
        
        // ---- m_T choice systematic contribution (percent) ----
        // Add separate energy header for m_T choice block
        std::ostringstream h_mt;
        h_mt << "[" << energies[ienergy] << " GeV]";
        csv_lines.push_back(h_mt.str());
        
        // compute mT center and neighbor differences in MeV
        double default_mt_mev = mtbin_centers[2] * 1000.0;
        double upper_diff_mev = (mtbin_centers[3] - mtbin_centers[2]) * 1000.0;
        double lower_diff_mev = (mtbin_centers[2] - mtbin_centers[1]) * 1000.0;
        double mt_pct_up = 0.0, mt_pct_dn = 0.0;
        TGraphAsymmErrors* gdef_mt = alpha_default[ienergy];
        if(gdef_mt && gdef_mt->GetN() > 2){
            double center = gdef_mt->GetY()[2];
            double neigh_up_sq = 0., neigh_dn_sq = 0.;
            for(int j : {1,3}){
                if(j<0 || j>=gdef_mt->GetN()) continue;
                double nb = gdef_mt->GetY()[j];
                double diff = nb - center;
                if(diff >= 0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
            }
            double neigh_up = sqrt(neigh_up_sq);
            double neigh_dn = sqrt(neigh_dn_sq);
            if(std::isfinite(center) && center != 0.){
                mt_pct_up = 100.0 * neigh_up / fabs(center);
                mt_pct_dn = 100.0 * neigh_dn / fabs(center);
            }
        }
        
        // Subheader with m_T details
        std::ostringstream mtsub;
        mtsub << "m_T choice (m_T=" << std::fixed << std::setprecision(0) << default_mt_mev
              << " MeV +" << std::fixed << std::setprecision(0) << upper_diff_mev
              << "/-" << std::fixed << std::setprecision(0) << lower_diff_mev << " MeV)";
        csv_lines.push_back(mtsub.str());
        
        // Data line with energy and percentages
        std::ostringstream mtline;
        mtline << energies[ienergy] << "," << std::fixed << std::setprecision(3)
               << " +" << mt_pct_up << " / -" << mt_pct_dn;
        csv_lines.push_back(mtline.str());
        
        // Store per-energy percent medians for later sqrt(sNN) plotting
        alpha_qlcms_pct_up_arr[ienergy] = alpha_qlcms_percent_up;
        alpha_qlcms_pct_dn_arr[ienergy] = alpha_qlcms_percent_dn;
        alpha_rhofit_pct_up_arr[ienergy] = alpha_rhofitmax_percent_up;
        alpha_rhofit_pct_dn_arr[ienergy] = alpha_rhofitmax_percent_dn;
        alpha_nevt_pct_up_arr[ienergy] = alpha_nevtavg_percent_up;
        alpha_nevt_pct_dn_arr[ienergy] = alpha_nevtavg_percent_dn;

        R_qlcms_pct_up_arr[ienergy] = R_qlcms_percent_up;
        R_qlcms_pct_dn_arr[ienergy] = R_qlcms_percent_dn;
        R_rhofit_pct_up_arr[ienergy] = R_rhofitmax_percent_up;
        R_rhofit_pct_dn_arr[ienergy] = R_rhofitmax_percent_dn;
        R_nevt_pct_up_arr[ienergy] = R_nevtavg_percent_up;
        R_nevt_pct_dn_arr[ienergy] = R_nevtavg_percent_dn;

        N_qlcms_pct_up_arr[ienergy] = N_qlcms_percent_up;
        N_qlcms_pct_dn_arr[ienergy] = N_qlcms_percent_dn;
        N_rhofit_pct_up_arr[ienergy] = N_rhofitmax_percent_up;
        N_rhofit_pct_dn_arr[ienergy] = N_rhofitmax_percent_dn;
        N_nevt_pct_up_arr[ienergy] = N_nevtavg_percent_up;
        N_nevt_pct_dn_arr[ienergy] = N_nevtavg_percent_dn;

        // mT-averaged defaults and ikt=2 values
        double sum_a=0., sum_R=0., sum_Nv=0.;
        int cnt=0;
        for(int ikt=0; ikt<NKT; ikt++){
            double aya = alpha_default[ienergy]->GetY()[ikt];
            double rya = R_default[ienergy]->GetY()[ikt];
            double nya = N_default[ienergy]->GetY()[ikt];
            if(!std::isfinite(aya) || !std::isfinite(rya) || !std::isfinite(nya)) continue;
            sum_a += aya; sum_R += rya; sum_Nv += nya; cnt++;
        }
        if(cnt>0){ avg_alpha[ienergy] = sum_a / cnt; avg_R[ienergy] = sum_R / cnt; avg_N[ienergy] = sum_Nv / cnt; }
        // store ikt=2 (3rd mT bin) default values if available
        if(NKT>2){ val_alpha_ikt2[ienergy] = alpha_default[ienergy]->GetY()[2]; val_R_ikt2[ienergy] = R_default[ienergy]->GetY()[2]; val_N_ikt2[ienergy] = N_default[ienergy]->GetY()[2]; }
        
         
    } // end energy loop

    // NOTE: From now on, even though collected, I won't deal with "statistical" uncertainties here (StdDev)
    
    ////////////////////////////////////////////////////
    /////////////  PLOTTING   /////////////////////////
    ////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////
    ////////  PLOTTING mT vs parameter with sys bands  
    ////////////////////////////////////////////////////
    // Use `energydouble[]` from header_for_all_emissionsource.h for numeric energies

    // Helper: Find the least-crowded quadrant for legend placement (NDC coords)
    // Uses `energydouble[]` from header_for_all_emissionsource.h
    auto findBestLegendPos = [](const std::vector<TGraphAsymmErrors*>& graphs, double xmin, double xmax, double ymin, double ymax,
                                             int iparam, double cutoff, int NENERGIES, int side) 
        -> std::tuple<double,double,double,double> {
        double xmid = 0.5*(xmin + xmax);
        double ymid = 0.5*(ymin + ymax);
        int count[4] = {0,0,0,0}; // TL,TR,BL,BR
        
        for(int ie=0; ie<NENERGIES; ie++) {
            double e = energydouble[ie];
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
        // Return an anchor box near the chosen quadrant; the caller will size it.
        // We return a default small box (caller will override width/height).
        double anchorX, anchorYlow; // lower-left corner
        if(minQuad == 0) { anchorX = 0.15; anchorYlow = 0.72; } // TL
        else if(minQuad == 1) { anchorX = 0.62; anchorYlow = 0.72; } // TR (slightly more left to allow wider legends)
        else if(minQuad == 2) { anchorX = 0.15; anchorYlow = 0.25; } // BL
        else { anchorX = 0.62; anchorYlow = 0.25; } // BR
        // default tiny box; caller will compute final x2,y2
        return std::make_tuple(anchorX, anchorYlow, anchorX+0.10, anchorYlow+0.10);
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
            double ymins[] = {1.3, 3, 0.97}; // {0.5, 0, 0.5};
            double ymaxs[] = {2.1, 8., 1.13}; // {2., 10., 1.2};
            double ymin = ymins[iparam];//global_ymin;
            double ymax = ymaxs[iparam];//global_ymax;
            double ypad = 0.12*(ymax - ymin);
            //ymin -= ypad; ymax += ypad;

            // Draw an empty frame with axis titles
            double xmin = mtbin_centers[0]-0.05;
            double xmax = mtbin_centers[NKT-1]+0.05;
            TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
            const char* ytitle = (iparam==0)?"#alpha":(iparam==1)?"R [fm]":"#lambda"; // N
            frame->GetXaxis()->SetTitle("m_{T} (GeV/c^{2})");
            frame->GetYaxis()->SetTitle(ytitle);

            // Legend: auto-position in least-crowded quadrant, with dynamic sizing
            auto [leg_x1, leg_y2, leg_x2, leg_y1] = findBestLegendPos(
                (iparam==0)? std::vector<TGraphAsymmErrors*>(alpha_default, alpha_default+NENERGIES) :
                (iparam==1)? std::vector<TGraphAsymmErrors*>(R_default, R_default+NENERGIES) :
                         std::vector<TGraphAsymmErrors*>(N_default, N_default+NENERGIES),
                xmin, xmax, ymin, ymax, iparam, cutoff, NENERGIES, side);
            // Count entries and choose columns
            int nEntries = 0;
            for(int ie=0; ie<NENERGIES; ie++){
                if(energy_to_plot>-1 && ie!=energy_to_plot) continue;
                double e = energydouble[ie];
                if((side==0 && e>cutoff) || (side==1 && e<=cutoff)) continue;
                TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
                TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
                if(!gdef || !gsys) continue;
                nEntries++;
            }
            int nCols = (nEntries>14)? 3 : (nEntries>7)? 2 : 1;
            int nRows = (nEntries==0)? 0 : ( (nEntries + nCols - 1) / nCols );
            // Derive legend height/width in NDC
            double lineH = 0.038; // per-entry height
            double vPad = 0.012;  // top/bottom padding
            double legH = std::max(0.06, nRows*lineH + 2*vPad);
            double colW = (nCols==1)? 0.24 : (nCols==2)? 0.42 : 0.58; // per total width for 1/2/3 cols
            double legW = colW;
            // Clamp to pad bounds
            double maxX2 = 0.98, maxY2 = 0.96;
            double x1 = leg_x1;
            double y1 = leg_y2; // lower y
            double x2 = x1 + legW;
            double y2 = y1 + legH;
            if(x2 > maxX2){ x1 = std::max(0.10, maxX2 - legW); x2 = x1 + legW; }
            if(y2 > maxY2){ y1 = std::max(0.12, maxY2 - legH); y2 = y1 + legH; }
            TLegend* leg = new TLegend(x1, y1, x2, y2);
            leg->SetNColumns(nCols);
            leg->SetFillColor(0);
            leg->SetBorderSize(1);
            leg->SetTextSize(0.028);
            leg->SetMargin(0.20);
            leg->SetEntrySeparation(0.005);
            leg->SetColumnSeparation(0.02);

            // Colors
            int colors[] = {kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+2,kOrange+7,kViolet+1,kCyan+1};

            int colorIdx=0;
            for(int ie=0; ie<NENERGIES; ie++)
            {
                if(energy_to_plot>-1 && ie!=energy_to_plot) continue; // plot only selected energy if specified
                double e = energydouble[ie];
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
                gsys->SetMarkerStyle(20+colorIdx); // differentiate points
                gsys->SetMarkerSize(2.0); // instead of 0.0, points needed for better visibility in case of mT vs param
                gsys->SetMarkerColor(col);
                gsys->Draw("LPX same");
                //gStyle->SetLegendTextSize(0.015);

                std::stringstream label;
                //label.precision(2);
                label << std::fixed << std::setprecision(1) << energydouble[ie];//energies[ie];
                // leg->AddEntry(gdef, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "l");
                leg->AddEntry(gsys, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "lpx"); // changed from "l", let the marker show in legend
                colorIdx++;
            }
            leg->Draw();
        }

        // Save canvas
        const char* energyplotted = (energy_to_plot>-1)? Form("_energy_%s", energies[energy_to_plot]) : "_allenergies";
        can->SetTitle(Form("m_{T} vs %s systematic uncertainties%s",
                            (iparam==0)?"#alpha":(iparam==1)?"R":"#lambda", energyplotted));
        //can->SaveAs(Form("figs/syserr/mT_vs_param_%d%s.png", iparam, energyplotted));
        can->SaveAs(Form("figs/syserr/mT_vs_param_%s%s.png", levy_params[iparam], energyplotted));
        delete can;
    }

    
    /////////////////////////////////////////////////////
    // SYSTEMATICS TO sqrt(s_NN) vs Lévy parameter /////
    ///////////////////////////////////////////////////
    
    // Prepare x-axis (collision energies as numbers)
    double x_energy[NENERGIES];
    for(int ie=0; ie<NENERGIES; ie++) x_energy[ie] = energydouble[ie];
    double xerr_low_s[NENERGIES] = {0};
    double xerr_high_s[NENERGIES] = {0};

    // For each parameter make a two-panel canvas: left = mT-averaged, right = ikt=2
    for(int iparam=0; iparam<3; iparam++){
        TCanvas* can = new TCanvas(Form("can_sqrtS_param_%d", iparam), "", 1600, 800);
        can->Divide(2,1);

        // If alpha analytic curve will be drawn, precompute its y-range over the plotted x-range
        double func_min = 1e99, func_max = -1e99;
        if(iparam==0){
            double fxmin = x_energy[0];
            double fxmax = x_energy[NENERGIES-1];
            int Nsample = 500;
            for(int is=0; is<=Nsample; ++is){
                double x = fxmin + (fxmax - fxmin) * (double)is / (double)Nsample;
                double y = 0.85 + pow(x, -0.14);
                if(std::isfinite(y)){
                    func_min = std::min(func_min, y);
                    func_max = std::max(func_max, y);
                }
            }
        }

        // Build arrays for left panel (mT-averaged)
        double yavg[NENERGIES];
        double yavg_errup[NENERGIES];
        double yavg_errdn[NENERGIES];
        for(int ie=0; ie<NENERGIES; ie++){
            double avg = (iparam==0)? avg_alpha[ie] : (iparam==1)? avg_R[ie] : avg_N[ie];
            yavg[ie] = avg;
            // compute sys contributions from stored percent medians
            double up_sq = 0., dn_sq = 0.;
            if(iparam==0){
                up_sq += pow(alpha_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            } else if(iparam==1){
                up_sq += pow(R_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            } else {
                up_sq += pow(N_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            }
            yavg_errup[ie] = sqrt(up_sq);
            yavg_errdn[ie] = sqrt(dn_sq);
        }

        // Left pad: mT-averaged
        can->cd(1);
        //gPad->SetLogx(1); // !!! FIXME uncomment if comparison with STAR not needed
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        double xmin = x_energy[0] - 1.0;
        double xmax = x_energy[NENERGIES-1] + 1.0;
        double ymin = 1e9, ymax = -1e9;
        for(int ie=0; ie<NENERGIES; ie++){
            if(!std::isfinite(yavg[ie])) continue;
            ymin = std::min(ymin, yavg[ie] - yavg_errdn[ie]);
            ymax = std::max(ymax, yavg[ie] + yavg_errup[ie]);
        }
        if(ymin>ymax){ ymin = 0.; ymax = 1.; }
        // Ensure analytic function (if computed) fits into the pad range
        if(iparam==0 && std::isfinite(func_min) && std::isfinite(func_max)){
            ymin = std::min(ymin, func_min);
            ymax = std::max(ymax, func_max);
            double yrange = ymax - ymin;
            if(yrange==0) yrange = fabs(ymax)>0? fabs(ymax)*0.1 : 1.0;
            double padMargin = 0.06 * yrange;
            ymin -= padMargin;
            ymax += padMargin;
        }
        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
        if(iparam==0)
        {
            //frame = gPad->DrawFrame(1., 0.1, 10000., 2.1); // !!! FIXME uncomment if comparison with STAR not needed
        }
        frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
        const char* ytitle = (iparam==0)?"#alpha":(iparam==1)?"R [fm]":"#lambda";
        frame->GetYaxis()->SetTitle(ytitle);

        // sys band graph
        TGraphAsymmErrors* gavg = new TGraphAsymmErrors(NENERGIES, x_energy, yavg, xerr_low_s, xerr_high_s, yavg_errdn, yavg_errup);
        int col = kBlue+2;
        int fillcol = TColor::GetColorTransparent(col, 0.35);
        gavg->SetFillColor(fillcol);
        gavg->SetFillStyle(1001);
        gavg->SetLineColor(col);
        gavg->SetLineWidth(3);
        gavg->Draw("3 same");
        gavg->SetMarkerStyle(21);
        gavg->SetMarkerSize(0.0); //0.8
        gavg->Draw("LPX same");

        // Right pad: ikt=2
        can->cd(2);
        //gPad->SetLogx(1); // !!! FIXME uncomment if comparison with STAR not needed
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        double ymin2 = 1e9, ymax2 = -1e9;
        double yvals[NENERGIES];
        double yerrup2[NENERGIES];
        double yerrdn2[NENERGIES];
        for(int ie=0; ie<NENERGIES; ie++){
            double val = (iparam==0)? val_alpha_ikt2[ie] : (iparam==1)? val_R_ikt2[ie] : val_N_ikt2[ie];
            yvals[ie] = val;
            // base sys from alpha_syserr (absolute) for ikt=2
            TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
            double base_up = 0., base_dn = 0.;
            if(gsys && gsys->GetN() > 2){ base_up = gsys->GetEYhigh()[2]; base_dn = gsys->GetEYlow()[2]; }
            // neighbor diffs
            double neigh_up_sq = 0., neigh_dn_sq = 0.;
            TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
            if(gdef && gdef->GetN()>2){
                double center = gdef->GetY()[2];
                for(int j : {1,3}){
                    if(j<0 || j>=gdef->GetN()) continue;
                    double nb = gdef->GetY()[j];
                    double diff = nb - center;
                    if(diff>=0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
                }
            }
            double neigh_up = sqrt(neigh_up_sq);
            double neigh_dn = sqrt(neigh_dn_sq);
            // Save m_T choice systematics for further use
            if(iparam==0)
            {
                alpha_mTchoice_syst_up[ie] = neigh_up; alpha_mTchoice_syst_dn[ie] = neigh_dn;
            }else if(iparam==1)
            {
                R_mTchoice_syst_up[ie] = neigh_up; R_mTchoice_syst_dn[ie] = neigh_dn;
            }else
            {
                N_mTchoice_syst_up[ie] = neigh_up; N_mTchoice_syst_dn[ie] = neigh_dn;
            }

            yerrup2[ie] = sqrt(base_up*base_up + neigh_up*neigh_up);
            yerrdn2[ie] = sqrt(base_dn*base_dn + neigh_dn*neigh_dn);
            if(std::isfinite(yvals[ie])){
                ymin2 = std::min(ymin2, yvals[ie] - yerrdn2[ie]);
                ymax2 = std::max(ymax2, yvals[ie] + yerrup2[ie]);
            }
        }
        if(ymin2>ymax2){ ymin2=0.; ymax2=1.; }
        // Ensure analytic function (if computed) fits into the right pad range as well
        if(iparam==0 && std::isfinite(func_min) && std::isfinite(func_max)){
            ymin2 = std::min(ymin2, func_min);
            ymax2 = std::max(ymax2, func_max);
            double yrange2 = ymax2 - ymin2;
            if(yrange2==0) yrange2 = fabs(ymax2)>0? fabs(ymax2)*0.1 : 1.0;
            double padMargin2 = 0.06 * yrange2;
            ymin2 -= padMargin2;
            ymax2 += padMargin2;
        }
        TH1F* frame2 = gPad->DrawFrame(xmin, ymin2, xmax, ymax2);
        if(iparam==0)
        {
            //frame2 = gPad->DrawFrame(1., 0.1, 10000., 2.1); // !!! FIXME uncomment if comparison with STAR not needed
        }
        frame2->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
        frame2->GetYaxis()->SetTitle(ytitle);
        TGraphAsymmErrors* gikt2 = new TGraphAsymmErrors(NENERGIES, x_energy, yvals, xerr_low_s, xerr_high_s, yerrdn2, yerrup2);
        int col2 = kRed+1;
        int fillcol2 = TColor::GetColorTransparent(col2, 0.35);
        gikt2->SetFillColor(fillcol2);
        gikt2->SetFillStyle(1001);
        gikt2->SetLineColor(col2);
        gikt2->SetLineWidth(3);
        gikt2->Draw("3 same");
        gikt2->SetMarkerStyle(21);
        gikt2->SetMarkerSize(0.0); //0.8
        gikt2->Draw("LPX same");

        // If alpha, draw analytic curve y = 0.85 + x^-0.14 (make thicker and draw as line on top)
        TF1* f = nullptr;
        if(iparam==0){
            double xminf = x_energy[0];
            double xmaxf = x_energy[NENERGIES-1];
            f = new TF1("alpha_curve", "myfunc(x)", xminf, xmaxf);
            f->SetLineColor(TColor::GetColorTransparent(kBlack,0.4));
            f->SetLineStyle(1);
            f->SetLineWidth(3);
            f->SetNpx(500);
            // draw on left pad
            can->cd(1);
            f->Draw("L same");
            // draw on right pad
            can->cd(2);
            f->Draw("L same");
        }

        // Smart legend placement helper (counts points by quadrant)
        auto chooseLegendPos = [](const double xs[], const double ys[], int n, double xmin, double xmax, double ymin, double ymax)
            -> std::tuple<double,double,double,double> {
            double xmid = 0.5*(xmin + xmax);
            double ymid = 0.5*(ymin + ymax);
            int count[4] = {0,0,0,0}; // TL,TR,BL,BR
            for(int i=0;i<n;i++){
                double x = xs[i];
                double y = ys[i];
                if(!std::isfinite(x) || !std::isfinite(y)) continue;
                int quad = (x < xmid ? 0 : 1) + (y > ymid ? 0 : 2); // 0=TL,1=TR,2=BL,3=BR
                count[quad]++;
            }
            int minQuad = 0, minCount = count[0];
            for(int q=1;q<4;q++) if(count[q] < minCount){ minQuad = q; minCount = count[q]; }
            double legW = 0.20, legH = 0.12;
            double x1,y2;
            if(minQuad==0){ x1=0.15; y2=0.78; }
            else if(minQuad==1){ x1=0.65; y2=0.78; }
            else if(minQuad==2){ x1=0.15; y2=0.35; }
            else { x1=0.65; y2=0.35; }
            return std::make_tuple(x1,y2,x1+legW,y2+legH);
        };

        // Put requested texts in pad title area instead of legends
        // Left pad: just "<m_T>"
        can->cd(1);
        {
            TLatex t; t.SetNDC(); t.SetTextAlign(22); t.SetTextSize(0.035);
            t.DrawLatex(0.5, 0.94, "<m_{T}>");
        }

        // Add auto-placed legends for alpha (both pads) so user can identify analytic curve vs UrQMD
        if(iparam==0){
            // Left pad legend: use mT-averaged positions
            {
                std::tuple<double,double,double,double> pos = chooseLegendPos(x_energy, yavg, NENERGIES, xmin, xmax, ymin, ymax);
                double lx = std::get<0>(pos);
                double ytop = std::get<1>(pos);
                double rx = std::get<2>(pos);
                double legH = 0.12; // must match helper
                double lly = ytop - legH;
                can->cd(1);
                TLegend* leg1 = new TLegend(lx, lly, rx, ytop);
                leg1->SetBorderSize(0);
                leg1->SetFillStyle(0);
                leg1->SetTextSize(0.032);
                leg1->AddEntry(gavg, "UrQMD (m_{T}-averaged)", "f");
                leg1->AddEntry(f, "#alpha=0.85 + #sqrt{s_{NN}}^{-0.14} (trend of STAR data)", "l");
                leg1->Draw("same");
            }
            // Right pad legend: use ikt=2 positions
            {
                std::tuple<double,double,double,double> pos2 = chooseLegendPos(x_energy, yvals, NENERGIES, xmin, xmax, ymin2, ymax2);
                double lx2 = std::get<0>(pos2);
                double ytop2 = std::get<1>(pos2);
                double rx2 = std::get<2>(pos2);
                double legH2 = 0.12;
                double lly2 = ytop2 - legH2;
                can->cd(2);
                TLegend* leg2 = new TLegend(lx2, lly2, rx2, ytop2);
                leg2->SetBorderSize(0);
                leg2->SetFillStyle(0);
                leg2->SetTextSize(0.030);
                leg2->AddEntry(gikt2, "UrQMD (m_{T} bin)", "f");
                leg2->AddEntry(f, "#alpha=0.85 + #sqrt{s_{NN}}^{-0.14} (trend of STAR data)", "l");
                leg2->Draw("same");
            }
        }

        // Right pad: show m_T bin center (3rd bin) and note about sys. unc. of m_T choice
        can->cd(2);
        {
            TLatex t; t.SetNDC(); t.SetTextAlign(22); t.SetTextSize(0.032);
            double bincenter_mev = mtbin_centers[2]*1000.0;
            t.DrawLatex(0.5, 0.955, Form("<m_{T}> = %.0f MeV", bincenter_mev));
            t.SetTextSize(0.030);
            t.DrawLatex(0.5, 0.925, "incl. sys. unc. of m_{T} choice");
        }

        // Save
        can->SetTitle(Form("sqrtS vs %s systematic uncertainties", (iparam==0)?"#alpha":(iparam==1)?"R":"N"));
        can->SaveAs(Form("figs/syserr/sqrtS_vs_%s.png", levy_params[iparam]));
        delete can;
    }

    /////////////////////////////////////////////////////
    // Single-panel overlays: each parameter         ////
    // showing both mT-averaged and ikt=2 together   ////
    /////////////////////////////////////////////////////
    
    for(int iparam=0; iparam<3; iparam++){
        TCanvas* can_overlay = new TCanvas(Form("can_overlay_param_%d", iparam), "", 1200, 800);
        //can_overlay->SetLogx(1); // !!! FIXME if not needed to compare with STAR data

        
        // Precompute y-range for both datasets
        double ymin_both = 1e9, ymax_both = -1e9;
        
        // mT-averaged range
        for(int ie=0; ie<NENERGIES; ie++){
            double avg = (iparam==0)? avg_alpha[ie] : (iparam==1)? avg_R[ie] : avg_N[ie];
            double up_sq = 0., dn_sq = 0.;
            if(iparam==0){
                up_sq += pow(alpha_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            } else if(iparam==1){
                up_sq += pow(R_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            } else {
                up_sq += pow(N_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
            }
            double errup = sqrt(up_sq);
            double errdn = sqrt(dn_sq);
            if(std::isfinite(avg)){
                ymin_both = std::min(ymin_both, avg - errdn);
                ymax_both = std::max(ymax_both, avg + errup);
            }
        }
        
        // ikt=2 range
        for(int ie=0; ie<NENERGIES; ie++){
            double val = (iparam==0)? val_alpha_ikt2[ie] : (iparam==1)? val_R_ikt2[ie] : val_N_ikt2[ie];
            TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
            double base_up = 0., base_dn = 0.;
            if(gsys && gsys->GetN() > 2){ base_up = gsys->GetEYhigh()[2]; base_dn = gsys->GetEYlow()[2]; }
            double neigh_up_sq = 0., neigh_dn_sq = 0.;
            TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
            if(gdef && gdef->GetN()>2){
                double center = gdef->GetY()[2];
                for(int j : {1,3}){
                    if(j<0 || j>=gdef->GetN()) continue;
                    double nb = gdef->GetY()[j];
                    double diff = nb - center;
                    if(diff>=0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
                }
            }
            double neigh_up = sqrt(neigh_up_sq);
            double neigh_dn = sqrt(neigh_dn_sq);
            double errup2 = sqrt(base_up*base_up + neigh_up*neigh_up);
            double errdn2 = sqrt(base_dn*base_dn + neigh_dn*neigh_dn);
            if(std::isfinite(val)){
                ymin_both = std::min(ymin_both, val - errdn2);
                ymax_both = std::max(ymax_both, val + errup2);
            }
        }
        
        // For alpha, also account for analytic curve range
        if(iparam==0){
            double fxmin = x_energy[0];
            double fxmax = x_energy[NENERGIES-1];
            int Nsample = 500;
            for(int is=0; is<=Nsample; ++is){
                double x = fxmin + (fxmax - fxmin) * (double)is / (double)Nsample;
                double y = 0.85 + pow(x, -0.14);
                if(std::isfinite(y)){
                    ymin_both = std::min(ymin_both, y);
                    ymax_both = std::max(ymax_both, y);
                }
            }
        }
        
        if(ymin_both > ymax_both){ ymin_both = 0.; ymax_both = 1.; }
        double ypadmargin = 0.06 * (ymax_both - ymin_both);
        ymin_both -= ypadmargin;
        ymax_both += ypadmargin;
        
        // Draw frame
        double xmin = x_energy[0] - 1.0;
        double xmax = x_energy[NENERGIES-1] + 1.0;
        TH1F* frame_overlay = gPad->DrawFrame(xmin, ymin_both, xmax, ymax_both);
        if(iparam==2)
        {
            frame_overlay = gPad->DrawFrame(xmin, 0.7, xmax, 1.08);
        }
        else if(iparam==0)
        {
            //frame_overlay = gPad->DrawFrame(1.0, 0.1, 10000., 2.1); // FIXME if not needed to compare with STAR data
        }
        frame_overlay->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
        const char* ytitle_overlay = (iparam==0)?"#alpha":(iparam==1)?"R [fm]":"#lambda"; // =N
        frame_overlay->GetYaxis()->SetTitle(ytitle_overlay);
        
        // Draw mT-averaged data
        TGraphAsymmErrors* g_avg_overlay;
        {
            double yavg_overlay[NENERGIES];
            double yavg_errup_overlay[NENERGIES];
            double yavg_errdn_overlay[NENERGIES];
            for(int ie=0; ie<NENERGIES; ie++){
                double avg = (iparam==0)? avg_alpha[ie] : (iparam==1)? avg_R[ie] : avg_N[ie];
                yavg_overlay[ie] = avg;
                double up_sq = 0., dn_sq = 0.;
                if(iparam==0){
                    up_sq += pow(alpha_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(alpha_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(alpha_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    // also add m_T choice sys uncert
                    up_sq += pow(alpha_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(alpha_mTchoice_syst_dn[ie], 2);
                } else if(iparam==1){
                    up_sq += pow(R_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(R_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(R_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(R_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(R_mTchoice_syst_dn[ie], 2);
                } else {
                    up_sq += pow(N_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(N_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(N_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(N_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(N_mTchoice_syst_dn[ie], 2);
                }
                yavg_errup_overlay[ie] = sqrt(up_sq);
                yavg_errdn_overlay[ie] = sqrt(dn_sq);
            }
            g_avg_overlay = new TGraphAsymmErrors(NENERGIES, x_energy, yavg_overlay, xerr_low_s, xerr_high_s, yavg_errdn_overlay, yavg_errup_overlay);
            int col_avg = kBlue+2;
            int fillcol_avg = TColor::GetColorTransparent(col_avg, 0.35);
            g_avg_overlay->SetFillColor(fillcol_avg);
            g_avg_overlay->SetFillStyle(1001);
            g_avg_overlay->SetLineColor(col_avg);
            g_avg_overlay->SetLineWidth(3);
            g_avg_overlay->Draw("3 same");
            g_avg_overlay->SetMarkerStyle(21);
            g_avg_overlay->SetMarkerSize(0.0);
            if(iparam==2)
            {
                //g_avg_overlay->GetYaxis()->SetRangeUser(0.9,1.1);
            }
            g_avg_overlay->Draw("LPX same");
        }
        
        // Draw ikt=2 data
        TGraphAsymmErrors* g_ikt2_overlay;
        {
            double yikt2_overlay[NENERGIES];
            double yikt2_errup_overlay[NENERGIES];
            double yikt2_errdn_overlay[NENERGIES];
            for(int ie=0; ie<NENERGIES; ie++){
                double val = (iparam==0)? val_alpha_ikt2[ie] : (iparam==1)? val_R_ikt2[ie] : val_N_ikt2[ie];
                yikt2_overlay[ie] = val;
                TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
                double base_up = 0., base_dn = 0.;
                if(gsys && gsys->GetN() > 2){ base_up = gsys->GetEYhigh()[2]; base_dn = gsys->GetEYlow()[2]; }
                double neigh_up_sq = 0., neigh_dn_sq = 0.;
                TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
                if(gdef && gdef->GetN()>2){
                    double center = gdef->GetY()[2];
                    for(int j : {1,3}){
                        if(j<0 || j>=gdef->GetN()) continue;
                        double nb = gdef->GetY()[j];
                        double diff = nb - center;
                        if(diff>=0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
                    }
                }
                double neigh_up = sqrt(neigh_up_sq);
                double neigh_dn = sqrt(neigh_dn_sq);
                yikt2_errup_overlay[ie] = sqrt(base_up*base_up + neigh_up*neigh_up);
                yikt2_errdn_overlay[ie] = sqrt(base_dn*base_dn + neigh_dn*neigh_dn);
            }
            g_ikt2_overlay = new TGraphAsymmErrors(NENERGIES, x_energy, yikt2_overlay, xerr_low_s, xerr_high_s, yikt2_errdn_overlay, yikt2_errup_overlay);
            int col_ikt2 = kRed+1;
            int fillcol_ikt2 = TColor::GetColorTransparent(col_ikt2, 0.35);
            g_ikt2_overlay->SetFillColor(fillcol_ikt2);
            g_ikt2_overlay->SetFillStyle(1001);
            g_ikt2_overlay->SetLineColor(col_ikt2);
            g_ikt2_overlay->SetLineWidth(3);
            g_ikt2_overlay->Draw("3 same");
            g_ikt2_overlay->SetMarkerStyle(21);
            g_ikt2_overlay->SetMarkerSize(0.0);
            if(iparam==2)
            {
                //g_ikt2_overlay->GetYaxis()->SetRangeUser(0.9,1.1);
            }
            g_ikt2_overlay->Draw("LPX same");
        }
        
        // For alpha, draw analytic curve
        TF1* f_analytic = nullptr;
        if(iparam==0){
            double fxmin = x_energy[0];
            double fxmax = x_energy[NENERGIES-1];
            f_analytic = new TF1("alpha_curve_overlay", "0.85 + pow(x, -0.14)", fxmin, fxmax);
            f_analytic->SetLineColor(kBlack);
            f_analytic->SetLineStyle(1);
            f_analytic->SetLineWidth(4);
            f_analytic->SetNpx(500);
            f_analytic->Draw("L same");
        }
        
        // Legend
        TLegend* leg_overlay = new TLegend(0.30, 0.71, 0.55, 0.91);
        leg_overlay->SetBorderSize(0);
        leg_overlay->SetFillStyle(0);
        leg_overlay->SetTextSize(0.025);
        leg_overlay->AddEntry(g_avg_overlay, "UrQMD <m_{T}>, incl. m_{T} choice sys.unc.", "f");
        leg_overlay->AddEntry(g_ikt2_overlay, "UrQMD single m_{T} = 331 MeV bin", "f");
        if(iparam==0){
            leg_overlay->AddEntry(f_analytic, "#alpha=0.85 + #sqrt{s_{NN}}^{-0.14} (trend of STAR data)", "l");
        }
        leg_overlay->Draw("L same");
        
        can_overlay->SaveAs(Form("figs/syserr/sqrtS_overlay_%s.png", levy_params[iparam]));
        delete can_overlay;
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

    // Write CSV percent breakdown collected during the energy loop
    if(!csv_lines.empty()){
        std::ofstream csvout("syserr_percent_breakdown.csv");
        if(csvout.is_open()){
            for(const auto &ln : csv_lines) csvout << ln << "\n";
            csvout.close();
        } else {
            cerr << "Warning: could not open syserr_percent_breakdown.csv for writing" << endl;
        }
    }

    return 0;
}
