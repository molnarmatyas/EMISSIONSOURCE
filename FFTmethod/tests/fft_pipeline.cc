// Integration and numerical regressions using small, deterministic ROOT fixtures.
#include "../../pair_source_tree.h"
#include "../../header_for_all_emissionsource.h"
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TF3.h>
#include <TParameter.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TROOT.h>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;
namespace {
struct TestPair {
    Int_t event;
    Short_t charge, kt;
    float rho;
    std::array<float,3> osl;
};
fs::path directory, programs;
void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }
void near(double actual, double expected, double tolerance, const std::string& message) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
            message + ": got " + std::to_string(actual) + ", expected " + std::to_string(expected));
}
std::string quote(const std::string& text) {
    std::string result = "'";
    for (char c : text) result += c == '\'' ? "'\\''" : std::string(1, c);
    return result + "'";
}
void run(const std::string& program, const std::vector<std::string>& args, int expected = 0) {
    std::string command = "cd " + quote(directory.string()) + " && " + quote(program);
    for (const auto& arg : args) command += " " + quote(arg);
    command += " > " + quote((directory / "last_command.log").string()) + " 2>&1";
    int status = std::system(command.c_str());
    if (status < 0 || !WIFEXITED(status) ||
        (WEXITSTATUS(status) != expected && !(expected == -1 && (WEXITSTATUS(status) == 0 || WEXITSTATUS(status) == 2)))) {
        std::ifstream log(directory / "last_command.log");
        std::cerr << log.rdbuf();
        throw std::runtime_error("Unexpected status for " + program + ": " + std::to_string(status));
    }
}
void transform(std::vector<std::string> args, int expected = 0) {
    run((programs / "transform_fft.exe").string(), args, expected);
}
void fit(std::vector<std::string> args, int expected = 0) {
    run((programs / "CqFit.exe").string(), args, expected);
}

std::vector<TestPair> readPairs(TTree* tree, bool signedCoordinates = false) {
    PairSourceData pair;
    require(setPairSourceBranches(tree, pair, signedCoordinates), "Cannot bind pair branches");
    std::vector<TestPair> rows;
    for(Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        TestPair row{pair.ievent, pair.ich, pair.iKT, pair.rho,
                     {{pair.rhoout, pair.rhoside, pair.rholong}}};
        if(signedCoordinates)
            for(int axis = 0; axis < 3; ++axis)
                if(pair.signs & (1 << axis)) row.osl[axis] = -row.osl[axis];
        rows.push_back(row);
    }
    tree->ResetBranchAddresses();
    return rows;
}
const std::vector<TestPair> fixture = {
    {0, 0, 0, std::sqrt(6.f), {{1, 2, -1}}},
    {0, 1, 0, std::sqrt(6.f), {{-1, 1, 2}}},
    {2, 0, 0, std::sqrt(14.f), {{3, -2, 1}}},
    {3, 1, 0, 10, {{10, 0, 0}}}, // overflow must contribute to normalization
    {2, 0, 1, 0, {{0, 0, 0}}}
};
void makePairs(const char* filename, bool signedCoordinates) {
    TFile file((directory / filename).c_str(), "RECREATE");
    TTree tree("pairTree", "test");
    TestPair p{};
    float out, side, longitudinal;
    UChar_t signs;
    tree.Branch("ievent", &p.event); tree.Branch("ich", &p.charge); tree.Branch("iKT", &p.kt);
    tree.Branch("rho", &p.rho); tree.Branch("rhoout", &out);
    tree.Branch("rhoside", &side); tree.Branch("rholong", &longitudinal);
    if (signedCoordinates) {
        tree.Branch("rho_osl_signs", &signs);
    }
    for (const auto& row : fixture) {
        p = row; out = std::abs(p.osl[0]); side = std::abs(p.osl[1]); longitudinal = std::abs(p.osl[2]);
        signs = (p.osl[0] < 0 ? 1 : 0) | (p.osl[1] < 0 ? 2 : 0) | (p.osl[2] < 0 ? 4 : 0);
        tree.Fill();
    }
    tree.Write(); TParameter<Long64_t>("source_events", 6).Write();
}
void checkDirect(const TH3& cq, const TH3& source) {
    const double normalization = source.Integral(0, source.GetNbinsX() + 1, 0, source.GetNbinsY() + 1, 0, source.GetNbinsZ() + 1);
    for (int x = 1; x <= cq.GetNbinsX(); ++x)
        for (int y = 1; y <= cq.GetNbinsY(); ++y)
            for (int z = 1; z <= cq.GetNbinsZ(); ++z) {
                const double qo = cq.GetXaxis()->GetBinCenter(x), qs = cq.GetYaxis()->GetBinCenter(y), ql = cq.GetZaxis()->GetBinCenter(z);
                double expected = 1.;
                for (int a = 1; a <= source.GetNbinsX(); ++a)
                    for (int b = 1; b <= source.GetNbinsY(); ++b)
                        for (int c = 1; c <= source.GetNbinsZ(); ++c)
                            expected += source.GetBinContent(a,b,c) * std::cos((qo * source.GetXaxis()->GetBinCenter(a) +
                                qs * source.GetYaxis()->GetBinCenter(b) + ql * source.GetZaxis()->GetBinCenter(c)) / fftHbarC) / normalization;
                near(cq.GetBinContent(x,y,z), expected, 2.e-12, "FFT vs direct cosine sum");
                near(cq.GetBinError(x,y,z), 1., 1.e-14, "Uniform bin error independent of cutoff");
            }
}
void pairTests() {
    makePairs("pairs.root", true); makePairs("legacy.root", false);
    std::array<std::array<double,4>,2> sums{};
    for (int f = 0; f < 2; ++f) {
        TFile input((directory / (f ? "pairs.root" : "legacy.root")).c_str());
        auto* tree = input.Get<TTree>("pairTree");
        for(const TestPair& p : readPairs(tree)) {
            sums[f][0] += p.rho;
            for (int i = 0; i < 3; ++i) sums[f][i+1] += std::abs(p.osl[i]);
        }
    }
    require(sums[0] == sums[1], "Legacy projection coordinates changed");
    transform({"pairs.root", "--bins", "9", "--range", "4.5", "--save-source"});
    {
        TFile result((directory / "pairs_Cq_FFT.root").c_str());
        auto* cq = result.Get<TH3D>("Cq_D_lcms_3d_kt0");
        auto* source = result.Get<TH3D>("D_lcms_3d_kt0");
        require(cq && source, "Missing double precision histograms");
        checkDirect(*cq, *source);
        near(cq->GetXaxis()->GetBinCenter(5), 0., 1.e-14, "q=0 centered");
        near(cq->GetBinContent(5,5,5), 1.75, 1.e-14, "Flow normalization");
        require(result.Get<TParameter<Long64_t>>("source_events")->GetVal() == 6, "Empty event provenance lost");
        auto* metadata = result.Get<TTree>("fft_sources");
        require(metadata && metadata->GetEntries() == 2, "kT grouping metadata");
        TTreeReader reader(metadata);
        TTreeReaderValue<double> retained(reader, "retained_fraction");
        require(reader.Next(), "Missing source summary"); near(*retained, .75, 1.e-14, "Retained fraction metadata");
    }
    transform({"pairs.root", "--bins", "9", "--range", "4.5", "--output", "selected_FFT.root",
               "--kt", "0", "--first-event", "2", "--events", "1", "--charge", "negative", "--save-source"});
    TFile selected((directory / "selected_FFT.root").c_str());
    auto* selectedSource = selected.Get<TH3D>("D_lcms_3d_kt0");
    near(selectedSource->Integral(), 1., 1.e-14, "Pair selection");
    checkDirect(*selected.Get<TH3D>("Cq_D_lcms_3d_kt0"), *selectedSource);
    require(!selected.Get("Cq_D_lcms_3d_kt1"), "kT filter ignored");
    transform({"legacy.root"}, 1);
    require(!fs::exists(directory / "legacy_Cq_FFT.root"), "Legacy unsigned tree accepted by FFT");
    transform({"pairs.root", "--bins", "8"}, 1);
    transform({"pairs.root", "--cutoff", "nan"}, 1);
    transform({"pairs.root", "--error", "0"}, 1);
    transform({"pairs.root", "--events", "0"}, 1);
    transform({"pairs.root", "--first-event", "50", "--output", "empty_FFT.root"}, 1);
    transform({"pairs.root", "--bins", "9", "--range", "4.5", "--output", "legacy.root"});
    require(fs::exists(directory / "legacy_FFT.root"), "Custom output did not get the FFT suffix");
    {
        TFile legacy((directory / "legacy.root").c_str());
        require(legacy.Get<TTree>("pairTree") && !legacy.Get("fft_complete"), "D(rho) file was overwritten");
    }
    transform({"pairs.root", "--bins", "9", "--range", "4.5", "--kt", "0"});
    {
        TFile replaced((directory / "pairs_Cq_FFT.root").c_str());
        require(replaced.Get<TH3D>("Cq_D_lcms_3d_kt0"), "Repeated FFT lost its result");
        require(!replaced.Get("Cq_D_lcms_3d_kt1") && !replaced.Get("D_lcms_3d_kt0"), "Repeated FFT kept stale output keys");
    }
    fs::copy_file(directory / "pairs.root", directory / "input_FFT.root");
    transform({"input_FFT.root", "--output", "input_FFT.root"}, 1);
    require(!fs::exists(directory / "empty_FFT.root"), "Empty selection published output");
    run("/tmp/emissionsource_transform_disabled.exe", {"pairs.root", "--output", "disabled_FFT.root"}, 1);
    require(!fs::exists(directory / "disabled_FFT.root"), "Disabled FFT switch was ignored");
    run("/tmp/emissionsource_projection_cache_tests.exe", {"legacy.root", "pairs.root"});
    std::cout << "PASS signed trees, old projections, selection, flows, FFT replacement and separate D(rho) names\n";
}
void histogramTests() {
    {
        TFile input((directory / "shifted.root").c_str(), "RECREATE");
        TH3D h("shifted", "", 9, -3.8, 5.2, 9, -9, 9, 9, 1, 10);
        h.SetBinContent(2,3,4,2.); h.SetBinContent(7,5,8,1.); h.SetBinContent(0,1,1,1.); h.Write();
        TH3D h2("smaller", "", 7,-3.5,3.5,7,-3.5,3.5,7,-3.5,3.5);
        h2.Fill(1,-1,0); h2.Write();
    }
    transform({"shifted.root", "--save-source"});
    TFile result((directory / "shifted_Cq_FFT.root").c_str());
    for (const char* name : {"shifted", "smaller"})
        checkDirect(*result.Get<TH3D>((std::string("Cq_") + name).c_str()), *result.Get<TH3D>(name));
    for (int test = 0; test < 4; ++test) {
        const std::string name = "invalid" + std::to_string(test) + ".root";
        {
            TFile input((directory / name).c_str(), "RECREATE");
            if (test == 0) { TH3D h("even", "", 4,-2,2,4,-2,2,4,-2,2); h.Fill(0.,0.,0.); h.Write(); }
            else if (test == 1) {
                double edges[] = {-3.,-2.,-1.,0.,1.,3.};
                TH3D h("variable", "",5,edges,5,edges,5,edges); h.Fill(0.,0.,0.); h.Write();
            } else {
                TH3D h("bad_counts", "",5,-2.5,2.5,5,-2.5,2.5,5,-2.5,2.5);
                h.SetBinContent(1,1,1,test == 2 ? -1. : std::numeric_limits<double>::quiet_NaN()); h.Write();
            }
        }
        transform({name}, 1);
        require(!fs::exists(directory / ("invalid" + std::to_string(test) + "_Cq_FFT.root")), "Invalid histogram published output");
    }
    std::cout << "PASS direct DFT including mixed signs, shifted origins, anisotropic axes and plan sizes\n";
}

double levy(double x, double y, double z, const std::array<double,5>& p) {
    return 1 + p[0] * std::exp(-std::pow((x*x*p[1]*p[1] + y*y*p[2]*p[2] + z*z*p[3]*p[3]) / (fftHbarC*fftHbarC), p[4]/2));
}
struct FitResult { std::array<double,5> p{}, raw{}, errors{}; std::array<double,25> covariance{}, covarianceRaw{}; double scale; Long64_t ndf; bool valid, atLimit; };
void checkParameterGraphs(TFile& input) {
    const char* names[] = {"lambda_vs_KT_FFT", "R_out_vs_KT_FFT", "R_side_vs_KT_FFT", "R_long_vs_KT_FFT", "alpha_vs_KT_FFT"};
    TGraphAsymmErrors* graphs[5];
    for(int p=0; p<5; ++p) {
        graphs[p]=input.Get<TGraphAsymmErrors>(names[p]);
        require(graphs[p], "Missing parameter graph: " + std::string(names[p]));
    }
    auto* tree=input.Get<TTree>("fit_results_FFT");
    tree->ResetBranchAddresses();
    int kt; bool valid; double parameters[5], errors[5];
    tree->SetBranchAddress("kt", &kt); tree->SetBranchAddress("valid", &valid);
    tree->SetBranchAddress("parameters", parameters); tree->SetBranchAddress("errors", errors);
    int point=0;
    // Compare graphs to valid tree rows in ascending kT order, also for sparse fits.
    for(int expectedKT=0; expectedKT<NKT; ++expectedKT) for(Long64_t i=0; i<tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if(!valid || kt!=expectedKT) continue;
        for(int p=0; p<5; ++p) {
            auto* graph=graphs[p];
            require(point<graph->GetN(), "Parameter graph lost a valid fit");
            near(graph->GetPointX(point),kT_center[kt],1.e-14,"Graph kT order");
            near(graph->GetPointY(point),parameters[p],1.e-14,"Graph parameter value");
            near(graph->GetErrorXlow(point),kT_center[kt]-ktbins[kt],1.e-14,"Graph lower kT edge");
            near(graph->GetErrorXhigh(point),ktbins[kt+1]-kT_center[kt],1.e-14,"Graph upper kT edge");
            near(graph->GetErrorYlow(point),errors[p],1.e-14,"Graph lower rescaled error");
            near(graph->GetErrorYhigh(point),errors[p],1.e-14,"Graph upper rescaled error");
        }
        ++point;
    }
    for(auto* graph : graphs) require(graph->GetN()==point, "Parameter graph contains an invalid/unlabelled fit or zero placeholder");
    tree->ResetBranchAddresses();
}
FitResult readFit(const std::string& name, bool requireValid = true) {
    TFile input((directory / name).c_str());
    auto* tree = input.Get<TTree>("fit_results_FFT");
    require(tree && tree->GetEntries() == 1, "Missing fit result row");
    FitResult result;
    tree->SetBranchAddress("parameters", result.p.data()); tree->SetBranchAddress("errors_raw", result.raw.data());
    tree->SetBranchAddress("errors", result.errors.data()); tree->SetBranchAddress("covariance", result.covariance.data());
    tree->SetBranchAddress("covariance_raw", result.covarianceRaw.data()); tree->SetBranchAddress("error_scale", &result.scale);
    tree->SetBranchAddress("ndf", &result.ndf); tree->SetBranchAddress("valid", &result.valid); tree->GetEntry(0);
    tree->SetBranchAddress("at_parameter_limit", &result.atLimit); tree->GetEntry(0);
    if (requireValid) require(result.valid, "Fit is not valid");
    auto* function = input.Get<TF3>("fit_Cq_model_FFT");
    if (function)
        near(function->Eval(.013,.021,-.019), levy(.013,.021,-.019,result.p), 1.e-13, "Saved fit function after reload");
    for (int p = 0; p < 5; ++p) near(result.errors[p], result.raw[p] * result.scale, 1.e-12, "Error rescaling");
    for (int p = 0; p < 25; ++p)
        near(result.covariance[p], result.covarianceRaw[p] * result.scale * result.scale, 1.e-12, "Covariance rescaling");
    checkParameterGraphs(input);
    return result;
}
void fitTests() {
    const std::array<double,5> truth = {{.72,4.2,5.3,6.4,1.35}};
    Long64_t points = 0;
    for (int sample = 0; sample < 2; ++sample) {
        TFile input((directory / ("model" + std::to_string(sample) + ".root")).c_str(), "RECREATE");
        TH3D h("Cq_model", "",41,-.205,.205,41,-.205,.205,41,-.205,.205);
        for (int x = 1; x <= 41; ++x) for (int y = 1; y <= 41; ++y) for (int z = 1; z <= 41; ++z) {
            const double qo=h.GetXaxis()->GetBinCenter(x), qs=h.GetYaxis()->GetBinCenter(y), ql=h.GetZaxis()->GetBinCenter(z);
            const bool inside = qo*qo + qs*qs + ql*ql <= .12*.12;
            const double value = inside ? levy(qo,qs,ql,truth) + .0002 * std::cos(71*qo+53*qs+31*ql) : 100.;
            h.SetBinContent(x,y,z,value); h.SetBinError(x,y,z,sample ? 5. : 1.);
            if (!sample && inside) ++points;
        }
        h.Write();
        std::string correlation = "Cq_model";
        int kt = 2;
        TTree metadata("fft_sources", "test provenance");
        metadata.Branch("correlation", &correlation); metadata.Branch("kt", &kt);
        metadata.Fill(); metadata.Write();
        TParameter<Long64_t>("source_events", 42).Write();
    }
    std::string limits;
    for(int kt = 0; kt < NKT; ++kt) limits += (kt ? "," : "") + std::string(kt == 2 ? "0.12" : "0.10");
    fit({"model0.root", "--qmax", "0.12"}); fit({"model1.root", "--qmax", limits});
    auto a = readFit("model0_fitresults_FFT.root"), b = readFit("model1_fitresults_FFT.root");
    require(a.ndf == points-5, "Spherical fit domain or NDF incorrect");
    for (int p = 0; p < 5; ++p) {
        near(a.p[p], truth[p], .01, "Recovery of known Levy parameters");
        near(a.p[p], b.p[p], 1.e-5, "Fit depends on arbitrary uniform error");
        near(a.errors[p], b.errors[p], .01 * a.errors[p], "Rescaled errors depend on arbitrary weight");
    }
    {
        TFile previous((directory / "model0_fitresults_FFT.root").c_str(), "UPDATE");
        TNamed("old_fit_output", "stale key").Write();
    }
    const fs::path plot=directory / "figs/fitting/lcms/FFT/model0_fitresults_FFT/Cq_model_FFT";
    std::ofstream(plot.string()+".png") << "old plot";
    std::ofstream(plot.string()+".pdf") << "old plot";
    fit({"model0.root", "--qmax", "0.12", "--output", "model0_fitresults.root"});
    const auto repeated=readFit("model0_fitresults_FFT.root");
    for(int p=0; p<5; ++p) near(repeated.p[p],a.p[p],1.e-10,"Repeated fit changed parameters");
    {
        TFile replaced((directory / "model0_fitresults_FFT.root").c_str());
        require(!replaced.Get("old_fit_output"), "Repeated fit kept stale output keys");
    }
    require(fs::file_size(plot.string()+".png")>1000 && fs::file_size(plot.string()+".pdf")>1000, "Repeated fit did not replace its figures");
    fs::copy_file(directory / "model0.root", directory / "input_Cq_FFT.root");
    fit({"input_Cq_FFT.root", "--output", "input_Cq_FFT.root"},1);
    fit({"model0.root", "--qmax", "0.3", "--output", "outside_FFT.root"}, 1);
    fit({"model0.root", "--qmax", "0.00001", "--output", "tiny_FFT.root"}, 1);
    require(!fs::exists(directory / "outside_FFT.root"), "Out-of-grid fit published output");
    std::cout << "PASS spherical Levy fit, recovery, covariance scaling, weight invariance and replacing results/figures\n";
}
void qmaxDefaultsTest() {
    const std::array<double,5> truth = {{.72,4.2,5.3,6.4,1.35}};
    {
        TFile input((directory / "default_ranges.root").c_str(), "RECREATE");
        std::string correlation;
        int kt;
        TTree metadata("fft_sources", "per-kT fit limits");
        metadata.Branch("correlation", &correlation); metadata.Branch("kt", &kt);
        // Reverse key order: selection must use the metadata kT index.
        for(kt=NKT-1; kt>=0; --kt) {
            correlation = "Cq_test_kt" + std::to_string(kt);
            TH3D h(correlation.c_str(), "",41,-.205,.205,41,-.205,.205,41,-.205,.205);
            for(int x=1; x<=41; ++x) for(int y=1; y<=41; ++y) for(int z=1; z<=41; ++z) {
                double qo=h.GetXaxis()->GetBinCenter(x), qs=h.GetYaxis()->GetBinCenter(y), ql=h.GetZaxis()->GetBinCenter(z);
                h.SetBinContent(x,y,z,levy(qo,qs,ql,truth) + .0002*std::cos(71*qo+53*qs+31*ql));
                h.SetBinError(x,y,z,1.);
            }
            h.Write(); metadata.Fill();
        }
        metadata.Write();
        TNamed("source_energy", "27").Write();
        TNamed("source_centrality", "0-10").Write();
        TParameter<int>("source_qlcms_cut",0).Write();
    }
    const std::string selected = "Cq_test_kt" + std::to_string(NKT-1);
    fit({"default_ranges.root"});
    fit({"default_ranges.root", "--hist", selected, "--output", "one_range_FFT.root"});
    for(const char* filename : {"default_ranges_fitresults_FFT.root", "one_range_FFT.root"}) {
        TFile input((directory / filename).c_str());
        auto* tree = input.Get<TTree>("fit_results_FFT");
        const bool single = std::string(filename)=="one_range_FFT.root";
        require(tree && tree->GetEntries()==(single ? 1 : NKT), "Default fit limits lost kT bins");
        int kt; double qmax, parameters[5]; bool valid;
        tree->SetBranchAddress("kt", &kt); tree->SetBranchAddress("qmax", &qmax);
        tree->SetBranchAddress("parameters", parameters); tree->SetBranchAddress("valid", &valid);
        std::vector<bool> seen(NKT,false);
        for(Long64_t i=0; i<tree->GetEntries(); ++i) {
            tree->GetEntry(i);
            require(kt>=0 && kt<NKT && !seen[kt] && valid, "Invalid or duplicate per-kT fit");
            seen[kt]=true;
            if(single) require(kt==NKT-1, "--hist used the wrong kT index");
            near(qmax,fftQmaxFit[kt],0.,"Default per-kT qmax");
            for(int p=0; p<5; ++p) near(parameters[p],truth[p],.01,"Default-range parameter recovery");
        }
        tree->ResetBranchAddresses();
        checkParameterGraphs(input);
    }
    fit({"shifted_Cq_FFT.root", "--hist", "Cq_shifted", "--output", "unknown_kt_FFT.root"},1);
    require(!fs::exists(directory / "unknown_kt_FFT.root"), "Default limits guessed an unknown kT index");
    std::cout << "PASS header qmax defaults, sorted parameter graphs with rescaled errors and --hist selection\n";
}
void plotTests() {
    TFile output((directory / "default_ranges_fitresults_FFT.root").c_str());
    const char* directions[]={"out","side","long"};
    for(int kt=0; kt<NKT; ++kt) {
        const std::string histogram="Cq_test_kt" + std::to_string(kt);
        const fs::path stem=directory / "figs/fitting/lcms/FFT/default_ranges_fitresults_FFT" / (histogram+"_FFT");
        require(fs::file_size(stem.string()+".png")>1000 && fs::file_size(stem.string()+".pdf")>1000, "Missing kT fit PNG/PDF");
        auto* canvas=output.Get<TCanvas>(("fitplot_"+histogram+"_FFT").c_str());
        require(canvas, "Missing saved fit canvas");
        auto* model=output.Get<TF3>(("fit_"+histogram+"_FFT").c_str());
        require(model, "Missing plotted model");
        for(int d=0; d<3; ++d) {
            const std::string direction=directions[d];
            auto* pad=(TPad*)canvas->GetPrimitive(("fitpad_"+direction+"_FFT").c_str());
            require(pad, "Missing out/side/long plot panel");
            auto* data=(TGraphErrors*)pad->GetPrimitive(("Cq_slice_"+direction+"_FFT").c_str());
            auto* curve=(TGraph*)pad->GetPrimitive(("Cq_fit_"+direction+"_FFT").c_str());
            auto* full=(TGraph*)pad->GetPrimitive(("Cq_fullfit_"+direction+"_FFT").c_str());
            require(data && curve && full && data->GetN()>0, "Missing slice data or model");
            near(curve->GetPointX(0),-fftQmaxFit[kt],1.e-14,"Solid curve lower fit edge");
            near(curve->GetPointX(curve->GetN()-1),fftQmaxFit[kt],1.e-14,"Solid curve upper fit edge");
            require(curve->GetLineStyle()==1 && full->GetLineStyle()==2, "Fit/extrapolation styles lost");
            TFile input((directory / "default_ranges.root").c_str());
            auto* source=input.Get<TH3D>(histogram.c_str());
            const int zero=source->GetXaxis()->FindFixBin(0.);
            TTreeReader reader(output.Get<TTree>("fit_results_FFT"));
            TTreeReaderValue<int> rowKT(reader,"kt");
            TTreeReaderValue<double> scale(reader,"error_scale");
            while(reader.Next() && *rowKT!=kt) {}
            require(*rowKT==kt,"Missing plot fit result");
            for(int point=0; point<data->GetN(); ++point) {
                int bins[]={zero,zero,zero}; bins[d]=source->GetXaxis()->FindFixBin(data->GetPointX(point));
                near(data->GetPointY(point),source->GetBinContent(bins[0],bins[1],bins[2]),1.e-14,"Plot uses axis slice");
                near(data->GetErrorY(point),source->GetBinError(bins[0],bins[1],bins[2])*(*scale),1.e-14,"Plot error rescaling");
            }
            for(int point=0; point<curve->GetN(); point+=100) {
                double q[]={0.,0.,0.}; q[d]=curve->GetPointX(point);
                near(curve->GetPointY(point),model->Eval(q[0],q[1],q[2]),1.e-14,"Plot uses full 3D fit model");
            }
        }
    }
    fit({"default_ranges.root", "--hist", "Cq_test_kt0", "--plot-dir", "custom figures", "--output", "custom_plot_FFT.root"});
    require(fs::exists(directory / "custom figures/custom_plot_FFT/Cq_test_kt0_FFT.png"), "Custom plot directory ignored");
    fit({"default_ranges.root", "--hist", "Cq_test_kt0", "--no-plots", "--output", "no_plot_FFT.root"});
    TFile noPlots((directory / "no_plot_FFT.root").c_str());
    require(!noPlots.Get("fitplot_Cq_test_kt0_FFT") && !fs::exists(directory / "figs/fitting/lcms/FFT/no_plot_FFT"), "--no-plots ignored");
    std::cout << "PASS fit figures for all kT bins, saved central slices/model, scaled errors and plot options\n";
}
void systematicOptionsTest() {
    const auto baseline=readFit("one_range_FFT.root");
    for(int range=1; range<=2; ++range) {
        fit({"default_ranges.root","--hist","Cq_test_kt9","--fit-range",std::to_string(range),"--no-plots"});
        const std::string name=std::string("default_ranges_fitresults")+fftFitRangeSuffix[range]+"_FFT.root";
        const auto varied=readFit(name);
        require(range==1 ? varied.ndf<baseline.ndf : varied.ndf>baseline.ndf, "Fit-range variation did not change selected bins");
        TFile input((directory/name).c_str());
        require(input.Get<TParameter<int>>("cqfit_range_variation")->GetVal()==range,"Missing fit-range provenance");
        near(input.Get<TParameter<double>>("cqfit_qmax_scale")->GetVal(),fftQmaxFitScale[range],0.,"Fit-range scale");
        auto* tree=input.Get<TTree>("fit_results_FFT");
        double qmax; tree->SetBranchAddress("qmax",&qmax); tree->GetEntry(0);
        near(qmax,fftQmaxFit[NKT-1]*fftQmaxFitScale[range],1.e-14,"Per-kT scaled range");
        tree->ResetBranchAddresses();
    }
    fit({"default_ranges.root","--fit-range","3","--no-plots"},1);
    transform({"pairs.root","--binning","1","--kt","0"});
    {
        TFile input((directory/"pairs_Cq_coarseFFT_FFT.root").c_str());
        auto* cq=input.Get<TH3D>("Cq_D_lcms_3d_kt0");
        require(cq && cq->GetNbinsX()==fftSourceBinsSyst[1],"Grid variation did not rebin the pair source");
        require(input.Get<TParameter<int>>("fft_binning_variation")->GetVal()==1,"Missing grid provenance");
        near(input.Get<TParameter<double>>("fft_source_half_range")->GetVal(),fftSourceHalfRange,0.,"Grid variation changed spatial extent");
        near(cq->GetBinContent(cq->FindBin(0.,0.,0.)),2.,1.e-14,"Coarse-grid normalization");
        near(cq->GetXaxis()->GetBinWidth(1),acos(-1.)*fftHbarC/fftSourceHalfRange,1.e-12,"Grid variation changed q spacing");
    }
    transform({"pairs.root","--binning","3"},1);
    transform({"pairs.root","--binning","2","--bins","9"},1);
    transform({"shifted.root","--binning","1"},1);
    std::cout << "PASS per-kT fit-range systematics, separate names, and actual source rebinning at fixed extent\n";
}
void gaussianTest() {
    const std::array<double,5> truth = {{1.,4.3,5.2,6.4,2.}};
    {
        TFile input((directory / "gaussian.root").c_str(), "RECREATE");
        TH3D h("gaussian", "",65,-100,100,65,-100,100,65,-100,100);
        for (int x = 1; x <= 65; ++x) for (int y = 1; y <= 65; ++y) for (int z = 1; z <= 65; ++z) {
            const double a=h.GetXaxis()->GetBinCenter(x), b=h.GetYaxis()->GetBinCenter(y), c=h.GetZaxis()->GetBinCenter(z);
            h.SetBinContent(x,y,z,std::exp(-.25*(a*a/(truth[1]*truth[1])+b*b/(truth[2]*truth[2])+c*c/(truth[3]*truth[3]))));
        }
        h.Write();
    }
    transform({"gaussian.root"});
    TFile input((directory / "gaussian_Cq_FFT.root").c_str());
    auto* cq = input.Get<TH3D>("Cq_gaussian");
    for (int x = 1; x <= 65; ++x) for (int y = 1; y <= 65; ++y) for (int z = 1; z <= 65; ++z) {
        const double a=cq->GetXaxis()->GetBinCenter(x), b=cq->GetYaxis()->GetBinCenter(y), c=cq->GetZaxis()->GetBinCenter(z);
        if (a*a+b*b+c*c <= .12*.12) near(cq->GetBinContent(x,y,z), levy(a,b,c,truth), 1.e-10, "Analytic Gaussian transform");
    }
    // At alpha=2 Minuit may return status 0 or an EDM warning; the boundary flag is required.
    fit({"gaussian_Cq_FFT.root", "--qmax", "0.12"}, -1);
    const auto fitted = readFit("gaussian_Cq_fitresults_FFT.root", false);
    require(fitted.atLimit, "Gaussian boundary fit was not flagged");
    for (int p = 0; p < 5; ++p) near(fitted.p[p], truth[p], .002, "Gaussian source to FFT to fit");
    std::cout << "PASS analytic anisotropic Gaussian transform, radius convention and full fit chain\n";
}
void analyzerTest() {
    {
        TFile input((directory / "AuAu_synthetic_tree.root").c_str(), "RECREATE");
        TTree tree("urqmd_tree", "test events");
        std::vector<double> px,py,pz,energy,rx,ry,rz,time;
        std::vector<int> pid;
        tree.Branch("pX",&px); tree.Branch("pY",&py); tree.Branch("pZ",&pz); tree.Branch("E",&energy);
        tree.Branch("rX",&rx); tree.Branch("rY",&ry); tree.Branch("rZ",&rz); tree.Branch("t",&time); tree.Branch("PID",&pid);
        for (int event = 0; event < 3; ++event) {
            px.clear();py.clear();pz.clear();energy.clear();rx.clear();ry.clear();rz.clear();time.clear();pid.clear();
            if (event != 1) for (int track = 0; track < 4; ++track) {
                px.push_back(.28 + .005*track); py.push_back(.01*track); pz.push_back(0.);
                energy.push_back(std::sqrt(px.back()*px.back()+py.back()*py.back()+Mass2_pi));
                rx.push_back(track); ry.push_back(-track); rz.push_back(track*2); time.push_back(0.);
                pid.push_back(track < 2 ? -211 : 211);
            }
            tree.Fill();
        }
        tree.Write();
    }
    run("/tmp/emissionsource_pairsource_legacy.exe", {"synthetic","0"});
    run("/tmp/emissionsource_pairsource_fft.exe", {"synthetic","0"});
    TFile old((directory / "UrQMD_3d_source_0-10cent_all_synthetic.root").c_str());
    TFile fft((directory / "UrQMD_3d_source_0-10cent_all_synthetic_FFT.root").c_str());
    std::vector<TestPair> a = readPairs(old.Get<TTree>("pairTree"));
    std::vector<TestPair> b = readPairs(fft.Get<TTree>("pairTree"));
    std::vector<TestPair> signedRows = readPairs(fft.Get<TTree>("pairTree"), true);
    require(a.size() == 4 && a.size() == b.size(), "Analyser accepted pairs changed");
    bool negative = false;
    for (std::size_t i=0;i<a.size();++i) {
        require(a[i].rho == b[i].rho && a[i].osl == b[i].osl, "Analyser projection output changed");
        for(int axis=0;axis<3;++axis) {
            require(b[i].osl[axis] == std::abs(signedRows[i].osl[axis]), "Signed/magnitude coordinates inconsistent");
            negative = negative || signedRows[i].osl[axis] < 0;
        }
    }
    require(negative, "Analyser lost separation signs");
    require(fft.Get<TParameter<Long64_t>>("source_events")->GetVal() == 3, "Analyser lost empty events");
    old.Close(); fft.Close();
    {
        TFile previous((directory / "UrQMD_3d_source_0-10cent_all_synthetic_FFT.root").c_str(), "UPDATE");
        TNamed("old_pair_output", "stale key").Write();
    }
    run("/tmp/emissionsource_pairsource_fft.exe", {"synthetic","0"});
    {
        TFile replaced((directory / "UrQMD_3d_source_0-10cent_all_synthetic_FFT.root").c_str());
        require(replaced.Get<TTree>("pairTree")->GetEntries()==4 && !replaced.Get("old_pair_output"), "Analyser did not replace its FFT output");
    }
    run("/tmp/emissionsource_pairsource_fft.exe", {"synthetic","0","1"});
    require(fs::exists(directory / "UrQMD_3d_source_0-10cent_all_synthetic_pos_FFT.root"), "Charge output names collide");
    std::cout << "PASS analyser signed branches, unchanged projections, empty events, charge names and FFT output replacement\n";
}

void ktRangeTest() {
    const double values[] = {ktbins[0]-.01, std::nextafter(ktbins[0], 0.), ktbins[0],
        kT_center[NKT/2], ktbins[NKT], std::nextafter(ktbins[NKT], 1.), ktbins[NKT]+.05};
    {
        TFile input((directory / "AuAu_ktlimits_tree.root").c_str(), "RECREATE");
        TTree tree("urqmd_tree", "kT range test");
        std::vector<double> px, py(2,0.), pz(2,0.), energy, rx={0.,1.}, ry={0.,2.}, rz={0.,3.}, time(2,0.);
        std::vector<int> pid(2,-211);
        tree.Branch("pX",&px); tree.Branch("pY",&py); tree.Branch("pZ",&pz); tree.Branch("E",&energy);
        tree.Branch("rX",&rx); tree.Branch("rY",&ry); tree.Branch("rZ",&rz); tree.Branch("t",&time); tree.Branch("PID",&pid);
        // Each event has one otherwise accepted pair with qLCMS=0.
        for(double kt : values) {
            px.assign(2,kt); energy.assign(2,std::sqrt(kt*kt+Mass2_pi)); tree.Fill();
        }
        tree.Write();
    }
    const char* suffix[] = {"", "_strictqLCMS", "_looseqLCMS"};
    for(int mode=0; mode<2; ++mode) for(int cut=0; cut<3; ++cut) {
        run(mode ? "/tmp/emissionsource_pairsource_fft.exe" : "/tmp/emissionsource_pairsource_legacy.exe",
            {"ktlimits",std::to_string(cut)});
        const std::string filename = std::string("UrQMD_3d_source_0-10cent_all_ktlimits") + suffix[cut] + (mode ? "_FFT.root" : ".root");
        TFile input((directory / filename).c_str());
        const auto rows = readPairs(input.Get<TTree>("pairTree"),mode==1);
        require(rows.size()==3, "Out-of-range kT pair entered an edge bin");
        require(rows[0].event==2 && rows[0].kt==0, "Lower kT endpoint was rejected");
        require(rows[1].event==3 && rows[1].kt==NKT/2, "Interior kT bin changed");
        require(rows[2].event==4 && rows[2].kt==NKT-1, "Upper kT endpoint was rejected");
        near(input.Get<TParameter<double>>("source_kt_min")->GetVal(),ktbins[0],0.,"Stored kT minimum");
        near(input.Get<TParameter<double>>("source_kt_max")->GetVal(),ktbins[NKT],0.,"Stored kT maximum");
        require(input.Get<TParameter<Long64_t>>("source_events")->GetVal()==7, "Empty events lost after kT cut");
        near(input.Get<TH1D>("ktdist")->GetEntries(),7.,0.,"Pre-cut kT diagnostic changed");
    }
    std::cout << "PASS strict kT range, exact endpoints and adjacent values for both methods/all qLCMS cuts\n";
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 2) return 1;
    gROOT->SetBatch(true);
    programs = fs::absolute(argv[1]);
    char temporary[] = "/tmp/emissionsource_fft_XXXXXX";
    const char* created = mkdtemp(temporary);
    if (!created) return 1;
    directory = created;
    try {
        pairTests(); histogramTests(); fitTests(); qmaxDefaultsTest(); plotTests(); systematicOptionsTest(); gaussianTest(); analyzerTest(); ktRangeTest();
        std::cout << "All FFT regressions passed. Fixtures: " << directory << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\nFixtures/log: " << directory << '\n'; return 1;
    }
}
