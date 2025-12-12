#ifndef _levy_reader_h_
#define _levy_reader_h_

#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <complex>
 
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TTree.h>
#include <TFile.h>
 
#include <TSystem.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TMath.h>
#include <TGaxis.h>
#include <TPad.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
 
#define SQR(x) ((x)*(x))

class Levy_reader
{
 public:
  Levy_reader(const char* tablefile_name);
  ~Levy_reader();
  double read_table_3d(const double alpha, const double x) const;
  double read_table_1d(const double alpha, const double x) const;
  double getValue_3d(const double alpha, const double x) const; // this is the proper one, with linear interpolation
  double getValue_1d(const double alpha, const double x) const; // this is the proper one, with linear interpolation
 private:
  float** levy_3d_array;
  float** levy_1d_array;
//  float** levy_cyl_array;
  double alpha_min;
  double alpha_max;
  double x_min;
  double x_max;
  int N_alpha;
  int N_x;
  double d_alpha;
  double d_x;

  void InitTable(const char* tablefile_name);
};

#endif // _levy_reader_h_
