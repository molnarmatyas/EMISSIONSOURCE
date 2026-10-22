#ifndef PAIR_SOURCE_TREE_H
#define PAIR_SOURCE_TREE_H
#include <TTree.h>
#include <iostream>

// Branch variables shared by the D(rho) fitter and FFT transformer.
struct PairSourceData {
  Int_t ievent;
  Short_t ich, iKT;
  Float_t rho, rhoout, rhoside, rholong;
  UChar_t signs; // bits 0,1,2: negative out, side, long
};

inline bool setPairSourceBranches(TTree* tree, PairSourceData& pair, bool needSigns = false)
{
  const char* names[] = {"ievent", "ich", "iKT", "rho", "rhoout", "rhoside", "rholong"};
  for(const char* name : names) {
    if(!tree->GetBranch(name)) {
      std::cerr << "Missing pairTree branch: " << name << std::endl;
      return false;
    }
  }
  if(needSigns && !tree->GetBranch("rho_osl_signs")) {
    std::cerr << "FFT needs the separation signs. Rerun pairsource_urqmd with FFT enabled." << std::endl;
    return false;
  }
  tree->SetBranchAddress("ievent", &pair.ievent);
  tree->SetBranchAddress("ich", &pair.ich);
  tree->SetBranchAddress("iKT", &pair.iKT);
  tree->SetBranchAddress("rho", &pair.rho);
  tree->SetBranchAddress("rhoout", &pair.rhoout);
  tree->SetBranchAddress("rhoside", &pair.rhoside);
  tree->SetBranchAddress("rholong", &pair.rholong);
  if(needSigns) tree->SetBranchAddress("rho_osl_signs", &pair.signs);
  return true;
}
#endif
