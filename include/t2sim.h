//////////////////////////////////////////////////////////
// Sun May 18 21:37:33 2014 by ROOT version 5.34/03
// data structure of the simulation input data
//////////////////////////////////////////////////////////

#ifndef t2sim_h
#define t2sim_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class t2sim {
public:
  TTree *fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Declaration of leaf types
  Char_t particle_name[10];
  Double_t time;
  Double_t energy;
  Double_t position[3];
  Double_t direction[3];
  Double_t polarization[3];
  //  bool multiple_particles;

  // List of branches
  TBranch *b_pname;        //!
  TBranch *b_time;         //!
  TBranch *b_energy;       //!
  TBranch *b_emiPosition;  //!
  TBranch *b_emiDirection; //!
  TBranch *b_polDirection; //!
                           // TBranch        *b_multiple_particles;   //!

  t2sim(TTree *tree = 0);
  virtual ~t2sim();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

#endif
