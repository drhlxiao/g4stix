
/***************************************************************
 * Class to define input data struction for simulation.
 * Author  : Hualin Xiao
 * Date    : May, 2014
 * Version : 1.10
 *
 ***************************************************************/

#include "t2sim.h"
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>

t2sim::t2sim(TTree *tree) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    tree = new TTree("t2sim", "t2sim");
    tree->Branch("particle_name", particle_name, "particle_name[10]/C");
    tree->Branch("time", &time, "time/D");
    tree->Branch("energy", &energy, "energy/D");
    tree->Branch("position", position, "position[3]/D");
    tree->Branch("direction", direction, "direction[3]/D");
    //    tree->Branch("polarization", polarization, "polarization[3]/D");
    //    tree->Branch("multiple_particles", &multiple_particles,
    //    "multiple_particles/o");
  }
  Init(tree);
}

t2sim::~t2sim() {
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t t2sim::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t t2sim::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void t2sim::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  // multiple_particles=false;

  fChain->SetBranchAddress("particle_name", particle_name, &b_pname);
  fChain->SetBranchAddress("time", &time, &b_time);
  fChain->SetBranchAddress("energy", &energy, &b_energy);
  fChain->SetBranchAddress("position", position, &b_emiPosition);
  fChain->SetBranchAddress("direction", direction, &b_emiDirection);
  fChain->SetBranchAddress("polarization", polarization, &b_polDirection);
  // fChain->SetBranchAddress("multiple_particles", &multiple_particles,
  // &b_multiple_particles);
  Notify();
}

Bool_t t2sim::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void t2sim::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t t2sim::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void t2sim::Loop() {
  //   In a ROOT session, you can do:
  //      Root > .L t2sim.C
  //      Root > t2sim t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  // by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0)
    return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}
