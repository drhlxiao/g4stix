//
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include <vector>

#include "globals.hh"

//#include "G4Step.hh"
//#include "G4Event.hh"
//#include "G4Run.hh"
#include "TString.h"
#include "G4ThreeVector.hh"

#define NUM_CHANNELS 12 * 32
class G4Run;
class G4Event;
class G4Step;


class TCanvas;
class TH1F;
class TH2F;
class TFile;
class TTree;
const int MAX_TRACKS = 30;

class AnalysisManager {
 public:
  static AnalysisManager *GetInstance();
  static void Dispose();

  void CloseROOT();
  void SetEventID(G4int eventid) { eventID = eventid; }
  void SetOutputFileName(TString filen) { outputFilename = filen; }

  void SetCommandLine(G4String s) { commandLine = s; }
  void InitEvent(const G4Event *event);
  void ProcessEvent(const G4Event *event);
  void ProcessStep(const G4Step *aStep);
  void InitRun(const G4Run *) { InitROOT(); };
  void InitROOT();
  void ProcessRun(const G4Run *);

  void AddEnergy(G4int detId, G4double edep);
  void AddCollectedEnergy(G4int detId, G4double edep);
  void CopyMacrosToROOT(TFile *f, TString &);
  G4double ComputeCollectedEnergy(G4ThreeVector &pos, G4double Ek);
  G4double GetEnergyResolution(G4double Ek);
  void SetMacroFileName(G4String &name) { macroFilename = name; }
  void KillTracksInGrids() {
    killTracksEnteringGrids = true;
    G4cout << "# Tracks entering Grids will be killed" << G4endl;
  }

  AnalysisManager();
  ~AnalysisManager();

 private:
  static AnalysisManager *fManager;

  TString outputFilename;
  TString macroFilename;

  TFile *rootFile;
  TTree *evtTree;
  long long totalHits;
  G4String commandLine;

  G4int itrack;

  G4double gunPosition[3];
  G4double gunDirection[3];
  G4double gunEnergy;
  TCanvas *c1;
  TH1F *hd[NUM_CHANNELS];
  TH1F *hEdepSum;
  TH1F *hz;
  TH1F *hcol;
  TH1F *hpc;
  TH1F *hdc;
  TH1F *hpat[33];
  TH1F *hpatsum[33];

  TH1F *hEdepSci[33];
  TH1F *hRealSci[33];
  TH1F *hEdep[33];
  TH1F *hReal[33];

  TH1F *hEdepSciSingleHit[33]; //deposited energies, stix science channel, only single-hit accepted, The 32nd histogram for sum spectrum
  TH1F *hRealSciSingleHit[33];//measured energies, stix science channel, only single-hit accepted, The 32nd histogram for sum spectrum
  TH1F *hEdepSingleHit[33];
  TH1F *hRealSingleHit[33];
  // energy spectrum with single hit only

  TH2F *h2xy;
  TTree *primTree;
  G4int numKilled;

  G4double sci[NUM_CHANNELS];
  G4double edepSum[NUM_CHANNELS];
  G4double collectedEdepSum[NUM_CHANNELS];
  G4double edepWithoutNoise[NUM_CHANNELS];
  G4double collectedEdepSumRealistic[NUM_CHANNELS];
  G4int nHits[32];
  G4int eventID;
  G4bool killTracksEnteringGrids;

  G4int parent[MAX_TRACKS];
  G4double energy[MAX_TRACKS];
  G4double hitx[MAX_TRACKS];
  G4double hity[MAX_TRACKS];
  G4double hitz[MAX_TRACKS];
  G4int hitPixelID[MAX_TRACKS];
  G4int pdg[MAX_TRACKS];
  G4double time[MAX_TRACKS];

  G4int numEventIn;
  G4int numEventOut;
};

#endif
