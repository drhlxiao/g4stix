//
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "globals.hh"

//#include "G4Step.hh"
//#include "G4Event.hh"
//#include "G4Run.hh"
#include "TString.h"

#define NUM_CHANNELS 12*32
class G4Run;
class G4Event;
class G4Step;

class TCanvas;
class TH1F;
class TH2F;
class TFile;
class TTree;
const int MAX_TRACKS=30;

class AnalysisManager {
	public:
		static AnalysisManager *GetInstance();
		static void Dispose();

		void OpenFile();
		void CreateTree();
		void flush();
		void SetEventID(G4int eventid) { eventID = eventid; }
		void SetOutputFileName(TString filen) { outputFilename = filen; }
		//void SetAttenuatorStatus(G4bool s){attenuatorIn=s;}

		void AddEnergy(G4int detId, G4double edep);
		void AddCollectedEnergy(G4int detId, G4double edep);

		void BeginOfEventAction(const G4Event *event);
		void EndOfEventAction(const G4Event *event);
		void SteppingAction(const G4Step *aStep);
		void BeginOfRunAction(const G4Run *);
		void EndOfRunAction(const G4Run *);
		void CopyMacrosToROOT(TFile *f, TString &);
		G4double ComputeCollectedEnergy(G4ThreeVector &pos, G4double Ek);
		G4double GetEnergyResolution(G4double Ek);
		void SetMacroFileName(G4String &name){macroFilename=name;}

		AnalysisManager();
		~AnalysisManager();

	private:
		static AnalysisManager *fManager;

		TString outputFilename;
		TString macroFilename;

		TFile *fTFile;
		TTree *fTTree;
		long long totalHits;


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
		TH1F *hpat[32];
		TH1F *hpatsum[32];
		TH2F *h2xy;
		TTree *fTSource;

		G4int eventID;

		G4double fEdepSum[NUM_CHANNELS];
		G4double fCollectedEdepSum[NUM_CHANNELS];
		G4double fEdepWithoutNoise[NUM_CHANNELS];
		G4double fCollectedEdepSumRealistic[NUM_CHANNELS];
		G4int fEventID;

		G4int parent[MAX_TRACKS];
		G4double energy[MAX_TRACKS];
		G4double hitx[MAX_TRACKS];
		G4double hity[MAX_TRACKS];
		G4double hitz[MAX_TRACKS];
		G4int hitPixelID[MAX_TRACKS];
		G4int pdg[MAX_TRACKS];
		G4double time[MAX_TRACKS];


		G4int fnEventIn;
		G4int fnEventOut;
};

#endif
