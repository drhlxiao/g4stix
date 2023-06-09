/***************************************************************
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/
#include "AnalysisManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4UnitsTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNamed.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"

const G4double highVoltage =
300;                        // CdTe HV is 300 during the nominal operations
const G4double ENOISE = 0.442;  // keV
const G4double threshold = 4;
/* At 30 keV, the energy resolution is 2%  * 30 keV= 0.6 keV
 * using fano factor, one could know the intrinsic resolution of CdTe is
 * rho=(pairs * 0.15)/pairs =0.0047 absolute resolution is rho*30 = 0.14 keV ADC
 * resolution is 0.5 ADC channel,this is equivalent  to 0.5 /2.3 =0.2 keV for
 * calibration spectrum electronics noise = sqrt(0.6*0.6 -0.14*0.14 - 0.2*0.2)
 */

const G4double CdTe_SURFACE_X =
13.362385;  // surface x-coordinates of CdTe detectors
const G4double PY_ORIGIN = 103.1;
const G4double PZ_ORIGIN = 127.5;
const G4double PAIR_CREATION_ENERGY = 4.43;
const G4double FANO_FACTOR = 0.15;
// CdTe fano factor is 0.15 according to
// https://www.researchgate.net/figure/Fano-factor-for-different-semiconductor-at-room-temperature_tbl5_343053397
double energyRanges[] = {0,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13,
	14, 15, 16, 18, 20, 22, 25, 28,  32,  36,  40,
	45, 50, 56, 63, 70, 76, 84, 100, 120, 150, 250};
int getScienceBin(Double_t energy) {
	if (energy < 4) {
		return 0;
	} else if (energy >= 150) {
		return 31;
	} else {
		for (int i = 0; i < 31; i++) {
			if (energy >= energyRanges[i] && energy < energyRanges[i + 1]) return i;
		}
	}
	return 31;
}
AnalysisManager *AnalysisManager::fManager = 0;
AnalysisManager *AnalysisManager::GetInstance() {
	if (!fManager) {
		fManager = new AnalysisManager();
	}
	return fManager;
}
AnalysisManager::AnalysisManager() {
	numKilled = 0;
	numEventIn = 0;
	numEventOut = 0;
	killTracksEnteringGrids = false;
}
void AnalysisManager::CopyMacrosToROOT(TFile *f, TString &macfilename) {
	G4cout << "Copying macros from file " << macfilename << " to root file"
		<< G4endl;
	if (macfilename == "" || !f) return;
	std::ifstream infile(macfilename.Data());
	if (!infile.good()) {
		G4cout << "can not open the macro file, existing..." << G4endl;
		return;
	}
	G4String macros = Form("\nCommand:%s\n", commandLine.data());
	macros += Form("\nMacro filename: %s\n ", macroFilename.Data());
	std::string line;
	G4cout << "Macros read from the macro file:" << G4endl;
	while (std::getline(infile, line)) {
		macros += line + "\n";
		G4cout << line << G4endl;
	}
	TNamed cmd;
	cmd.SetTitle(macros);
	f->cd();
	cmd.Write("metadata");
	infile.close();
}

void AnalysisManager::InitROOT() {
	rootFile = new TFile(outputFilename.Data(), "recreate");
	evtTree = new TTree("events", "events");
	evtTree->Branch("edep", edepSum, Form("edep[%d]/D", NUM_CHANNELS));
	evtTree->Branch("sci", sci, Form("sci[%d]/D", NUM_CHANNELS));
	evtTree->Branch("collected", collectedEdepSum,
			Form("collected[%d]/D", NUM_CHANNELS));
	evtTree->Branch("charge", edepWithoutNoise,
			Form("charge[%d]/D", NUM_CHANNELS));
	evtTree->Branch("charge2", collectedEdepSumRealistic,
			Form("charge2[%d]/D", NUM_CHANNELS));
	evtTree->Branch("eventID", &eventID, "eventID/I");
	evtTree->Branch("E0", &gunEnergy, "E0/D");
	evtTree->Branch("gunPos", gunPosition, Form("gunPos[%d]/D", 3));
	evtTree->Branch("gunVec", gunDirection, Form("gunVec[%d]/D", 3));
	evtTree->Branch("numTracks", &itrack, Form("numTracks/I"));
	evtTree->Branch("nHits", nHits, Form("nHits[32]/I"));
	evtTree->Branch("hitx", hitx, Form("hitx[%d]/D", MAX_TRACKS));
	evtTree->Branch("hity", hity, Form("hity[%d]/D", MAX_TRACKS));
	evtTree->Branch("hitz", hitz, Form("hitz[%d]/D", MAX_TRACKS));
	evtTree->Branch("parent", parent, Form("parent[%d]/I", MAX_TRACKS));
	evtTree->Branch("pdg", pdg, Form("pdg[%d]/I", MAX_TRACKS));
	evtTree->Branch("pixel", hitPixelID, Form("pixel[%d]/I", MAX_TRACKS));
	evtTree->Branch("energy", energy, Form("energy[%d]/D", MAX_TRACKS));
	evtTree->Branch("time", time, Form("time[%d]/D", MAX_TRACKS));
	primTree = new TTree("source", "source");
	primTree->Branch("pos", gunPosition, "pos[3]/D");
	primTree->Branch("vec", gunDirection, "vec[3]/D");
	primTree->Branch("E0", &gunEnergy, "E0/D");

	c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
/*	for (int i = 0; i < NUM_CHANNELS; i++) {
		hd[i] = new TH1F(Form("hd%d", i),
				Form("Spectrum of pixel %d energy depositions; Energy "
					"deposition (keV); Counts",
					i),
				200, 0, 500);
	}
	*/
	for (int i = 0; i < 33; i++) {
		hpat[i] = new TH1F(Form("hDetCntPat_%d", i),
				"Detector count pattern; Pixel #; counts;", 12, 0, 12);
		hpatsum[i] = new TH1F(
				Form("hDetCntSummedPat_%d", i),
				"Pattern of Big+small summed counts; Pixel column ; Total counts;", 4,
				0, 4);

		hRealSci[i] = new TH1F(Form("hRealSci%d", i),
				Form("Recorded energy spectrum - Rebinned to SCI "
					"channels (D%d); Energy (keV)",
					i),
				32, energyRanges);
		hEdepSci[i] = new TH1F(Form("hEdepSci%d", i),
				Form("Deposited energy spectrum - Rebinned to SCI "
					"channels (D%d); Energy (keV)",
					i),
				32, energyRanges);
		hReal[i] = new TH1F(
				Form("hReal%d", i),
				Form("Recorded energy spectrum  (D%d); Energy (keV)", i), 300, 0, 150);
		hEdep[i] = new TH1F(
				Form("hEdep%d", i),
				Form("Deposited energy spectrum (D%d); Energy (keV)", i), 300, 0, 150);

		hRealSciSingleHit[i] =
			new TH1F(Form("hRealSciSingleHit%d", i),
					Form("Recorded energy spectrum - Rebinned to SCI channels "
						"(D%d); Energy (keV)",
						i),
					32, energyRanges);
		hEdepSciSingleHit[i] =
			new TH1F(Form("hEdepSciSingleHit%d", i),
					Form("Deposited energy spectrum - Rebinned to SCI channels "
						"(D%d); Energy (keV)",
						i),
					32, energyRanges);
		hRealSingleHit[i] = new TH1F(
				Form("hRealSingleHit%d", i),
				Form("Recorded energy spectrum  (D%d); Energy (keV)", i), 300, 0, 150);
		hEdepSingleHit[i] = new TH1F(
				Form("hEdepSingleHit%d", i),
				Form("Deposited energy spectrum (D%d); Energy (keV)", i), 300, 0, 150);
	}
	h2xy = new TH2F("h2xy", "Locations of hits; X (mm); Y(mm)", 1800, -90, 90,
			1800, -90, 90);
	hz = new TH1F("h1depth", "Energy deposition depth; Depth (mm); Counts;", 100,
			0, 1);
	hEdepSum = new TH1F("hEdepSum",
			"Detector summed energy spectrum; Energy (keV); Counts;",
			200, 0, 100);
	hdc = new TH1F("h1DetTotCnts", "Detector total counts; Detector ID; counts;",
			32, 0, 32);
	hpc = new TH1F("h1PixelTotCnts", "Pixel total counts; Pixel ID; counts;", 384,
			0, 384);
	hcol = new TH1F(
			"h1ChargeColEff",
			"Distribution of Charge collection efficiency ; Efficiency; Counts ;",
			200, 0, 1);

	hEdepSum->SetCanExtend(TH1::kXaxis);

	// for ROOT version >6.0
}

//////////////////////////////////////////////////////////////////////////

void AnalysisManager::ProcessRun(const G4Run *run) {
	CloseROOT();
	G4cout << "Events entered detectors:" << numEventIn << G4endl;
	G4cout << "Events escaped from detectors:" << numEventOut << G4endl;
	G4cout << "Events captured:" << numEventIn - numEventOut << G4endl;
}

/// EventAction
void AnalysisManager::InitEvent(const G4Event *event) {
	for (G4int i = 0; i < NUM_CHANNELS; i++) {
		edepSum[i] = 0.0;
		collectedEdepSum[i] = 0;
		edepWithoutNoise[i] = 0;
		collectedEdepSumRealistic[i] = 0;
		sci[i] = -1;
	}
	for (int i = 0; i < 32; i++) nHits[i] = 0;

	for (int i = 0; i < MAX_TRACKS; i++) {
		pdg[i] = 0;
		energy[i] = 0;
		hitPixelID[i] = 0;
		parent[i] = 0;
		time[i] = 0;
		hitx[i] = 0;
		hity[i] = 0;
		hitz[i] = 0;
	}
	itrack = 0;
}
void AnalysisManager::ProcessEvent(const G4Event *event) {
	eventID = event->GetEventID();
	G4bool toFill = false;
	for (int i = 0; i < NUM_CHANNELS; i++) {
		if (edepSum[i] > 0) {
			hEdepSum->Fill(edepSum[i]);
			//hd[i]->Fill(edepSum[i]);
			toFill = true;
			G4double sigma =
				collectedEdepSum[i] * GetEnergyResolution(collectedEdepSum[i]);
			// std. Deviation of charge, in units of eV

			edepWithoutNoise[i] = gRandom->Gaus(collectedEdepSum[i], sigma) *
				PAIR_CREATION_ENERGY /
				1000;  // convert back to keV, we randomize it to
					   // smear the energy resolution
					   //
			collectedEdepSumRealistic[i] = gRandom->Gaus(edepWithoutNoise[i], ENOISE);
			// we asssue the electronics noise

			int detIdx = i / 12;
			int pixIdx = i % 12;
			hpat[detIdx]->Fill(pixIdx);
			hpatsum[detIdx]->Fill(pixIdx / 4);
			hpc->Fill(i);
			hdc->Fill(detIdx);

			if (collectedEdepSumRealistic[i] > threshold) {
				nHits[detIdx]++;
			}
			sci[i] = getScienceBin(collectedEdepSumRealistic[i]);

			hEdep[detIdx]->Fill(edepSum[i]);
			hReal[detIdx]->Fill(collectedEdepSumRealistic[i]);
			//not binned to stix 
			hEdep[32]->Fill(edepSum[i]);
			hReal[32]->Fill(collectedEdepSumRealistic[i]);
			//histograms to store spectra of events of all detectors

			hEdepSci[detIdx]->Fill(edepSum[i]);
			hRealSci[detIdx]->Fill(collectedEdepSumRealistic[i]);
			//
			hEdepSci[32]->Fill(edepSum[i]);
			hRealSci[32]->Fill(collectedEdepSumRealistic[i]);
			// summed spectrum
		}
	}
	//iterate over detectors
	for (int i = 0; i < 32; i++) {
		if (nHits[i] != 1) continue;
		// only accepts single hit events
		for (int j = 0; j < 12; j++) {
			int ch = i * 12 + j;
			if (edepSum[ch] < 0) continue;

			hEdepSingleHit[i]->Fill(edepSum[ch]);
			hRealSingleHit[i]->Fill(collectedEdepSumRealistic[ch]);
			hEdepSciSingleHit[i]->Fill(edepSum[ch]);
			hRealSciSingleHit[i]->Fill(collectedEdepSumRealistic[ch]);
			//stix energy bins

			hEdepSingleHit[32]->Fill(edepSum[ch]);
			hRealSingleHit[32]->Fill(collectedEdepSumRealistic[ch]);
			hEdepSciSingleHit[32]->Fill(edepSum[ch]);
			hRealSciSingleHit[32]->Fill(collectedEdepSumRealistic[ch]);
			//sum spectrum 

		}
	}

	///	hdc->Fill(detectorID);
	// toFill=true;

	// write particle information to the tree
	G4RunManager *runManager = G4RunManager::GetRunManager();
	PrimaryGeneratorAction *primaryAction =
		(PrimaryGeneratorAction *)runManager->GetUserPrimaryGeneratorAction();
	G4ThreeVector position, direction;
	G4double energy;
	primaryAction->GetGPS(position, direction, energy);

	// G4cout<<"#energy :"<<energy<<G4endl;
	//    G4cout<<"position:"<<position[0]<<" "<<position[1]<<"
	//    "<<position[2]<<G4endl;
	gunPosition[0] = position.getX() / mm;
	gunPosition[1] = position.getY() / mm;
	gunPosition[2] = position.getZ() / mm;
	gunDirection[0] = direction.getX();
	gunDirection[1] = direction.getY();
	gunDirection[2] = direction.getZ();
	gunEnergy = energy;
	if (eventID < 1e5) primTree->Fill();
	if (toFill) evtTree->Fill();
}
// Stepping Action

G4double AnalysisManager::ComputeCollectedEnergy(G4ThreeVector &pos,
		G4double edep) {
	// Hecht equation, see "Recent Progress in CdTe and CdZnTe Detectors"
	// Tadayuki Takahashi and Shin Watanabe
	// Oliver's paper
	// Energy must be in units of eV
	// Spectral signature of near-surface damage in CdTe X-ray detectors
	//
	G4double xx = pos.x() / mm;

	// G4cout<<xx<<G4endl;

	G4double depth = xx - CdTe_SURFACE_X;  // depth from cathode
	hz->Fill(depth);
	G4double freePathElectron = 1100 * 100 * 3e-6 * highVoltage;
	G4double freePathHoles = 100 * 100 * 2e-6 * highVoltage;

	G4double nPairs = (edep / 4.6);

	// see Oliver's paper
	G4double detThickness = 1;
	// G4double charge = nPairs;
	G4double factor =
		(1 - exp(-(detThickness - depth) / freePathElectron)) * freePathElectron /
		detThickness +
		(freePathHoles / detThickness) * (1.0 - exp(-depth / freePathHoles));
	// G4double charge=(1-exp((depth-detThickness)/freePathElectron));
	hcol->Fill(factor);
	// G4cout<<freePathElectron<<" "<<depth<<" "<<detThickness<<"
	// "<<charge<<G4endl;
	return factor * nPairs;
}
G4double AnalysisManager::GetEnergyResolution(G4double edep) {
	//	Recent Progress in CdTe and CdZnTe Detectors
	// Tadayuki Takahashi and Shin Watanabe
	// energy in units of keV
	G4double rho = 0;
	if (edep > 0) {
		G4double numPairs = edep / PAIR_CREATION_ENERGY * 1000;
		// keV to ev
		//  see oliver's paper
		G4double numSigma2 = FANO_FACTOR * numPairs;
		// sigma^2 = fano * n  according to the text-book
		rho = sqrt(numSigma2) / numPairs;
	}
	return rho;
}

void AnalysisManager::ProcessStep(const G4Step *aStep) {
	const G4Track *track = aStep->GetTrack();
	G4String volName;
	if (track->GetVolume()) volName = track->GetVolume()->GetName();

	G4int detIdx = -1;
	G4int detectorID = -1;
	G4int pixelID = -1;
	if (killTracksEnteringGrids && volName == "gridStrip") {
		aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		numKilled++;
	}

	if (volName == "pixel") {
		detIdx = track->GetTouchable()->GetCopyNumber(3);
		pixelID = track->GetTouchable()->GetCopyNumber(0);
		G4double edep = aStep->GetTotalEnergyDeposit() / keV;

		// G4cout<<track->GetTouchable()->GetCopyNumber(0)<<","
		//  <<track->GetTouchable()->GetCopyNumber(1)<<","
		//  <<track->GetTouchable()->GetCopyNumber(3)<<" EDEP: "<<edep<<G4endl;
		//	*/

		detectorID = detIdx * 12 + pixelID;

		// TString detName=volName.data();

		if (edep > 0.0) {
			AddEnergy(detectorID, edep);
		G4ThreeVector postPos= aStep->GetPostStepPoint()->GetPosition();
			G4double charge = ComputeCollectedEnergy(postPos, edep * 1000);
			AddCollectedEnergy(detectorID, charge);
		}
		G4StepPoint *preStep = aStep->GetPreStepPoint();
		G4ThreeVector prePos= aStep->GetPostStepPoint()->GetPosition();
		G4double px = prePos.x() / mm;
		G4double py = prePos.y() / mm;
		G4double pz = prePos.z() / mm;
		if (preStep->GetStepStatus() == fGeomBoundary) {
			h2xy->Fill(py - PY_ORIGIN, pz - PZ_ORIGIN);
		}

		if (itrack >= MAX_TRACKS) return;

		if (preStep->GetStepStatus() == fGeomBoundary) {
			hitx[itrack] = px;
			hity[itrack] = py;
			hitz[itrack] = pz;

			G4int parentID = track->GetParentID();  //
			pdg[itrack] = track->GetDefinition()->GetPDGEncoding();
			energy[itrack] = preStep->GetKineticEnergy() / keV;
			hitPixelID[itrack] = detectorID;
			parent[itrack] = parentID;
			time[itrack] = preStep->GetGlobalTime() / ns;
			itrack++;
		}
	}
}

////////////////////////////////////////////////////////////////////

inline void AnalysisManager::AddEnergy(G4int detId, G4double edep) {
	if (detId > NUM_CHANNELS || detId < 0) {
		G4cout << "invalid index" << G4endl;
		return;
	}
	edepSum[detId] += edep;
}
inline void AnalysisManager::AddCollectedEnergy(G4int detId, G4double dep) {
	if (detId > NUM_CHANNELS || detId < 0) {
		G4cout << "invalid index" << G4endl;
		return;
	}
	collectedEdepSum[detId] += dep;
}

AnalysisManager::~AnalysisManager() {
	if (rootFile) delete rootFile;
	if (fManager) delete fManager;
	fManager = 0;
}

void AnalysisManager::CloseROOT() {
	rootFile->cd();
	evtTree->Write();
	primTree->Write();

	rootFile->cd();
	TDirectory *cdhist = rootFile->mkdir("hist");
	cdhist->cd();
	c1->cd();
	c1->Divide(2, 5);
//	for (int i = 0; i < NUM_CHANNELS; i++) {
//		hd[i]->Write();
//		c1->cd(i + 1);
//		hd[i]->Draw();
//	}
	for (int i = 0; i < 33; i++) {
		hRealSci[i]->Write();
		hEdepSci[i]->Write();
		hReal[i]->Write();
		hEdep[i]->Write();

		hRealSciSingleHit[i]->Write();
		hEdepSciSingleHit[i]->Write();
		hRealSingleHit[i]->Write();
		hEdepSingleHit[i]->Write();




		if (i >= 32) continue;
		hpat[i]->Write();
		hpatsum[i]->Write();
	}
	hEdepSum->Write();
	h2xy->Write();
	hdc->Write();
	hpc->Write();
	c1->Write();
	hz->Write();
	hcol->Write();
	G4cout << ">> Number of event recorded:" << evtTree->GetEntries() << G4endl;
	G4cout << ">> Number of track killed:" << numKilled << G4endl;
	CopyMacrosToROOT(rootFile, macroFilename);

	rootFile->Close();
}
