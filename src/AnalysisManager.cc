/***************************************************************
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TString.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TNamed.h"

#include "G4Event.hh"
#include "G4TrackStatus.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"

#include "PrimaryGeneratorAction.hh"

#include "AnalysisManager.hh"

const G4double highVoltage = 200; // 200 V
const G4double ENOISE = 0.7;      // keV
const G4double CdTe_SURFACE_X= 13.362385;     //surface x-coordinates of CdTe detectors
const G4double PY_ORIGIN=103.1;
const G4double PZ_ORIGIN=127.5;
const G4double PAIR_CREATION_ENERGY=4.43;
const G4double FANO_FACTOR=0.1;


AnalysisManager *AnalysisManager::fManager = 0;
AnalysisManager *AnalysisManager::GetInstance() {
	if (!fManager) {
		fManager = new AnalysisManager();
	}
	return fManager;
}
AnalysisManager::AnalysisManager() {

	numKilled=0;
	fnEventIn = 0;
	fnEventOut = 0;
	killTracksEnteringGrids=false;
}
void AnalysisManager::CopyMacrosToROOT(TFile *f, TString &macfilename)
{
	G4cout<<"Copying macros from file "<<macfilename<<" to root file"<<G4endl;
	if(macfilename=="" ||!f)return;
	std::ifstream infile(macfilename.Data()); 
	if(!infile.good())
	{ 
		G4cout<<"can not open the macro file, existing..."<<G4endl; 
		return ;
	} 
	G4String macros=Form("Macro filename: %s\n ",macroFilename.Data());
	std::string line;
	G4cout<<"Macros read from the macro file:"<<G4endl;
	while(std::getline(infile,line))
	{
		macros+=line+"\n";
		G4cout<<line<<G4endl;
	}
	TNamed cmd;
	cmd.SetTitle(macros);
	f->cd();
	cmd.Write("metadata");
	infile.close(); 
}

void AnalysisManager::CreateTree() {
	fTTree = new TTree("events", "events");
	fTTree->Branch("edep", fEdepSum, Form("edep[%d]/D", NUM_CHANNELS));
	fTTree->Branch("collected", fCollectedEdepSum,
			Form("collected[%d]/D", NUM_CHANNELS));
	fTTree->Branch("charge", fEdepWithoutNoise,
			Form("charge[%d]/D", NUM_CHANNELS));
	fTTree->Branch("charge2", fCollectedEdepSumRealistic,
			Form("charge2[%d]/D", NUM_CHANNELS));
	fTTree->Branch("eventID", &fEventID, "eventID/I");
	fTTree->Branch("gunPos", gunPosition,
			Form("gunPos[%d]/D", 3));
	fTTree->Branch("gunVec", gunDirection,
			Form("gunVec[%d]/D", 3));
	fTTree->Branch("E0", &gunEnergy,
			Form("E0/D"));
	fTTree->Branch("numTracks", &itrack, Form("numTracks/I"));
	fTTree->Branch("hitx", hitx, 	Form("hitx[%d]/D", MAX_TRACKS));
	fTTree->Branch("hity", hity, 	Form("hity[%d]/D", MAX_TRACKS));
	fTTree->Branch("hitz", hitz, 	Form("hitz[%d]/D", MAX_TRACKS));
	fTTree->Branch("parent", parent, 	Form("parent[%d]/I", MAX_TRACKS));
	fTTree->Branch("pdg", pdg, 	Form("pdg[%d]/I", MAX_TRACKS));
	fTTree->Branch("pixel", hitPixelID, 	Form("pixel[%d]/I", MAX_TRACKS));
	fTTree->Branch("energy", energy, 	Form("energy[%d]/D", MAX_TRACKS));
	fTTree->Branch("time", time, 	Form("time[%d]/D", MAX_TRACKS));
	fTSource=new TTree("source","source");
	fTSource->Branch("pos",gunPosition,"pos[3]/D");    
	fTSource->Branch("vec",gunDirection,"vec[3]/D");    
	fTSource->Branch("E0",&gunEnergy,"E0/D");    

	c1 = new TCanvas("c1", "c1", 10, 10, 800, 800);
	for (int i = 0; i < NUM_CHANNELS; i++) {
		hd[i] = new TH1F(Form("hd%d", i),
				Form("Spectrum of pixel %d energy depositions; Energy deposition (keV); Counts", i),
				200, 0, 500);
	}
	for(int i=0;i<32;i++){
		hpat[i]=new TH1F(Form("hDetCntPat_%d",i), "Detector count pattern; Pixel #; counts;", 12, 0, 12);
		hpatsum[i]=new TH1F(Form("hDetCntSummedPat_%d",i), "Pattern of Big+small summed counts; Pixel column ; Total counts;", 4, 0, 4);
	}
	h2xy=new TH2F("h2xy","Locations of hits; X (mm); Y(mm)", 1800,-90,90, 1800,-90,90);
	hz = new TH1F("h1depth", "Energy deposition depth; Depth (mm); Counts;", 100, 0, 1);
	hEdepSum = new TH1F("hEdepSum", "Detector summed energy spectrum; Energy (keV); Counts;", 200, 0, 100);
	hdc = new TH1F("h1DetTotCnts", "Detector total counts; Detector ID; counts;", 32, 0, 32);
	hpc = new TH1F("h1PixelTotCnts", "Pixel total counts; Pixel ID; counts;", 384, 0, 384);
	hcol = new TH1F("h1ChargeColEff", "Distribution of Charge collection efficiency ; Efficiency; Counts ;", 200, 0, 1);

	hEdepSum->SetCanExtend(TH1::kXaxis);
	//for ROOT version >6.0
}

//////////////////////////////////////////////////////////////////////////

void AnalysisManager::initRun(const G4Run *run) {

	OpenFile();
	CreateTree();
}

void AnalysisManager::processRun(const G4Run *run) {
	flush();
	G4cout << "Events entered detectors:" << fnEventIn << G4endl;
	G4cout << "Events escaped from detectors:" << fnEventOut << G4endl;
	G4cout << "Events captured:" << fnEventIn - fnEventOut << G4endl;
}

/// EventAction
void AnalysisManager::initEvent(const G4Event *event) {
	for (G4int i = 0; i < NUM_CHANNELS; i++) {
		fEdepSum[i] = 0.0;
		fCollectedEdepSum[i] = 0;
		fEdepWithoutNoise[i] = 0;
		fCollectedEdepSumRealistic[i] = 0;
	}

	for(int i=0;i<MAX_TRACKS;i++){
		pdg[i]=0;
		energy[i]=0;
		hitPixelID[i]=0;
		parent[i]=0;
		time[i]=0;
		hitx[i]=0;
		hity[i]=0;
		hitz[i]=0;
	}
	itrack=0;
}
void AnalysisManager::processEvent(const G4Event *event) {
	fEventID = event->GetEventID();
	G4bool toFill = false;
	for (int i = 0; i < NUM_CHANNELS; i++) {
		if (fEdepSum[i] > 0) {
			hEdepSum->Fill(fEdepSum[i]);
			hd[i]->Fill(fEdepSum[i]);
			toFill = true;
			G4double sigma =
				fCollectedEdepSum[i] * GetEnergyResolution(fCollectedEdepSum[i]);
			G4double real = gRandom->Gaus(fCollectedEdepSum[i], sigma) * PAIR_CREATION_ENERGY /
				1000; // convert back to keV
			fEdepWithoutNoise[i] = real;
			fCollectedEdepSumRealistic[i] = gRandom->Gaus(real, ENOISE);
			int detIdx=i/12;
			int pixIdx=i%12;
			hpat[detIdx]->Fill(pixIdx);
			hpatsum[detIdx]->Fill(pixIdx/4);
			hpc->Fill(i);
			hdc->Fill(detIdx);
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

	//    G4cout<<"#energy :"<<energy<<G4endl;
	//   G4cout<<"position:"<<position[0]<<" "<<position[1]<<"
	//   "<<position[2]<<G4endl;
	gunPosition[0] = position.getX()/mm;
	gunPosition[1] = position.getY()/mm;
	gunPosition[2] = position.getZ()/mm;
	gunDirection[0] = direction.getX();
	gunDirection[1] = direction.getY();
	gunDirection[2] = direction.getZ();
	gunEnergy = energy;
	if(fEventID<1e5)fTSource->Fill();
	if (toFill)fTTree->Fill();
}
// Stepping Action

G4double AnalysisManager::ComputeCollectedEnergy(G4ThreeVector &pos,
		G4double energy) {
	// Hecht equation, see "Recent Progress in CdTe and CdZnTe Detectors"
	// Tadayuki Takahashi and Shin Watanabe
	// Oliver's paper
	// Energy must be in units of eV
	// Spectral signature of near-surface damage in CdTe X-ray detectors
	//
	G4double xx = pos.x() / mm;

	//G4cout<<xx<<G4endl;

	G4double depth = xx - CdTe_SURFACE_X; // depth from cathode
	hz->Fill(depth);
	G4double freePathElectron = 1100 * 100 * 3e-6 * highVoltage;
	G4double freePathHoles = 100 * 100 * 2e-6 * highVoltage;

	G4double nPairs = (energy / 4.6);

	// see Oliver's paper
	G4double detThickness = 1;
	G4double charge = nPairs;
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
G4double AnalysisManager::GetEnergyResolution(G4double energy) {
	//	Recent Progress in CdTe and CdZnTe Detectors
	// Tadayuki Takahashi and Shin Watanabe
	// energy in units of keV
	G4double rho = 0;
	if (energy > 0) {
		G4double numPairs = energy / PAIR_CREATION_ENERGY * 1000;
		// see oliver's paper
		G4double numSigma2 = FANO_FACTOR* numPairs; 
		rho = sqrt(numSigma2) / numPairs;           
	}
	return rho;
}

void AnalysisManager::processStep(const G4Step *aStep) {

	const G4Track *track = aStep->GetTrack();
	G4String volName;
	if (track->GetVolume())volName = track->GetVolume()->GetName();

	G4int detIdx= -1;
	G4int detectorID= -1;
	G4int pixelID= -1;
	if(killTracksEnteringGrids&& volName=="gridStrip" ){
		aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		numKilled++;
	}

	if (volName=="pixel")
	{
		detIdx= track->GetTouchable()->GetCopyNumber(2);
		pixelID= track->GetTouchable()->GetCopyNumber(0);
		G4double edep = aStep->GetTotalEnergyDeposit() / keV;

		/*G4cout<<track->GetTouchable()->GetCopyNumber(0)<<","
		  <<track->GetTouchable()->GetCopyNumber(1)<<","
		  <<track->GetTouchable()->GetCopyNumber(2)<<" EDEP: "<<edep<<G4endl;
		  */
		//	*/

		detectorID=detIdx*12+pixelID;

		// TString detName=volName.data();

		if (edep > 0.0) {

			AddEnergy(detectorID, edep);
			G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
			G4double charge = ComputeCollectedEnergy(position, edep * 1000);
			AddCollectedEnergy(detectorID, charge);
		}
		G4StepPoint* preStep= aStep->GetPreStepPoint();
		G4ThreeVector	position = aStep->GetPostStepPoint()-> GetPosition();
		G4double px=position.x()/mm;
		G4double py=position.y()/mm;
		G4double pz=position.z()/mm;
		if (preStep->GetStepStatus() == fGeomBoundary) 
		{
			h2xy->Fill(py- PY_ORIGIN,pz -PZ_ORIGIN);
		}



		if(itrack>=MAX_TRACKS)return;

		if (preStep->GetStepStatus() == fGeomBoundary) 
		{
			hitx[itrack]=px;
			hity[itrack]=py;
			hitz[itrack]=pz;

			G4int parentID=track->GetParentID(); //
			pdg[itrack]=track->GetDefinition()->GetPDGEncoding();
			energy[itrack]=preStep->GetKineticEnergy()/keV;
			hitPixelID[itrack]=detectorID;
			parent[itrack]=parentID;
			time[itrack]=preStep->GetGlobalTime()/ns;
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
	fEdepSum[detId] += edep;
}
inline void AnalysisManager::AddCollectedEnergy(G4int detId, G4double dep) {
	if (detId > NUM_CHANNELS || detId < 0) {
		G4cout << "invalid index" << G4endl;
		return;
	}
	fCollectedEdepSum[detId] += dep;
}

AnalysisManager::~AnalysisManager() {
	if(fTFile)delete fTFile;
	if (fManager)delete fManager;
	fManager = 0;
}
void AnalysisManager::OpenFile() {
	fTFile = new TFile(outputFilename.Data(), "recreate");
}

void AnalysisManager::flush() {
	fTFile->cd();
	fTTree->Write();
	fTSource->Write();

	fTFile->cd();
	TDirectory *cdhist=fTFile->mkdir("hist");
	cdhist->cd();
	c1->cd();
	c1->Divide(2, 5);
	for (int i = 0; i < NUM_CHANNELS; i++) {
		hd[i]->Write();
		c1->cd(i + 1);
		hd[i]->Draw();
	}
	for(int i=0;i<32;i++){
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
	G4cout<<">> Number of event recorded:"<<fTTree->GetEntries()<<G4endl;
	G4cout<<">> Number of track killed:"<<numKilled<<G4endl;
	CopyMacrosToROOT(fTFile, macroFilename);

	fTFile->Close();
}
