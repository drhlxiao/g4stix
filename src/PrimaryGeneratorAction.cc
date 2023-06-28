/***************************************************************
 * descriptions of partcile gun
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/

#include "PrimaryGeneratorAction.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "t2sim.h"
// For Random Generator
#include <CLHEP/Random/RandFlat.h>

#include <fstream>

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
//Ba133  energies
#define NUM_GAMMAS 29

//ba133 emission lines
//taken from decay.exe
const G4double prob[NUM_GAMMAS] = {
8.94,
62.05,
18.33,
7.164,
0.45,
0.645,
34.06,
2.62,
2.199,
0.74,
3.58,
0.123,
11.6,
5.99,
64.5,
34.9,
0.004,
0.22,
0.15,
0.54,
1.19,
0.048,
0.93,
0.56,
3.8,
6,
0.66,
0.11,
0.24};
const G4double gammaEnergy[NUM_GAMMAS] = {
383.85,
356.02,
302.85,
276.4,
223.23,
160.61,
81,
79.61,
53.16,
35.907,
35.818,
35.252,
34.987,
34.92,
30.973,
30.625,
30.27,
5.553,
5.542,
5.281,
4.934,
4.781,
4.717,
4.649,
4.62,
4.286,
4.272,
4.142,
3.795
};

PrimaryGeneratorAction::PrimaryGeneratorAction()
	: sourceType(0), fTree(NULL), fFile(NULL), nEntries(0), iEntry(0) {
		// fParticleGunMessenger = new PrimaryGeneratorMessenger(this);
		fHistoEnergy = NULL;

		// particle source
		particleSource = "";
		G4int n_particle = 1;
		fParticleGun = new G4ParticleGun(n_particle);
		// default particle kinematic
		particleTable = G4ParticleTable::GetParticleTable();
		G4String particleName;
		G4ParticleDefinition *particle =
			particleTable->FindParticle(particleName = "gamma");
		fParticleGun->SetParticleDefinition(particle);
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
		fParticleGun->SetParticleEnergy(30. * keV);
		//also defined in the DetectorConstruction
		fParticleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, 15 * cm));
		// particle source
		fParticleSource = new G4GeneralParticleSource();
		// fParticleSource->SetParticleEnergy(50*keV);
		histoEnergy = NULL;
	}
PrimaryGeneratorAction::~PrimaryGeneratorAction() {
	delete fParticleSource;
	// delete fParticleGunMessenger;
	delete fParticleGun;
	if (fHistoEnergy) fHistoEnergy->Close();
}

void PrimaryGeneratorAction::InitParticleSpectrumFromROOT(G4String val) {
	// don't use this method, use /gps/hist/file instead June 7, 2023
	// because using gps with setparticleenergy doesn't work
	/*
	   particleSourceFile=val;
	   particleSource="fromROOT";
	   G4cout<<"Read spectrum file ROOT:"<<particleSourceFile<<G4endl;
	//file generated in the folder:
	/home/xiaohl/FHNW/STIX/SolarFlareAnalysis/response_matrix fHistoEnergy=new
	TFile(val.data()); histoEnergy=(TH1F*)fHistoEnergy->Get("hspec");
	if(!histoEnergy){
	G4cout<<"Can not read hspec from ROOT file"<<G4endl;
	}

*/
}
void PrimaryGeneratorAction::GetGPS(G4ThreeVector &position,
		G4ThreeVector &direction,
		G4double &energy) {
	//noted that this dosen't work with older versions of geant4
	if (particleSource == "") {
		position = fParticleSource->GetParticlePosition();
		direction = fParticleSource->GetParticleMomentumDirection();
		energy = fParticleSource->GetParticleEnergy() / keV;
	} else {
		position = fParticleGun->GetParticlePosition();
		direction = fParticleGun->GetParticleMomentumDirection();
		energy = fParticleGun->GetParticleEnergy() / keV;
	}
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
	G4double energy;
	G4double rnd;
	G4double phi, theta,px,py,pz, radius;
	G4double sourcePlaneRadius=157/2;
	//

	//particleTable->FindParticle("gamma");
	//fParticleGun->SetParticleDefinition(particle);

	if (particleSource == "Ba133") {
		G4int numGamma = 0;

		phi= 2.0 * 3.14159 * (G4UniformRand());
		radius= sourcePlaneRadius* sqrt( G4UniformRand() )*mm;
		G4double shiftY=radius*cos(phi);
		G4double shiftZ=radius*sin(phi);
		G4double extraMargin=-2*mm;
		//G4ThreeVector calFoilKaptonCenter(9.251*mm+ extraMargin-1*mm, 104.074*mm, 126.709*mm);
		G4ThreeVector sourcePos(9.251*mm+ extraMargin-1*mm, 104.074*mm+shiftY, 126.709*mm+shiftZ);
		//generate source on the surface of Kapton foil
		fParticleGun->SetParticlePosition(sourcePos);

		while (numGamma == 0) {
			for (G4int i = 0; i < NUM_GAMMAS; i++) {
				rnd = 100 * G4UniformRand();
				if (rnd < prob[i]) {
					phi = 2.0 * 3.14159 * (G4UniformRand());
					theta = acos(2 * G4UniformRand() - 1);
					px = cos(phi) * sin(theta);
					py = sin(phi) * sin(theta);
					pz = cos(theta);

					fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
					energy = gammaEnergy[i];
					fParticleGun->SetParticleEnergy(energy * keV);
					fParticleGun->GeneratePrimaryVertex(anEvent);
					numGamma++;
				}
			}
		}

		// G4cout<<numGamma<<G4endl;

		iEntry++;
		if (iEntry % 10000 == 0) {
			G4cout << "ba133 source" << G4endl;
		}
	}

		else {
			fParticleSource->GeneratePrimaryVertex(anEvent);
		}
}

G4bool PrimaryGeneratorAction::InitFile() {}
