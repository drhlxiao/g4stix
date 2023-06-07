/***************************************************************
 * descriptions of partcile gun
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "t2sim.h"

#include "PrimaryGeneratorAction.hh"
// For Random Generator
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include <CLHEP/Random/RandFlat.h>
#include <fstream>

PrimaryGeneratorAction::PrimaryGeneratorAction()
	: sourceType(0), fTree(NULL), fFile(NULL), nEntries(0), iEntry(0) {
		// fParticleGunMessenger = new PrimaryGeneratorMessenger(this);
		fHistoEnergy=NULL;

		// particle source
		particleSource="";
		G4int n_particle = 1;
		fParticleGun = new G4ParticleGun(n_particle);
		// default particle kinematic
		particleTable = G4ParticleTable::GetParticleTable();
		G4String particleName;
		G4ParticleDefinition *particle =
			particleTable->FindParticle(particleName = "proton");
		fParticleGun->SetParticleDefinition(particle);
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
		fParticleGun->SetParticleEnergy(30. * MeV);
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
	 if(fHistoEnergy)fHistoEnergy->Close();
}

void PrimaryGeneratorAction::InitParticleSpectrumFromROOT(G4String val) {
	particleSourceFile=val;
	particleSource="fromROOT";
	G4cout<<"Read spectrum file ROOT:"<<particleSourceFile<<G4endl;
	//file generated in the folder: /home/xiaohl/FHNW/STIX/SolarFlareAnalysis/response_matrix
	fHistoEnergy=new TFile(val.data());
	histoEnergy=(TH1F*)fHistoEnergy->Get("hspec");
	if(!histoEnergy){
		G4cout<<"Can not read hspec from ROOT file"<<G4endl;
	}

	

}
void PrimaryGeneratorAction::GetGPS(G4ThreeVector &position,
		G4ThreeVector &direction,
		G4double &energy) {
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
	if (particleSource == "Ba133") {
		G4int numGamma = 0;
		G4double prob[12] = {34.9, 64.5,  5.99, 11.6,  3.58,  2.2,
			2.62, 34.06, 7.1,  18.85, 62.05, 8.94};
		// http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=560133
		G4double gammaEnergy[12] = {30.6, 30.9,  34.9,  34.98, 35.8, 53.1,
			79.6, 80.99, 276.3, 303,   356,  383};
		G4double rnd;
		while (numGamma == 0) {
			for (G4int i = 0; i < 12; i++) {
				rnd = 100 * G4UniformRand();
				if (rnd < prob[i]) {
					G4double phi = 2.0 * 3.14159 * (G4UniformRand());
					G4double theta = acos(2 * G4UniformRand() - 1);
					G4double px = cos(phi) * sin(theta);
					G4double py = sin(phi) * sin(theta);
					G4double pz = cos(theta);

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
	else if (particleSource=="fromROOT"){
		if(histoEnergy==NULL){
			G4cout<<"Histogram not initialized..."<<G4endl;
			return;
		}

		while(energy<8){
			energy=histoEnergy->GetRandom();
		}




		
		fParticleGun->GeneratePrimaryVertex(anEvent);

	}

	else {

		fParticleSource->GeneratePrimaryVertex(anEvent);
	}
}

G4bool PrimaryGeneratorAction::InitFile() {

}
