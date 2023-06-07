//
/// \file g4
/// \brief Main program 
//   author: Hualin Xiao (hualin.xiao@fhnw.ch)
//   History:
//    Jun 25, 2022, first release
//
//

#define G4VIS_USE 1

#include "G4RunManager.hh"
//#include "G4MTRunManager.hh"

#include "FTFP_BERT.hh"
#include "G4UImanager.hh"
#include "G4VUserPhysicsList.hh"
#include "QGSP_BIC.hh"
#include "XrayFluoPhysicsList.hh"
// FTFP_BERT.hh"
#include "time.h"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4UIExecutive.hh"

void Help() {
	G4cout << "g4STIX GEANT4 Simulation package" << G4endl;
	G4cout << "Author: Hualin Xiao (hualin.xiao@fhnw.ch)" << G4endl;
	G4cout << "Usage:" << G4endl << "./g4STIX [OPTIONS] -m run.mac  -o OUTPUT"
		<< G4endl;
	G4cout << "Options:" << G4endl << " -i                  <input.root> "
		<< G4endl
		<< "                     Read incident particle information from "
		"input.root."
		<< G4endl
		<< "                     The structure of the root file is defined in "
		"t2sim.h ."
		<< G4endl << " -k                  kill track in grids " << G4endl
		<< G4endl << " -s                 spec.csv.   The CSV file defined the energy spectrum" << G4endl
		<< G4endl << " -qgsp               Use QGSP_EMX model" << G4endl
		<< " -h                  print help information" << G4endl;
}

int main(int argc, char **argv) {

	G4String outputFilename = "defaultOuput.root";
	G4String macFilename;
	G4String particleSourceType = "";
	G4String particleSourceFile= "";

	if (argc == 1)
		Help();

	int s = 0;
	bool killTrackInGrids = false;
	bool useQGSP=false;
	G4String sel;
	G4String commandLine;

	for(int i=0;i<argc;i++){
		commandLine+=" ";
		commandLine+=argv[i];
	}
	while (s < argc - 1) {


		sel = argv[++s];
		if (sel == "-h" || sel == "--help" || sel == "--h") {
			Help();
			return 0;

		} else if (sel == "-o") {
			outputFilename = argv[++s];
			if (!outputFilename.contains(".root")) {
				Help();
				return 0;
			}
		} 
		else if (sel == "-m") {
			macFilename = argv[++s];
			if (!macFilename.contains(".mac")) {
				Help();
				return 0;
			}
		} 
		else if (sel == "-s") {
			particleSourceFile= argv[++s];
			if (!particleSourceFile.contains(".root")) {
				Help();
				return 0;
			}
		} 


		else if (sel == "-k") {
			killTrackInGrids = true;
		} else if (sel == "-qgsp") {
			useQGSP= true;
		} else if (sel == "-Ba133") {
			particleSourceType = "Ba133";
			/*if(!particleSourceType.contains(".root"))
			  {
			  Help();
			  return 0;
			  }
			  */
		}
	}

	time_t DateTime = time(NULL);
	CLHEP::HepRandom::setTheSeed(DateTime);

	AnalysisManager *analysisManager = AnalysisManager::GetInstance();
	// CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	analysisManager->SetOutputFileName(outputFilename);

	if (killTrackInGrids)
		analysisManager->KillTracksInGrids();
	analysisManager->SetCommandLine(commandLine);

	//#ifdef G4MULTITHREADED
	/*
	   G4MTRunManager * runManager = new G4MTRunManager();
	   G4int nThreads = G4Threading::G4GetNumberOfCores()-2;
	   if(nThreads<1)nThreads=1;
	   runManager->SetNumberOfThreads(nThreads);
	   */
	//#else
	G4RunManager *runManager = new G4RunManager();
	//#endif

	DetectorConstruction *detConstruction = new DetectorConstruction();

	G4cout << "Initializing detector" << G4endl;
	runManager->SetUserInitialization(detConstruction);
	G4cout << "Initializing physicslist" << G4endl;
	//    G4VUserPhysicsList* physics = new QGSP_BIC;
	// choose of physics model:
	// http://ivana.home.cern.ch/ivana/ED-Geant4-2014/presentations/X-1-physics_lists.pdf
	// G4VUserPhysicsList* physics = new FTFP_BERT;
	//   runManager->SetUserInitialization(physics);

	if(useQGSP)
	{
		G4String plname = "QGSP_BIC_EMY";  // set whatever you like ...
		G4PhysListFactory factory;
		G4VModularPhysicsList* physlist = factory.GetReferencePhysList(plname);
		runManager->SetUserInitialization(physlist);
	}
	else{
		runManager->SetUserInitialization(new XrayFluoPhysicsList());
	}

	G4cout << "Initializing primary generation" << G4endl;
	PrimaryGeneratorAction *primarygen = new PrimaryGeneratorAction();
	runManager->SetUserAction(primarygen);
	if (particleSourceType != "") {
		primarygen->SetParticleSource(particleSourceType);
	}
	if (particleSourceFile!= "") {
		primarygen->InitParticleSpectrumFromROOT(particleSourceFile);
		//primarygen->SetParticleSource(particleSourceType);
	}
	G4cout << "Set particle type:" << particleSourceType << G4endl;

	G4cout << "Initializing run action" << G4endl;
	runManager->SetUserAction(new RunAction());
	EventAction *evtAction = new EventAction();
	G4cout << "Initializing event action" << G4endl;
	runManager->SetUserAction(evtAction);
	G4cout << "Initializing stepping action" << G4endl;
	runManager->SetUserAction(new SteppingAction());

	//#ifdef G4VIS_USE
	G4VisManager *visManager = new G4VisExecutive;
	visManager->Initialize();
	//#endif
	if (macFilename == "") {
		//#ifdef G4UI_USE
		G4UIExecutive *ui = new G4UIExecutive(argc, argv, "qt");
		ui->SessionStart();
		delete ui;
		//#endif
	} else {
		G4UIExecutive *ui = new G4UIExecutive(argc, argv, "tcsh");

		G4UImanager *UImanager = G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		G4String fileName = argv[2];
		G4cout << "Applying command: " << command + fileName << G4endl;
		analysisManager->SetMacroFileName(fileName);
		UImanager->ApplyCommand(command + fileName);
		delete ui;
	}

#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;
	G4cout << "Output filename: " << outputFilename << G4endl;

	return 0;
}
