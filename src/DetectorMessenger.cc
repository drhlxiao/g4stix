/***************************************************************
 * Class to define detector messeger.
 * Author  : Hualin Xiao
 * Date    : Jan, 2015
 * Version : 1.10
 *
 ***************************************************************/

#include "globals.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction *theDet)
    : fDetector(theDet) {

  fSetAttenStatusCmd= new G4UIcmdWithABool("/stix/geo/att", this);
  fSetAttenStatusCmd->SetGuidance("Attenuator is in, true or false?, it will load WorldAttOut.gdml if att is out else WorldAttIn.gdml ");
  fSetAttenStatusCmd->SetParameterName("Set ATT status", false);
  fSetAttenStatusCmd->AvailableForStates(G4State_PreInit);

  fSetGridStatusCmd= new G4UIcmdWithABool("/stix/geo/grids", this);
  fSetGridStatusCmd->SetGuidance("add grids or not");
  fSetGridStatusCmd->SetParameterName("Set grid status", false);
  fSetGridStatusCmd->AvailableForStates(G4State_PreInit);


  fSetGdmlCmd= new G4UIcmdWithAString("/stix/geo/gdml", this);
  fSetGdmlCmd->SetGuidance("Import mass model from GDML file?");
  fSetGdmlCmd->SetParameterName("Import gdml file", false);
  fSetGdmlCmd->AvailableForStates(G4State_PreInit);

  fDetectorSelectionCmd= new G4UIcmdWithAnInteger("/stix/geo/det", this);
  fDetectorSelectionCmd->SetGuidance("Only construct a single detector if 0<= det < 32 else all detectors");
  fDetectorSelectionCmd->SetParameterName("detector number", false);
  fDetectorSelectionCmd->AvailableForStates(G4State_PreInit);
}

DetectorMessenger::~DetectorMessenger() {
  delete fSetAttenStatusCmd;
  delete fSetGdmlCmd;
  delete fDetectorSelectionCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fSetAttenStatusCmd) {
    		fDetector->SetAttenuatorStatus(fSetAttenStatusCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSetGdmlCmd) {
    		fDetector->SetGdmlFile(newValue);
  }
  else if (command == fSetGridStatusCmd) {
    		fDetector->SetGridsStatus(fSetGridStatusCmd->GetNewBoolValue(newValue));
  }
  else if (command == fDetectorSelectionCmd) {
    		fDetector->SetActivatedDetectorFlag(fDetectorSelectionCmd->GetNewIntValue(newValue));
  }
  
}
