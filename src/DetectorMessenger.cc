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
  fSetAttenStatusCmd->SetGuidance("Attenuator is in, true or false?");
  fSetAttenStatusCmd->SetParameterName("Set ATT status", false);
  fSetAttenStatusCmd->AvailableForStates(G4State_PreInit);

  fSetGridStatusCmd= new G4UIcmdWithABool("/stix/geo/grids", this);
  fSetGridStatusCmd->SetGuidance("add grids or not");
  fSetGridStatusCmd->SetParameterName("Set grid status", false);
  fSetGridStatusCmd->AvailableForStates(G4State_PreInit);


  fImportCADCmd= new G4UIcmdWithABool("/stix/geo/cad", this);
  fImportCADCmd->SetGuidance("Import cad modules, true or false?");
  fImportCADCmd->SetParameterName("Import cad modules", false);
  fImportCADCmd->AvailableForStates(G4State_PreInit);

  //fSetCADTypeCommand = new G4UIcmdWithAString("/MCP/det/setCADType", this);
}

DetectorMessenger::~DetectorMessenger() {
  delete fSetAttenStatusCmd;
  delete fImportCADCmd;
//  delete fSetCADTypeCommand;
  //delete fDetectorDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue) {
  if (command == fSetAttenStatusCmd) {
    		fDetector->SetAttenuatorStatus(fSetAttenStatusCmd->GetNewBoolValue(newValue));
  }
  else if (command == fImportCADCmd) {
    		fDetector->SetImportCADFlag(fImportCADCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSetGridStatusCmd) {
    		fDetector->SetGridsStatus(fSetGridStatusCmd->GetNewBoolValue(newValue));
  }
}
