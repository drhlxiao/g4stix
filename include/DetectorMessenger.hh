//
/// \file /include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithABool;

class DetectorMessenger : public G4UImessenger {
public:
  DetectorMessenger(DetectorConstruction *);
  ~DetectorMessenger();

  void SetNewValue(G4UIcommand *, G4String);

private:
 DetectorConstruction *fDetector;
  //G4UIdirectory *fDetectorDir;
  G4UIcmdWithABool *fSetAttenStatusCmd;
  G4UIcmdWithABool *fSetGridStatusCmd;
  G4UIcmdWithABool *fImportCADCmd;
 // G4UIcmdWithAString *fSetCADTypeCommand;
};

#endif
