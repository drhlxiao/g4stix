

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

// STL //
#include <string>

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
#include "DetectorMessenger.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume *Construct();

  void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green,
                    G4double blue, G4double alpha, G4bool wireFrame, G4bool solid);
  void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green,
                    G4double blue, G4double alpha);
  void RandomizeColor();
  void SetAttenuatorStatus(G4bool att);
  void ConstructCFL();
  void ConstructGrids();
  void ConstructBKG();
  void SetImportCADFlag(G4bool v){ importCADFlag=v; }

private:

  G4bool importCADFlag;

  G4bool attenuatorIn;
  G4String fWorldFile;
  G4LogicalVolume *worldLogical;
  G4Material *CdTe;
  G4Material *Tungsten, *Alum, *Iron, *Vacuum,*Air, *Alu25, *goldLayerMaterial, *Platinum, *Copper;
  G4LogicalVolume *ConstructCdTeDetector();
  G4RotationMatrix  rotMatrix;

  G4VPhysicalVolume *worldPhysical;
  G4LogicalVolume *ConstructDEMBackCover();

  DetectorMessenger *detMsg;
};

#endif
