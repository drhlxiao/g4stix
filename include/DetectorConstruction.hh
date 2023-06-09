

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

// STL //
#include <string>

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4AssemblyVolume;
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
                    G4double blue, G4double alpha, G4bool wireFrame,
                    G4bool solid);
  void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green,
                    G4double blue, G4double alpha);
  void SetVisColors();
  void SetAttenuatorStatus(G4bool att);
  void SetGridsStatus(G4bool v) { gridsEnabled = v; };
  void ConstructCFL();
  void ConstructGrids();
  void ConstructSpacecraft();

  void ConstructBKG();
  void SetGdmlFile(G4String v) { fWorldFile = v; }
  void SetActivatedDetectorFlag(G4int v) { activatedDetectorFlag = v; }

 private:
  G4bool gridsEnabled;
  G4int activatedDetectorFlag;

  G4bool attenuatorIn;
  G4String fWorldFile;
  G4LogicalVolume *worldLogical;
  G4LogicalVolume *padsLogical, *cdTeLogical, *calisteBaseLogical;
  G4Material *CdTe;
  G4Material *Epoxy, *FR4, *Resin, *SilverEpoxy;
  G4Material *Tungsten, *Alum, *Alum7075, *Iron, *Vacuum, *Air, *Alu25, *Gold,*Kapton,
      *Nickle, *Siliver, *LeadPadMat, *goldLayerMaterial, *Platinum, *Copper,
      *SiO2, *padStackMaterial;

  G4LogicalVolume *ConstructCaliste();
  G4LogicalVolume *ConstructCalisteBase();
  G4LogicalVolume *ConstructCdTe();
  G4AssemblyVolume *ConstructPads();
  bool checkOverlaps;

  // caliste

  G4RotationMatrix rotMatrix;

  G4VPhysicalVolume *worldPhysical;
  DetectorMessenger *detMsg;
};

#endif
