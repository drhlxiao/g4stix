/// description of stix instrument
// Author: Hualin Xiao(hualin.xiao@fhnw.ch)
// Date: Fri Jun 10, 2022
#include "DetectorConstruction.hh"

#include <TFile.h>
#include <TTree.h>

#include <G4AssemblyVolume.hh>
#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Cons.hh>
#include <G4Element.hh>
#include <G4EllipticalTube.hh>
#include <G4GDMLParser.hh>
#include <G4GenericTrap.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4Polycone.hh>
#include <G4Polyhedra.hh>
#include <G4RotationMatrix.hh>
#include <G4RunManager.hh>
#include <G4SolidStore.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4Tubs.hh>
#include <G4TwoVector.hh>
#include <G4UImanager.hh>
#include <G4UnionSolid.hh>
#include <G4VisAttributes.hh>
#include <vector>

#include "AnalysisManager.hh"
#include "G4ExtrudedSolid.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "GridParameters.hh"
#include "globals.hh"
const G4double pi = CLHEP::pi;
const G4ThreeVector singleDetectorPosition(0, 0, -0.8 * mm);
const G4double bigW = 2.15 * mm;
const G4double bigH = 4.55 * mm;
const G4double bigSW = 1.1 * mm;
const G4double bigSH = 0.455 * mm;
const G4double smallW = 1.05 * mm;
const G4double smallH = 0.86 * mm;
const G4double pixel0CenterX = -3.3 * mm;
const G4double pixel0CenterY = 2.3 * mm;
const G4double pixel4CenterX = -3.3 * mm;
const G4double pixel4CenterY = -2.3 * mm;
const G4double deltaW = 2.2 * mm;
const G4double pixel8CenterX = -3.85 * mm;
const G4double pixel8CenterY = 0 * mm;

const G4double cdteThickness = 1 * mm;
const G4double anodeThickness = (15 + 15 + 100) * 1e-9 * m; // 130 nm
const G4double cathodeThickness = 15 * 1e-9 * m;            // 15 nm

const G4double calisteWidth = 11 * mm;
const G4double calisteLength = 12 * mm;
const G4double calisteBaseThickness = 14.4 * mm;
const G4double CdTeCalisteOffsetY =
    0.8 * mm; //  (top margion) 0.2 + 10 + 1.8 (Bottom margin)  = 12

const G4double bondingLandZoneOnModuleLength = 0.4 * mm;
const G4double bondingLandZoneOnModuleWidth = 2 * mm;
const G4double bondingLandZoneOnModuleThickness = 0.25 * mm;
const G4double pinsThickness = 1.2 * mm; //
G4double padCdTeBondingWidth = 0.4 * mm;
G4double padCdTeBondingLength = 1 * mm;
G4double padCdTeBondingThickness = 0.250 * mm;
G4double silverThickness = 0.2 * mm;
G4double platingThickness = (2 + 20 + 2.5 + 2) * 1e-3 * mm;
G4double padTotalThickness =
    silverThickness + platingThickness; // pads underneath the  CdTe

G4RotationMatrix *noRotation = new G4RotationMatrix(0., 0., 0.);
G4double CdTeTotalThickness = anodeThickness + cdteThickness + cathodeThickness;

const G4double calisteTotalHight =
    calisteBaseThickness + cdteThickness + anodeThickness + cathodeThickness +
    padTotalThickness + pinsThickness + padCdTeBondingThickness;
// 16.6
const G4double BKG_BIGEST_APER_X[80] = {0.5,
                                        0.498378654067105,
                                        0.4935251313189564,
                                        0.485470908713026,
                                        0.47426822097357274,
                                        0.4599897218294121,
                                        0.44272801282660496,
                                        0.42259504277189736,
                                        0.3997213817017506,
                                        0.3742553740855506,
                                        0.34636217675479974,
                                        0.31622268779768864,
                                        0.28403237336557796,
                                        0.25000000000000006,
                                        0.21434628070152711,
                                        0.1773024435212679,
                                        0.13910873195822637,
                                        0.10001284688802223,
                                        0.060268340127661614,
                                        0.02013297005470762,
                                        -0.02013297005470745,
                                        -0.06026834012766144,
                                        -0.10001284688802217,
                                        -0.1391087319582262,
                                        -0.17730244352126773,
                                        -0.21434628070152695,
                                        -0.2499999999999999,
                                        -0.28403237336557785,
                                        -0.3162226877976886,
                                        -0.3463621767547997,
                                        -0.37425537408555043,
                                        -0.3997213817017505,
                                        -0.4225950427718973,
                                        -0.4427280128266049,
                                        -0.4599897218294121,
                                        -0.4742682209735727,
                                        -0.48547090871302595,
                                        -0.4935251313189564,
                                        -0.49837865406710496,
                                        -0.5,
                                        -0.5,
                                        -0.498378654067105,
                                        -0.4935251313189564,
                                        -0.48547090871302606,
                                        -0.4742682209735728,
                                        -0.4599897218294121,
                                        -0.44272801282660496,
                                        -0.4225950427718974,
                                        -0.3997213817017507,
                                        -0.37425537408555065,
                                        -0.3463621767547998,
                                        -0.31622268779768886,
                                        -0.28403237336557796,
                                        -0.2500000000000002,
                                        -0.21434628070152706,
                                        -0.17730244352126795,
                                        -0.13910873195822665,
                                        -0.10001284688802228,
                                        -0.06026834012766179,
                                        -0.020132970054707572,
                                        0.020132970054707388,
                                        0.06026834012766116,
                                        0.1000128468880221,
                                        0.13910873195822604,
                                        0.17730244352126778,
                                        0.2143462807015269,
                                        0.24999999999999967,
                                        0.2840323733655778,
                                        0.3162226877976887,
                                        0.34636217675479963,
                                        0.3742553740855504,
                                        0.3997213817017503,
                                        0.4225950427718973,
                                        0.442728012826605,
                                        0.45998972182941206,
                                        0.4742682209735727,
                                        0.4854709087130259,
                                        0.4935251313189564,
                                        0.498378654067105,
                                        0.5};
const G4double BKG_BIGEST_APER_Y[80] = {
    0.8036504591506379,  0.8438837435090009,  0.883856099579518,
    0.9233082912944168,  0.9619844560513742,  0.9996337640806754,
    1.036012045172522,   1.0708833722145386,  1.1040215912696274,
    1.1352117882710355,  1.164251682822545,   1.1909529400644652,
    1.2151423920974662,  1.2366631610428573,  1.255375676455829,
    1.2711585804933452,  1.283909514966324,   1.2935457851717718,
    1.3000048961996649,  1.3032449582364227,  1.3032449582364227,
    1.3000048961996649,  1.2935457851717718,  1.2839095149663242,
    1.2711585804933454,  1.255375676455829,   1.2366631610428573,
    1.2151423920974662,  1.1909529400644652,  1.1642516828225453,
    1.1352117882710357,  1.1040215912696274,  1.0708833722145386,
    1.0360120451725223,  0.9996337640806755,  0.9619844560513744,
    0.9233082912944169,  0.8838560995795182,  0.843883743509001,
    0.803650459150638,   -0.8036504591506378, -0.8438837435090009,
    -0.883856099579518,  -0.9233082912944166, -0.9619844560513741,
    -0.9996337640806754, -1.036012045172522,  -1.0708833722145386,
    -1.1040215912696272, -1.1352117882710355, -1.164251682822545,
    -1.190952940064465,  -1.2151423920974662, -1.2366631610428571,
    -1.255375676455829,  -1.2711585804933452, -1.283909514966324,
    -1.2935457851717718, -1.3000048961996649, -1.3032449582364227,
    -1.3032449582364227, -1.300004896199665,  -1.2935457851717718,
    -1.2839095149663242, -1.2711585804933454, -1.255375676455829,
    -1.2366631610428573, -1.2151423920974662, -1.1909529400644652,
    -1.1642516828225453, -1.1352117882710357, -1.1040215912696278,
    -1.0708833722145386, -1.036012045172522,  -0.9996337640806756,
    -0.9619844560513744, -0.9233082912944173, -0.8838560995795183,
    -0.8438837435090009, -0.803650459150638};
// bigest aperture in BKG detector, width=0.5 mm, height ~ 1.6

DetectorConstruction::DetectorConstruction() {
  // for (int i = 0; i < 32; i++)G4cout <<Grid::getDetectorCenterCoordsCAD(i) <<
  // G4endl;
  fWorldFile = "";
  usingCADModels = false;
  gridsEnabled = true;
  detectorEnabled = true;
  activatedDetectorFlag = 100;
  cflConstructed = false;
  bkgConstructed = false;

  isSingleDetector = false;

  // all detectors will be constructed  if it is not between 0 -- 31

  G4RotationMatrix rotY;
  rotY.rotateY(-90 * deg);
  G4RotationMatrix rotX;
  rotX.rotateX(90 * deg);
  rotMatrix = rotX * rotY;
  checkOverlaps = true;

  detMsg = new DetectorMessenger(this);
}
void DetectorConstruction::ConstructSpacecraft() {
  // definition moved to gdml
  /*
     G4double scWidth= 2.5 *m;
     G4double scThickness= 2.7 *m;
     G4double scLength= 3.1*m;
  //taken from wikipedia
  G4double scWallThickness=1*mm;

  G4Box *spacecraftBox= new G4Box("spacecaft",
  scLength/2,
  scThickness/2,
  scWidth/2
  );//in stix coordinate frame
    //stix_x
    //
    G4Box *spacecraftBoxInner= new G4Box("spacecaft",
    scLength/2,
    scThickness/2-1,
    scWidth/2-1
    );//in stix coordinate frame
          //

          G4SubtractionSolid *spacecraftHollowBox= new
          G4SubtractionSolid("SpaceLab", SL_123, SL3hole_cons); G4LogicalVolume
   *spacecraftLog=new G4LogicalVolume(spacecraftBox, Alum, "spacecraftLog", 0,
   0, 0); G4Box *spacecraftBox= new G4Box("spacecaft", scLength/2,
   scThickness/2, scWidth/2
   );//in stix coordinate frame

*/
}

void DetectorConstruction::ConstructGrids() {
  // don't construct if disabled by macro

  // created by grid_data_creator.py in
  // /home/xiaohl/FHNW/STIX/SolarFlareAnalysis/stix_simulator/pySimulator/export_grid_data.inpy
  G4String fname = "./grid_data/stix_grid_parameters.root";
  G4cout << "Loading grid parameters from ROOT file: " << fname << G4endl;
  TFile f(fname);
  TTree *grids = (TTree *)f.Get("grids");
  Int_t is_front;
  Int_t det_idx;
  Int_t i_strip;
  Float_t polygon_x[5];
  Float_t polygon_y[5];
  Int_t is_nominal_parameters;

  // Set branch addresses.
  grids->SetBranchAddress("is_front", &is_front);
  grids->SetBranchAddress("det_idx", &det_idx);
  grids->SetBranchAddress("i_strip", &i_strip);
  grids->SetBranchAddress("polygon_x", polygon_x);
  grids->SetBranchAddress("polygon_y", polygon_y);
  grids->SetBranchAddress("is_nominal_parameters", &is_nominal_parameters);

  Long64_t nentries = grids->GetEntries();

  G4ThreeVector pos;

  G4bool nominalOrRealFlag[32] = {0};

  G4Box *frontGridBox =
      new G4Box("frontGridBox", 11.05 * mm, 10.05 * mm, 0.205 * mm);
  G4Box *rearGridBox =
      new G4Box("rearGridBox", 6.505 * mm, 6.505 * mm, 0.205 * mm);

  G4double frameThickness = 0.5 * mm / 2;
  G4Box *frontGridBoxOuter =
      new G4Box("frontGridBox", 11.05 * mm + frameThickness,
                10.05 * mm + frameThickness, 0.205 * mm);
  G4Box *rearGridBoxOuter =
      new G4Box("rearGridBox", 6.505 * mm + frameThickness,
                6.505 * mm + frameThickness, 0.205 * mm);

  // dummy boxes
  // added extra 50 um to avoid overlapping

  G4LogicalVolume *frontGridContainerLog[32];
  G4LogicalVolume *rearGridContainerLog[32];

  G4LogicalVolume *frontGridContainerOuterLog[32];
  G4LogicalVolume *rearGridContainerOuterLog[32];
  for (int i = 0; i < 32; i++) {

    frontGridContainerOuterLog[i] = new G4LogicalVolume(
        frontGridBoxOuter, Tungsten, "frontGridContainerOuter", 0, 0, 0);
    rearGridContainerOuterLog[i] = new G4LogicalVolume(
        rearGridBoxOuter, Tungsten, "rearGridContainerOuter", 0, 0, 0);

    frontGridContainerLog[i] = new G4LogicalVolume(
        frontGridBox, Vacuum, "frontGridContainer", 0, 0, 0);
    rearGridContainerLog[i] =
        new G4LogicalVolume(rearGridBox, Vacuum, "rearGridContainer", 0, 0, 0);
  }

  for (int i = 0; i < nentries; i++) {
    grids->GetEntry(i);
    nominalOrRealFlag[det_idx] = is_nominal_parameters;
    std::vector<G4TwoVector> vertexCoords;
    for (int j = 0; j < 4; j++) {
      vertexCoords.push_back(G4TwoVector(polygon_x[j] * mm, polygon_y[j] * mm));
    }
    // Left wedge Solid and logical Volume
    G4ExtrudedSolid *strip =
        new G4ExtrudedSolid("strip", vertexCoords,
                            0.4 * 0.5 * mm, // 0.4 mm thick, halfz
                            G4TwoVector(), 1, G4TwoVector(), 1);
    G4LogicalVolume *gridStripLog =
        new G4LogicalVolume(strip, Tungsten, "gridStripLog", 0, 0, 0);
    G4cout << "creating grids:" << i << " " << pos << G4endl;
    if (is_front == 1) {
      new G4PVPlacement(0, G4ThreeVector(), gridStripLog, "gridStrip",
                        frontGridContainerLog[det_idx], false, i, false);
    } else {
      new G4PVPlacement(0, G4ThreeVector(), gridStripLog, "gridStrip",
                        rearGridContainerLog[det_idx], false, i, false);
    }
  }

  for (int i = 0; i < 32; i++) {
    new G4PVPlacement(0, G4ThreeVector(), frontGridContainerLog[i],
                      "frontGridContainerOuter", frontGridContainerOuterLog[i],
                      false, i, false);
    new G4PVPlacement(0, G4ThreeVector(), rearGridContainerLog[i],
                      "rearGridContainerOuter", rearGridContainerOuterLog[i],
                      false, i, false);
  }

  for (int i = 0; i < 32; i++) {
    if (i == 8 || i == 9)
      continue;
    // don't construct grids for CFL and BKG
    G4String gtype = nominalOrRealFlag[i] == 1 ? "NOMINAL" : "REAL";
    G4cout << "Grid #" << i << " parameter type: " << gtype << G4endl;

    if (isSingleDetector) {
      // if it is a single detector model, only reconstruct grids for this
      // detector
      if (i != activatedDetectorFlag)
        continue;
    }

    pos = isSingleDetector ? G4ThreeVector(-21 * mm, 0, 0)
                           : Grid::getGridCenterCAD(i, 1);
    new G4PVPlacement(G4Transform3D(rotMatrix, pos),
                      frontGridContainerOuterLog[i], "frontGridWindow",
                      worldLogical, false, i, false);

    pos = isSingleDetector ? G4ThreeVector(-19 * mm, 0, 0)
                           : Grid::getGridCenterCAD(i, 0);
    new G4PVPlacement(G4Transform3D(rotMatrix, pos),
                      rearGridContainerOuterLog[i], "rearGridWindow",
                      worldLogical, false, i, false);
  }

  // fluorescence test

  G4cout << "Grid construction finished." << G4endl;
  f.Close();
}
void DetectorConstruction::ConstructCalibrationFoil() {

  //
  //
  // add calibration foils, the one converted from cad has issues, absorbs all
  // photons below 6 keV, there might be overlapping
  //  2023, June 9th, Hualin
  //
  G4cout << "Constructing Calibration Foil..." << G4endl;
  std::vector<G4TwoVector> calibrationFoilVertexCoords;
  G4double extraMargin = -2 * mm;
  G4ThreeVector calFoilAluCenter(9.251 * mm + extraMargin, 104.074 * mm,
                                 126.709 * mm);
  G4ThreeVector calFoilKaptonCenter(9.251 * mm + extraMargin - 1 * mm,
                                    104.074 * mm, 126.709 * mm);
  // measured from cad, extraMargin added manually to avoid overlapping
  G4double calCoords[8][2] = {{-44.0, 78.5},  {44.0, 78.5},  {85.0, 39.0},
                              {85.0, -39.0},  {44.0, -78.5}, {-44.0, -78.5},
                              {-85.0, -39.0}, {-85.0, 39.0}};
  for (int i = 0; i < 8; i++) {
    calibrationFoilVertexCoords.push_back(
        G4TwoVector(calCoords[i][1] * mm, calCoords[i][0] * mm));
  }
  G4ExtrudedSolid *calFoilKapton =
      new G4ExtrudedSolid("calFoilKapton", calibrationFoilVertexCoords,
                          0.2032 * 0.5 * mm, // 0.2032 mm thick, halfz
                          G4TwoVector(), 1, G4TwoVector(), 1);
  G4ExtrudedSolid *calFoilAlum =
      new G4ExtrudedSolid("calFoilAlum", calibrationFoilVertexCoords,
                          4000 * 1e-10 * 0.5 * m, // 4000 angstom thick, halfz
                          G4TwoVector(), 1, G4TwoVector(), 1);

  G4LogicalVolume *calFoilKaptonLog =
      new G4LogicalVolume(calFoilKapton, Kapton, "calKaptonLog", 0, 0, 0);
  G4LogicalVolume *calFoilAluLog =
      new G4LogicalVolume(calFoilAlum, Alum, "calAlumLog", 0, 0, 0);

  new G4PVPlacement(G4Transform3D(rotMatrix, calFoilKaptonCenter),
                    calFoilKaptonLog, "calFoilKaptonPhy", worldLogical, false,
                    0, false);

  new G4PVPlacement(G4Transform3D(rotMatrix, calFoilAluCenter), calFoilAluLog,
                    "calFoilAluLogPhy", worldLogical, false, 0, false);
}

G4LogicalVolume *DetectorConstruction::ConstructCdTe() {
  // cdte detector with pixels
  // the orgin is the top of surface

  G4Box *CdTeBox =
      new G4Box("CdTeBox", 5 * mm, 5 * mm, CdTeTotalThickness / 2.);
  G4LogicalVolume *CdTeLog =
      new G4LogicalVolume(CdTeBox, Vacuum, "CdTeLog", 0, 0, 0);
  // it doesn't matter what materials it is

  // thickness negligible  compared to the thickness  uncertainty
  checkOverlaps = true;
  ////
  ///// calisteResin base with different coating layers
  // caliste base material, is unknown
  // can be considered filled resin

  /////construct CdTe and pixels

  G4Box *CdTeSDBox = new G4Box("CdTeDetSD", 5 * mm, 5 * mm, cdteThickness / 2.);
  G4Box *CdTeAnodeBox =
      new G4Box("CdTeAnode", 5 * mm, 5 * mm, anodeThickness / 2.);
  G4Box *CdTeCathodeBox =
      new G4Box("CdTeCathode", 5 * mm, 5 * mm, cathodeThickness / 2.);

  // Top and Small pixels
  std::vector<G4TwoVector> vertexCoordsTop;
  vertexCoordsTop.push_back(G4TwoVector(-bigW / 2, bigH / 2));
  vertexCoordsTop.push_back(G4TwoVector(bigW / 2, bigH / 2));
  vertexCoordsTop.push_back(G4TwoVector(bigW / 2, -bigH / 2));
  vertexCoordsTop.push_back(G4TwoVector(-bigW / 2 + bigSW, -bigH / 2));
  vertexCoordsTop.push_back(G4TwoVector(-bigW / 2 + bigSW, -bigH / 2 + bigSH));
  vertexCoordsTop.push_back(G4TwoVector(-bigW / 2, -bigH / 2 + bigSH));

  G4ExtrudedSolid *bigPixelTopGeo =
      new G4ExtrudedSolid("topBigPixel", vertexCoordsTop, cdteThickness / 2.,
                          G4TwoVector(), 1, G4TwoVector(), 1);
  std::vector<G4TwoVector> vertexCoordsBott;
  vertexCoordsBott.push_back(G4TwoVector(-bigW / 2, bigH / 2 - bigSH));
  vertexCoordsBott.push_back(G4TwoVector(-bigW / 2 + bigSW, bigH / 2 - bigSH));
  vertexCoordsBott.push_back(G4TwoVector(-bigW / 2 + bigSW, bigH / 2));
  vertexCoordsBott.push_back(G4TwoVector(bigW / 2, bigH / 2));
  vertexCoordsBott.push_back(G4TwoVector(bigW / 2, -bigH / 2));
  vertexCoordsBott.push_back(G4TwoVector(-bigW / 2, -bigH / 2));

  G4ExtrudedSolid *bigPixelBottomGeo =
      new G4ExtrudedSolid("bottBigPixel", vertexCoordsBott, cdteThickness / 2.,
                          G4TwoVector(), 1, G4TwoVector(), 1);
  G4Box *smallPixelGeo =
      new G4Box("smallPixel", smallW / 2, smallH / 2, cdteThickness / 2.);

  G4LogicalVolume *CdTeDetSDLog =
      new G4LogicalVolume(CdTeSDBox, Vacuum, "CdTeDetSDLog", 0, 0, 0);

  G4LogicalVolume *CdTeAnodeLog =
      new G4LogicalVolume(CdTeAnodeBox, Vacuum, "CdTeAnodeLog", 0, 0, 0);
  G4LogicalVolume *CdTeCathodeLog =
      new G4LogicalVolume(CdTeCathodeBox, Vacuum, "CdTeCathodeLog", 0, 0, 0);
  G4LogicalVolume *bigPixelTopLog =
      new G4LogicalVolume(bigPixelTopGeo, CdTe, "bigPixelTopLog", 0, 0, 0);
  G4LogicalVolume *bigPixelBottomLog = new G4LogicalVolume(
      bigPixelBottomGeo, blackHole, "bigPixelBottomLog", 0, 0, 0);
  G4LogicalVolume *smallPixelLog =
      new G4LogicalVolume(smallPixelGeo, blackHole, "smallPixelLog", 0, 0, 0);
  G4int copyNb;
  G4double detZ = 0;
  // place electrode
  G4String name = "pixel";
  for (int i = 0; i < 4; i++) {
    // it becomes pixel4 after rotation
    copyNb = i;
    G4ThreeVector posBigPixelTop(pixel0CenterX + deltaW * i, pixel0CenterY,
                                 detZ);
    new G4PVPlacement(0, posBigPixelTop, bigPixelTopLog, name, CdTeDetSDLog,
                      false, copyNb, checkOverlaps);

    // it becomes pixel4 after rotation
    copyNb = i + 4;
    G4ThreeVector posBigPixelBottom(pixel4CenterX + deltaW * i, pixel4CenterY,
                                    detZ);
    new G4PVPlacement(0, posBigPixelBottom, bigPixelBottomLog, name,
                      CdTeDetSDLog, false, copyNb, checkOverlaps);

    // small pixels
    copyNb = i + 8;
    G4ThreeVector posSmallPixel(pixel8CenterX + deltaW * i, pixel8CenterY,
                                detZ);

    new G4PVPlacement(0, posSmallPixel, smallPixelLog, name, CdTeDetSDLog,
                      false, copyNb, checkOverlaps);
  }

  G4ThreeVector TmCdTeCathod(0, 0,
                             CdTeTotalThickness / 2. - cathodeThickness / 2);
  G4ThreeVector TmDet(
      0, 0, CdTeTotalThickness / 2. - cathodeThickness - cdteThickness / 2);
  G4ThreeVector TmAnode(0, 0,
                        CdTeTotalThickness / 2. - cathodeThickness -
                            cdteThickness - anodeThickness / 2);

  new G4PVPlacement(0, TmDet, CdTeDetSDLog, "CdTeDetSD", CdTeLog, false, 0,
                    checkOverlaps);
  // new G4PVPlacement(0, TmCdTeCathod, CdTeCathodeLog, "CdTeCathodPhys",
  // CdTeLog, 		false, 0, checkOverlaps); new G4PVPlacement(0, TmAnode,
  // CdTeAnodeLog, "CdTeAnodePhys", CdTeLog, false, 		0, checkOverlaps);

  return CdTeLog;
  // done
}
G4LogicalVolume *DetectorConstruction::ConstructCalisteBase() {
  G4Box *calisteBaseOuter =
      new G4Box("CdTeModuleBaseOuter", calisteWidth / 2, calisteLength / 2,
                calisteBaseThickness / 2);

  G4LogicalVolume *calisteBaseOuterLog =
      new G4LogicalVolume(calisteBaseOuter, Gold, "calisteBaseOuter", 0, 0, 0);

  // Au+Ni+Cu+2u
  ////////Plating from inside to outside:
  //+ All surfaces
  //+ Ni 2u + Cu 20u + Ni 2.5u + Au 2 u
  //
  G4double coatingThickness[4] = {2e-3 * mm, 2.5e-3 * mm, 20e-3 * mm,
                                  2e-3 * mm};
  // we need to reverse the order
  G4Material *coatingMateris[4] = {Nickle, Copper, Nickle, Resin};
  // outest layer is gold
  //  different layer of the coating material is  considered as mixture

  G4Box *calisteBaseInner[4];
  G4String calisteBaseInnerLayerName[4] = {"nickle", "copper", "nickle",
                                           "resin"};

  G4LogicalVolume *calisteBaseInnerLog[4];
  G4double totalCoatingThickness = 0;

  G4LogicalVolume *coatingMotherLog;

  for (int i = 0; i < 4; i++) {
    totalCoatingThickness += coatingThickness[i];
    calisteBaseInner[i] = new G4Box(
        G4String("CdTeModuleBaseInner_") + calisteBaseInnerLayerName[i],
        calisteWidth / 2 - totalCoatingThickness,
        calisteLength / 2 - totalCoatingThickness, calisteBaseThickness / 2);
    calisteBaseInnerLog[i] = new G4LogicalVolume(
        calisteBaseInner[i], coatingMateris[i],
        G4String("log_CalisteBaseInner_") + calisteBaseInnerLayerName[i], 0, 0,
        0);
    if (i == 0) {
      coatingMotherLog = calisteBaseOuterLog;
    } else {
      coatingMotherLog = calisteBaseInnerLog[i - 1];
    }

    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), calisteBaseInnerLog[i],
                      calisteBaseInnerLayerName[i] + "_phys", coatingMotherLog,
                      false, 0, checkOverlaps);
  }
  return calisteBaseOuterLog;
}
G4AssemblyVolume *DetectorConstruction::ConstructPads() {
  //////construct pads
  //////
  G4AssemblyVolume *padAssembly = new G4AssemblyVolume();
  G4double padZ = 0 * mm;

  G4ThreeVector TmStackInPadMother(
      0, 0, -padTotalThickness / 2 + platingThickness / 2);

  // round pad
  G4Tubs *roundPadGeo = new G4Tubs("roundPadGeo", 0, 0.3 * mm,
                                   padTotalThickness / 2., 0, 360. * deg);
  G4Tubs *roundPadStackGeo = new G4Tubs("roundPadStackGeo", 0, 0.3 * mm,
                                        platingThickness / 2., 0, 360. * deg);

  G4LogicalVolume *roundPadLog =
      new G4LogicalVolume(roundPadGeo, SilverEpoxy, "roundPadLog", 0, 0, 0);
  G4LogicalVolume *roundPadStackLog = new G4LogicalVolume(
      roundPadStackGeo, padStackMaterial, "roundPadStackLog", 0, 0, 0);
  new G4PVPlacement(0, TmStackInPadMother, roundPadStackLog,
                    "roundPadStackPhys", roundPadLog, false, 0, checkOverlaps);

  // elliptical pads
  G4EllipticalTube *ellipticalPadGeo = new G4EllipticalTube(
      "ellipticalGeo", 0.3, 0.4 * mm, padTotalThickness / 2);
  G4EllipticalTube *ellipticalPadStackGeo = new G4EllipticalTube(
      "ellipticalStackGeo", 0.3, 0.4 * mm, platingThickness / 2);
  G4LogicalVolume *ellipticalPadLog = new G4LogicalVolume(
      ellipticalPadGeo, SilverEpoxy, "ellipticalPadLog", 0, 0, 0);
  G4LogicalVolume *ellipticalPadStackLog =
      new G4LogicalVolume(ellipticalPadStackGeo, padStackMaterial,
                          "ellipticalPadStackLog", 0, 0, 0);
  new G4PVPlacement(0, TmStackInPadMother, ellipticalPadStackLog,
                    "ellipticalPadStackLog", ellipticalPadLog, false, 0,
                    checkOverlaps);

  G4Box *rectPadGeo =
      new G4Box("rectPadGeo", 0.3 * mm, 4.45 / 2 * mm, padTotalThickness / 2.);
  G4LogicalVolume *rectPadLog =
      new G4LogicalVolume(rectPadGeo, SilverEpoxy, "rectPadLog", 0, 0, 0);

  G4Box *rectPadStackGeo = new G4Box("rectPadStckGeo", 0.3 * mm, 4.45 / 2 * mm,
                                     platingThickness / 2.);
  G4LogicalVolume *rectPadStackLog = new G4LogicalVolume(
      rectPadStackGeo, padStackMaterial, "rectPadStackLog", 0, 0, 0);
  new G4PVPlacement(0, TmStackInPadMother, rectPadStackLog, "rectPadStackPhys",
                    rectPadLog, false, 0, checkOverlaps);

  G4double roundPadX[] = {
      -3.744, -2.832, -1.541, -0.646, 0.714,  1.489,  2.849,  3.744,  -3.727,
      -2.832, -1.541, -0.628, 0.680,  1.558,  2.832,  3.744,  -3.744, -2.849,
      -1.523, -0.646, 0.646,  1.541,  2.832,  3.727,  -3.865, 2.746,  -3.710,
      -2.832, -1.558, -0.628, 1.541,  2.832,  3.744,  3.727,  2.815,  1.523,
      0.611,  -0.663, -1.523, -2.832, -3.710, 3.744,  2.849,  1.523,  0.628,
      -0.663, -1.523, -2.832, -3.744, 0.559,  -1.627, 0.628};
  G4double roundPadY[] = {
      -4.678, -4.678, -4.695, -4.678, -4.678, -4.695, -4.644, -4.661, -3.768,
      -3.768, -3.785, -3.768, -3.768, -3.751, -3.768, -3.768, -2.858, -2.876,
      -2.910, -2.893, -2.893, -2.876, -2.876, -2.893, -0.747, -0.798, 1.296,
      1.313,  1.313,  1.313,  1.348,  1.348,  1.365,  2.258,  2.240,  2.240,
      2.258,  2.223,  2.240,  2.240,  2.240,  3.150,  3.150,  3.150,  3.150,
      3.133,  3.116,  3.133,  3.133,  -0.764, -0.764, 1.365};
  G4double ellipitalPadX[] = {-2.746, -2.763, -0.525, -0.542,
                              1.644,  1.644,  3.830,  3.847};
  G4double ellipitalPadY[] = {-1.433, -0.060, -1.451, -0.077,
                              -1.451, -0.060, -1.433, -0.077};
  G4double rectPadX[] = {4.71, -4.75};
  G4double rectPadY[] = {-0.755, -0.755};
  // positions extracted from the graph and verified using pyplot
  //	cdteThickness - anodeThickness - padTotalThickness/2.;

  // placing pads
  //
  for (int i = 0; i < 52; i++) {
    G4ThreeVector Tm(roundPadX[i] * mm, -roundPadY[i] * mm, padZ);
    padAssembly->AddPlacedVolume(roundPadLog, Tm, noRotation);
  }
  for (int i = 0; i < 8; i++) {
    G4ThreeVector Tm(ellipitalPadX[i] * mm, -ellipitalPadY[i] * mm, padZ);
    padAssembly->AddPlacedVolume(ellipticalPadLog, Tm, noRotation);
  }
  for (int i = 0; i < 2; i++) {
    G4ThreeVector Tm(rectPadX[i] * mm, -rectPadY[i] * mm, padZ);
    padAssembly->AddPlacedVolume(rectPadLog, Tm, noRotation);
  }
  // y reversed
  //  see email from olivier, sent on Feb. 09 2023: CAD model Caliste-SO
  //
  // Place the gold bonding rod for HV
  G4Box *bondingLandZoneOnModuleGeo = new G4Box(
      "bondingLandZoneOnModuleGeo", bondingLandZoneOnModuleWidth / 2,
      bondingLandZoneOnModuleLength / 2, bondingLandZoneOnModuleThickness / 2);

  G4LogicalVolume *bondingLandZoneOnModuleLog = new G4LogicalVolume(
      bondingLandZoneOnModuleGeo, Gold, "bondingLandZoneOnModuleLog", 0, 0, 0);
  G4Box *padCdTeBondingGeo =
      new G4Box("padCdTeBondingGeo", padCdTeBondingWidth / 2,
                padCdTeBondingLength / 2, padCdTeBondingThickness / 2);
  G4LogicalVolume *padCdTeBondingLog = new G4LogicalVolume(
      padCdTeBondingGeo, Gold, "padCdTeBondingLog", 0, 0, 0);
  G4ThreeVector TmHVPad(
      calisteWidth / 2 - bondingLandZoneOnModuleWidth / 2,
      -calisteLength / 2. + bondingLandZoneOnModuleLength / 2.,
      padZ + (bondingLandZoneOnModuleThickness - padTotalThickness) / 2.);

  padAssembly->AddPlacedVolume(bondingLandZoneOnModuleLog, TmHVPad, noRotation);
  // placed at the corner corner
  G4ThreeVector TmCdTePadInAssembly(4.8 * mm, -3.5 * mm,
                                    padZ + padTotalThickness / 2. +
                                        CdTeTotalThickness +
                                        padCdTeBondingThickness / 2);
  padAssembly->AddPlacedVolume(padCdTeBondingLog, TmCdTePadInAssembly,
                               noRotation);
  // center is at the pad center
  return padAssembly;
}

G4LogicalVolume *DetectorConstruction::ConstructCaliste() {
  G4Box *calisteWorld = new G4Box("CdTeModuleWorld", calisteWidth / 2,
                                  calisteLength / 2, calisteTotalHight / 2);
  G4LogicalVolume *calisteLog =
      new G4LogicalVolume(calisteWorld, Vacuum, "calisteWorld", 0, 0, 0);
  // CdTe detector module world
  G4LogicalVolume *CdTeLog = ConstructCdTe();
  G4LogicalVolume *calisteBaseLog = ConstructCalisteBase();
  G4ThreeVector TmHVBondingPad(
      0, 0,
      calisteTotalHight / 2 - (padCdTeBondingThickness + CdTeTotalThickness) -
          padTotalThickness / 2);

  // G4AssemblyVolume *padAssembly = ConstructPads();
  // padAssembly->MakeImprint(calisteLog, TmHVBondingPad, noRotation);

  G4ThreeVector TmCdTe(0, CdTeCalisteOffsetY,
                       calisteTotalHight / 2 - padCdTeBondingThickness -
                           CdTeTotalThickness / 2);
  // shifted by 0.8 mm
  G4ThreeVector TmCalisteBase(0, 0,
                              calisteTotalHight / 2 - padCdTeBondingThickness -
                                  CdTeTotalThickness - padTotalThickness -
                                  calisteBaseThickness / 2);

  new G4PVPlacement(0, TmCdTe, CdTeLog, "CdTeLog", calisteLog, false, 0,
                    checkOverlaps);
  // new G4PVPlacement(0, TmCalisteBase, calisteBaseLog, "calisteBasePhys",
  //		calisteLog, false, 0, checkOverlaps);
  return calisteLog;
}

DetectorConstruction::~DetectorConstruction() { delete detMsg; }

G4VPhysicalVolume *DetectorConstruction::Construct() {
  // AnalysisManager->SetAttenuatorStatus(attenuatorIn);

  G4NistManager *nist_manager = G4NistManager::Instance();
  Alum = nist_manager->FindOrBuildMaterial("G4_Al");
  Vacuum = nist_manager->FindOrBuildMaterial("G4_Galactic");
  Tungsten = nist_manager->FindOrBuildMaterial("G4_W");
  Platinum = nist_manager->FindOrBuildMaterial("G4_Pt");
  Copper = nist_manager->FindOrBuildMaterial("G4_Cu");

  Siliver = nist_manager->FindOrBuildMaterial("G4_Ag");
  Gold = nist_manager->FindOrBuildMaterial("G4_Au");
  Nickle = nist_manager->FindOrBuildMaterial("G4_Ni");

  G4Element *elAl = nist_manager->FindOrBuildElement("Al");
  G4Element *elH = nist_manager->FindOrBuildElement("H");
  G4Element *elC = nist_manager->FindOrBuildElement("C");
  G4Element *elTi = nist_manager->FindOrBuildElement("Ti");
  G4Element *elCd = nist_manager->FindOrBuildElement("Cd");
  G4Element *elTe = nist_manager->FindOrBuildElement("Te");
  G4Element *elAu = nist_manager->FindOrBuildElement("Au");
  G4Element *elAg = nist_manager->FindOrBuildElement("Ag");

  G4Element *elCr = nist_manager->FindOrBuildElement("Cr");
  G4Element *elZn = nist_manager->FindOrBuildElement("Zn");
  G4Element *elMn = nist_manager->FindOrBuildElement("Mn");
  G4Element *elCu = nist_manager->FindOrBuildElement("Cu");
  G4Element *elFe = nist_manager->FindOrBuildElement("Fe");
  G4Element *elO = nist_manager->FindOrBuildElement("O");
  G4Element *elSi = nist_manager->FindOrBuildElement("Si");
  G4Element *elMg = nist_manager->FindOrBuildElement("Mg");
  G4Element *elNi = nist_manager->FindOrBuildElement("Ni");
  G4Element *elN = nist_manager->FindOrBuildElement("N");

  G4double density = 5.85 * g / cm3; // STIX-DS-0017-PSI
  G4int nelements, natoms;
  G4double fractionmass;
  CdTe = new G4Material("CdTe", density, nelements = 2);
  CdTe->AddElement(elCd, natoms = 1);
  CdTe->AddElement(elTe, natoms = 1);

  goldLayerMaterial =
      new G4Material("goldLayerMaterial", density, nelements = 3);
  goldLayerMaterial->AddElement(elAu, 0.94700687); // mass fraction
  goldLayerMaterial->AddElement(elTi, 0.033120707);
  goldLayerMaterial->AddElement(elAl, 0.019872424);

  // Same as the body: Ni 2u + Cu 20u + Ni 2.5u + Au 2 u

  density = 9.0 * g / cm3;
  padStackMaterial = new G4Material("padStackMaterial", density, nelements = 3);
  padStackMaterial->AddElement(elCu, 0.695); // mass fraction
  padStackMaterial->AddElement(elAu, 0.15);
  padStackMaterial->AddElement(elNi, 0.155);
  // it is considered as a mixture, maybe it is OK for X rays

  SiO2 = new G4Material("SiO2", density = 2.2 * g / cm3, nelements = 2);
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);

  blackHole =
      new G4Material("blackHole", density = 4e14 * g / cm3, nelements = 1);
  blackHole->AddElement(elH, 1);
  //

  Alum7075 = new G4Material("Alum7075", density = 2.8 * g / cm3, nelements = 7);
  Alum7075->AddElement(elAl, 0.876);
  Alum7075->AddElement(elZn, 0.056);
  Alum7075->AddElement(elCr, 0.023);
  Alum7075->AddElement(elMg, 0.025);
  Alum7075->AddElement(elCu, 0.016);
  Alum7075->AddElement(elFe, 0.0025);
  Alum7075->AddElement(elMn, 0.0015);

  density = 1.43 * g / cm3;
  Kapton = new G4Material("Kapton", density, nelements = 4);
  Kapton->AddElement(elH, 0.026362);
  Kapton->AddElement(elC, 0.691133);
  Kapton->AddElement(elN, 0.07327);
  Kapton->AddElement(elO, 0.209235);

  density = 1.2 * g / cm3;

  Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, nelements = 2);
  Epoxy->AddElement(elH, natoms = 2);
  Epoxy->AddElement(elC, natoms = 2);

  density = 1.86 * g / cm3;
  FR4 = new G4Material("FR4", density, nelements = 2);
  FR4->AddMaterial(SiO2, 0.528);
  FR4->AddMaterial(Epoxy, 0.472);

  // 2/ Molding resin
  // Use the following mixture
  //+ density of the mixture 1.77
  //+ 73% resin (C2H2)
  //+ 27% SiO2

  density = 1.77 * g / cm3;
  Resin = new G4Material("Resin", density, nelements = 2);
  Resin->AddMaterial(SiO2, fractionmass = 0.27);
  Resin->AddMaterial(Epoxy, fractionmass = 0.73);

  SilverEpoxy =
      new G4Material("SilverEpoxy", density = 2.566 * g / cm3, nelements = 2);
  // not 3.3 g/cm3, should be 2.566 measured by Olivier
  SilverEpoxy->AddMaterial(Epoxy, fractionmass = 0.7);
  SilverEpoxy->AddMaterial(Siliver, fractionmass = 0.3);

  LeadPadMat =
      new G4Material("LeadPadMat", density = 11 * g / cm3, nelements = 2);
  LeadPadMat->AddMaterial(Nickle, fractionmass = 0.587);
  LeadPadMat->AddMaterial(Gold, fractionmass = 0.413);

  // construct world
  G4GDMLParser parser;
  if (fWorldFile != "") {
    G4cout << "==== === Loading mass model from gdml files " << fWorldFile
           << "..." << G4endl;
    parser.Read("gdml/" + fWorldFile + ".gdml");
    worldPhysical = parser.GetWorldVolume();
    worldLogical = worldPhysical->GetLogicalVolume();
  } else {
    G4cout << "## CAD models will not be imported!" << G4endl;
    G4Box *worldSolid = new G4Box("worldSolid", 500 * cm, 500 * cm, 500 * cm);
    worldLogical =
        new G4LogicalVolume(worldSolid, Vacuum, "worldLogical", 0, 0, 0);
    worldPhysical = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), worldLogical,
                                      "worldPhysical", 0, false, 0);
  }

  if (detectorEnabled) {
    // default it is 100
    // don't create detectors is it is -1
    // only one detector
    if (activatedDetectorFlag >= 0 && activatedDetectorFlag < 32)
      isSingleDetector = true;
    G4cout << "Constructing Caliste..." << G4endl;
    G4LogicalVolume *CalisteLog = ConstructCaliste();
    G4cout << "Placing detector Caliste..." << G4endl;
    for (int i = 0; i < 32; i++) {
      // placing detector
      G4ThreeVector pos = Grid::getCalisteCenterCoordsCAD(i);
      if (isSingleDetector) {
        // single detector only, for testing
        if (i != activatedDetectorFlag)
          continue;
        else {
          pos.setY(0);
          pos.setZ(-0.8 * mm);
          // move the detector to the center
        }
      }
      G4cout << ">> detector:  " << i << "  , position: " << pos << G4endl;
      new G4PVPlacement(G4Transform3D(rotMatrix, pos), CalisteLog, "Caliste",
                        worldLogical, false, i, true);
      // break;
      //
    }
    /*
    if (activatedDetectorFlag == 8||!isSingleDetector ){
            if(gridsEnabled)ConstructCFLAperture();
    }
    if (activatedDetectorFlag == 9||!isSingleDetector){
            if(gridsEnabled)ConstructBKGAperture();
    }
    if(!isSingleDetector   && usingCADModels){
            ConstructCalibrationFoil();
    }
    */
    CalisteLog->SetVisAttributes(G4VisAttributes::Invisible);
  }

  // X-ray window

  if (gridsEnabled) {
    ConstructGrids();
    ConstructCFLAperture();
    ConstructBKGAperture();
  }

  SetVisColors();

  // ConstructSpaceCraft();
  // can be activated if needed

  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4cout << "World construction completed" << G4endl;
  return worldPhysical;
}

void DetectorConstruction::ConstructCFLAperture() {

  if (cflConstructed)
    return;

  G4cout << "Constructing CFL Aperture..." << G4endl;
  G4ThreeVector pos = Grid::getGridCenterCAD(8, 1);
  if (activatedDetectorFlag >= 0 && activatedDetectorFlag < 32) {
    // construct single detector
    if (activatedDetectorFlag != 8)
      return;
    else {
      pos.setY(0);
      pos.setZ(0);
      // move it the center, keep the distance
      // for testing only
    }
  }

  // construct CFL
  G4Box *CFLBox = new G4Box("CFLBox", 12 * mm, 11 * mm, 0.2 * mm);
  G4LogicalVolume *CFLLog =
      new G4LogicalVolume(CFLBox, Tungsten, "CFLLog", 0, 0, 0);

  G4Box *CFLBigHole =
      new G4Box("CFLBigHole", 8.8 * mm / 2, 6.55 * mm / 2, 0.2 * mm);
  G4LogicalVolume *CFLBigHoleLog =
      new G4LogicalVolume(CFLBigHole, Vacuum, "CFLBigHoleLog", 0, 0, 0);

  G4Box *CFLSmallHole =
      new G4Box("CFLSmallHole", 1.1 * mm / 2, 6.55 * mm / 2, 0.2 * mm);
  G4LogicalVolume *CFLSmallHoleLog =
      new G4LogicalVolume(CFLSmallHole, Vacuum, "CFLSmallHoleLog", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 6.725 * mm, 0.2 * mm), CFLBigHoleLog,
                    "CFL_BIG_TOP_Aperture", CFLLog, false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(0, -6.725 * mm, 0.2 * mm), CFLBigHoleLog,
                    "CFL_BIG_BOTTOM_Aperture", CFLLog, false, 0, false);

  new G4PVPlacement(0, G4ThreeVector(-10.45 * mm, 6.725 * mm, 0.2 * mm),
                    CFLSmallHoleLog, "CFL_SMALL_TOP_LFET_Aperture", CFLLog,
                    false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(-10.45 * mm, -6.725 * mm, 0.2 * mm),
                    CFLSmallHoleLog, "CFL_SMALL_BOTTOM_LFET_Aperture", CFLLog,
                    false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(10.45 * mm, -6.725 * mm, 0.2 * mm),
                    CFLSmallHoleLog, "CFL_SMALL_BOTTOM_RIGHT_Aperture", CFLLog,
                    false, 0, false);
  new G4PVPlacement(0, G4ThreeVector(10.45 * mm, 6.725 * mm, 0.2 * mm),
                    CFLSmallHoleLog, "CFL_SMALL_TOP_RIGHT", CFLLog, false, 0,
                    false);
  new G4PVPlacement(G4Transform3D(rotMatrix, pos), CFLLog, "CFLLogical",
                    worldLogical, false, 0, false);
}

void DetectorConstruction::ConstructBKGAperture() {

  if (bkgConstructed)
    return;
  G4cout << "Constructing BKG Aperture..." << G4endl;

  G4ThreeVector pos = Grid::getGridCenterCAD(9, 0);

  if (activatedDetectorFlag >= 0 && activatedDetectorFlag < 32) {
    // construct single detector
    if (activatedDetectorFlag != 9)
      return;
    else {
      pos.setY(0);
      pos.setZ(0);
      // center collimator
      // for testing only
    }
  }

  G4Box *BKGBox = new G4Box("BKGBox", 12 * mm, 11 * mm, 0.2 * mm);
  G4LogicalVolume *BKGLog =
      new G4LogicalVolume(BKGBox, Tungsten, "BKGLog", 0, 0, 0);

  std::vector<G4TwoVector> vertexCoordsBKGBigAperture;
  for (int j = 0; j < 80; j++) {
    vertexCoordsBKGBigAperture.push_back(
        G4TwoVector(BKG_BIGEST_APER_X[j] * mm, BKG_BIGEST_APER_Y[j] * mm));
  }
  G4ExtrudedSolid *BKGBigAperture =
      new G4ExtrudedSolid("BKGBigAperture", vertexCoordsBKGBigAperture,
                          0.4 * 0.5 * mm, // 0.4 mm thick, halfz
                          G4TwoVector(), 1, G4TwoVector(), 1);

  G4LogicalVolume *BKGBigBigRectApertureLog =
      new G4LogicalVolume(BKGBigAperture, Vacuum, "BKGBigHoleLog", 0, 0, 0);

  G4Tubs *BKGBigRoundHole =
      new G4Tubs("BKGBigRoundHole", 0, 0.1784 * mm, 0.2 * mm, 0, 360. * deg);
  G4Tubs *BKGSmallRoundHole =
      new G4Tubs("BKGBigSmallHole", 0, 0.0564 * mm, 0.2 * mm, 0, 360 * deg);

  G4LogicalVolume *BKGBigRoundHoleLog = new G4LogicalVolume(
      BKGBigRoundHole, Vacuum, "BKGBigRoundHoleLog", 0, 0, 0);
  G4LogicalVolume *BKGSmallRoundHoleLog = new G4LogicalVolume(
      BKGSmallRoundHole, Vacuum, "BKGSmallRoundHoleLog", 0, 0, 0);

  // estimated by looking at the picture in the STIX instrument paper, area,
  // 0.01, 0.1 and 1 mm

  new G4PVPlacement(0, G4ThreeVector(-3.5 * mm, 3 * mm, 0.2 * mm),
                    BKGBigBigRectApertureLog, "PIXEL_0_Aperture", BKGLog, false,
                    0, false);
  new G4PVPlacement(0, G4ThreeVector(-1.5 * mm, 3 * mm, 0.2 * mm),
                    BKGSmallRoundHoleLog, "PIXEL_1_Aperture", BKGLog, false, 0,
                    false);
  new G4PVPlacement(0, G4ThreeVector(1.5 * mm, 3 * mm, 0.2 * mm),
                    BKGBigRoundHoleLog, "PIXEL_2_Aperture", BKGLog, false, 0,
                    false);
  new G4PVPlacement(0, G4ThreeVector(-1.5 * mm, -3 * mm, 0.2 * mm),
                    BKGBigRoundHoleLog, "PIXEL_6_Aperture", BKGLog, false, 0,
                    false);
  new G4PVPlacement(0, G4ThreeVector(1.5 * mm, -3 * mm, 0.2 * mm),
                    BKGSmallRoundHoleLog, "PIXEL_7_Aperture", BKGLog, false, 0,
                    false);
  new G4PVPlacement(0, G4ThreeVector(3.5 * mm, -3 * mm, 0.2 * mm),
                    BKGBigBigRectApertureLog, "PIXEL_8_Aperture", BKGLog, false,
                    0, false);
  // estimated by looking at the picture in the STIX instrument paper

  new G4PVPlacement(G4Transform3D(rotMatrix, pos), BKGLog, "BKG", worldLogical,
                    false, 0, false);
}

void DetectorConstruction::SetVisAttrib(G4LogicalVolume *log, G4double red,
                                        G4double green, G4double blue,
                                        G4double alpha, G4bool wireFrame,
                                        G4bool solid) {
  G4VisAttributes *visAttrib =
      new G4VisAttributes(G4Colour(red, green, blue, alpha));
  visAttrib->SetForceWireframe(wireFrame);
  visAttrib->SetForceSolid(solid);
  log->SetVisAttributes(visAttrib);
}
void DetectorConstruction::SetVisAttrib(G4LogicalVolume *log, G4double red,
                                        G4double green, G4double blue,
                                        G4double alpha) {
  SetVisAttrib(log, red, green, blue, alpha, true, true);
}

void DetectorConstruction::SetAttenuatorStatus(G4bool att) {
  // change attenuator status
  if (!att) {
    fWorldFile = "WorldAttOut";
    usingCADModels = true;
    attenuatorIn = false;
    G4cout << "Attenuator is out" << G4endl;
  } else {
    fWorldFile = "WorldAttIn";
    usingCADModels = true;
    attenuatorIn = true;
    G4cout << "Attenuator is inserted!" << G4endl;
  }
}

void DetectorConstruction::SetVisColors() {
  G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
  std::vector<G4LogicalVolume *>::const_iterator lvciter;
  for (lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++) {
    G4String volumeName = (*lvciter)->GetName();
    G4double red = G4UniformRand() * 0.7 + 0.15;
    G4double green = G4UniformRand() * 0.7 + 0.15;
    G4double blue = G4UniformRand() * 0.7 + 0.15;
    G4double alpha = 0;

    if (volumeName == "CdTeAnodeLog" || volumeName == "CdTeCathodeLog") {
      red = 1;
      green = 0;
      blue = 0;
      alpha = 0.05;
      SetVisAttrib(*lvciter, red, green, blue, alpha, true, false);
    } else if (volumeName.contains("frontGrid") ||
               volumeName.contains("rearGrid")) {
      (*lvciter)->SetVisAttributes(G4VisAttributes::Invisible);
    } else if (volumeName.contains("grid")) {
      red = 0.28;
      green = 0.28;
      blue = 0.28;
      alpha = 0.1;
      SetVisAttrib(*lvciter, red, green, blue, alpha, true, true);
    } else {
      SetVisAttrib(*lvciter, red, green, blue, alpha);
    }
    // randomize color

    double mass = (*lvciter)->GetMass() / g;
    G4cout << "~~~ The MASS of " << volumeName << " is " << mass << " g. ~~~"
           << G4endl;
  }
}
