//
/// \brief Implementation of the OPPP tracker class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Sphere.hh"
#include "G4GenericTrap.hh"
#include "G4Torus.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4GeometryManager.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"
#include "LXSetUp.hh"
#include "LxAux.hh"
#include "LxDetector.hh"
#include "LxTrackerOPPP.hh"


////////////////////////////////////////////////////////////////////////
//// LxDetectorOPPP

LxTrackerOPPP::LxTrackerOPPP(DetectorConstruction *detc) : LxDetector(detc), 
  fCoolingPipeConnectorAssembly(0), fCableTerm1Assembly(0), fCableTerm2Assembly(0),
  fOutServFlatCblAssmblyZpos(0.0),
  fInServCoolPipeAssmblyZpos(0.0),
  fInServFlatCbl0AssmblyZpos(0.0),
  fInServFlatCbl1AssmblyZpos(0.0),
  fInServFlatCbl2AssmblyZpos(0.0)
{}


void LxTrackerOPPP::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* opppContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);

  G4double dxdetc = lxs->OPPPDetTopPlateX;
  G4double dydetc = lxs->OPPPStaveLiftY + lxs->OPPPColdPlateY;
  G4double dzdetc = lxs->OPPPDetTopPlateZ;

  G4Box *solidOPPPTrackerL = new G4Box("solidOPPPTrackerL", dxdetc/2.0, dydetc/2.0, dzdetc/2.0);
  G4LogicalVolume *logicOPPPTrackerL = new G4LogicalVolume(solidOPPPTrackerL, opppContainerMaterial, "logicOPPPTrackerL");
  G4Box *solidOPPPTrackerR = new G4Box("solidOPPPTrackerR", dxdetc/2.0, dydetc/2.0, dzdetc/2.0);
  G4LogicalVolume *logicOPPPTrackerR = new G4LogicalVolume(solidOPPPTrackerR, opppContainerMaterial, "logicOPPPTrackerR");

  G4AssemblyVolume* staveAssembly = ConstructStaveAssembly();
  G4AssemblyVolume* holderAssembly = ConstructStaveHolderAssembly();
  G4LogicalVolume* staveFlexVolume = ConstructSensorFlex();
  G4AssemblyVolume* innerStaveServiceLinesAssembly = ConstructInnerServiceLinesAssembly();
  G4AssemblyVolume* outerStaveServiceLinesAssembly = ConstructOuterServiceLinesAssembly();
  G4AssemblyVolume* innerTermAssembly = ConstructInnerServiceLinesTermAssembly();

  G4double opppdetmodz = lxs->OPPPDetOffsetZ;
  G4double staveextraz = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness + lxs->OPPPStaveSideZpos;
  G4double sasxpos = (lxs->OPPPStaveX-dxdetc)/2.0;
  G4double sasypos = (lxs->OPPPColdPlateY-dydetc)/2.0 + lxs->OPPPStaveLiftY;
  G4double saszpos = (lxs->OPPPTrussZ-dzdetc)/2.0 + staveextraz + opppdetmodz;
  G4ThreeVector trmin(sasxpos, sasypos, saszpos + lxs->OPPPStaveInOutZ);
  G4double osasxpos = sasxpos + lxs->OPPPTrackerActiveX - lxs->OPPPStaveInOutX;
  G4ThreeVector trmout(osasxpos, sasypos, saszpos);

  G4double sensstavez = lxs->OPPPSensorZ + lxs->OPPPSesorFPCLayer1Z + lxs->OPPPSesorFPCLayer2Z;
  G4double flexxpos = (lxs->OPPPStaveX-dxdetc)/2.0;
  G4double flexoxpos = flexxpos + lxs->OPPPTrackerActiveX - lxs->OPPPStaveInOutX;
  G4double flexzpos = -dzdetc/2.0 + lxs->OPPPStaveSideZpos  + opppdetmodz - sensstavez/2.0;
  G4ThreeVector trflex(flexxpos, (dydetc-lxs->OPPPSensorY)/2.0, flexzpos + lxs->OPPPStaveInOutZ);
  G4ThreeVector trflexout(flexoxpos, (dydetc-lxs->OPPPSensorY)/2.0, flexzpos);

  for (G4int il = 0; il < lxs->OPPPTrackerNLayers; ++il) {
    staveAssembly->MakeImprint(logicOPPPTrackerL, trmin, 0, 0, lxs->OverlapTest);
    holderAssembly->MakeImprint(logicOPPPTrackerL, trmin, 0, 0, lxs->OverlapTest);
    staveAssembly->MakeImprint(logicOPPPTrackerL, trmout, 0, 0, lxs->OverlapTest);
    holderAssembly->MakeImprint(logicOPPPTrackerL, trmout, 0, 0, lxs->OverlapTest);

    innerStaveServiceLinesAssembly->MakeImprint(logicOPPPTrackerL, trmin, 0, 0, lxs->OverlapTest);
    outerStaveServiceLinesAssembly->MakeImprint(logicOPPPTrackerL, trmout, 0, 0, lxs->OverlapTest);
    innerTermAssembly->MakeImprint(logicOPPPTrackerL, trmin, 0, 0, lxs->OverlapTest);

    G4ThreeVector trminr(-trmin.x(), trmin.y(), trmin.z());
    G4ThreeVector trmoutr(-trmout.x(), trmout.y(), trmout.z());
    staveAssembly->MakeImprint(logicOPPPTrackerR, trminr, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI),
                               0, lxs->OverlapTest);
    holderAssembly->MakeImprint(logicOPPPTrackerR, trminr, 0, 0, lxs->OverlapTest);
    staveAssembly->MakeImprint(logicOPPPTrackerR, trmoutr, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI),
                               0, lxs->OverlapTest);
    holderAssembly->MakeImprint(logicOPPPTrackerR, trmoutr, 0, 0, lxs->OverlapTest);

    innerStaveServiceLinesAssembly->MakeImprint(logicOPPPTrackerR, trminr,
                                  new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), 0, lxs->OverlapTest);
    outerStaveServiceLinesAssembly->MakeImprint(logicOPPPTrackerR, trmoutr,
                                  new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), 0, lxs->OverlapTest);
    G4ThreeVector trmincterm(trmin.x()-lxs->OPPPDetTopPlateX+lxs->OPPPServiceCablTermBase1X, trmin.y(), trmin.z());
    innerTermAssembly->MakeImprint(logicOPPPTrackerR, trmincterm, 0, 0, lxs->OverlapTest);

    trmin.setZ(trmin.z() + lxs->OPPPTrackerInterLayerZ);
    trmout.setZ(trmout.z() + lxs->OPPPTrackerInterLayerZ);

    new G4PVPlacement (0, trflex, staveFlexVolume, "OPPPSensorFlex", logicOPPPTrackerL, false, 2*il, lxs->OverlapTest);
    new G4PVPlacement (0, trflexout, staveFlexVolume, "OPPPSensorFlex", logicOPPPTrackerL, false, 2*il+1, lxs->OverlapTest);

    G4int layerid = 2*lxs->OPPPTrackerNLayers + 2*il;
    G4ThreeVector trflexr(-trflex.x(), trflex.y(), trflex.z());
    G4ThreeVector trflexoutr(-trflexout.x(), trflexout.y(), trflexout.z());
    new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), trflexr, staveFlexVolume, "OPPPSensorFlex",
                     logicOPPPTrackerR, false, layerid, lxs->OverlapTest);
    new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), trflexoutr, staveFlexVolume, "OPPPSensorFlex",
                     logicOPPPTrackerR, false, layerid+1, lxs->OverlapTest);

    trflex.setZ(trflex.z() + lxs->OPPPTrackerInterLayerZ);
    trflexout.setZ(trflexout.z() + lxs->OPPPTrackerInterLayerZ);
  }

  G4double detxpos = dxdetc/2.0 + lxs->OPPPDetXPos;
  G4double detypos = (lxs->OPPPColdPlateY-dydetc)/2.0;
  G4double detzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet
                     + lxs->OPPPDetTopPlateZ/2.0 - opppdetmodz;
  new G4PVPlacement (0, G4ThreeVector(detxpos, detypos, detzpos), logicOPPPTrackerL, "OPPPTrackerL",
                     fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-detxpos, detypos - lxs->OPPPDetYLowPos, detzpos), logicOPPPTrackerR, "OPPPTrackerR",
                     fLogicWorld, false, 0, lxs->OverlapTest);

  //MakeImprint support assembly
  G4double supporthight;
  G4AssemblyVolume* supportAssembly = ConstructSupportAssembly(supporthight);
  G4double supypos = (lxs->OPPPColdPlateY - lxs->OPPPDetTopPlateY)/2.0 - dydetc;
  G4ThreeVector trsupport(detxpos, supypos, detzpos);
  supportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);
  trsupport.setX(-detxpos);
  trsupport.setY(trsupport.y() - lxs->OPPPDetYLowPos);
  supportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);

  //MakeImprint service support assembly
  G4double servicesupportwidth;
  G4AssemblyVolume* serviceSupportAssembly = ConstructServiceSupportAssembly(servicesupportwidth);
  G4ThreeVector trservsupport(dxdetc + lxs->OPPPDetXPos + servicesupportwidth/2.0, supypos, detzpos);
  serviceSupportAssembly->MakeImprint(fLogicWorld, trservsupport, 0, 0, lxs->OverlapTest);
  trservsupport.setX(-trservsupport.x());
  trservsupport.setY(trservsupport.y() - lxs->OPPPDetYLowPos);
  serviceSupportAssembly->MakeImprint(fLogicWorld, trservsupport,
                                      new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), 0, lxs->OverlapTest);

  ConstructSideServiceLines(fLogicWorld);

//   //MakeImprint hexapod (or something similar) assembly
  G4double hexhightl = 0.0;

  G4AssemblyVolume* hexAssemblyL = LxAux::BuildHexapod("TrackerL", hexhightl);
  G4double hexzpos = detzpos + (lxs->OPPPDetBottomSupportZ - lxs->OPPPDetTopPlateZ)/2.0;
  G4ThreeVector hextrsupport(detxpos, supypos + lxs->OPPPDetTopPlateY/2.0 - supporthight - hexhightl/2.0, hexzpos);
  hexAssemblyL->MakeImprint(fLogicWorld, hextrsupport, 0, 0, lxs->OverlapTest);

  G4double hexhightr = lxs->OPPPDetYLowPos;
  G4AssemblyVolume* hexAssemblyR = LxAux::BuildHexapod("TrackerR", hexhightr);
  hextrsupport.setX(-detxpos);
  hextrsupport.setY(hextrsupport.y() + (hexhightl - hexhightr)/2.0 - lxs->OPPPDetYLowPos);
  hexAssemblyR->MakeImprint(fLogicWorld, hextrsupport, 0, 0, lxs->OverlapTest);

  supporthight += hexhightl;

  G4double baseypos = supypos - supporthight + (lxs->OPPPColdPlateY - lxs->OPPPBasePlateY)/2.0;
  G4double basezpos = detzpos + (lxs->OPPPDetBottomSupportZ - lxs->OPPPDetTopPlateZ)/2.0;

// Table support
  G4double ypestal = 1.5*m;
  G4double ylevel = baseypos + lxs->HICSBasePlateY/2.0;
  G4double tblhight = ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4double tblx = lxs->OPPPBasePlateX;
  G4double tblz = lxs->OPPPBasePlateZ;
  G4AssemblyVolume *tablesupport = LxAux::BuildTable("TrackerEcalTable", tblx, tblhight, tblz, 4);
  G4ThreeVector trsupporttbl(0.0, ylevel, basezpos);
  tablesupport->MakeImprint(fLogicWorld, trsupporttbl, 0, 0, lxs->OverlapTest);

  G4LogicalVolume *pedstal = LxAux::BuildPedestal("TrackerEcal", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.5*ypestal + lxs->FloorSurfaceYpos, basezpos), pedstal,
                      "TrackerEcalPedestal", fLogicWorld, false, 0, lxs->OverlapTest);

  AddSegmentation();

  ConstructElectronicsRack();
}



void LxTrackerOPPP::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4double sensx = lxs->OPPPSensorNCellX * lxs->OPPPSensorPixelX;
  G4double sensy = lxs->OPPPSensorNCellY * lxs->OPPPSensorPixelY;

  fDetector->AddSensorSegmentation("OPPPSensitive", sensx, sensy,
                                                    lxs->OPPPSensorNCellX, lxs->OPPPSensorNCellY);
}



G4AssemblyVolume* LxTrackerOPPP::ConstructStaveAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* opppColdPlateMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPColdPlateMaterial);
  G4Material* opppSpaceFrameMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPFrameMaterial);
  G4Material* opppCoolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingPipeMaterial);
  G4Material* opppCoolingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* opppStaveSideMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPStaveSideMaterial);

  G4int NTruss = lxs->OPPPNTruss;
  G4double dxbot = lxs->OPPPTrussX;
  G4double dybot = lxs->OPPPColdPlateY;
  G4double dxtop = lxs->OPPPTrussTopXY;
  G4double dytop = lxs->OPPPTrussTopXY;
  G4double dz = lxs->OPPPTrussZ;
  G4double dtop = lxs->OPPPTrussTopThickness;
  G4double dside = 0.5 *mm;

  G4Trd *solidTruss1 = new G4Trd("solidTruss1", dxbot/2.0, dxtop/2.0, dybot/2.0, dytop/2.0, dz/2.0);
  G4Trd *solidTrussCut1 = new G4Trd("solidTrussCut1", dxbot, dxbot, dybot/2.0, dytop/2.0, dz/2.0);
  G4Trd *solidTrussCut2 = new G4Trd("solidTrussCut2", (dxbot-dside)/2.0, (dxtop-dside)/2.0, dybot, dybot, dz/2.0);
  G4Transform3D transform2(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, -dtop));
  G4SubtractionSolid* solidTrussS1 = new G4SubtractionSolid("solidTrussS1", solidTruss1, solidTrussCut1, transform2);
  G4SubtractionSolid* solidTruss = new G4SubtractionSolid("solidTruss", solidTrussS1, solidTrussCut2, transform2);
  G4LogicalVolume *logicTruss = new G4LogicalVolume(solidTruss, opppSpaceFrameMaterial, "logicTruss");

  G4AssemblyVolume *staveAssembly = new G4AssemblyVolume();
  G4ThreeVector asstr(0.0, 0.0, 0.0);
  for (G4int ii = 0; ii < NTruss; ++ii) {
    asstr.setX(dxbot*(ii+0.5-0.5*NTruss));
    staveAssembly->AddPlacedVolume(logicTruss, asstr, 0);
  }

  G4Box *solidTrussBar = new G4Box("solidTrussBar", (NTruss*dxbot)/2.0, dytop/2.0, dtop/2.0);
  G4LogicalVolume *logicTrussBar = new G4LogicalVolume(solidTrussBar, opppSpaceFrameMaterial, "logicTrussBar");
  asstr = G4ThreeVector(0.0, 0.0, (dz-dtop)/2.0-dtop);
  staveAssembly->AddPlacedVolume(logicTrussBar, asstr, 0);

  G4double coolpiper = lxs->OPPPCoolingPipeR;
  G4double coolwaterr = lxs->OPPPCoolingPipeR - lxs->OPPPCoolingPipeWall;
  G4double coolpipeypos = 3.0 *mm;

  G4Tubs *solidTrkCoolPipe = new G4Tubs("solidTrkCoolPipe", 0.0, coolpiper, (NTruss*dxbot)/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicTrkCoolPipe = new G4LogicalVolume(solidTrkCoolPipe, opppCoolPipeMaterial, "logicTrkCoolPipe");
  G4Tubs *solidTrkCoolWater = new G4Tubs("solidTrkCoolWater", 0.0, coolwaterr, (NTruss*dxbot)/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicTrkCoolWater = new G4LogicalVolume(solidTrkCoolWater, opppCoolingMaterial, "logicTrkCoolWater");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                    logicTrkCoolWater, "TrkCoolWater", logicTrkCoolPipe, false, 0, lxs->OverlapTest);
  asstr = G4ThreeVector(0.0, coolpipeypos, -dz/2.0+coolpiper);
  staveAssembly->AddPlacedVolume(logicTrkCoolPipe, asstr,
                                   new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0));
  asstr = G4ThreeVector(0.0, -coolpipeypos, -dz/2.0+coolpiper);
  staveAssembly->AddPlacedVolume(logicTrkCoolPipe, asstr,
                                   new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0));

  G4double dxcplate = lxs->OPPPColdPlateX;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4Box *solidColdPlate = new G4Box("solidColdPlate", dxcplate/2.0, dybot/2.0, dzcplate/2.0);
  G4LogicalVolume *logicColdPlate = new G4LogicalVolume(solidColdPlate, opppColdPlateMaterial, "logicColdPlate");
  asstr = G4ThreeVector(0.0, 0.0, -(dzcplate+dz)/2.0);
  staveAssembly->AddPlacedVolume(logicColdPlate, asstr, 0);

  G4double dxendcap = 0.5 * (dxcplate - NTruss*dxbot);
  G4Trd *solidTrussEndCap = new G4Trd("solidTrussEndCap", dxendcap/2.0, dxendcap/2.0, dybot/2.0, dytop/2.0, dz/2.0);
  G4LogicalVolume *logicTrussEndCap = new G4LogicalVolume(solidTrussEndCap, opppStaveSideMaterial, "logicTrussEndCap");
  asstr = G4ThreeVector(0.5*(dxendcap+NTruss*dxbot), 0.0, 0.0);
  staveAssembly->AddPlacedVolume(logicTrussEndCap, asstr, 0);
  asstr = G4ThreeVector(-0.5*(dxendcap+NTruss*dxbot), 0.0, 0.0);
  staveAssembly->AddPlacedVolume(logicTrussEndCap, asstr, 0);

  G4double stavex = lxs->OPPPStaveX;
  G4double stavesidey = 10.0 *mm;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double stavesidepincutr = 0.6 *mm;
  G4double stavesidescrewcutr = 1.0 *mm;
  G4double stavesidepincutx = 9.0 *mm;
  G4double stavesidescrewcutx = 5.0 *mm;

  G4double dxstaveside = 0.5*(stavex - dxcplate);
  G4Box *solidStaveSide1 = new G4Box("solidStaveSide1", dxstaveside/2.0, stavesidey/2.0, stavesidez/2.0);
  G4Tubs *solidStvSidePinCut = new G4Tubs("solidStvSidePinCut", 0.0, stavesidepincutr, stavesidez, 0.0, 2.0*M_PI);
  G4Tubs *solidStvSideScrewCut = new G4Tubs("solidStvSideScrewCut", 0.0, stavesidescrewcutr, stavesidez, 0.0, 2.0*M_PI);
  G4Transform3D transform3(G4RotationMatrix(), G4ThreeVector(stavesidepincutx-dxstaveside/2.0, 0.0, 0.0));

  G4SubtractionSolid* solidStaveSide2 = new G4SubtractionSolid("solidStaveSide2", solidStaveSide1,
                                                               solidStvSidePinCut, transform3);
  G4Transform3D transform4(G4RotationMatrix(), G4ThreeVector(stavesidescrewcutx-dxstaveside/2.0, 0.0, 0.0));
  G4SubtractionSolid* solidStaveSideIn = new G4SubtractionSolid("solidStaveSideIn", solidStaveSide2,
                                                                    solidStvSideScrewCut, transform4);
  G4SubtractionSolid* solidStaveSideOut = solidStaveSide2;
  G4LogicalVolume *logicStaveSideIn = new G4LogicalVolume(solidStaveSideIn, opppStaveSideMaterial, "logicStaveSideIn");
  G4LogicalVolume *logicStaveSideOut = new G4LogicalVolume(solidStaveSideOut, opppStaveSideMaterial, "logicStaveSideOut");

  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;
  asstr = G4ThreeVector(-0.5*(dxcplate + dxstaveside), 0.0, ecdz);
  staveAssembly->AddPlacedVolume(logicStaveSideIn, asstr, 0);
  asstr = G4ThreeVector(0.5*(dxcplate + dxstaveside), 0.0, ecdz);
  staveAssembly->AddPlacedVolume(logicStaveSideOut, asstr, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));

  return staveAssembly;
}



G4AssemblyVolume* LxTrackerOPPP::ConstructSupportAssembly(G4double &sphight)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppDetSupportMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPDetSupportMaterial);

  G4Box *solidOPPPBottomSupport = new G4Box("solidOPPPBottomSupport", lxs->OPPPDetBottomSupportX/2.0,
                                            lxs->OPPPDetBottomSupportY/2.0, lxs->OPPPDetBottomSupportZ/2.0);
  G4LogicalVolume *logicOPPPBottomSupport = new G4LogicalVolume(solidOPPPBottomSupport, opppDetSupportMaterial,
                                                                "logicOPPPBottomSupport");

  G4Box *solidOPPPTopPlate = new G4Box("solidOPPPTopPlate", lxs->OPPPDetTopPlateX/2.0,
                                            lxs->OPPPDetTopPlateY/2.0, lxs->OPPPDetTopPlateZ/2.0);
  G4LogicalVolume *logicOPPPTopPlate = new G4LogicalVolume(solidOPPPTopPlate, opppDetSupportMaterial,
                                                                "logicOPPPTopPlate");

  G4Sphere *solidOPPPDetSupportBall = new G4Sphere("solidOPPPDetSupportBall", 0.0, lxs->OPPPDetSupportBallH/2.0,
                                                   0.0, 2.0*M_PI, 0.0, M_PI);
  G4LogicalVolume *logicOPPPDetSupportBall = new G4LogicalVolume(solidOPPPDetSupportBall, opppDetSupportMaterial,
                                                                "logicOPPPDetSupportBall");

  G4AssemblyVolume *oppSupportAssembly = new G4AssemblyVolume();

  G4ThreeVector asstr(0.0, 0.0, 0.0);
  oppSupportAssembly->AddPlacedVolume(logicOPPPTopPlate, asstr, 0);
  asstr.setY(-(lxs->OPPPDetBottomSupportY + lxs->OPPPDetTopPlateY)/2.0 - lxs->OPPPDetSupportBallH);
  asstr.setZ((lxs->OPPPDetBottomSupportZ - lxs->OPPPDetTopPlateZ)/2.0);
  oppSupportAssembly->AddPlacedVolume(logicOPPPBottomSupport, asstr, 0);

  G4double dbz = 0.5*lxs->OPPPDetTopPlateZ - 1.05*lxs->OPPPDetSupportBallH;
  G4double dbx = 0.5*lxs->OPPPDetTopPlateX - 1.05*lxs->OPPPDetSupportBallH;
  G4ThreeVector bpos(dbx, asstr.y() + (lxs->OPPPDetBottomSupportY+lxs->OPPPDetSupportBallH)/2.0 , dbz);
  oppSupportAssembly->AddPlacedVolume(logicOPPPDetSupportBall, bpos, 0);
  bpos.setX(-dbx);
  oppSupportAssembly->AddPlacedVolume(logicOPPPDetSupportBall, bpos, 0);
  bpos.setX(0.0);
  bpos.setZ(-dbz);
  oppSupportAssembly->AddPlacedVolume(logicOPPPDetSupportBall, bpos, 0);

  sphight = lxs->OPPPDetBottomSupportY + lxs->OPPPDetTopPlateY + lxs->OPPPDetSupportBallH;
  return oppSupportAssembly;
}



G4AssemblyVolume* LxTrackerOPPP::ConstructStaveHolderAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppStaveHolderMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPStaveHolderMaterial);
//Base
  G4Box *solidAngleHold1 = new G4Box("solidAngleHold1", lxs->OPPPAngleHolderBaseX/2.0,
                                            lxs->OPPPAngleHolderBaseY/2.0, lxs->OPPPAngleHolderZ/2.0);
  G4Box *solidAngleHoldCut = new G4Box("solidAngleHoldCut", lxs->OPPPAngleHolderBaseX,
                                            lxs->OPPPAngleHolderBaseY/2.0, lxs->OPPPAngleHolderZ/2.0);
  G4Transform3D transform1(G4RotationMatrix(), G4ThreeVector(0.0, lxs->OPPPAngleHolderH, lxs->OPPPAngleHolderBaseT));
  G4SubtractionSolid* solidAngleHold2 = new G4SubtractionSolid("solidAngleHold2", solidAngleHold1,
                                                              solidAngleHoldCut, transform1);

  G4Box *solidAngleHoldLCut = new G4Box("solidAngleHoldLCut", lxs->OPPPAngleHolderLCutX/2.0,
                                            lxs->OPPPAngleHolderLCutY/2.0, lxs->OPPPAngleHolderBaseT);
  G4double trlcuty = 0.5*(lxs->OPPPAngleHolderLCutY-lxs->OPPPAngleHolderBaseY)
                     + lxs->OPPPAngleHolderH + lxs->OPPPAngleHolderRCutR;
  G4double trlcutz = -0.5 * (lxs->OPPPAngleHolderZ - lxs->OPPPAngleHolderBaseT);
  G4Transform3D trlcut(G4RotationMatrix(), G4ThreeVector(0.0, trlcuty, trlcutz));
  G4SubtractionSolid* solidAngleHold3 = new G4SubtractionSolid("solidAngleHold3", solidAngleHold2,
                                                              solidAngleHoldLCut, trlcut);

  G4Tubs *solidAngleHoldRCut = new G4Tubs("solidAngleHoldRCut", 0.0, lxs->OPPPAngleHolderRCutR,
                                                 lxs->OPPPAngleHolderBaseT, 0.0, 2.0*M_PI);
  G4Transform3D trrcut0(G4RotationMatrix(), G4ThreeVector(0.0, trlcuty-lxs->OPPPAngleHolderLCutY/2.0, trlcutz));
  G4SubtractionSolid* solidAngleHold = new G4SubtractionSolid("solidAngleHold4", solidAngleHold3,
                                                              solidAngleHoldRCut, trrcut0);
// //Tip without filet
//   G4double atipy = lxs->OPPPAngleHolderY - lxs->OPPPAngleHolderBaseY;
//   G4Box *solidAngleHoldTip1 = new G4Box("solidAngleHoldTip1", lxs->OPPPAngleHolderX/2.0, atipy/2.0,
//                                            lxs->OPPPAngleHolderT/2.0);
//
//   G4double trtiplcuty = trlcuty - lxs->OPPPAngleHolderBaseY;
//   G4Transform3D trtiplcut(G4RotationMatrix(), G4ThreeVector(0.0, trtiplcuty, 0.0));
//   G4SubtractionSolid* solidAngleHoldTip2 = new G4SubtractionSolid("solidAngleHoldTip2", solidAngleHoldTip1,
//                                                               solidAngleHoldLCut, trtiplcut);
//   G4Transform3D trrcut1(G4RotationMatrix(), G4ThreeVector(0.0, trtiplcuty+lxs->OPPPAngleHolderLCutY/2.0, 0.0));
//   G4SubtractionSolid* solidAngleHoldTip3 = new G4SubtractionSolid("solidAngleHoldTip3", solidAngleHoldTip2,
//                                                               solidAngleHoldRCut, trrcut1);
//   G4double topcutx = sqrt(2) * lxs->OPPPAngleHolderTopCutXY;
//   G4Box *solidAngleHoldTopCut = new G4Box("solidAngleHoldTopCut", topcutx/2.0, topcutx, lxs->OPPPAngleHolderT);
//   G4Transform3D trtopcut0(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI/4.0),
//                           G4ThreeVector(lxs->OPPPAngleHolderX/2.0, atipy/2.0, 0.0));
//   G4SubtractionSolid* solidAngleHoldTip4 = new G4SubtractionSolid("solidAngleHoldTip4", solidAngleHoldTip3,
//                                                                 solidAngleHoldTopCut, trtopcut0);
//   G4Transform3D trtopcut1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, -1.0), M_PI/4.0),
//                           G4ThreeVector(-lxs->OPPPAngleHolderX/2.0, atipy/2.0, 0.0));
//   G4SubtractionSolid* solidAngleHoldTip = new G4SubtractionSolid("solidAngleHoldTip", solidAngleHoldTip4,
//                                                                 solidAngleHoldTopCut, trtopcut1);

//Tip with filet
  G4double atipy = lxs->OPPPAngleHolderY - lxs->OPPPAngleHolderBaseY;
  G4Box *solidAngleHoldTip0 = new G4Box("solidAngleHoldTip0", lxs->OPPPAngleHolderX/2.0, atipy/2.0,
                                        lxs->OPPPAngleHolderBaseT/2.0);
  G4Box *solidAngleHoldTipFCut = new G4Box("solidAngleHoldTipFCut", lxs->OPPPAngleHolderX, atipy/2.0,
                                        lxs->OPPPAngleHolderBaseT/2.0);
  G4double trtipfcutr = lxs->OPPPAngleHolderBaseT - lxs->OPPPAngleHolderT;
  G4double trtipfcutz = lxs->OPPPAngleHolderT;
  G4Transform3D trtipfcut(G4RotationMatrix(), G4ThreeVector(0.0, trtipfcutr, trtipfcutz));
  G4SubtractionSolid* solidAngleHoldTip01 = new G4SubtractionSolid("solidAngleHoldTip01", solidAngleHoldTip0,
                                                                   solidAngleHoldTipFCut, trtipfcut);
  G4Tubs *solidAngleHoldRFCut = new G4Tubs("solidAngleHoldRFCut", 0.0, trtipfcutr,
                                                 lxs->OPPPAngleHolderX, 0.0, 2.0*M_PI);
  G4Transform3D trtipfrcut(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0),
                           G4ThreeVector(0.0, trtipfcutr-atipy/2.0, lxs->OPPPAngleHolderBaseT/2.0));
  G4SubtractionSolid* solidAngleHoldTip1 = new G4SubtractionSolid("solidAngleHoldTip1", solidAngleHoldTip01,
                                                                   solidAngleHoldRFCut, trtipfrcut);
  G4double trtiplcuty = trlcuty - lxs->OPPPAngleHolderBaseY;
  G4Transform3D trtiplcut(G4RotationMatrix(), G4ThreeVector(0.0, trtiplcuty, 0.0));
  G4SubtractionSolid* solidAngleHoldTip2 = new G4SubtractionSolid("solidAngleHoldTip2", solidAngleHoldTip1,
                                                              solidAngleHoldLCut, trtiplcut);
  G4Transform3D trrcut1(G4RotationMatrix(), G4ThreeVector(0.0, trtiplcuty+lxs->OPPPAngleHolderLCutY/2.0, 0.0));
  G4SubtractionSolid* solidAngleHoldTip3 = new G4SubtractionSolid("solidAngleHoldTip3", solidAngleHoldTip2,
                                                              solidAngleHoldRCut, trrcut1);
  G4double topcutx = sqrt(2) * lxs->OPPPAngleHolderTopCutXY;
  G4Box *solidAngleHoldTopCut = new G4Box("solidAngleHoldTopCut", topcutx/2.0, topcutx, lxs->OPPPAngleHolderBaseT);
  G4Transform3D trtopcut0(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI/4.0),
                          G4ThreeVector(lxs->OPPPAngleHolderX/2.0, atipy/2.0, 0.0));
  G4SubtractionSolid* solidAngleHoldTip4 = new G4SubtractionSolid("solidAngleHoldTip4", solidAngleHoldTip3,
                                                                solidAngleHoldTopCut, trtopcut0);
  G4Transform3D trtopcut1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, -1.0), M_PI/4.0),
                          G4ThreeVector(-lxs->OPPPAngleHolderX/2.0, atipy/2.0, 0.0));
  G4SubtractionSolid* solidAngleHoldTip = new G4SubtractionSolid("solidAngleHoldTip", solidAngleHoldTip4,
                                                                solidAngleHoldTopCut, trtopcut1);

  G4LogicalVolume *logicAngleHold = new G4LogicalVolume(solidAngleHold, opppStaveHolderMaterial, "logicAngleHold");
  G4LogicalVolume *logicAngleHoldTip = new G4LogicalVolume(solidAngleHoldTip, opppStaveHolderMaterial, "logicAngleHoldTip");

  G4double dybot = lxs->OPPPColdPlateY;
  G4double dz = lxs->OPPPTrussZ;
  G4double stavex = lxs->OPPPStaveX;
  G4double dxcplate = lxs->OPPPColdPlateX;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double stavesidez = lxs->OPPPStaveSideZ;

  G4double dxstaveside = 0.5*(stavex - dxcplate);
  G4double ahxpos = 0.5*(dxcplate + dxstaveside);
//   G4double ahypos = (lxs->OPPPAngleHolderY - dybot)/2.0 - lxs->OPPPStaveLiftY;
  G4double ahypos = (lxs->OPPPAngleHolderBaseY - dybot)/2.0 - lxs->OPPPStaveLiftY;
  G4double ahzpos = 0.5*(lxs->OPPPAngleHolderZ-dz) - dzcplate - stavesidezpos + stavesidez;

  G4AssemblyVolume *holderAssembly = new G4AssemblyVolume();
  G4ThreeVector asstr(ahxpos, ahypos, ahzpos);
  holderAssembly->AddPlacedVolume(logicAngleHold, asstr, 0);
  asstr.setX(-0.5*(dxcplate + dxstaveside));
  holderAssembly->AddPlacedVolume(logicAngleHold, asstr, 0);

  G4double ahtipy = ahypos + 0.5*(lxs->OPPPAngleHolderBaseY + atipy);
//   G4double ahtipz = ahzpos + 0.5*(lxs->OPPPAngleHolderT - lxs->OPPPAngleHolderZ);
  G4double ahtipz = ahzpos + 0.5*(lxs->OPPPAngleHolderBaseT - lxs->OPPPAngleHolderZ);
  G4ThreeVector tratip(ahxpos, ahtipy, ahtipz);
  holderAssembly->AddPlacedVolume(logicAngleHoldTip, tratip, 0);
  tratip.setX(-0.5*(dxcplate + dxstaveside));
  holderAssembly->AddPlacedVolume(logicAngleHoldTip, tratip, 0);

  return holderAssembly;
}



G4LogicalVolume* LxTrackerOPPP::ConstructSensorFlex()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* sensorContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* opppFPCMaterial1 = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPSesorFPCMaterial1);
  G4Material* opppFPCMaterial2 = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPSesorFPCMaterial2);
  G4Material* opppSensorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPTrackerMaterial);

  G4double sensstavex = lxs->OPPPNStaveSensors * (lxs->OPPPStaveSensorsGapX + lxs->OPPPSensorX) - lxs->OPPPStaveSensorsGapX;
  G4double sensstavez = lxs->OPPPSensorZ + lxs->OPPPSesorFPCLayer1Z + lxs->OPPPSesorFPCLayer2Z;
  
  G4Box *solidStaveSensorContainer = new G4Box("solidStaveSensorContainer", sensstavex/2.0, 
                                            lxs->OPPPSensorY/2.0, sensstavez/2.0);
  G4LogicalVolume *logicStaveSensorContainer = new G4LogicalVolume(solidStaveSensorContainer, sensorContainerMaterial,
                                                                "logicStaveSensorContainer");
  G4Box *solidFPCLayer1 = new G4Box("solidFPCLayer1", sensstavex/2.0, 
                                            lxs->OPPPSensorY/2.0, lxs->OPPPSesorFPCLayer1Z/2.0);
  G4LogicalVolume *logicFPCLayer1 = new G4LogicalVolume(solidFPCLayer1, opppFPCMaterial1, "logicFPCLayer1");
  G4Box *solidFPCLayer2 = new G4Box("solidFPCLayer2", sensstavex/2.0, 
                                            lxs->OPPPSensorY/2.0, lxs->OPPPSesorFPCLayer2Z/2.0);
  G4LogicalVolume *logicFPCLayer2 = new G4LogicalVolume(solidFPCLayer2, opppFPCMaterial2, "logicFPCLayer2");
  
  G4double lzpos = (sensstavez-lxs->OPPPSesorFPCLayer1Z)/2.0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lzpos),
                    logicFPCLayer1, "FPCLayer1", logicStaveSensorContainer, false, 0, lxs->OverlapTest);
  lzpos -= (lxs->OPPPSesorFPCLayer1Z + lxs->OPPPSesorFPCLayer2Z)/2.0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lzpos),
                    logicFPCLayer2, "FPCLayer2", logicStaveSensorContainer, false, 0, lxs->OverlapTest);

  G4Box *solidOPPPSensor1 = new G4Box("solidOPPPSensor1", lxs->OPPPSensorX/2.0, lxs->OPPPSensorY/2.0, lxs->OPPPSensorZ/2.0);
  G4double sensx = lxs->OPPPSensorNCellX * lxs->OPPPSensorPixelX;
  G4double sensy = lxs->OPPPSensorNCellY * lxs->OPPPSensorPixelY;
  G4Box *solidOPPPSensorSensitive = new G4Box("solidOPPPSensorSensitive", sensx/2.0, sensy/2.0, lxs->OPPPSensorPixelZ/2.0);
  G4Transform3D trnssens(G4RotationMatrix(), 
                           G4ThreeVector(0.0, (lxs->OPPPSensorY-sensy)/2.0, (lxs->OPPPSensorPixelZ-lxs->OPPPSensorZ)/2.0));
  G4SubtractionSolid* solidOPPPSensor = new G4SubtractionSolid("solidOPPPSensor", solidOPPPSensor1, 
                                                                solidOPPPSensorSensitive, trnssens);

  G4LogicalVolume *logicOPPPSensor = new G4LogicalVolume(solidOPPPSensor, opppSensorMaterial, "logicOPPPSensor");
  G4LogicalVolume *logicOPPPSensitive = new G4LogicalVolume(solidOPPPSensorSensitive,
                                                            opppSensorMaterial, "logicOPPPSensitive");

  G4double detxpos = (lxs->OPPPSensorX - sensstavex)/2.0;
  lzpos -= (lxs->OPPPSesorFPCLayer2Z + lxs->OPPPSensorZ)/2.0;
  for (G4int idet = 0; idet < lxs->OPPPNStaveSensors; ++idet) {
    new G4PVPlacement(0, G4ThreeVector(detxpos, 0.0, lzpos),
                    logicOPPPSensor, "OPPPSensor", logicStaveSensorContainer, false, idet, lxs->OverlapTest);

    new G4PVPlacement(0, G4ThreeVector(detxpos, (lxs->OPPPSensorY-sensy)/2.0,
                                       lzpos - (lxs->OPPPSensorZ-lxs->OPPPSensorPixelZ)/2.0),
                    logicOPPPSensitive, "OPPPSensitive", logicStaveSensorContainer, false, idet, lxs->OverlapTest);
    detxpos += lxs->OPPPSensorX + lxs->OPPPStaveSensorsGapX;
  }

  return logicStaveSensorContainer;
}



G4AssemblyVolume* LxTrackerOPPP::ConstructServiceSupportAssembly(G4double &width)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppServiceSupMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPDetSupportMaterial);
  G4Material* opppSSTMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceSupportTubeMaterial);

  G4Box *solidServiceSupport = new G4Box("solidServiceSupport", lxs->OPPPServiceSupportX/2.0,
                                          lxs->OPPPServiceSupportY/2.0, lxs->OPPPServiceSupportZ/2.0);
  G4LogicalVolume *logicOPPPServiceSupport = new G4LogicalVolume(solidServiceSupport, opppServiceSupMaterial,
                                                                 "logicOPPPServiceSupport");

  G4double rout = lxs->OPPPServiceSupportTubeRin + lxs->OPPPServiceSupportTubeThickness;
  G4Tubs *solidOPPPServiceSupportTube = new G4Tubs("solidOPPPServiceSupportTube", lxs->OPPPServiceSupportTubeRin,
                                                    rout, lxs->OPPPServiceSupportTubeX/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPServiceSupportTube = new G4LogicalVolume(solidOPPPServiceSupportTube,
                                                            opppSSTMaterial, "logicOPPPServiceSupportTube");

  G4AssemblyVolume *serviceSupportAssembly = new G4AssemblyVolume();
  G4ThreeVector asstr(lxs->OPPPServiceSupportTubeX/2.0, 0.0, 0.0);
  serviceSupportAssembly->AddPlacedVolume(logicOPPPServiceSupport, asstr, 0);
  asstr.setX(-lxs->OPPPServiceSupportX/2.0);
  asstr.setZ(lxs->OPPPServiceSupportZ/4.0);
  serviceSupportAssembly->AddPlacedVolume(logicOPPPServiceSupportTube, asstr,
                                  new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0));
  asstr.setZ(-lxs->OPPPServiceSupportZ/4.0);
  serviceSupportAssembly->AddPlacedVolume(logicOPPPServiceSupportTube, asstr,
                                  new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0));
  width = lxs->OPPPServiceSupportX + lxs->OPPPServiceSupportTubeX;

  return serviceSupportAssembly;

}



G4AssemblyVolume* LxTrackerOPPP::ConstructCoolingPipeConnectorAssembly()
{
  if (fCoolingPipeConnectorAssembly) return fCoolingPipeConnectorAssembly;

  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* coolPipeCnctrMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeCnctrMaterial);
  G4Material* opppCoolantMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* coolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeMaterial);

// Connector for cooling pipe
  G4Tubs *solidInServCoolPipeCnctr = new G4Tubs("solidInServCoolPipeCnctr", lxs->OPPPServiceCoolPipeCnctrRin,
                lxs->OPPPServiceCoolPipeCnctrRout, lxs->OPPPServiceCoolPipeCnctrL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipeCnctr = new G4LogicalVolume(solidInServCoolPipeCnctr, coolPipeCnctrMaterial,
                                                  "logicInServCoolPipeCnctr");
  G4Tubs *solidInServCoolPipeCnctrW = new G4Tubs("solidInServCoolPipeCnctrW", 0.0,
                lxs->OPPPServiceCoolPipeCnctrRin, lxs->OPPPServiceCoolPipeCnctrL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipeCnctrW = new G4LogicalVolume(solidInServCoolPipeCnctrW, opppCoolantMaterial,
                                                  "logicInServCoolPipeCnctrW");
  G4Tubs *solidInServCoolPipeCnctrS = new G4Tubs("solidInServCoolPipeCnctrS", lxs->OPPPServiceCoolPipeCnctrSRin,
                lxs->OPPPServiceCoolPipeCnctrSRout, lxs->OPPPServiceCoolPipeCnctrSL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipeCnctrS = new G4LogicalVolume(solidInServCoolPipeCnctrS, coolPipeCnctrMaterial,
                                                  "logicInServCoolPipeCnctrS");
  G4Tubs *solidInServCoolPipeCnctrSW = new G4Tubs("solidInServCoolPipeCnctrSW", 0.0,
                lxs->OPPPServiceCoolPipeCnctrSRin, lxs->OPPPServiceCoolPipeCnctrSL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipeCnctrSW = new G4LogicalVolume(solidInServCoolPipeCnctrSW, opppCoolantMaterial,
                                                  "logicInServCoolPipeCnctrSW");
  G4Tubs *solidInServCoolPipeCnctrSP = new G4Tubs("solidInServCoolPipeCnctrSP", lxs->OPPPServiceCoolPipeCnctrSRout,
                lxs->OPPPServiceCoolPipeRout, lxs->OPPPServiceCoolPipeCnctrSL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipeCnctrSP = new G4LogicalVolume(solidInServCoolPipeCnctrSP, coolPipeMaterial,
                                                  "logicInServCoolPipeCnctrSP");
// Assembling cooling pipe with connector and water
  G4AssemblyVolume *coolingPipeConnectorAssembly = new G4AssemblyVolume();

  G4ThreeVector trcntr(lxs->OPPPServiceCoolPipeCnctrL/2.0, 0.0, 0.0);
  G4RotationMatrix *rotcntr = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  coolingPipeConnectorAssembly->AddPlacedVolume(logicInServCoolPipeCnctr, trcntr, rotcntr);
  coolingPipeConnectorAssembly->AddPlacedVolume(logicInServCoolPipeCnctrW, trcntr, rotcntr);
  trcntr.setX(trcntr.x() + (lxs->OPPPServiceCoolPipeCnctrL + lxs->OPPPServiceCoolPipeCnctrSL)/2.0);
  coolingPipeConnectorAssembly->AddPlacedVolume(logicInServCoolPipeCnctrS, trcntr, rotcntr);
  coolingPipeConnectorAssembly->AddPlacedVolume(logicInServCoolPipeCnctrSW, trcntr, rotcntr);
  coolingPipeConnectorAssembly->AddPlacedVolume(logicInServCoolPipeCnctrSP, trcntr, rotcntr);

  fCoolingPipeConnectorAssembly = coolingPipeConnectorAssembly;
  return coolingPipeConnectorAssembly;
}


G4AssemblyVolume* LxTrackerOPPP::ConstructInnerServiceLinesAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCoolantMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* coolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeMaterial);
  G4Material* cableMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablMaterial);

  G4double coolprin = lxs->OPPPServiceCoolPipeRout - lxs->OPPPServiceCoolPipeThick;
  G4double inStaveCoolPipeX = lxs->OPPPDetTopPlateX - lxs->OPPPStaveX 
                            - lxs->OPPPServiceCoolPipeCnctrL - lxs->OPPPServiceCoolPipeCnctrSL;

// Bent cooling pipe
  G4double pbphi = acos(1.0 - lxs->OPPPInnerServiceCoolPipeDZ / (2.0 * lxs->OPPPInnerServiceCoolPipeBentR));
  G4double bentx = 2.0 * lxs->OPPPInnerServiceCoolPipeBentR * sin(pbphi);
  G4double inStaveCoolPipeL0 = 0.2 * (inStaveCoolPipeX - bentx);
  G4double inStaveCoolPipeL1 = inStaveCoolPipeX - bentx - inStaveCoolPipeL0;

  G4Tubs *solidInServCoolPipe0 = new G4Tubs("solidInServCoolPipe0", coolprin, lxs->OPPPServiceCoolPipeRout,
                                           inStaveCoolPipeL0/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipe0 = new G4LogicalVolume(solidInServCoolPipe0, coolPipeMaterial,
                                                  "logicInServCoolPipe0");
  G4Tubs *solidInServCoolPipe0W = new G4Tubs("solidInServCoolPipe0W", 0.0, coolprin, inStaveCoolPipeL0/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipe0W = new G4LogicalVolume(solidInServCoolPipe0W, opppCoolantMaterial,
                                                  "logicInServCoolPipe0W");
  G4Tubs *solidInServCoolPipe1 = new G4Tubs("solidInServCoolPipe1", coolprin, lxs->OPPPServiceCoolPipeRout,
                                           inStaveCoolPipeL1/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipe1 = new G4LogicalVolume(solidInServCoolPipe1, coolPipeMaterial,
                                                  "logicInServCoolPipe1");
  G4Tubs *solidInServCoolPipe1W = new G4Tubs("solidInServCoolPipe1W", 0.0, coolprin, inStaveCoolPipeL1/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicInServCoolPipe1W = new G4LogicalVolume(solidInServCoolPipe1W, opppCoolantMaterial,
                                                  "logicInServCoolPipe1W");
  G4Torus *solidInServCoolPipeR = new G4Torus("solidInServCoolPipeR", coolprin, lxs->OPPPServiceCoolPipeRout,
                                     lxs->OPPPInnerServiceCoolPipeBentR, 0.0, pbphi);
  G4LogicalVolume *logicInServCoolPipeR = new G4LogicalVolume(solidInServCoolPipeR, coolPipeMaterial,
                                                  "logicInServCoolPipeR");
  G4Torus *solidInServCoolPipeRW = new G4Torus("solidInServCoolPipeRW", 0.0, coolprin,
                                               lxs->OPPPInnerServiceCoolPipeBentR, 0.0, pbphi);
  G4LogicalVolume *logicInServCoolPipeRW = new G4LogicalVolume(solidInServCoolPipeRW, opppCoolantMaterial,
                                                  "logicInServCoolPipeRW");
// Assembling cooling pipe with connector and water
  G4AssemblyVolume *coolingPipeConnectorAssembly = ConstructCoolingPipeConnectorAssembly();
  G4AssemblyVolume *coolingPipeAssembly = new G4AssemblyVolume();

  G4ThreeVector trcntr(lxs->OPPPStaveX/2.0, 0.0, 0.0);
  G4RotationMatrix *rotcntr = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
//   coolingPipeAssembly->AddPlacedAssembly(coolingPipeConnectorAssembly, trcntr, 0);
  LxAux::AddAssmblyVolumes(coolingPipeAssembly, coolingPipeConnectorAssembly, trcntr, 0);
  trcntr.setX(trcntr.x() + lxs->OPPPServiceCoolPipeCnctrL + lxs->OPPPServiceCoolPipeCnctrSL/2.0);
  trcntr.setX(trcntr.x() + (lxs->OPPPServiceCoolPipeCnctrSL + inStaveCoolPipeL0)/2.0);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipe0, trcntr, rotcntr);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipe0W, trcntr, rotcntr);
  trcntr.setX(trcntr.x() + bentx + (inStaveCoolPipeL0 + inStaveCoolPipeL1)/2.0);
  trcntr.setZ(trcntr.z() + lxs->OPPPInnerServiceCoolPipeDZ);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipe1, trcntr, rotcntr);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipe1W, trcntr, rotcntr);
  fInServCoolPipeAssmblyZpos = trcntr.z();

  G4double trpiperx = lxs->OPPPStaveX/2.0 + lxs->OPPPServiceCoolPipeCnctrL
                    + lxs->OPPPServiceCoolPipeCnctrSL + inStaveCoolPipeL0;
  G4ThreeVector trpiper(trpiperx, 0.0, lxs->OPPPInnerServiceCoolPipeBentR);
  G4RotationMatrix *rotpiper = new G4RotationMatrix();
  rotpiper->rotateX(M_PI/2.0);
  rotpiper->rotateY(M_PI/2.0);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipeR, trpiper, rotpiper);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipeRW, trpiper, rotpiper);
  rotpiper->rotateY(M_PI);
  trpiper.setX(trpiperx + bentx);
  trpiper.setZ(lxs->OPPPInnerServiceCoolPipeDZ-lxs->OPPPInnerServiceCoolPipeBentR);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipeR, trpiper, rotpiper);
  coolingPipeAssembly->AddPlacedVolume(logicInServCoolPipeRW, trpiper, rotpiper);

// Flat cbales assembly
  G4double sensstavex = lxs->OPPPNStaveSensors * (lxs->OPPPStaveSensorsGapX + lxs->OPPPSensorX) - lxs->OPPPStaveSensorsGapX;
  G4double flatcabledx = lxs->OPPPDetTopPlateX - 0.5*(lxs->OPPPStaveX + sensstavex);

  G4double thetsec0 = atan2(lxs->OPPPInnerServiceCablSec0Z, lxs->OPPPInnerServiceCablSec0X);
  G4double sec0dx = lxs->OPPPInnerServiceCablSec0Thick / cos(thetsec0);
  G4Para *solidInServFlatCblSec0 = new G4Para("solidInServFlatCblSec0", sec0dx/2.0, lxs->OPPPInnerServiceCablSec0Y/2.0,
                                              lxs->OPPPInnerServiceCablSec0X/2.0, 0.0, thetsec0, 0.0);
  G4LogicalVolume *logicInServFlatCblSec0 = new G4LogicalVolume(solidInServFlatCblSec0, cableMaterial,
                                                  "logicInServFlatCblSec0");

  G4Box *solidInServFlatCblSec1 = new G4Box("solidInServFlatCblSec1", lxs->OPPPInnerServiceCablSec1X/2.0,
                                            lxs->OPPPInnerServiceCablSec1Y/2.0, lxs->OPPPInnerServiceCablSec1Thick/2.0);
  G4LogicalVolume *logicInServFlatCblSec1 = new G4LogicalVolume(solidInServFlatCblSec1, cableMaterial,
                                                  "logicInServFlatCblSec1");

  G4double sec20dx = flatcabledx - lxs->OPPPInnerServiceCablSec0X - lxs->OPPPInnerServiceCablSec1X;
  G4Box *solidInServFlatCblSec20 = new G4Box("solidInServFlatCblSec20", sec20dx/2.0, lxs->OPPPInnerServiceCablSec20Y/2.0,
                                              lxs->OPPPInnerServiceCablSec20Thick/2.0);
  G4LogicalVolume *logicInServFlatCblSec20 = new G4LogicalVolume(solidInServFlatCblSec20, cableMaterial,
                                                  "logicInServFlatCblSec20");

  G4double thetsec21 = atan2(lxs->OPPPInnerServiceCablSec21Z, lxs->OPPPInnerServiceCablSec21X);
  G4double sec21dx = lxs->OPPPInnerServiceCablSec21Thick / cos(thetsec21);
  G4Para *solidInServFlatCblSec21 = new G4Para("solidInServFlatCblSec21", sec21dx/2.0, lxs->OPPPInnerServiceCablSec21Y/2.0,
                                              lxs->OPPPInnerServiceCablSec21X/2.0, 0.0, thetsec21, M_PI);
  G4LogicalVolume *logicInServFlatCblSec21 = new G4LogicalVolume(solidInServFlatCblSec21, cableMaterial,
                                                  "logicInServFlatCblSec21");

  G4double thetsec22 = atan2(lxs->OPPPInnerServiceCablSec22Z, lxs->OPPPInnerServiceCablSec22X);
  G4double sec22dx = lxs->OPPPInnerServiceCablSec22Thick / cos(thetsec22);
  G4Para *solidInServFlatCblSec22 = new G4Para("solidInServFlatCblSec22", sec22dx/2.0, lxs->OPPPInnerServiceCablSec22Y/2.0,
                                              lxs->OPPPInnerServiceCablSec22X/2.0, 0.0, thetsec22, M_PI);
  G4LogicalVolume *logicInServFlatCblSec22 = new G4LogicalVolume(solidInServFlatCblSec22, cableMaterial,
                                                  "logicInServFlatCblSec22");

  G4double sec31dx = sec20dx - lxs->OPPPInnerServiceCablSec21X;
  G4Box *solidInServFlatCblSec31 = new G4Box("solidInServFlatCblSec31", sec31dx/2.0, lxs->OPPPInnerServiceCablSec31Y/2.0,
                                              lxs->OPPPInnerServiceCablSec31Thick/2.0);
  G4LogicalVolume *logicInServFlatCblSec31 = new G4LogicalVolume(solidInServFlatCblSec31, cableMaterial,
                                                  "logicInServFlatCblSec31");

  G4double sec32dx = sec20dx - lxs->OPPPInnerServiceCablSec22X;
  G4Box *solidInServFlatCblSec32 = new G4Box("solidInServFlatCblSec32", sec32dx/2.0, lxs->OPPPInnerServiceCablSec32Y/2.0,
                                              lxs->OPPPInnerServiceCablSec32Thick/2.0);
  G4LogicalVolume *logicInServFlatCblSec32 = new G4LogicalVolume(solidInServFlatCblSec32, cableMaterial,
                                                  "logicInServFlatCblSec32");
// Assembling differnet sections of the flat cables
  G4AssemblyVolume *innerServiceCablesAssembly = new G4AssemblyVolume();

  G4double trcblsec0x = 0.5*(sensstavex + lxs->OPPPInnerServiceCablSec0X);
  G4double trcblsec0z = lxs->OPPPStaveSideZpos -lxs->OPPPStaveSideZ/2.0 - 0.5*(lxs->OPPPInnerServiceCablSec0Z + sec0dx);
  G4ThreeVector trcblsec0(trcblsec0x, 0.0, trcblsec0z);
  G4RotationMatrix *rotcblsec0 = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec0, trcblsec0, rotcblsec0);
  G4ThreeVector trcblsec1(trcblsec0x + 0.5*(lxs->OPPPInnerServiceCablSec0X+lxs->OPPPInnerServiceCablSec1X), 0.0,
               trcblsec0z - 0.5*(lxs->OPPPInnerServiceCablSec0Z + sec0dx - lxs->OPPPInnerServiceCablSec1Thick));
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec1, trcblsec1, 0);
  G4ThreeVector trcblsec20(trcblsec1.x() + 0.5*(lxs->OPPPInnerServiceCablSec1X + sec20dx), 0.0,
                 trcblsec1.z() + 0.5*(lxs->OPPPInnerServiceCablSec20Thick - lxs->OPPPInnerServiceCablSec1Thick));
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec20, trcblsec20, 0);
  fInServFlatCbl0AssmblyZpos = trcblsec20.z() - lxs->OPPPInnerServiceCablSec20Thick/2.0;

  G4double trsec21x = trcblsec1.x() + 0.5*(lxs->OPPPInnerServiceCablSec21X + lxs->OPPPInnerServiceCablSec1X);
  G4double trsec21z = trcblsec20.z() + 0.5*(lxs->OPPPInnerServiceCablSec21Z + sec21dx + lxs->OPPPInnerServiceCablSec20Thick);
  G4ThreeVector trcblsec21(trsec21x, 0.0, trsec21z);
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec21, trcblsec21, rotcblsec0);

  G4double trsec22x = trcblsec1.x() + 0.5*(lxs->OPPPInnerServiceCablSec22X + lxs->OPPPInnerServiceCablSec1X);
  G4double trsec22z = trsec21z + 0.5*(lxs->OPPPInnerServiceCablSec22Z + sec22dx + sec21dx - lxs->OPPPInnerServiceCablSec21Z);
  G4ThreeVector trcblsec22(trsec22x, 0.0, trsec22z);
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec22, trcblsec22, rotcblsec0);

  G4double trsec31x = trsec21x + 0.5*(lxs->OPPPInnerServiceCablSec21X + sec31dx);
  G4double trsec31z = trsec21z + 0.5*(lxs->OPPPInnerServiceCablSec21Z + sec21dx - lxs->OPPPInnerServiceCablSec31Thick);
  G4ThreeVector trcblsec31(trsec31x, 0.0, trsec31z);
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec31, trcblsec31, 0);
  fInServFlatCbl1AssmblyZpos = trcblsec31.z() - lxs->OPPPInnerServiceCablSec31Thick/2.0;

  G4double trsec32x = trsec22x + 0.5*(lxs->OPPPInnerServiceCablSec22X + sec32dx);
  G4double trsec32z = trsec22z + 0.5*(lxs->OPPPInnerServiceCablSec22Z + sec22dx - lxs->OPPPInnerServiceCablSec32Thick);
  G4ThreeVector trcblsec32(trsec32x, 0.0, trsec32z);
  innerServiceCablesAssembly->AddPlacedVolume(logicInServFlatCblSec32, trcblsec32, 0);
  fInServFlatCbl2AssmblyZpos = trcblsec32.z() - lxs->OPPPInnerServiceCablSec32Thick/2.0;

// Assembling service components
  G4AssemblyVolume *innerServiceLinesAssembly = new G4AssemblyVolume();

  G4double dz = lxs->OPPPTrussZ;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;
  G4ThreeVector trcollpipeassembly(0.0, lxs->OPPPServiceCoolPipeYdist/2.0, ecdz);
//   innerServiceLinesAssembly->AddPlacedAssembly(coolingPipeAssembly, trcollpipeassembly, 0);
  LxAux::AddAssmblyVolumes(innerServiceLinesAssembly, coolingPipeAssembly, trcollpipeassembly, 0);
  trcollpipeassembly.setY(-lxs->OPPPServiceCoolPipeYdist/2.0);
//   innerServiceLinesAssembly->AddPlacedAssembly(coolingPipeAssembly, trcollpipeassembly, 0);
  LxAux::AddAssmblyVolumes(innerServiceLinesAssembly, coolingPipeAssembly, trcollpipeassembly, 0);

  G4ThreeVector trcablesec0assembly(0.0, 0.0, ecdz);
//   innerServiceLinesAssembly->AddPlacedAssembly(innerServiceCablesAssembly, trcablesec0assembly, 0);
  LxAux::AddAssmblyVolumes(innerServiceLinesAssembly, innerServiceCablesAssembly, trcablesec0assembly, 0);

  return innerServiceLinesAssembly;
}



G4AssemblyVolume* LxTrackerOPPP::ConstructInnerServiceLinesTermAssembly()
{
// Flat cable terminators
  G4AssemblyVolume *innerFlatCableTermAssembly = new G4AssemblyVolume();

  G4AssemblyVolume* cableTerm1Assembly = ConstructFlatCbleTerminator1Assembly();

  LXSetUp *lxs = LXSetUp::Instance();
  G4double dz = lxs->OPPPTrussZ;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;

  G4double trterm1x = lxs->OPPPDetTopPlateX - lxs->OPPPStaveX/2.0 - lxs->OPPPServiceCablTermBase1X;
  G4double trterm1y = -lxs->OPPPStaveLiftY - lxs->OPPPColdPlateY/2.0;
  G4double trterm1z = fInServFlatCbl0AssmblyZpos + lxs->OPPPInnerServiceCablSec20Thick + ecdz;
  G4ThreeVector trterm1(trterm1x, trterm1y, trterm1z);
//   innerFlatCableTermAssembly->AddPlacedAssembly(cableTerm1Assembly, trterm1, 0);
  LxAux::AddAssmblyVolumes(innerFlatCableTermAssembly, cableTerm1Assembly, trterm1, 0);

  return innerFlatCableTermAssembly;
}


G4AssemblyVolume* LxTrackerOPPP::ConstructOuterServiceLinesAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCoolantMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* coolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeMaterial);
  G4Material* cableMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablMaterial);

// Cooling pipes for outer stave
  G4double coolprin = lxs->OPPPServiceCoolPipeRout - lxs->OPPPServiceCoolPipeThick;
  G4double outStaveCoolPipeX = lxs->OPPPDetTopPlateX - lxs->OPPPStaveX - lxs->OPPPTrackerActiveX + lxs->OPPPStaveInOutX
                            - lxs->OPPPServiceCoolPipeCnctrL - lxs->OPPPServiceCoolPipeCnctrSL;
  G4double outStaveCoolPipeL0 = outStaveCoolPipeX;

  G4Tubs *solidOutServCoolPipe0 = new G4Tubs("solidOutServCoolPipe0", coolprin, lxs->OPPPServiceCoolPipeRout,
                                           outStaveCoolPipeL0/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOutServCoolPipe0 = new G4LogicalVolume(solidOutServCoolPipe0, coolPipeMaterial,
                                                  "logicOutServCoolPipe0");
  G4Tubs *solidOutServCoolPipe0W = new G4Tubs("solidOutServCoolPipe0W", 0.0, coolprin, outStaveCoolPipeL0/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOutServCoolPipe0W = new G4LogicalVolume(solidOutServCoolPipe0W, opppCoolantMaterial,
                                                  "logicOutServCoolPipe0W");
// Assembling cooling pipe with connector and water
  G4AssemblyVolume *coolingPipeConnectorAssembly = ConstructCoolingPipeConnectorAssembly();
  G4AssemblyVolume *outerCoolingPipeAssembly = new G4AssemblyVolume();

  G4ThreeVector trcntr(lxs->OPPPStaveX/2.0, 0.0, 0.0);
//   outerCoolingPipeAssembly->AddPlacedAssembly(coolingPipeConnectorAssembly, trcntr, 0);
  LxAux::AddAssmblyVolumes(outerCoolingPipeAssembly, coolingPipeConnectorAssembly, trcntr, 0);
  G4RotationMatrix *rotcntr = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  trcntr.setX(trcntr.x() + lxs->OPPPServiceCoolPipeCnctrL + lxs->OPPPServiceCoolPipeCnctrSL/2.0);
  trcntr.setX(trcntr.x() + (lxs->OPPPServiceCoolPipeCnctrSL + outStaveCoolPipeL0)/2.0);
  outerCoolingPipeAssembly->AddPlacedVolume(logicOutServCoolPipe0, trcntr, rotcntr);
  outerCoolingPipeAssembly->AddPlacedVolume(logicOutServCoolPipe0W, trcntr, rotcntr);

// Flat cbales assembly
  G4double flatcabledx = lxs->OPPPDetTopPlateX - 0.5*(lxs->OPPPStaveX + lxs->OPPPTrackerActiveX)
                         - lxs->OPPPTrackerActiveX + lxs->OPPPStaveInOutX;

  G4double thetsec0 = atan2(lxs->OPPPOuterServiceCablSec0Z, lxs->OPPPOuterServiceCablSec0X);
  G4double sec0dx = lxs->OPPPOuterServiceCablSec0Thick / cos(thetsec0);
  G4Para *solidOutServFlatCblSec0 = new G4Para("solidOutServFlatCblSec0", sec0dx/2.0, lxs->OPPPOuterServiceCablSec0Y/2.0,
                                              lxs->OPPPOuterServiceCablSec0X/2.0, 0.0, thetsec0, 0);
  G4LogicalVolume *logicOutServFlatCblSec0 = new G4LogicalVolume(solidOutServFlatCblSec0, cableMaterial,
                                                  "logicOutServFlatCblSec0");

  G4double sec1dx = flatcabledx - lxs->OPPPOuterServiceCablSec0X;
  G4Box *solidOutServFlatCblSec1 = new G4Box("solidOutServFlatCblSec1", sec1dx/2.0,
                                            lxs->OPPPOuterServiceCablSec1Y/2.0, lxs->OPPPOuterServiceCablSec1Thick/2.0);
  G4LogicalVolume *logicOutServFlatCblSec1 = new G4LogicalVolume(solidOutServFlatCblSec1, cableMaterial,
                                                  "logicOutServFlatCblSec1");

// Assembling differnet sections of the flat cables
  G4AssemblyVolume *outerServiceCablesAssembly = new G4AssemblyVolume();

  G4double trcblsec0x = 0.5*(lxs->OPPPTrackerActiveX + lxs->OPPPOuterServiceCablSec0X);
  G4double trcblsec0z = lxs->OPPPStaveSideZpos - lxs->OPPPStaveSideZ/2.0 - 0.5*(lxs->OPPPOuterServiceCablSec0Z + sec0dx);
  G4ThreeVector trcblsec0(trcblsec0x, 0.0, trcblsec0z);
  G4RotationMatrix *rotcblsec0 = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  outerServiceCablesAssembly->AddPlacedVolume(logicOutServFlatCblSec0, trcblsec0, rotcblsec0);
  G4ThreeVector trcblsec1(trcblsec0x + 0.5*(lxs->OPPPOuterServiceCablSec0X + sec1dx), 0.0,
                    trcblsec0z - 0.5*(lxs->OPPPOuterServiceCablSec0Z + sec0dx - lxs->OPPPOuterServiceCablSec1Thick));
  outerServiceCablesAssembly->AddPlacedVolume(logicOutServFlatCblSec1, trcblsec1, 0);
  fOutServFlatCblAssmblyZpos = trcblsec1.z() - lxs->OPPPOuterServiceCablSec1Thick/2.0;

// Assembling service components
  G4AssemblyVolume *outerServiceLinesAssembly = new G4AssemblyVolume();

  G4double dz = lxs->OPPPTrussZ;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;
  G4ThreeVector trcollpipeassembly(0.0, lxs->OPPPServiceCoolPipeYdist/2.0, ecdz);
//   outerServiceLinesAssembly->AddPlacedAssembly(outerCoolingPipeAssembly, trcollpipeassembly, 0);
  LxAux::AddAssmblyVolumes(outerServiceLinesAssembly, outerCoolingPipeAssembly, trcollpipeassembly, 0);
  trcollpipeassembly.setY(-lxs->OPPPServiceCoolPipeYdist/2.0);
//   outerServiceLinesAssembly->AddPlacedAssembly(outerCoolingPipeAssembly, trcollpipeassembly, 0);
  LxAux::AddAssmblyVolumes(outerServiceLinesAssembly, outerCoolingPipeAssembly, trcollpipeassembly, 0);

  G4ThreeVector trcablesec0assembly(0.0, 0.0, ecdz);
//   outerServiceLinesAssembly->AddPlacedAssembly(outerServiceCablesAssembly, trcablesec0assembly, 0);
  LxAux::AddAssmblyVolumes(outerServiceLinesAssembly, outerServiceCablesAssembly, trcablesec0assembly, 0);

  return outerServiceLinesAssembly;

}



void LxTrackerOPPP::ConstructSideServiceLines(G4LogicalVolume  *worldVol)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppServContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);

  G4double sservconx = lxs->OPPPServiceSupportX + lxs->OPPPServiceSupportTubeX;
  G4double sservcony = lxs->OPPPStaveLiftY + lxs->OPPPColdPlateY;
  G4Box *solidOPPPServSideContainer = new G4Box("solidOPPPServSideContainer", sservconx/2.0,
                                            sservcony/2.0, lxs->OPPPServiceSupportZ/2.0);
  //Cut container as it slightly overlaps with container of vacuum chamber
  G4double ssccutx = 0.2*sservconx;
  G4double ssccutz = 0.5*(lxs->OPPPServiceSupportZ - lxs->OPPPDetTopPlateZ - lxs->VacChambertoOPPPDetZGap);
  G4Box *solidOPPPServSideContainerCut = new G4Box("solidOPPPServSideContainerCut", ssccutx, sservcony, ssccutz);
  G4Transform3D trssccutl(G4RotationMatrix(), G4ThreeVector(-sservconx/2.0, 0.0, -lxs->OPPPServiceSupportZ/2.0));
  G4SubtractionSolid* solidOPPPServSideContainerL = new G4SubtractionSolid("solidOPPPServSideContainerL",
                         solidOPPPServSideContainer, solidOPPPServSideContainerCut, trssccutl);
  G4Transform3D trssccutr(G4RotationMatrix(), G4ThreeVector(sservconx/2.0, 0.0, -lxs->OPPPServiceSupportZ/2.0));
  G4SubtractionSolid* solidOPPPServSideContainerR = new G4SubtractionSolid("solidOPPPServSideContainerR",
                         solidOPPPServSideContainer, solidOPPPServSideContainerCut, trssccutr);

  G4LogicalVolume *logicOPPPServSideContainerL = new G4LogicalVolume(solidOPPPServSideContainerL,
                                                     opppServContainerMaterial, "logicOPPPServSideContainerL");
  G4LogicalVolume *logicOPPPServSideContainerR = new G4LogicalVolume(solidOPPPServSideContainerR,
                                                     opppServContainerMaterial, "logicOPPPServSideContainerR");

//   G4AssemblyVolume* sidesOutServAssembly = ConstructSideServiceLinesAssembly();
  G4AssemblyVolume* sidesOutServAssembly = 0;
  G4AssemblyVolume* sidesOutServCableTermAssemblyL = 0;
  G4AssemblyVolume* sidesOutServCableTermAssemblyR = 0;
  ConstructSideServiceLinesOutAssemblies(sidesOutServAssembly, sidesOutServCableTermAssemblyL, sidesOutServCableTermAssemblyR);

  G4AssemblyVolume* sidesInServAssembly = 0;
  G4AssemblyVolume* sidesInServCableTermAssemblyL = 0;
  G4AssemblyVolume* sidesInServCableTermAssemblyR = 0;
  ConstructSideServiceLinesInAssemblies(sidesInServAssembly, sidesInServCableTermAssemblyL, sidesInServCableTermAssemblyR);

  G4double opppdetmodz = lxs->OPPPDetOffsetZ;
  G4double staveextraz = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness + lxs->OPPPStaveSideZpos;
  G4double sasxpos = -sservconx/2.0;
  G4double sasypos = lxs->OPPPStaveLiftY/2.0;
  G4double saszpos = (lxs->OPPPTrussZ-lxs->OPPPDetTopPlateZ)/2.0 + staveextraz + opppdetmodz;
  G4ThreeVector trmin(sasxpos, sasypos, saszpos + lxs->OPPPStaveInOutZ);
  G4ThreeVector trmout(sasxpos, sasypos, saszpos);

  for (G4int il = 0; il < lxs->OPPPTrackerNLayers; ++il) {
    sidesInServAssembly->MakeImprint(logicOPPPServSideContainerL, trmin, 0, 0, lxs->OverlapTest);
    sidesInServCableTermAssemblyL->MakeImprint(logicOPPPServSideContainerL, trmin, 0, 0, lxs->OverlapTest);

    sidesOutServAssembly->MakeImprint(logicOPPPServSideContainerL, trmout, 0, 0, lxs->OverlapTest);
    sidesOutServCableTermAssemblyL->MakeImprint(logicOPPPServSideContainerL, trmout, 0, 0, lxs->OverlapTest);

    G4ThreeVector trminr(-trmin.x(), trmin.y(), trmin.z());
    G4ThreeVector trmoutr(-trmout.x(), trmout.y(), trmout.z());

    sidesInServAssembly->MakeImprint(logicOPPPServSideContainerR, trminr,
                          new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), 0, lxs->OverlapTest);
    sidesInServCableTermAssemblyR->MakeImprint(logicOPPPServSideContainerR, trminr, 0, 0, lxs->OverlapTest);

    sidesOutServAssembly->MakeImprint(logicOPPPServSideContainerR, trmoutr,
                          new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), 0, lxs->OverlapTest);
    sidesOutServCableTermAssemblyR->MakeImprint(logicOPPPServSideContainerR, trmoutr, 0, 0, lxs->OverlapTest);

    trmin.setZ(trmin.z() + lxs->OPPPTrackerInterLayerZ);
    trmout.setZ(trmout.z() + lxs->OPPPTrackerInterLayerZ);
  }

  G4double detxpos = lxs->OPPPDetTopPlateX/2.0 + lxs->OPPPDetXPos;
  G4double detypos = -lxs->OPPPStaveLiftY/2.0;
  G4double detzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet
                     + lxs->OPPPDetTopPlateZ/2.0 - lxs->OPPPDetOffsetZ;
  G4double scposx = detxpos + 0.5*(lxs->OPPPDetTopPlateX + sservconx);

  new G4PVPlacement (0, G4ThreeVector(scposx, detypos, detzpos), logicOPPPServSideContainerL,
                      "OPPPServSideContainerL", worldVol, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-scposx, detypos-lxs->OPPPDetYLowPos, detzpos), logicOPPPServSideContainerR,
                      "OPPPServSideContainerR", worldVol, false, 0, lxs->OverlapTest);
}



void LxTrackerOPPP::ConstructSideServiceLinesInAssemblies(G4AssemblyVolume* &sideServiceAssembly,
                       G4AssemblyVolume* &sideFlatCableTermAssemblyL, G4AssemblyVolume* &sideFlatCableTermAssemblyR)
{
  if (sideServiceAssembly) delete sideServiceAssembly;
  if (sideFlatCableTermAssemblyL) delete sideFlatCableTermAssemblyL;
  if (sideFlatCableTermAssemblyR) delete sideFlatCableTermAssemblyR;

  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCoolantMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* coolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeMaterial);
  G4Material* cableMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablMaterial);

// Cooling pipe with coolant
  G4double coolprin = lxs->OPPPServiceCoolPipeRout - lxs->OPPPServiceCoolPipeThick;
  G4double sideInCoolPipeL = lxs->OPPPServiceSupportX + lxs->OPPPServiceSupportTubeX;
  G4Tubs *solidSideInServCoolPipe = new G4Tubs("solidSideInServCoolPipe", 0.0, lxs->OPPPServiceCoolPipeRout,
                                           sideInCoolPipeL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSideInServCoolPipe = new G4LogicalVolume(solidSideInServCoolPipe, coolPipeMaterial,
                                                  "logicSideInServCoolPipe");
  G4Tubs *solidSideInServCoolPipeW = new G4Tubs("solidSideInServCoolPipeW", 0.0, coolprin,
                                                 sideInCoolPipeL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSideInServCoolPipeW = new G4LogicalVolume(solidSideInServCoolPipeW, opppCoolantMaterial,
                                                  "logicSideInServCoolPipeW");
  new G4PVPlacement (0, G4ThreeVector(), logicSideInServCoolPipeW,
                      "SideInServCoolPipeCoolant", logicSideInServCoolPipe, false, 0, lxs->OverlapTest);
// Flat cbales assembly
  G4Box *solidSdInServFlatCblSec01 = new G4Box("solidSdInServFlatCblSec01", lxs->OPPPSideInServiceCablSec01X/2.0,
        lxs->OPPPSideInServiceCablSec01Y/2.0, lxs->OPPPSideInServiceCablSec01Thick/2.0);
  G4LogicalVolume *logicSdInServFlatCblSec01 = new G4LogicalVolume(solidSdInServFlatCblSec01, cableMaterial,
                                                  "logicSdInServFlatCblSec01");

  G4Box *solidSdInServFlatCblSec02 = new G4Box("solidSdInServFlatCblSec02", lxs->OPPPSideInServiceCablSec02X/2.0,
        lxs->OPPPSideInServiceCablSec02Y/2.0, lxs->OPPPSideInServiceCablSec02Thick/2.0);
  G4LogicalVolume *logicSdInServFlatCblSec02 = new G4LogicalVolume(solidSdInServFlatCblSec02, cableMaterial,
                                                  "logicSdInServFlatCblSec02");
// Assembling differnet sections of the flat cables
  G4AssemblyVolume *outerSideServiceCablesAssembly = new G4AssemblyVolume();

  G4double trsec01x = 0.5 * lxs->OPPPSideInServiceCablSec01X;
  G4double trsec01z = fInServFlatCbl1AssmblyZpos + 0.5*lxs->OPPPSideInServiceCablSec01Thick;
  G4ThreeVector trcblsec01(trsec01x, 0.0, trsec01z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdInServFlatCblSec01, trcblsec01, 0);

  G4double trsec02x = 0.5 * lxs->OPPPSideInServiceCablSec02X;
  G4double trsec02z = fInServFlatCbl2AssmblyZpos + 0.5*lxs->OPPPSideInServiceCablSec02Thick;
  G4ThreeVector trcblsec02(trsec02x, 0.0, trsec02z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdInServFlatCblSec02, trcblsec02, 0);

// Assembling service components
  sideServiceAssembly = new G4AssemblyVolume();

  G4double dz = lxs->OPPPTrussZ;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;

// Assembling cooling pipes
  G4ThreeVector trcpipe(sideInCoolPipeL/2.0, lxs->OPPPServiceCoolPipeYdist/2.0, fInServCoolPipeAssmblyZpos + ecdz);
  G4RotationMatrix *rotcpipe = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  sideServiceAssembly->AddPlacedVolume(logicSideInServCoolPipe, trcpipe, rotcpipe);
  trcpipe.setY(-lxs->OPPPServiceCoolPipeYdist/2.0);
  sideServiceAssembly->AddPlacedVolume(logicSideInServCoolPipe, trcpipe, rotcpipe);

// Assembling flat cables
  G4ThreeVector trcablesassembly(0.0, 0.0, ecdz);
//   sideServiceAssembly->AddPlacedAssembly(outerSideServiceCablesAssembly, trcablesassembly, 0);
  LxAux::AddAssmblyVolumes(sideServiceAssembly, outerSideServiceCablesAssembly, trcablesassembly, 0);

// Flat cables terminators
  sideFlatCableTermAssemblyL = new G4AssemblyVolume();
  sideFlatCableTermAssemblyR = new G4AssemblyVolume();

  G4AssemblyVolume* cableTerm2Assembly = ConstructFlatCbleTerminator2Assembly();
  G4double trterm2x = trsec01x + lxs->OPPPSideInServiceCablSec01X/2.0;
  G4double trterm2y = -lxs->OPPPStaveLiftY - lxs->OPPPColdPlateY/2.0;
  G4double trterm2z = (trsec01z + trsec02z)/2.0 + ecdz;
  G4ThreeVector trterm2l(trterm2x, trterm2y, trterm2z);
//   sideFlatCableTermAssemblyL->AddPlacedAssembly(cableTerm2Assembly, trterm2l, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyL, cableTerm2Assembly, trterm2l, 0);
  G4ThreeVector trterm2r(-trterm2x - lxs->OPPPServiceCablTermBase2X, trterm2y, trterm2z);
//   sideFlatCableTermAssemblyR->AddPlacedAssembly(cableTerm2Assembly, trterm2r, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyR, cableTerm2Assembly, trterm2r, 0);
}



void LxTrackerOPPP::ConstructSideServiceLinesOutAssemblies(G4AssemblyVolume* &sideServiceAssembly,
                       G4AssemblyVolume* &sideFlatCableTermAssemblyL, G4AssemblyVolume* &sideFlatCableTermAssemblyR)
{
  if (sideServiceAssembly) delete sideServiceAssembly;
  if (sideFlatCableTermAssemblyL) delete sideFlatCableTermAssemblyL;
  if (sideFlatCableTermAssemblyR) delete sideFlatCableTermAssemblyR;

  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCoolantMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCoolingMaterial);
  G4Material* coolPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCoolPipeMaterial);
  G4Material* cableMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablMaterial);

// Cooling pipe with coolant
  G4double coolprin = lxs->OPPPServiceCoolPipeRout - lxs->OPPPServiceCoolPipeThick;
  G4double sideInCoolPipeL = lxs->OPPPServiceSupportX + lxs->OPPPServiceSupportTubeX;
  G4Tubs *solidSideServCoolPipe = new G4Tubs("solidSideServCoolPipe", 0.0, lxs->OPPPServiceCoolPipeRout,
                                           sideInCoolPipeL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSideServCoolPipe = new G4LogicalVolume(solidSideServCoolPipe, coolPipeMaterial,
                                                  "logicSideServCoolPipe");
  G4Tubs *solidSideServCoolPipeW = new G4Tubs("solidSideServCoolPipeW", 0.0, coolprin, sideInCoolPipeL/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSideServCoolPipeW = new G4LogicalVolume(solidSideServCoolPipeW, opppCoolantMaterial,
                                                  "logicSideServCoolPipeW");
  new G4PVPlacement (0, G4ThreeVector(), logicSideServCoolPipeW,
                      "SideServCoolPipeCoolant", logicSideServCoolPipe, false, 0, lxs->OverlapTest);

// Flat cbales assembly
  G4double thetsec0 = atan2(lxs->OPPPSideOutServiceCablSec0Z, lxs->OPPPSideOutServiceCablSec0X);
  G4double sec0dx = lxs->OPPPSideOutServiceCablSec0Thick / cos(thetsec0);
  G4Para *solidSdOutServFlatCblSec0 = new G4Para("solidSdOutServFlatCblSec0", sec0dx/2.0,
                    lxs->OPPPSideOutServiceCablSec0Y/2.0, lxs->OPPPSideOutServiceCablSec0X/2.0, 0.0, thetsec0, 0.0);
  G4LogicalVolume *logicSdOutServFlatCblSec0 = new G4LogicalVolume(solidSdOutServFlatCblSec0, cableMaterial,
                                                  "logicSdOutServFlatCblSec0");

  G4Box *solidSdOutServFlatCblSec1 = new G4Box("solidSdOutServFlatCblSec1", lxs->OPPPSideOutServiceCablSec1X/2.0,
                        lxs->OPPPSideOutServiceCablSec1Y/2.0, lxs->OPPPSideOutServiceCablSec1Thick/2.0);
  G4LogicalVolume *logicSdOutServFlatCblSec1 = new G4LogicalVolume(solidSdOutServFlatCblSec1, cableMaterial,
                                                  "logicSdOutServFlatCblSec1");

  G4Box *solidSdOutServFlatCblSec20 = new G4Box("solidSdOutServFlatCblSec20", lxs->OPPPSideOutServiceCablSec20X/2.0,
                            lxs->OPPPSideOutServiceCablSec20Y/2.0, lxs->OPPPSideOutServiceCablSec20Thick/2.0);
  G4LogicalVolume *logicSdOutServFlatCblSec20 = new G4LogicalVolume(solidSdOutServFlatCblSec20, cableMaterial,
                                                  "logicSdOutServFlatCblSec20");

  G4double thetsec21 = atan2(lxs->OPPPSideOutServiceCablSec21Z, lxs->OPPPSideOutServiceCablSec21X);
  G4double sec21dx = lxs->OPPPSideOutServiceCablSec21Thick / cos(thetsec21);
  G4Para *solidSdOutServFlatCblSec21 = new G4Para("solidSdOutServFlatCblSec21", sec21dx/2.0,
                    lxs->OPPPSideOutServiceCablSec21Y/2.0, lxs->OPPPSideOutServiceCablSec21X/2.0, 0.0, thetsec21, 0);
  G4LogicalVolume *logicSdOutServFlatCblSec21 = new G4LogicalVolume(solidSdOutServFlatCblSec21, cableMaterial,
                                                  "logicSdOutServFlatCblSec21");

  G4double thetsec22 = atan2(lxs->OPPPSideOutServiceCablSec22Z, lxs->OPPPSideOutServiceCablSec22X);
  G4double sec22dx = lxs->OPPPSideOutServiceCablSec22Thick / cos(thetsec22);
  G4Para *solidSdOutServFlatCblSec22 = new G4Para("solidSdOutServFlatCblSec22", sec22dx/2.0,
           lxs->OPPPSideOutServiceCablSec22Y/2.0, lxs->OPPPSideOutServiceCablSec22X/2.0, 0.0, thetsec22, 0.0);
  G4LogicalVolume *logicSdOutServFlatCblSec22 = new G4LogicalVolume(solidSdOutServFlatCblSec22, cableMaterial,
                                                  "logicSdOutServFlatCblSec22");

  G4Box *solidSdOutServFlatCblSec31 = new G4Box("solidSdOutServFlatCblSec31", lxs->OPPPSideOutServiceCablSec31X/2.0,
        lxs->OPPPSideOutServiceCablSec31Y/2.0, lxs->OPPPSideOutServiceCablSec31Thick/2.0);
  G4LogicalVolume *logicSdOutServFlatCblSec31 = new G4LogicalVolume(solidSdOutServFlatCblSec31, cableMaterial,
                                                  "logicSdOutServFlatCblSec31");

  G4Box *solidSdOutServFlatCblSec32 = new G4Box("solidSdOutServFlatCblSec32", lxs->OPPPSideOutServiceCablSec32X/2.0,
        lxs->OPPPSideOutServiceCablSec32Y/2.0, lxs->OPPPSideOutServiceCablSec32Thick/2.0);
  G4LogicalVolume *logicSdOutServFlatCblSec32 = new G4LogicalVolume(solidSdOutServFlatCblSec32, cableMaterial,
                                                  "logicSdOutServFlatCblSec32");
// Assembling differnet sections of the flat cables
  G4AssemblyVolume *outerSideServiceCablesAssembly = new G4AssemblyVolume();

  G4double trcblsec0x = 0.5*lxs->OPPPSideOutServiceCablSec0X;
  G4double trcblsec0z = fOutServFlatCblAssmblyZpos + 0.5*(sec0dx - lxs->OPPPSideOutServiceCablSec0Z);
  G4ThreeVector trcblsec0(trcblsec0x, 0.0, trcblsec0z);
  G4RotationMatrix *rotcblsec0 = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec0, trcblsec0, rotcblsec0);

  G4ThreeVector trcblsec1(trcblsec0x + 0.5*(lxs->OPPPSideOutServiceCablSec0X+lxs->OPPPSideOutServiceCablSec1X), 0.0,
               trcblsec0z - 0.5*(lxs->OPPPSideOutServiceCablSec0Z + sec0dx - lxs->OPPPSideOutServiceCablSec1Thick));
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec1, trcblsec1, 0);

  G4double trcblsec20x = trcblsec1.x() + 0.5*(lxs->OPPPSideOutServiceCablSec1X + lxs->OPPPSideOutServiceCablSec20X);
  G4double trcblsec20z = trcblsec1.z() - 0.5*(lxs->OPPPSideOutServiceCablSec20Thick - lxs->OPPPSideOutServiceCablSec1Thick);
  G4ThreeVector trcblsec20(trcblsec20x, 0.0, trcblsec20z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec20, trcblsec20, 0);

  G4double trsec21x = trcblsec1.x() + 0.5*(lxs->OPPPSideOutServiceCablSec1X + lxs->OPPPSideOutServiceCablSec21X);
  G4double trsec21z = trcblsec20.z() - 0.5*(lxs->OPPPSideOutServiceCablSec21Z + sec21dx
                    + lxs->OPPPSideOutServiceCablSec20Thick);
  G4ThreeVector trcblsec21(trsec21x, 0.0, trsec21z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec21, trcblsec21, rotcblsec0);

  G4double trsec22x = trcblsec1.x() + 0.5*(lxs->OPPPSideOutServiceCablSec1X + lxs->OPPPSideOutServiceCablSec22X);
  G4double trsec22z = trsec21z - 0.5*(lxs->OPPPSideOutServiceCablSec22Z + sec22dx + sec21dx
                                       - lxs->OPPPSideOutServiceCablSec21Z);
  G4ThreeVector trcblsec22(trsec22x, 0.0, trsec22z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec22, trcblsec22, rotcblsec0);

  G4double trsec31x = trsec21x + 0.5*(lxs->OPPPSideOutServiceCablSec21X + lxs->OPPPSideOutServiceCablSec31X);
  G4double trsec31z = trsec21z - 0.5*(lxs->OPPPSideOutServiceCablSec21Z + sec21dx - lxs->OPPPSideOutServiceCablSec31Thick);
  G4ThreeVector trcblsec31(trsec31x, 0.0, trsec31z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec31, trcblsec31, 0);

  G4double trsec32x = trsec21x + 0.5*(lxs->OPPPSideOutServiceCablSec22X + lxs->OPPPSideOutServiceCablSec32X);
  G4double trsec32z = trsec22z - 0.5*(lxs->OPPPSideOutServiceCablSec22Z + sec22dx - lxs->OPPPSideOutServiceCablSec32Thick);
  G4ThreeVector trcblsec32(trsec32x, 0.0, trsec32z);
  outerSideServiceCablesAssembly->AddPlacedVolume(logicSdOutServFlatCblSec32, trcblsec32, 0);

// Assembling service components
  sideServiceAssembly = new G4AssemblyVolume();

  G4double dz = lxs->OPPPTrussZ;
  G4double stavesidez = lxs->OPPPStaveSideZ;
  G4double dzcplate = lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness;
  G4double stavesidezpos = lxs->OPPPStaveSideZpos;
  G4double ecdz = 0.5*(stavesidez-dz) - dzcplate - stavesidezpos;

// Assembling cooling pipes
  G4ThreeVector trcpipe(sideInCoolPipeL/2.0, lxs->OPPPServiceCoolPipeYdist/2.0, ecdz);
  G4RotationMatrix *rotcpipe = new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0);
  sideServiceAssembly->AddPlacedVolume(logicSideServCoolPipe, trcpipe, rotcpipe);
  trcpipe.setY(-lxs->OPPPServiceCoolPipeYdist/2.0);
  sideServiceAssembly->AddPlacedVolume(logicSideServCoolPipe, trcpipe, rotcpipe);

// Assembling flat cables
  G4ThreeVector trcablesassembly(0.0, 0.0, ecdz);
//   sideServiceAssembly->AddPlacedAssembly(outerSideServiceCablesAssembly, trcablesassembly, 0);
  LxAux::AddAssmblyVolumes(sideServiceAssembly, outerSideServiceCablesAssembly, trcablesassembly, 0);

// Flat cables terminators
  sideFlatCableTermAssemblyL = new G4AssemblyVolume();
  sideFlatCableTermAssemblyR = new G4AssemblyVolume();

  G4AssemblyVolume* cableTerm1Assembly = ConstructFlatCbleTerminator1Assembly();
  G4double trterm1x = trsec32x + lxs->OPPPSideOutServiceCablSec32X/2.0 - lxs->OPPPServiceCablTermBase1X;
  G4double trterm1y = -lxs->OPPPStaveLiftY - lxs->OPPPColdPlateY/2.0;
  G4double trterm1z = trsec32z + lxs->OPPPSideOutServiceCablSec32Thick/2.0 + ecdz;
  G4ThreeVector trterm1l(trterm1x, trterm1y, trterm1z);
//   sideFlatCableTermAssemblyL->AddPlacedAssembly(cableTerm1Assembly, trterm1l, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyL, cableTerm1Assembly, trterm1l, 0);
  G4ThreeVector trterm1r(-trterm1x - lxs->OPPPServiceCablTermBase1X, trterm1y, trterm1z);
//   sideFlatCableTermAssemblyR->AddPlacedAssembly(cableTerm1Assembly, trterm1r, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyR, cableTerm1Assembly, trterm1r, 0);

  G4AssemblyVolume* cableTerm2Assembly = ConstructFlatCbleTerminator2Assembly();
  G4double trterm2x = trsec31x + lxs->OPPPSideOutServiceCablSec31X/2.0;
  G4double trterm2y = -lxs->OPPPStaveLiftY - lxs->OPPPColdPlateY/2.0;
  G4double trterm2z = (trsec31z + trcblsec20z)/2.0 + ecdz;
  G4ThreeVector trterm2l(trterm2x, trterm2y, trterm2z);
//   sideFlatCableTermAssemblyL->AddPlacedAssembly(cableTerm2Assembly, trterm2l, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyL, cableTerm2Assembly, trterm2l, 0);
  G4ThreeVector trterm2r(-trterm2x - lxs->OPPPServiceCablTermBase2X, trterm2y, trterm2z);
//   sideFlatCableTermAssemblyR->AddPlacedAssembly(cableTerm2Assembly, trterm2r, 0);
  LxAux::AddAssmblyVolumes(sideFlatCableTermAssemblyR, cableTerm2Assembly, trterm2r, 0);

}



G4AssemblyVolume* LxTrackerOPPP::ConstructFlatCbleTerminator1Assembly()
{
  if (fCableTerm1Assembly) return fCableTerm1Assembly;

  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCTermMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablTermMaterial);

  G4double ctbasey = lxs->OPPPStaveLiftY;
  G4Box *solidOPPPCableTermBase1 = new G4Box("solidOPPPCableTermBase1", lxs->OPPPServiceCablTermBase1X/2.0,
                                                            ctbasey/2.0, lxs->OPPPServiceCablTermBase1Z/2.0);
  G4LogicalVolume *logicOPPPCableTermBase1 = new G4LogicalVolume(solidOPPPCableTermBase1, opppCTermMaterial,
                                                  "logicOPPPCableTermBase1");
  G4Box *solidOPPPCableTermMid1 = new G4Box("solidOPPPCableTermMid1", lxs->OPPPServiceCablTermMid1X/2.0,
                                     lxs->OPPPServiceCablTermMid1Y/2.0, lxs->OPPPServiceCablTermMid1Z/2.0);
  G4LogicalVolume *logicOPPPCableTermMid1 = new G4LogicalVolume(solidOPPPCableTermMid1, opppCTermMaterial,
                                                  "logicOPPPCableTermMid1");

  G4AssemblyVolume *serviceCablTermAssembly = new G4AssemblyVolume();
  G4ThreeVector trbase(lxs->OPPPServiceCablTermBase1X/2.0, ctbasey/2.0,
                             lxs->OPPPServiceCablTermMid1Z - lxs->OPPPServiceCablTermBase1Z/2.0);
  serviceCablTermAssembly->AddPlacedVolume(logicOPPPCableTermBase1, trbase, 0);
  trbase.setY(trbase.y() + 0.5*(ctbasey + lxs->OPPPServiceCablTermMid1Y));
  trbase.setZ(lxs->OPPPServiceCablTermMid1Z/2.0);
  serviceCablTermAssembly->AddPlacedVolume(logicOPPPCableTermMid1, trbase, 0);

  fCableTerm1Assembly = serviceCablTermAssembly;
  return serviceCablTermAssembly;
}



G4AssemblyVolume* LxTrackerOPPP::ConstructFlatCbleTerminator2Assembly()
{
  if (fCableTerm2Assembly) return fCableTerm2Assembly;
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppCTermMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPServiceCablTermMaterial);

  G4double ctbasey = lxs->OPPPStaveLiftY;
  G4Box *solidOPPPCableTermBase2 = new G4Box("solidOPPPCableTermBase2", lxs->OPPPServiceCablTermBase2X/2.0,
                                                            ctbasey/2.0, lxs->OPPPServiceCablTermBase2Z/2.0);
  G4LogicalVolume *logicOPPPCableTermBase2 = new G4LogicalVolume(solidOPPPCableTermBase2, opppCTermMaterial,
                                                  "logicOPPPCableTermBase2");
  G4Box *solidOPPPCableTermMid2 = new G4Box("solidOPPPCableTermMid2", lxs->OPPPServiceCablTermMid2X/2.0,
                                     lxs->OPPPServiceCablTermMid2Y/2.0, lxs->OPPPServiceCablTermMid2Z/2.0);
  G4LogicalVolume *logicOPPPCableTermMid2 = new G4LogicalVolume(solidOPPPCableTermMid2, opppCTermMaterial,
                                                  "logicOPPPCableTermMid2");

  G4AssemblyVolume *serviceCablTermAssembly = new G4AssemblyVolume();
  G4ThreeVector trbase(lxs->OPPPServiceCablTermBase2X/2.0, ctbasey/2.0, 0.0);
  serviceCablTermAssembly->AddPlacedVolume(logicOPPPCableTermBase2, trbase, 0);
  trbase.setY(trbase.y() + 0.5*(ctbasey + lxs->OPPPServiceCablTermMid2Y));
  serviceCablTermAssembly->AddPlacedVolume(logicOPPPCableTermMid2, trbase, 0);

  fCableTerm2Assembly = serviceCablTermAssembly;
  return serviceCablTermAssembly;
}



void LxTrackerOPPP::ConstructElectronicsRack()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

// Electroncis rack
  G4Material* trackerRackContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* trackerRackMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TrackerElectronicsRackMaterial);
  G4Material* trackerPCBMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Tracker_FR4");

  G4Box *solidTrackerElectronicsRackContainer = new G4Box("solidTrackerElectronicsRackContainer", lxs->TrackerElectronicsRackX/2.0,
                                            lxs->TrackerElectronicsRackY/2.0, lxs->TrackerElectronicsRackZ/2.0);
  G4LogicalVolume *logicTrackerElectronicsRackContainer = new G4LogicalVolume(solidTrackerElectronicsRackContainer,
                                          trackerRackContainerMaterial, "logicTrackerElectronicsRackContainer");
  //Rack box with a box cut for internal components
  std::vector<G4double> racksize{lxs->TrackerElectronicsRackX, lxs->TrackerElectronicsRackY, lxs->TrackerElectronicsRackZ};
  G4Box *solidTrackerElectronicsRack1 = new G4Box("solidTrackerElectronicsRack1", racksize[0]/2.0,
                                            racksize[1]/2.0, racksize[2]/2.0);
  std::for_each(racksize.begin(), racksize.end(), [&](G4double &x) {x -= 2.0*lxs->TrackerElectronicsRackThick; });
  G4Box *solidTrackerElectronicsRackCut = new G4Box("solidTrackerElectronicsRackCut", racksize[0]/2.0,
                                            racksize[1]/2.0, racksize[2]/2.0);
  G4Transform3D transform1(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid* solidTrackerElectronicsRack = new G4SubtractionSolid("solidTrackerElectronicsRack",
                                                        solidTrackerElectronicsRack1, solidTrackerElectronicsRackCut, transform1);

  G4LogicalVolume *logicTrackerElectronicsRack = new G4LogicalVolume(solidTrackerElectronicsRack, trackerRackMaterial,
                                                                     "logicTrackerElectronicsRack");
  //Rack internal PCBs
  G4Box *solidXYPCB = new G4Box("solidTrackerElectronicsXYPCB", racksize[0]/2.0, racksize[1]/2.0, lxs->TrackerElectronicsPCBThick/2.0);
  G4LogicalVolume *logicXYPCB = new G4LogicalVolume(solidXYPCB, trackerPCBMaterial, "logicTrackerElectronicsXYPCB");

  G4double yzpcbz = 0.5*(racksize[2] - lxs->TrackerElectronicsPCBThick);
  G4Box *solidYZPCB = new G4Box("solidTrackerElectronicsYZPCB", lxs->TrackerElectronicsPCBThick/2.0, racksize[1]/2.0, yzpcbz/2.0);
  G4LogicalVolume *logicYZPCB = new G4LogicalVolume(solidYZPCB, trackerPCBMaterial, "logicTrackerElectronicsYZPCB");

  G4double xzpcbx = 0.5*(racksize[0] - lxs->TrackerElectronicsPCBThick);
  G4Box *solidXZPCB = new G4Box("solidTrackerElectronicsXZPCB", xzpcbx/2.0, lxs->TrackerElectronicsPCBThick/2.0, yzpcbz/2.0);
  G4LogicalVolume *logicXZPCB = new G4LogicalVolume(solidXZPCB, trackerPCBMaterial, "logicTrackerElectronicsXZPCB");

  new G4PVPlacement(0, G4ThreeVector(), logicTrackerElectronicsRack, "TrackerElectronicsRackWalls",
                    logicTrackerElectronicsRackContainer, false, 0, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector(), logicXYPCB, "TrackerElectronicsXYPCB",
                    logicTrackerElectronicsRackContainer, false, 0, lxs->OverlapTest);

  G4double yzpcbzpos = 0.5 * (yzpcbz + lxs->TrackerElectronicsPCBThick);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, yzpcbzpos), logicYZPCB, "TrackerElectronicsYZPCB",
                    logicTrackerElectronicsRackContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -yzpcbzpos), logicYZPCB, "TrackerElectronicsYZPCB",
                    logicTrackerElectronicsRackContainer, false, 1, lxs->OverlapTest);

  G4double xzpcbxpos = 0.5 * (xzpcbx + lxs->TrackerElectronicsPCBThick);
  std::vector<G4double> xsign{1.0, 1.0, -1.0, -1.0};
  std::vector<G4double> zsign{1.0, -1.0, 1.0, -1.0};
  for (G4int ii = 0; ii < 4; ++ii) {
    new G4PVPlacement(0, G4ThreeVector(xzpcbxpos * xsign[ii], 0.0, yzpcbzpos) * zsign[ii], logicXZPCB, "TrackerElectronicsXZPCB",
                      logicTrackerElectronicsRackContainer, false, ii, lxs->OverlapTest);
  }

  G4double rackxpos = lxs->OPPPBasePlateX/2.0 + lxs->TrackerElectronicsRackXpos + lxs->TrackerElectronicsRackX/2.0;
  G4double opppdetmodz = lxs->OPPPDetOffsetZ;
  G4double detzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet
                     + lxs->OPPPDetTopPlateZ/2.0 - opppdetmodz;
  G4double basezpos = detzpos + (lxs->OPPPDetBottomSupportZ - lxs->OPPPDetTopPlateZ)/2.0;

  new G4PVPlacement(0, G4ThreeVector(rackxpos, 0.5*lxs->TrackerElectronicsRackY + lxs->FloorSurfaceYpos, basezpos),
                    logicTrackerElectronicsRackContainer, "TrackerElectronicsRack", fLogicWorld, false, 0, lxs->OverlapTest);

}



void LxTrackerOPPP::CreateMaterial()
{
  G4double density, fractionmass, z, a;
  G4int ncomponents, natoms;
  G4String symbol;

  LXSetUp *lxs = LXSetUp::Instance();

  G4NistManager* mman = G4NistManager::Instance();

  G4Element* H = mman->FindOrBuildElement("H");
  G4Element* C = mman->FindOrBuildElement("C");
  G4Element* O = mman->FindOrBuildElement("O");
  G4Material* Epoxy = new G4Material("EpoxyTracker", density= 1.3*g/cm3, ncomponents=3);
  Epoxy->AddElement(H, fractionmass=0.1310);
  Epoxy->AddElement(C, fractionmass=0.5357);
  Epoxy->AddElement(O, fractionmass=0.3333);

  G4Material* Cmat = mman->FindOrBuildMaterial("G4_C");
  G4Material* CarbonFiber = new G4Material("CarbonFiber", 1.6*g/cm3, ncomponents = 2);
  CarbonFiber->AddMaterial(Cmat, fractionmass=50.0*perCent);
  CarbonFiber->AddMaterial(Epoxy, fractionmass=50.0*perCent);

  G4double cpmtot = CarbonFiber->GetDensity()*lxs->OPPPColdPlateCarbonFiberThikness
                + Epoxy->GetDensity()*lxs->OPPPColdPlateEpoxyThikness;
  fractionmass = Epoxy->GetDensity()*lxs->OPPPColdPlateEpoxyThikness / cpmtot;

  density = cpmtot/(lxs->OPPPColdPlateEpoxyThikness + lxs->OPPPColdPlateCarbonFiberThikness);
  G4Material* ColdPlateMaterial = new G4Material(lxs->OPPPColdPlateMaterial, density, ncomponents = 2);
  ColdPlateMaterial->AddMaterial(Epoxy, fractionmass);
  ColdPlateMaterial->AddMaterial(CarbonFiber, 1.0-fractionmass);

  // C19_H12_O3 Polyetheretherketon (PEEK)
  G4Material* peek = new G4Material("PEEK", density=1.32*g/cm3, ncomponents=3);
  peek->AddElement(C, natoms=19);
  peek->AddElement(H, natoms=12);
  peek->AddElement(O, natoms=3);

  G4Material* staveEndCapMaterial = new G4Material(lxs->OPPPStaveSideMaterial, density=1.4*g/cm3, ncomponents=2);
  staveEndCapMaterial->AddMaterial(peek, fractionmass=70.0*perCent);
  staveEndCapMaterial->AddMaterial(Cmat, fractionmass=30.0*perCent);

  // fiber glass
  G4Element* Si = new G4Element("Silcone"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material* fiberglass = new G4Material("TrackerFiberglass", density=2.61*g/cm3, ncomponents=2);
  fiberglass->AddElement(Si, natoms=1);
  fiberglass->AddElement(O, natoms=2);
  // PCBoard material FR4
  G4Material* Copper = mman->FindOrBuildMaterial("G4_Cu");
  G4Material *FR4 = new G4Material("Tracker_FR4",density=1.93*g/cm3,ncomponents=3);
  FR4->AddMaterial( Epoxy, fractionmass=0.34);
  FR4->AddMaterial( fiberglass, fractionmass=0.56);
  FR4->AddMaterial( Copper, fractionmass=0.1);
}


