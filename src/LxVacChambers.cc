//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

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
#include "LxVacChambers.hh"


//////////////////////////////////////////////////////////////////////////////////////////////
// PipeChamberAssembly

PipeChamberAssembly::PipeChamberAssembly(const G4String mtypename) :
  fChamberType(mtypename), fChamberAssembly(0) 
{
  ConstructVacuumChamber();  
}



void PipeChamberAssembly::ConstructVacuumChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* bpipeWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeWindowMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  G4double bpipel = lxs->OPPPDetZtoMagnet;
//   G4double detxpos = lxs->OPPPDetX/2.0 + 4.0*lxs->BPipeR;
//   G4double dumpMagnetZ = lxs->TypMBFieldLength;
//   G4double thetad = atan2(detxpos, bpipel);

  G4Tubs *solidBPipeDContainer = new G4Tubs("solidBPipeDContainer", 0.0, lxs->BPipeR, bpipel/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBPipeDContainer = new G4LogicalVolume(solidBPipeDContainer, vacuumMaterial, "logicBPipeDContainer");

  // Beam pipe top and bottom part
  G4Tubs *solidBPipeD = new G4Tubs("solidBPipeD", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                    bpipel/2.0, lxs->DMBPipeWindow/2.0, M_PI-lxs->DMBPipeWindow);
  G4LogicalVolume *logicBPipeD = new G4LogicalVolume(solidBPipeD, beamPipeMaterial, "logicBPipeD");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicBPipeD, "BPipeD", logicBPipeDContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicBPipeD, "BPipeD", logicBPipeDContainer, false, 1, lxs->OverlapTest);

  // Beam pipe windows left and right sides
  G4Tubs *solidBPipeWindowD = new G4Tubs("solidBPipeWindowD", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                          bpipel/2.0, -lxs->DMBPipeWindow/2.0, lxs->DMBPipeWindow);
  G4LogicalVolume *logicBPipeWindowD = new G4LogicalVolume(solidBPipeWindowD, bpipeWindowMaterial, "logicBPipeWindowD");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicBPipeWindowD, "BPipeWindowD", logicBPipeDContainer,
                     false, 0, lxs->OverlapTest);
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicBPipeWindowD, "BPipeWindowD", logicBPipeDContainer, false, 1, lxs->OverlapTest);

  fChamberAssembly = new G4AssemblyVolume();

  G4double gaptomag = 0.5 * (lxs->TypMBMagnetCoreZ - lxs->TypMBFieldLength);
  G4ThreeVector pipetrans(0.0, 0.0, (lxs->VacChambertoOPPPDetZGap-gaptomag)/2.0);
  fChamberAssembly->AddPlacedVolume(logicBPipeDContainer, pipetrans, 0);
}


//////////////////////////////////////////////////////////////////////////////////////////////
// DetChamberAssembly

DetChamberAssembly::DetChamberAssembly(const G4String mtypename) :
  fChamberType(mtypename), fChamberAssembly(0) 
{
  ConstructVacuumChamber();
}



void DetChamberAssembly::ConstructVacuumChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* vacChamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->VacChamberMaterial);
  G4Material* vacChamberSideMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->VacChamberSideMaterial);
  G4Material* vcWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->VacChamberWindowMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);

//   G4double gaptomag = 0.5 * (lxs->TypMBMagnetCoreZ - lxs->TypMBFieldLength);
  G4double gaptomag = 0.5 * (lxs->FlashMagnetCoilZ - lxs->FlashMagnetZ);
  G4double vclength = lxs->OPPPDetZtoMagnet - lxs->VacChambertoOPPPDetZGap - gaptomag;
  G4double magnetx = 0.25 * lxs->VacChamberMagXZWidth;
  G4double detx = 0.5 * lxs->VacChamberDetXZWidth + lxs->VacChamberSideWallThickness;
  G4double dwall = 0.5 * lxs->VacChamberXZWallThickness;
  G4double vctheta = atan2(detx-magnetx, vclength + gaptomag);
  G4double vcminx = magnetx + gaptomag * tan(vctheta);

  G4Trap *solidVCContainer = new G4Trap("solidVCContainer", vclength/2.0, vctheta, 0.0, lxs->VacChamberHight/2.0,
                                    vcminx, vcminx, 0.0, lxs->VacChamberHight/2.0, detx, detx, 0.0);
  G4LogicalVolume *logicVCContainer = new G4LogicalVolume(solidVCContainer, vacuumMaterial, "logicVCContainer");

//Top and bottom wals of the chamber
  G4double thetaside = atan2(2.0*(detx-magnetx), (vclength + gaptomag));
  G4double fronttotz = lxs->VacChamberFrontThickness + lxs->VacChamberWindowPanelZ;
  G4double vcwalll = vclength - fronttotz;
  G4double vcwallx = detx - 0.5*fronttotz * tan(thetaside);

  G4Trap *solidVCXZWall = new G4Trap("solidVCXZWall", vcwalll/2.0, vctheta, 0.0, dwall, vcminx, vcminx, 0.0,
                                  dwall, vcwallx, vcwallx, 0.0);
  G4LogicalVolume *logicVCXZWall = new G4LogicalVolume(solidVCXZWall, vacChamberMaterial, "logicVCXZWall");

  G4ThreeVector wtrans((vcwallx-detx)/2.0, -(lxs->VacChamberHight/2.0 - dwall), -fronttotz/2.0);
  new G4PVPlacement(0, wtrans, logicVCXZWall, "VCXZWall", logicVCContainer, false, 0, lxs->OverlapTest);
  wtrans.setY(lxs->VacChamberHight/2.0 - dwall);
  new G4PVPlacement(0, wtrans, logicVCXZWall, "VCXZWall", logicVCContainer, false, 1, lxs->OverlapTest);

//Side wall of the chamber
  G4double  fpnldx = lxs->VacChamberSideWallThickness / cos(thetaside);
  G4double  fpnldy = lxs->VacChamberHight - 2.0 * lxs->VacChamberXZWallThickness;
  G4Para *solidSideWall = new G4Para("solidSideWall", fpnldx/2.0, fpnldy/2.0, vcwalll/2.0, 0.0, thetaside, 0.0);
  G4LogicalVolume *logicSideWall = new G4LogicalVolume(solidSideWall, vacChamberSideMaterial, "logicSideWall");

//   G4ThreeVector swtrans((vcwallx-detx-fpnldx)/2.0 + 0.5*(vcwallx + magnetx), 0.0, -fronttotz/2.0);
  G4ThreeVector swtrans((vcwallx-detx-fpnldx)/2.0 + 0.5*(vcwallx + vcminx), 0.0, -fronttotz/2.0);
  new G4PVPlacement(0, swtrans, logicSideWall, "VCSideWall", logicVCContainer, false, 0, lxs->OverlapTest);

//Front panel of the chamber
  G4double frontpnlz = lxs->VacChamberFrontThickness;
  G4Box *solidVCFront1 = new G4Box("solidVCFront1", vcwallx, lxs->VacChamberHight/2.0, frontpnlz/2.0);
  G4Tubs *solidCVBeamPipeHole = new G4Tubs("solidCVBeamPipeHole", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                           frontpnlz, 0.0, 2.0*M_PI);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(-vcwallx, 0.0, 0.0));
  G4SubtractionSolid* solidVCFront2 = new G4SubtractionSolid("solidVCFront2", solidVCFront1,
                                                                        solidCVBeamPipeHole, transform);
  G4Tubs *solidCVFrontCutR = new G4Tubs("solidCVFrontCutR", 0.0, lxs->VacChamberFrontCutY/2.0, frontpnlz, 0.0, 2.0*M_PI);
  G4Box *solidVCFrontCut1 = new G4Box("solidVCFrontCut1", lxs->VacChamberFrontCutX/2.0,
                                         lxs->VacChamberFrontCutY/2.0, frontpnlz);
  G4double reccutposx = lxs->VacChamberFrontCutXpos + lxs->VacChamberFrontCutX/2.0 - vcwallx;
  G4Transform3D trfrreccut(G4RotationMatrix(), G4ThreeVector(reccutposx, 0.0, 0.0));
  G4SubtractionSolid* solidVCFront3 = new G4SubtractionSolid("solidVCFront3", solidVCFront2, solidVCFrontCut1, trfrreccut);

  G4Transform3D trfrcutr1(G4RotationMatrix(), G4ThreeVector(reccutposx - lxs->VacChamberFrontCutX/2.0, 0.0, 0.0));
  G4SubtractionSolid* solidVCFront4 = new G4SubtractionSolid("solidVCFront4", solidVCFront3, solidCVFrontCutR, trfrcutr1);
  G4Transform3D trfrcutr2(G4RotationMatrix(), G4ThreeVector(reccutposx + lxs->VacChamberFrontCutX/2.0, 0.0, 0.0));
  G4SubtractionSolid* solidVCFront = new G4SubtractionSolid("solidVCFront", solidVCFront4, solidCVFrontCutR, trfrcutr2);

  G4LogicalVolume *logicVCFront = new G4LogicalVolume(solidVCFront, vacChamberMaterial, "logicVCFront");
  G4ThreeVector wndtrans(vcwallx - (vcminx+detx)/2.0, 0.0, (vcwalll-fronttotz+frontpnlz)/2.0);
  new G4PVPlacement(0, wndtrans, logicVCFront, "VCFront", logicVCContainer, false, 0, lxs->OverlapTest);

//Reinforcement for connecting the beam pipe
  G4double reinfy = lxs->VacChamberHight - 2.0*lxs->VacChamberXZWallThickness;
  G4Box *solidVCReinforceBar1 = new G4Box("solidVCReinforceBar1", lxs->VacChamberReinforceBarX/2.0,
                                         reinfy/2.0, lxs->VacChamberReinforceBarZ/2.0);
  G4Transform3D trreinfpipecut(G4RotationMatrix(), G4ThreeVector(-lxs->VacChamberReinforceBarX/2.0, 0.0, 0.0));
  G4SubtractionSolid* solidVCReinforceBar = new G4SubtractionSolid("solidVCReinforceBar", solidVCReinforceBar1,
                                                                        solidCVBeamPipeHole, trreinfpipecut);
  G4LogicalVolume *logicVCReinforceBar = new G4LogicalVolume(solidVCReinforceBar, vacChamberMaterial,
                                                             "logicVCReinforceBar");
  G4ThreeVector trreinbar(lxs->VacChamberReinforceBarX/2.0 - (vcminx+detx)/2.0, 0.0,
                          (vcwalll-fronttotz-lxs->VacChamberReinforceBarZ)/2.0);
  new G4PVPlacement(0, trreinbar, logicVCReinforceBar, "VCReinforceBar", logicVCContainer, false, 0, lxs->OverlapTest);

//Fillet for reinforcement part
  G4Box *solidVCReinforceFillet1 = new G4Box("solidVCReinforceFillet1", lxs->VacChamberReinforceFilletXZ/2.0,
                                         reinfy/2.0, lxs->VacChamberReinforceFilletXZ/2.0);
  G4Tubs *solidCVReinforceFilletCut = new G4Tubs("solidCVReinforceFilletCut", 0.0, lxs->VacChamberReinforceFilletXZ,
                                                 reinfy, 0.0, 2.0*M_PI);
  G4Transform3D trfiletcut(G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0),
                           G4ThreeVector(lxs->VacChamberReinforceFilletXZ/2.0, 0.0, -lxs->VacChamberReinforceFilletXZ/2.0));
  G4SubtractionSolid* solidVCReinforceFillet2 = new G4SubtractionSolid("solidVCReinforceFillet2", solidVCReinforceFillet1,
                                                           solidCVReinforceFilletCut, trfiletcut);

  G4Box *solidVCReiFiletCut = new G4Box("solidVCReiFiletCut", lxs->VacChamberReinforceFilletXZ/2.0,
                                         lxs->VacChamberFrontCutY/2.0, lxs->VacChamberReinforceFilletXZ);
  G4double fltcutposx = lxs->VacChamberFrontCutXpos - lxs->VacChamberReinforceBarX;
  G4Transform3D trfltcut(G4RotationMatrix(), G4ThreeVector(fltcutposx, 0.0, 0.0));
  G4SubtractionSolid* solidVCReinforceFillet3 = new G4SubtractionSolid("solidVCReinforceFillet3", solidVCReinforceFillet2,
                                                                       solidVCReiFiletCut, trfltcut);
  G4Transform3D fltcutrx(G4RotationMatrix(), G4ThreeVector(fltcutposx-lxs->VacChamberReinforceFilletXZ/2.0, 0.0, 0.0));
  G4SubtractionSolid* solidVCReinforceFillet = new G4SubtractionSolid("solidVCReinforceFillet",
                                                        solidVCReinforceFillet3, solidCVFrontCutR, fltcutrx);

  G4LogicalVolume *logicVCReinforceFillet = new G4LogicalVolume(solidVCReinforceFillet, vacChamberMaterial,
                                                             "logicVCReinforceFillet");
  G4ThreeVector trfilet(lxs->VacChamberReinforceBarX + lxs->VacChamberReinforceFilletXZ/2.0 - (vcminx+detx)/2.0, 0.0,
                          (vcwalll-fronttotz-lxs->VacChamberReinforceFilletXZ)/2.0);
  new G4PVPlacement(0, trfilet, logicVCReinforceFillet, "VCReinforceFillet", logicVCContainer, false, 0, lxs->OverlapTest);

//Window of the chamber
  G4Box *solidVCWindowPanel1 = new G4Box("solidVCWindowPanel1", lxs->VacChamberWindowPanelX/2.0,
                                         lxs->VacChamberWindowPanelY/2.0, lxs->VacChamberWindowPanelZ/2.0);
  G4Box *solidVCWindowCut = new G4Box("solidVCWindowCut", lxs->VacChamberWindowPanelCutX/2.0,
                                         lxs->VacChamberWindowPanelCutY/2.0, lxs->VacChamberWindowPanelZ/2.0);
  G4double wincutposx = lxs->VacChamberWindowPanelCutXpos
                        + (lxs->VacChamberWindowPanelCutX-lxs->VacChamberWindowPanelX)/2.0;
  G4Transform3D trwindowcut(G4RotationMatrix(), G4ThreeVector(wincutposx, 0.0, -lxs->VacChamberWindowThickness));
  G4SubtractionSolid* solidVCWindowPanel2 = new G4SubtractionSolid("solidVCWindowPanel2", solidVCWindowPanel1,
                                                                  solidVCWindowCut, trwindowcut);
  G4Tubs *solidCVWBeamPipeCut = new G4Tubs("solidCVWBeamPipeCut", 0.0, lxs->BPipeR,
                                           lxs->VacChamberWindowPanelZ, 0.0, 2.0*M_PI);
  G4double winpnlpipecutx = -lxs->VacChamberWindowPanelX/2.0 - lxs->VacChamberWindowPanelXpos;
  G4Transform3D trwindpipecut(G4RotationMatrix(), G4ThreeVector(winpnlpipecutx, 0.0, 0.0));
  G4SubtractionSolid* solidVCWindowPanel = new G4SubtractionSolid("solidVCWindowPanel", solidVCWindowPanel2,
                                                                        solidCVWBeamPipeCut, trwindpipecut);

  G4LogicalVolume *logicVCWindowPanel = new G4LogicalVolume(solidVCWindowPanel, vcWindowMaterial, "logicVCWindowPanel");
  G4double winpnlposx = lxs->VacChamberWindowPanelXpos + lxs->VacChamberWindowPanelX/2.0 - (vcminx+detx)/2.0;
  G4ThreeVector wndpnltrans(winpnlposx, 0.0, frontpnlz + (vcwalll - fronttotz + lxs->VacChamberWindowPanelZ)/2.0);
  new G4PVPlacement(0, wndpnltrans, logicVCWindowPanel, "VCWindowPanel", logicVCContainer, false, 0, lxs->OverlapTest);

//Short piece of beampipe between widow panels
  G4Tubs *solidVCWindBeamPipe = new G4Tubs("solidVCWindBeamPipe", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                         lxs->VacChamberWindowPanelZ/2.0, 0.0, M_PI);
  G4LogicalVolume *logicVCWindowBeamPipe = new G4LogicalVolume(solidVCWindBeamPipe, beamPipeMaterial,
                                                               "logicVCWindowBeamPipe");
  G4ThreeVector wndpipetrans(-(vcminx+detx)/2.0, 0.0, wndpnltrans.z());
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI/2.0),  wndpipetrans,
                    logicVCWindowBeamPipe, "VCWindowBeamPipe", logicVCContainer, false, 0, lxs->OverlapTest);
  
//Short piece of beampipe to join vacuum chamber with the pipe toward photon detectors
  G4double lbpipe = lxs->VacChambertoOPPPDetZGap;
  G4Tubs *solidBeamPipeVCG = new G4Tubs("solidBeamPipeVCG", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                         0.5*lbpipe, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeVCG = new G4LogicalVolume(solidBeamPipeVCG, beamPipeMaterial, "logicBeamPipeVCG");
  G4Tubs *solidBeamPipeVCGVac = new G4Tubs("solidBeamPipeVCGVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                           0.5*lbpipe, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeVCGVac = new G4LogicalVolume(solidBeamPipeVCGVac, vacuumMaterial, "logicBeamPipeVCGVac");

//Section to join vacuum chambers with the beam pipe in the field area
  G4double jclength = 0.8*gaptomag;
  G4double maptheta = atan2(2.0*(vcminx-magnetx), jclength);
  G4double hx1 = 2.0*magnetx;
  G4double hx2 = 2.0*magnetx + jclength * tan(maptheta);
  G4double jvch1 = lxs->FlashMagneteffY;
  G4double jvch2 = lxs->VacChamberHight;
  G4Trd *solidVCMagFieldJoin = new G4Trd("solidVCMagFieldJoin", hx1, hx2, jvch1/2.0, jvch2/2.0, jclength/2.0);
  G4LogicalVolume *logicVCMagFieldJoin = new G4LogicalVolume(solidVCMagFieldJoin, beamPipeMaterial, "logicVCMagFieldJoin");
  hx1 -= lxs->BPipeThickness;
  hx2 -= lxs->VacChamberSideWallThickness/cos(thetaside);
  jvch1 -= 2.0 * lxs->BPipeThickness;
  jvch2 -= 2.0 * lxs->VacChamberXZWallThickness;
  G4Trd *solidVCMagFieldJoinVac = new G4Trd("solidVCMagFieldJoinVac", hx1, hx2, jvch1/2.0, jvch2/2.0, jclength/2.0);
  G4LogicalVolume *logicVCMagFieldJoinVac = new G4LogicalVolume(solidVCMagFieldJoinVac, vacuumMaterial,
                                                                "logicVCMagFieldJoinVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicVCMagFieldJoinVac, "VCMagFieldJoinVac",
                                                   logicVCMagFieldJoin, false, 0, lxs->OverlapTest);
  G4double hx3 = 2.0*magnetx;
  G4double jclength1 = gaptomag - jclength;
  G4Box *solidVCMagFieldJoin1 = new G4Box("solidVCMagFieldJoin1", hx3, lxs->FlashMagneteffY/2.0, jclength1/2.0);
  G4LogicalVolume *logicVCMagFieldJoin1 = new G4LogicalVolume(solidVCMagFieldJoin1, beamPipeMaterial,
                                                              "logicVCMagFieldJoin1");
  hx3 -= lxs->BPipeThickness;
  G4Box *solidVCMagFieldJoin1Vac = new G4Box("solidVCMagFieldJoin1Vac", hx3,
                                             lxs->FlashMagneteffY/2.0-lxs->BPipeThickness, jclength1/2.0);
  G4LogicalVolume *logicVCMagFieldJoin1Vac = new G4LogicalVolume(solidVCMagFieldJoin1Vac, vacuumMaterial,
                                                                 "logicVCMagFieldJoin1Vac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicVCMagFieldJoin1Vac, "VCMagFieldJoin1Vac",
                                                   logicVCMagFieldJoin1, false, 0, lxs->OverlapTest);

  fVCLength = vclength;

  fChamberAssembly = new G4AssemblyVolume();
  G4ThreeVector pipetrans(0.5*vclength*tan(vctheta) + vcminx, 0.0, 0.0);
  fChamberAssembly->AddPlacedVolume(logicVCContainer, pipetrans, 0);
  pipetrans.setX(-pipetrans.getX());
  fChamberAssembly->AddPlacedVolume(logicVCContainer, pipetrans, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));
  G4ThreeVector vcbpg(0.0, 0.0, (vclength + lbpipe)/2.0);
  fChamberAssembly->AddPlacedVolume(logicBeamPipeVCG, vcbpg, 0);
  fChamberAssembly->AddPlacedVolume(logicBeamPipeVCGVac, vcbpg, 0);
  G4ThreeVector vcj(0.0, 0.0, -(vclength + jclength)/2.0);
  fChamberAssembly->AddPlacedVolume(logicVCMagFieldJoin, vcj, 0);
  G4ThreeVector vcj1(0.0, 0.0, -jclength - (vclength + jclength1)/2.0);
  fChamberAssembly->AddPlacedVolume(logicVCMagFieldJoin1, vcj1, 0);

  G4AssemblyVolume*  frameassembly = ConstructFrameAssembly();
  G4ThreeVector trframe(0.0, lxs->VacChamberHight/2.0, -lxs->VacChamberWindowPanelZ/2.0);
//   fChamberAssembly->AddPlacedAssembly(frameassembly, trframe, 0);
  LxAux::AddAssmblyVolumes(fChamberAssembly, frameassembly, trframe, 0);
  trframe.setY(-trframe.y());
//   fChamberAssembly->AddPlacedAssembly(frameassembly, trframe, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));
  LxAux::AddAssmblyVolumes(fChamberAssembly, frameassembly, trframe, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));

}



G4AssemblyVolume*  DetChamberAssembly::ConstructFrameAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* vacChamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->VacChamberMaterial);

  G4double frameLZ = fVCLength - lxs->VacChamberWindowPanelZ;
  G4Box *solidVCFrameLBar = new G4Box("solidVCFrameLBar", lxs->VacChamberFrameLX/2.0,
                                         lxs->VacChamberFrameLY/2.0, frameLZ/2.0);
  G4LogicalVolume *logicVCFrameLBar = new G4LogicalVolume(solidVCFrameLBar, vacChamberMaterial,
                                                                 "logicVCFrameLBar");

  G4Box *solidVCFrameHBar = new G4Box("solidVCFrameHBar", lxs->VacChamberFrameHX/2.0,
                                         lxs->VacChamberFrameHY/2.0, lxs->VacChamberFrameHZ/2.0);
  G4LogicalVolume *logicVCFrameHBar = new G4LogicalVolume(solidVCFrameHBar, vacChamberMaterial,
                                                                 "logicVCFrameHBar");

  G4AssemblyVolume *frameAssembly = new G4AssemblyVolume();
  G4ThreeVector trlbar((lxs->VacChamberFrameHX + lxs->VacChamberFrameLX)/2.0, lxs->VacChamberFrameLY/2.0, 0.0);
  frameAssembly->AddPlacedVolume(logicVCFrameLBar, trlbar, 0);
  trlbar.setX(-trlbar.x());
  frameAssembly->AddPlacedVolume(logicVCFrameLBar, trlbar, 0);

  G4ThreeVector trhbar(0.0, lxs->VacChamberFrameHY/2.0, 0.0);
  frameAssembly->AddPlacedVolume(logicVCFrameHBar, trhbar, 0);

  return frameAssembly;
}
