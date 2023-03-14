//
/// \brief Implementation of the LxECal class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4GenericTrap.hh"
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
#include "LxHICSBeam.hh"


////////////////////////////////////////////////////////////////////////
//// LxHICSBeam

void LxHICSBeam::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* hicsDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSDumpMaterial);
  G4Material* HDInsertAlMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSDumpAlInsertMaterial);
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  fHICSDumpAngle = atan2(lxs->HICSDumpFrontXPos, lxs->HICSDumpFrontZPos - lxs->IPMagnetZpos);

  G4Tubs *solidHICSDumpGammaBeamHoleCut = new G4Tubs("solidHICSDumpGammaBeamHoleCut", 0.0, lxs->CollimatorRin,
                                             lxs->HICSDumpZ, 0.0, 2.0*M_PI);

  G4VSolid *solidHICSDump = 0;
  G4VSolid *solidHICSDumpAlInsert = 0;
  G4VSolid *solidHICSDumpAirInsert = 0;

  if ( lxs->HICSDumpR > lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle) ) {
    G4Tubs *solidHICSDump0 = new G4Tubs("solidHICSDump0", 0.0, lxs->HICSDumpR, lxs->HICSDumpZ/2.0, 0.0, 2.0*M_PI);
    G4double dx = tan(fHICSDumpAngle) * lxs->HICSDumpZ/2.0 + lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle);
    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), fHICSDumpAngle),
                            G4ThreeVector(dx, 0.0, 0.0));
    solidHICSDump = new G4SubtractionSolid("solidHICSDump", solidHICSDump0, solidHICSDumpGammaBeamHoleCut, transform);
  } else {
    solidHICSDump = new G4Tubs("solidHICSDump", 0.0, lxs->HICSDumpR, lxs->HICSDumpZ/2.0, 0.0, 2.0*M_PI);
  }

  if ( lxs->HICSDumpAlInsertR > lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle) ) {
    G4Tubs *solidHICSDumpAlInsert0 = new G4Tubs("solidHICSDumpAlInsert0", 0.0, lxs->HICSDumpAlInsertR,
                                               lxs->HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);
    G4double dx = tan(fHICSDumpAngle) * lxs->HICSDumpAlInsertZ/2.0 *  + lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle);
    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), fHICSDumpAngle),
                            G4ThreeVector(dx, 0.0, 0.0));
    solidHICSDumpAlInsert = new G4SubtractionSolid("solidHICSDumpAlInsert", solidHICSDumpAlInsert0,
                                solidHICSDumpGammaBeamHoleCut, transform);
  } else {
    solidHICSDumpAlInsert = new G4Tubs("solidHICSDumpAlInsert", 0.0, lxs->HICSDumpAlInsertR,
                                               lxs->HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);
  }

  if ( lxs->HICSDumpAirInsertR > lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle) ) {
    G4Tubs *solidHICSDumpAirInsert0 = new G4Tubs("solidHICSDumpAirInsert0", 0.0, lxs->HICSDumpAirInsertR,
                                             lxs->HICSDumpAirInsertZ/2.0, 0.0, 2.0*M_PI);
    G4double dx = tan(fHICSDumpAngle) * lxs->HICSDumpAirInsertZ/2.0 *  + lxs->HICSDumpFrontXPos/cos(fHICSDumpAngle);
    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), fHICSDumpAngle),
                            G4ThreeVector(dx, 0.0, 0.0));
    solidHICSDumpAirInsert = new G4SubtractionSolid("solidHICSDumpAirInsert", solidHICSDumpAirInsert0,
                                solidHICSDumpGammaBeamHoleCut, transform);
  } else {
    solidHICSDumpAirInsert = new G4Tubs("solidHICSDumpAirInsert", 0.0, lxs->HICSDumpAirInsertR,
                                             lxs->HICSDumpAirInsertZ/2.0, 0.0, 2.0*M_PI);
  }

  G4LogicalVolume *logicHICSDump = new G4LogicalVolume(solidHICSDump, hicsDumpMaterial, "logicHICSDump");
  G4LogicalVolume *logicHICSDumpAlInsert = new G4LogicalVolume(solidHICSDumpAlInsert, HDInsertAlMaterial,
                                                               "logicHICSDumpAlInsert");
  G4LogicalVolume *logicHICSDumpAirInsert = new G4LogicalVolume(solidHICSDumpAirInsert, environmentMaterial,
                                                       "logicHICSDumpAirInsert");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->HICSDumpAirInsertZ-lxs->HICSDumpAlInsertZ)/2.0),
                     logicHICSDumpAirInsert, "HICSDumpAirInsert", logicHICSDumpAlInsert, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->HICSDumpAlInsertZ-lxs->HICSDumpZ)/2.0),
                     logicHICSDumpAlInsert, "HICSDumpAlInsert", logicHICSDump, false, 0, lxs->OverlapTest);

  G4double dupm2tposx = lxs->HICSDump2TPosX;
  new G4PVPlacement (0, G4ThreeVector(-dupm2tposx, 0.0, (lxs->HICSDumpAlInsertZ-lxs->HICSDumpZ)/2.0),
                     logicHICSDumpAlInsert, "HICSDumpAlInsert", logicHICSDump, false, 1, lxs->OverlapTest);

  G4double hicsdx = lxs->HICSDumpFrontXPos + 0.5*lxs->HICSDumpZ * sin(fHICSDumpAngle);
  G4double hicsdz = lxs->HICSDumpFrontZPos + 0.5*lxs->HICSDumpZ * cos(fHICSDumpAngle);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), fHICSDumpAngle),
                    G4ThreeVector(-hicsdx, 0.0, hicsdz),
                               logicHICSDump, "HICSDumpAssembly", fLogicWorld, false, 0, lxs->OverlapTest);

  //vacuum in the hole  for gamma beam
  G4double frontZpos = lxs->HICSDumpFrontZPos  + lxs->HICSDumpFrontXPos*tan(fHICSDumpAngle);
  G4double lholeVac = lxs->HICSDumpZ/cos(fHICSDumpAngle);
  G4double lpipe_pos_1st = frontZpos + lholeVac/2.0;
  G4CutTubs *solidBeamPipeGammaHICSdumpVac = new G4CutTubs("solidBeamPipeGammaHICSdumpVac", 0.0, lxs->CollimatorRin,
                         lholeVac/2.0, 0.0, 2.0*M_PI, G4ThreeVector(sin(fHICSDumpAngle), 0.0, -cos(fHICSDumpAngle)), G4ThreeVector(-sin(fHICSDumpAngle), 0.0, cos(fHICSDumpAngle)));

  G4LogicalVolume *logicBeamPipeGammaHICSdumpVac = new G4LogicalVolume(solidBeamPipeGammaHICSdumpVac, vacuumMaterial, "logicBeamPipeGammaHICSdumpVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_1st),
                    logicBeamPipeGammaHICSdumpVac, "BeamPipeGammaHICSdumpVac", fLogicWorld, false, 0, lxs->OverlapTest);

// 2T dump component
  G4double z2t = lxs->HICSDumpZ;
  G4CutTubs *solidHICSDump2T0 = new G4CutTubs( "solidHICSDump2T0", 0.0, lxs->HICSDumpR, z2t/2.0, 0.0, 2.0*M_PI,
          G4ThreeVector(-sin(fHICSDumpAngle), 0.0, -cos(fHICSDumpAngle)), G4ThreeVector(0.0, 0.0, 1.0));
  G4Tubs *solidHICSDump2TCut = new G4Tubs("solidHICSDump2TCut", 0.0, lxs->HICSDumpR, lxs->HICSDumpZ, 0.0, 2.0*M_PI);
  G4Transform3D transf2t(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), 0.0),
                            G4ThreeVector(dupm2tposx*cos(fHICSDumpAngle), 0.0, 0.0));
  G4SubtractionSolid *solidHICSDump2T = new G4SubtractionSolid("solidHICSDump2T", solidHICSDump2T0,
                                solidHICSDump2TCut, transf2t);
  G4LogicalVolume *logicHICSDump2T = new G4LogicalVolume(solidHICSDump2T, hicsDumpMaterial, "logicHICSDump2T");
  G4double zpos2t = hicsdz - dupm2tposx*sin(fHICSDumpAngle);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), fHICSDumpAngle),
                    G4ThreeVector(-hicsdx-dupm2tposx*cos(fHICSDumpAngle), 0.0, zpos2t),
                               logicHICSDump2T, "HICSDump2T", fLogicWorld, false, 0, lxs->OverlapTest);

  ConstructElectronShielding(fLogicWorld);
  ConstructNeutronAbsorber(fLogicWorld);

  G4double supy;
  G4AssemblyVolume *supportAssembly = ConstructSupportAssembly(supy);
  G4double supxpos = lxs->HICSShieldingXPos;
  G4double supypos = -(lxs->HICSShieldingY + supy)/2.0;
  G4double supzpos = lxs->HICSDumpFrontZPos - lxs->HICSDumpR*sin(fHICSDumpAngle) - lxs->HICSShieldingSideZ;
  G4ThreeVector hicssuptr(supxpos, supypos, supzpos);
  supportAssembly->MakeImprint(fLogicWorld, hicssuptr, 0, 1, lxs->OverlapTest);

//Probably not the best place for constructing detector support, but it is HICS related, so can be here
  G4double supporthight;
  G4AssemblyVolume* hicsDetSupport = ConstructHICSElDetSupportAssembly(supporthight);
  if (hicsDetSupport) {
    G4double detxpos = lxs->HICSDetBottomSupportX/2.0;// lxs->OPPPDetTopPlateX/2.0 + lxs->OPPPDetXPos;
    G4double detypos = -lxs->ComptonElectronBeamtoStageY - lxs->ScintFrameThickness;
    // G4double detzpos = lxs->IPMagnetZpos + lxs->TypMBFieldLength/2.0 + lxs->OPPPDetZtoMagnet + lxs->OPPPDetTopPlateZ/2.0;
//   G4double detzpos = lxs->IPMagnetZpos + lxs->TypMBFieldLength/2.0  + lxs->ComptonElBackshift;
  G4double detzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0  + lxs->ComptonElBackshift;
                  
    G4ThreeVector trsupport(-detxpos, detypos, detzpos);
    hicsDetSupport->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);
  }

}



void LxHICSBeam::ConstructElectronShielding(G4LogicalVolume  *logicWorld)
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* hicsShielSMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSShieldingSideMaterial);
  G4Material* hicsShielMMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSShieldingMiddleMaterial);

//Shielding front plane
  G4Box *solidHICSShieldingSide1 = new G4Box("solidHICSShieldingSide1", lxs->HICSShieldingX/2.0,
                                             lxs->HICSShieldingY/2.0, lxs->HICSShieldingSideZ/2.0);
  G4Tubs *solidHICSShieldSideCut = new G4Tubs("solidHICSShieldSideCut", 0.0, lxs->HICSDumpAlInsertR,
                                               lxs->HICSDumpAlInsertZ/2.0, 0.0, 2.0*M_PI);
  G4double dx = -(lxs->HICSShieldingXPos + lxs->HICSDumpFrontXPos);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(dx, 0.0, 0.0));
  G4SubtractionSolid *solidHICSShieldingSide0 = new G4SubtractionSolid("solidHICSShieldingSide0",
                                                solidHICSShieldingSide1, solidHICSShieldSideCut, transform);

//   G4double xscut = lxs->HICSDumpFrontXPos + 4.0*lxs->CollimatorRin;
//   G4Box *solidHICSShieldingSideCut1 = new G4Box("solidHICSShieldingSideCut1", xscut/2.0,
//                                              2.0*lxs->CollimatorRin, lxs->HICSShieldingSideZ);
//   G4Transform3D transscut(G4RotationMatrix(), G4ThreeVector(dx+xscut/2.0, 0.0, 0.0));
//   G4SubtractionSolid *solidHICSShieldingSide2 = new G4SubtractionSolid("solidHICSShieldingSide2",
//                                                 solidHICSShieldingSide0, solidHICSShieldingSideCut1, transscut);

  G4Tubs *solidHICSShieldSide4Bpipe = new G4Tubs("solidHICSShieldSideCut4Bpipe", 0.0, lxs->BPipeRLG,
                                                 lxs->HICSShieldingSideZ, 0.0, 2.0*M_PI);
  G4Transform3D transform2(G4RotationMatrix(), G4ThreeVector(-lxs->HICSShieldingXPos, 0.0, 0.0));
  G4SubtractionSolid* solidHICSShieldingSide2 = new G4SubtractionSolid("solidHICSShieldingSide2", solidHICSShieldingSide0,
                                      solidHICSShieldSide4Bpipe, transform2);

  G4double xgapcut = lxs->HICSShieldingGapX;
  G4Box *solidHICSShieldingSideGapCut = new G4Box("solidHICSShieldingSideGapCut", xgapcut/2.0,
                                             lxs->HICSShieldingGapY/2.0, lxs->HICSShieldingSideZ/2.0);
  G4Transform3D trnsgapcut(G4RotationMatrix(), G4ThreeVector(dx-xgapcut/2.0, 0.0, -lxs->HICSShieldingSideZ/2.0));
  G4SubtractionSolid *solidHICSShieldingSide3 = new G4SubtractionSolid("solidHICSShieldingSide3",
                                     solidHICSShieldingSide2, solidHICSShieldingSideGapCut, trnsgapcut);

  G4double dupm2tposx = lxs->HICSDump2TPosX;
  G4double dx2t = -(lxs->HICSShieldingXPos + lxs->HICSDumpFrontXPos + dupm2tposx);
  G4Transform3D transform2t(G4RotationMatrix(), G4ThreeVector(dx2t, 0.0, 0.0));
  G4SubtractionSolid *solidHICSShieldingSide = new G4SubtractionSolid("solidHICSShieldingSide",
                                     solidHICSShieldingSide3, solidHICSShieldSideCut, transform2t);

  G4LogicalVolume *logicHICSShieldingSide = new G4LogicalVolume(solidHICSShieldingSide, hicsShielSMaterial,
                                                       "logicHICSShieldingSide");

//Shielding back plane
  G4Box *solidHICSShieldingMiddle0 = new G4Box("solidHICSShieldingMiddle0", lxs->HICSShieldingX/2.0,
                                             lxs->HICSShieldingY/2.0, lxs->HICSShieldingMiddleZ/2.0);
  G4Tubs *solidHICSSieldMiddleCut = new G4Tubs("solidHICSSieldMiddleCut", 0.0, lxs->HICSDumpR,
                                               lxs->HICSShieldingMiddleZ, 0.0, 2.0*M_PI);
  G4double dxm = tan(fHICSDumpAngle)*(0.5*lxs->HICSShieldingMiddleZ - lxs->HICSDumpR*sin(fHICSDumpAngle));
  G4Transform3D transm(G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), fHICSDumpAngle),
                          G4ThreeVector(-dxm-lxs->HICSShieldingXPos-lxs->HICSDumpFrontXPos, 0.0, 0.0));
  G4SubtractionSolid *solidHICSShieldingMiddle1 = new G4SubtractionSolid("solidHICSShieldingMiddle1",
                                                solidHICSShieldingMiddle0, solidHICSSieldMiddleCut, transm);
  G4Transform3D transm1(G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), fHICSDumpAngle),
      G4ThreeVector(-dxm-lxs->HICSShieldingXPos-lxs->HICSDumpFrontXPos-dupm2tposx/cos(fHICSDumpAngle), 0.0, 0.0));
  G4SubtractionSolid *solidHICSShieldingMiddle = new G4SubtractionSolid("solidHICSShieldingMiddle",
                                                solidHICSShieldingMiddle1, solidHICSSieldMiddleCut, transm1);

  G4LogicalVolume *logicHICSShieldingMiddle = new G4LogicalVolume(solidHICSShieldingMiddle, hicsShielMMaterial,
                                                       "logicHICSShieldingMiddle");

  G4double zpos = lxs->HICSDumpFrontZPos - lxs->HICSDumpR*sin(fHICSDumpAngle) - lxs->HICSShieldingSideZ/2.0;
  new G4PVPlacement(0, G4ThreeVector(lxs->HICSShieldingXPos, 0.0, zpos),
                               logicHICSShieldingSide, "HICSShieldingSide", logicWorld, false, 0, lxs->OverlapTest);
  zpos += lxs->HICSShieldingSideZ/2.0 + lxs->HICSShieldingMiddleZ/2.0;
  new G4PVPlacement(0, G4ThreeVector(lxs->HICSShieldingXPos, 0.0, zpos),
                    logicHICSShieldingMiddle, "HICSShieldingMiddle", logicWorld, false, 0, lxs->OverlapTest);
}



void LxHICSBeam::ConstructNeutronAbsorber(G4LogicalVolume  *logicWorld)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* hicsNeutronAbsorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSNeutronAbsorberMaterial);

//Neutron absorber front plane
  G4Box *solidHICSNeutronAbsorberSide1 = new G4Box("solidHICSNeutronAbsorberSide1", lxs->HICSNeutronAbsorberX/2.0,
                                             lxs->HICSNeutronAbsorberY/2.0, lxs->HICSNeutronAbsorberZ/2.0);
  G4Tubs *solidNeutronAbsorberSideCut = new G4Tubs("solidNeutronAbsorberSideCut", 0.0, lxs->HICSDumpAlInsertR,
                                               lxs->HICSNeutronAbsorberZ, 0.0, 2.0*M_PI);
  G4double dx = -(lxs->HICSShieldingXPos + lxs->HICSDumpFrontXPos);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(dx, 0.0, 0.0));
  G4SubtractionSolid *solidNeutronAbsorberSide0 = new G4SubtractionSolid("solidNeutronAbsorberSide0",
                                                solidHICSNeutronAbsorberSide1, solidNeutronAbsorberSideCut, transform);

  G4Tubs *solidNeutronAbsorber4Bpipe = new G4Tubs("solidNeutronAbsorber4Bpipe", 0.0, lxs->BPipeRLG,
                                                 lxs->HICSShieldingSideZ, 0.0, 2.0*M_PI);
  G4Transform3D transform2(G4RotationMatrix(), G4ThreeVector(-lxs->HICSShieldingXPos, 0.0, 0.0));
  G4SubtractionSolid* solidNeutronAbsorberSide2 = new G4SubtractionSolid("solidNeutronAbsorberSide2", solidNeutronAbsorberSide0,
                                      solidNeutronAbsorber4Bpipe, transform2);

  G4double xgapcut = lxs->HICSShieldingGapX;
  G4Box *solidNeutronAbsorberGapCut = new G4Box("solidNeutronAbsorberGapCut", xgapcut/2.0,
                                             lxs->HICSShieldingGapY/2.0, lxs->HICSNeutronAbsorberZ);
  G4Transform3D trnsgapcut(G4RotationMatrix(), G4ThreeVector(dx-xgapcut/2.0, 0.0, 0.0));
  G4SubtractionSolid *solidNeutronAbsorberSide3 = new G4SubtractionSolid("solidNeutronAbsorberSide3",
                                     solidNeutronAbsorberSide2, solidNeutronAbsorberGapCut, trnsgapcut);

  G4double dupm2tposx = lxs->HICSDump2TPosX;
  G4double dx2t = -(lxs->HICSShieldingXPos + lxs->HICSDumpFrontXPos + dupm2tposx);
  G4Transform3D transform2t(G4RotationMatrix(), G4ThreeVector(dx2t, 0.0, 0.0));
  G4SubtractionSolid *solidNeutronAbsorberSide = new G4SubtractionSolid("solidNeutronAbsorberSide",
                                     solidNeutronAbsorberSide3, solidNeutronAbsorberSideCut, transform2t);

  G4LogicalVolume *logicNeutronAbsorberSide = new G4LogicalVolume(solidNeutronAbsorberSide, hicsNeutronAbsorberMaterial,
                                                       "logicNeutronAbsorberSide");
// Neutron absorber top part
  G4double supy = -lxs->FloorSurfaceYpos - lxs->HICSNeutronAbsorberY/2.0;
  G4double suptopx = 1.1 * lxs->HICSNeutronAbsorberX;
  G4double suptopy = lxs->HICSShieldingSupportY - supy;
  G4double suptopz = lxs->HICSNeutronAbsorberZ;

  G4Box *solidHICSNeutronAbsorberTop1 = new G4Box("solidHICSNeutronAbsorberTop1", suptopx/2.0, suptopy/2.0, suptopz/2.0);
  G4Box *solidHICSNeutronAbsorberTopCut1 = new G4Box("solidHICSNeutronAbsorberTopCut1", lxs->HICSNeutronAbsorberX/2.0,
                                             lxs->HICSNeutronAbsorberY, lxs->HICSShieldingSideZ);
  G4Transform3D transcut1(G4RotationMatrix(), G4ThreeVector(0.0, -suptopy/2.0, 0.0));
  G4SubtractionSolid *soliHICSNeutronAbsorberTop = new G4SubtractionSolid("soliHICSNeutronAbsorberTop", solidHICSNeutronAbsorberTop1,
                                solidHICSNeutronAbsorberTopCut1, transcut1);
  G4LogicalVolume *logicHICSNeutronAbsorberTop = new G4LogicalVolume(soliHICSNeutronAbsorberTop, hicsNeutronAbsorberMaterial,
                                                       "logicHICSNeutronAbsorberTop");
// Neutron absorber bottom part
  G4Box *solidHICSNeutronAbsorberBottom = new G4Box("solidHICSNeutronAbsorberBottom", suptopx/2.0, supy/2.0, suptopz/2.0);
  G4LogicalVolume *logicHICSNeutronAbsorberBottom = new G4LogicalVolume(solidHICSNeutronAbsorberBottom, hicsNeutronAbsorberMaterial,
                                                       "logicHICSNeutronAbsorberBottom");
// Placement of the parts
  G4double zpos = lxs->HICSDumpFrontZPos - lxs->HICSDumpR*sin(fHICSDumpAngle) - lxs->HICSShieldingSideZ - lxs->HICSNeutronAbsorberZ/2.0;
  new G4PVPlacement(0, G4ThreeVector(lxs->HICSShieldingXPos, 0.0, zpos),
                               logicNeutronAbsorberSide, "HICSNeutronAbsorberSide", logicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(lxs->HICSShieldingXPos, 0.5*(suptopy-lxs->HICSNeutronAbsorberY), zpos),
                               logicHICSNeutronAbsorberTop, "HICSNeutronAbsorberTop", logicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(lxs->HICSShieldingXPos, lxs->FloorSurfaceYpos + 0.5*supy, zpos),
                    logicHICSNeutronAbsorberBottom, "HICSNeutronAbsorberBottom", logicWorld, false, 0, lxs->OverlapTest);

}



G4AssemblyVolume* LxHICSBeam::ConstructSupportAssembly(G4double &supporthight)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* hicsShieldSupportMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSShieldingSupportMaterial);
  G4Material* hicsShieldTopMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSShieldingTopMaterial);

  G4double supx = 1.1 * lxs->HICSShieldingX;
  G4double supy = -lxs->FloorSurfaceYpos - lxs->HICSShieldingY/2.0;
  G4double supz = lxs->HICSShieldingSupportZ;
  G4LogicalVolume *pedstal = LxAux::BuildPedestal("HICSShield", supx, supy, supz);

// Top concrete part
  G4double suptopx = supx;
  G4double suptopy = lxs->HICSShieldingSupportY - supy;
  G4double suptopz = supz;
  G4Box *solidHICSSSupportTop1 = new G4Box("solidHICSSSupportTop1", suptopx/2.0, suptopy/2.0, suptopz/2.0);
  G4double suptopcutz = lxs->HICSShieldingSideZ + lxs->HICSShieldingMiddleZ;
  G4Box *solidHICSSSupportTopCut1 = new G4Box("solidHICSSSupportTopCut1", lxs->HICSShieldingX/2.0,
                                             lxs->HICSShieldingY, suptopcutz);
  G4Transform3D transcut1(G4RotationMatrix(), G4ThreeVector(0.0, -suptopy/2.0, -suptopz/2.0));
  G4SubtractionSolid *soliHICSSSupportTop2 = new G4SubtractionSolid("soliHICSSSupportTop2", solidHICSSSupportTop1,
                                solidHICSSSupportTopCut1, transcut1);

  G4double topdumpcutx = 1.2 * (2.0 * lxs->HICSDumpR + lxs->HICSDump2TPosX);
  G4double topdumpcuty = lxs->HICSShieldingY;
  G4double topdumpcutz = 1.1 * lxs->HICSDumpZ;
  G4Box *solidHICSSSupportTopCut2 = new G4Box("solidHICSSSupportTopCut2", topdumpcutx/2.0, topdumpcuty, topdumpcutz/2.0);
  G4double trdumpcutx = -lxs->HICSShieldingXPos - lxs->HICSDumpFrontXPos - 0.5*lxs->HICSDump2TPosX;
  G4double trdumpcuty = -suptopy/2.0;
  G4double trdumpcutz = -0.5*(suptopz - topdumpcutz) + lxs->HICSShieldingSideZ;
  G4Transform3D transcut2(G4RotationMatrix(), G4ThreeVector(trdumpcutx, trdumpcuty, trdumpcutz));
  G4SubtractionSolid *soliHICSSSupportTop3 = new G4SubtractionSolid("soliHICSSSupportTop3", soliHICSSSupportTop2,
                                solidHICSSSupportTopCut2, transcut2);

  G4Tubs *solidHICSSSupportTopCut3 = new G4Tubs("solidHICSSSupportTopCut3", 0.0, lxs->CollimatorRout,
                                                supz, 0.0, 2.0*M_PI);
  G4double trcolimcutx = -lxs->HICSShieldingXPos;
  G4double trcolimcuty = -0.5*(suptopy-lxs->HICSShieldingY);
  G4Transform3D transcut3(G4RotationMatrix(), G4ThreeVector(trcolimcutx, trcolimcuty, 0.0));
  G4SubtractionSolid *soliHICSSSupportTop = new G4SubtractionSolid("soliHICSSSupportTop4", soliHICSSSupportTop3,
                                solidHICSSSupportTopCut3, transcut3);

  G4LogicalVolume *logicHICSSSupportTop = new G4LogicalVolume(soliHICSSSupportTop, hicsShieldTopMaterial,
                                                       "logicHICSSSupportTop");
// Add volumes to the assembly
  G4AssemblyVolume *hicsShieldSupportAssembly = new G4AssemblyVolume();
  G4ThreeVector hicssuptr(0.0, 0.0, lxs->HICSShieldingSupportZ/2.0);
  hicsShieldSupportAssembly->AddPlacedVolume(pedstal, hicssuptr, 0);
  G4ThreeVector hicssuptoptr(0.0, 0.5*(supy+suptopy), lxs->HICSShieldingSupportZ/2.0);
  hicsShieldSupportAssembly->AddPlacedVolume(logicHICSSSupportTop, hicssuptoptr, 0);

  supporthight = supy;
  return hicsShieldSupportAssembly;
}



G4AssemblyVolume* LxHICSBeam::ConstructHICSElDetSupportAssembly(G4double &sphight)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppDetSupportMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->HICSDetSupportMaterial);

  G4Box *soliHICSScintBottomSupport = new G4Box("soliHICSScintBottomSupport", lxs->HICSDetBottomSupportX/2.0, 
                                            lxs->HICSDetBottomSupportY/2.0, lxs->HICSDetBottomSupportZ/2.0);
  G4LogicalVolume *logicHICSScintBottomSupport = new G4LogicalVolume(soliHICSScintBottomSupport, 
                                   opppDetSupportMaterial, "logicHICSScintBottomSupport");

//   G4double detzpos = lxs->IPMagnetZpos + lxs->TypMBFieldLength/2.0  + lxs->ComptonElBackshift;
  G4double detzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0  + lxs->ComptonElBackshift;
  
  G4int ccl = lxs->CerenkovChannelLayers;
  G4double cszf = lxs->CerenkovStrawZFrequency;
  G4double sir = lxs->CerenkovStrawInnerRadius;
  G4double sglt = lxs->CerenkovStrawGraphiteLayerThickness;
  G4double salt = lxs->CerenkovStrawAlLayerThickness;
  G4double sklt = lxs->CerenkovStrawKaptonLayerThickness;
  G4double splt = lxs->CerenkovStrawPolyurethaneLayerThickness;
  G4double cspb = lxs->CerenkovStrawPlateBuffer;
  G4double cbt = lxs->CerenkovBoxThickness;
  G4double cspt = lxs->CerenkovShieldingPlateThickness;

  G4double scmotherdz = 2*lxs->ComptonElectronBeamtoStageY + 3*lxs->ScintFrameThickness;
  G4double cherenkovdz = 2.0 * ((ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);
  G4double basezpos = detzpos + lxs->HICSDetBottomSupportZ/2. - scmotherdz/2.;

  G4AssemblyVolume *hicsSupportAssembly = new G4AssemblyVolume();

  G4ThreeVector asstr(0.0, 0.0, basezpos - detzpos);
  hicsSupportAssembly->AddPlacedVolume(logicHICSScintBottomSupport, asstr, 0);

  //Hexapod
  G4double hexhight = 0.0;
  G4AssemblyVolume* hexAssembly = LxAux::BuildHexapod("HICS", hexhight);
  G4ThreeVector hexpos(0.0, asstr.y() - (lxs->OPPPDetBottomSupportY+hexhight)/2.0,  basezpos - detzpos);
//   hicsSupportAssembly->AddPlacedAssembly(hexAssembly, hexpos, 0);
  LxAux::AddAssmblyVolumes(hicsSupportAssembly, hexAssembly, hexpos, 0);

  G4double baseypos = -lxs->ComptonElectronBeamtoStageY - hexhight - lxs->HICSDetBottomSupportY - (lxs->HICSBasePlateY)/2.0;
  sphight = (lxs->HICSDetBottomSupportY + hexhight)/2.0 - hexpos.y();

  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

// Table support
  G4double ypestal = 1.5*m;
  G4double ylevel = baseypos + lxs->HICSBasePlateY/2.0;
  G4double tblhight = ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4double tblx = lxs->HICSBasePlateX;
  G4double tblz = lxs->HICSBasePlateZ;
  G4AssemblyVolume *tablesupport = LxAux::BuildTable("HICSDetectorTable", tblx, tblhight, tblz, 4);
  G4ThreeVector trsupport(-0.5*tblx, ylevel, basezpos);
  tablesupport->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);

  G4LogicalVolume *pedstal = LxAux::BuildPedestal("HICSDetector", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(-0.5*tblx, 0.5*ypestal + lxs->FloorSurfaceYpos, basezpos), pedstal,
                      "HICSDetectorPedestal", fLogicWorld, false, 0, lxs->OverlapTest);

  return hicsSupportAssembly;
}



//Neutron absorber material B1PE50 (Journal of Nuclear Materials 452 (2014) 205–211)
void LxHICSBeam::CreateMaterial()
{
  G4double z, a, density, fractionmass;
  G4int ncomponents, natoms;
  G4String symbol;

  G4Element* H  = new G4Element("Hydrogen", symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  , symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  , symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silcone" , symbol="Si", z= 14., a= 28.09*g/mole);
  G4Element* Na = new G4Element("Sodium"  ,  symbol="Na", z=11,   a= 22.99*g/mole);
  G4Element* Mg = new G4Element("Magnesium", symbol="Mg", z=12.,  a= 24.305*g/mole);
  G4Element* Al = new G4Element("Aluminum",  symbol="Al", z= 13.,  a= 26.9815386*g/mole);
  G4Element* K  = new G4Element("Potassium", symbol="K",  z= 19.,  a= 39.098*g/mole);
  G4Element* Ca = new G4Element("Calcium"  , symbol="Ca", z= 20.,  a= 40.078*g/mole);
  G4Element* Fe = new G4Element("Iron"     , symbol="Fe", z= 26,   a= 55.845*g/mole);
  G4Element* B  = new G4Element("Boron"    , symbol="B" , z= 5,    a= 10.81*g/mole);

  G4Material *concrete = new G4Material("RegularConcrete", density= 2.336*g/cm3, ncomponents=9);
  concrete->AddElement(H, fractionmass=0.5568*perCent);
  concrete->AddElement(O, fractionmass=49.892*perCent);
  concrete->AddElement(Na, fractionmass=1.71623*perCent);
  concrete->AddElement(Mg, fractionmass=0.2592*perCent);
  concrete->AddElement(Al, fractionmass=4.58488*perCent);
  concrete->AddElement(Si, fractionmass=31.55*perCent);
  concrete->AddElement(K, fractionmass=1.91807*perCent);
  concrete->AddElement(Ca, fractionmass=8.29198*perCent);
  concrete->AddElement(Fe, fractionmass=1.23084*perCent);

  G4Material *boronCarbide = new G4Material("BoronCarbide", density= 2.55*g/cm3, ncomponents=2);
  boronCarbide->AddElement(B, fractionmass=79.99*perCent);
  boronCarbide->AddElement(C, fractionmass=20.01*perCent);

  G4Material *polyethylene = new G4Material("Polyethylene", density= 0.93*g/cm3, ncomponents=2);
  polyethylene->AddElement(H, fractionmass=14.37*perCent);
  polyethylene->AddElement(C, fractionmass=85.63*perCent);

  G4Material *boratedPolyethyleneConcrete = new G4Material("BoratedPolyethyleneConcrete", density= 3.5*g/cm3, ncomponents=3);
  boratedPolyethyleneConcrete->AddMaterial(concrete, fractionmass=49.0*perCent);
  boratedPolyethyleneConcrete->AddMaterial(polyethylene, fractionmass=50.0*perCent);
  boratedPolyethyleneConcrete->AddMaterial(boronCarbide, fractionmass=1.0*perCent);

}
