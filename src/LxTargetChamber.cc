//
/// \brief Implementation of the BeamProfiler tracker class
//

#include <numeric>
#include <cmath>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"
#include "LXSetUp.hh"
#include "LxAux.hh"
#include "LxDetector.hh"
#include "LxTargetChamber.hh"



void LxTargetChamber::ConstructTargetChamber()
{
//   CreateMaterial();
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* chamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TargetChamberMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  G4double tcshiftx = lxs->TargetChamberXPos;

  G4Box* solidTargetChamberContainer1 = new G4Box("solidTargetChamberContainer1", lxs->TargetChamberX/2.0,
                                                lxs->TargetChamberY/2.0, lxs->TargetChamberZ/2.0);
  G4Box* solidChamberTargetVolume = new G4Box("solidChamberTargetVolume", lxs->ChamberTargetVolX/2.0,
                                                lxs->ChamberTargetVolY/2.0, lxs->ChamberTargetVolZ/2.0);
  G4Transform3D trtarvolcut(G4RotationMatrix(), G4ThreeVector(-tcshiftx, 0.0, 0.0));
  G4SubtractionSolid* solidTargetChamberContainer = new G4SubtractionSolid("solidTargetChamberContainer",
                                        solidTargetChamberContainer1, solidChamberTargetVolume, trtarvolcut);

  G4LogicalVolume* logicTargetChamberContainer = new G4LogicalVolume(solidTargetChamberContainer, vacuumMaterial,
                                                                   "logicTargetChamberContainer");

  G4Box* solidTargetChamberBoxOut = new G4Box("solidTargetChamberBoxOut", lxs->TargetChamberX/2.0,
                                                lxs->TargetChamberY/2.0, lxs->TargetChamberZ/2.0);
  G4double tcwallthick = lxs->TargetChamberThickness;
  G4Box* solidTargetChamberBoxIn = new G4Box("solidTargetChamberBoxIn", lxs->TargetChamberX/2.0 - tcwallthick,
                               lxs->TargetChamberY/2.0 - tcwallthick, lxs->TargetChamberZ/2.0 - tcwallthick);
  G4Tubs *solidBeamPipeCut = new G4Tubs("solidBeamPipeCut", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                           tcwallthick, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidTargetChamberBox1 = new G4SubtractionSolid("solidTargetChamberBox1",
                                                solidTargetChamberBoxOut, solidTargetChamberBoxIn);

  G4Transform3D transform1(G4RotationMatrix(), G4ThreeVector(-tcshiftx, 0.0,  lxs->TargetChamberZ/2.0));
  G4SubtractionSolid* solidTargetChamberBox2 = new G4SubtractionSolid("solidTargetChamberBox2", solidTargetChamberBox1,
                                                                      solidBeamPipeCut, transform1);
  G4Transform3D transform2(G4RotationMatrix(), G4ThreeVector(-tcshiftx, 0.0, -lxs->TargetChamberZ/2.0));
  G4SubtractionSolid* solidTargetChamberBox = new G4SubtractionSolid("solidTargetChamberBox",
                                                  solidTargetChamberBox2, solidBeamPipeCut, transform2);
  G4LogicalVolume* logicTargetChamberBox = new G4LogicalVolume(solidTargetChamberBox, chamberMaterial, "logicTargetChamberBox");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicTargetChamberBox, "TargetChamberBox",
                                                  logicTargetChamberContainer, false, 0, lxs->OverlapTest);

  G4AssemblyVolume* targetSupportAssembly = ConstructTargetSupportAssembly();
  G4ThreeVector trsupport(-tcshiftx, 0.0, 0.0);
  targetSupportAssembly->MakeImprint(logicTargetChamberContainer, trsupport, 0, 0, lxs->OverlapTest);
  flogicTargetChamber = logicTargetChamberContainer;

}



G4LogicalVolume* LxTargetChamber::GetTargetVolume(const G4String volname)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4String soltvolname = G4String("solid") + volname + "Container";
  G4String logtvolname = G4String("logic") + volname + "Container";
  G4Box* solidChamberTargetVolume = new G4Box(soltvolname, lxs->ChamberTargetVolX/2.0,
                                                lxs->ChamberTargetVolY/2.0, lxs->ChamberTargetVolZ/2.0);
  G4LogicalVolume* logicChamberTargetVolume = new G4LogicalVolume(solidChamberTargetVolume, vacuumMaterial,
                                                                   logtvolname);
  return logicChamberTargetVolume;
}



G4AssemblyVolume* LxTargetChamber::ConstructTargetSupportAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* targetMotorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TargetChamberMotorMaterial);
  G4Material* targetFrameMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TargetFrameMaterial);

  G4Box *solidTCMotorH = new G4Box("solidTCMotorH", lxs->TargetChamberMotorX/2.0,
                                        lxs->TargetChamberMotorY/2.0, lxs->TargetChamberMotorZ/2.0);
  G4LogicalVolume *logicTCMotorH = new G4LogicalVolume(solidTCMotorH, targetMotorMaterial, "logicTCMotorH");

  G4double framex = lxs->ChamberTargetVolX + lxs->TargetFrameD;
  G4double framey = lxs->ChamberTargetVolY + 2.0*lxs->TargetFrameD;

  G4Box *solidTargetFrame1 = new G4Box("solidTargetFrame1", framex/2.0, framey/2.0, lxs->TargetFrameZ/2.0);
  G4Box *solidTargetFrameCut = new G4Box("solidTargetFrameCut", framex/2.0, lxs->ChamberTargetVolY/2.0,
                                                lxs->TargetFrameZ);

  G4Transform3D trframecut(G4RotationMatrix(), G4ThreeVector(-lxs->TargetFrameD, 0.0, 0.0));
  G4SubtractionSolid* solidTargetFrame = new G4SubtractionSolid("solidTargetFrame", solidTargetFrame1,
                                                                      solidTargetFrameCut, trframecut);
  G4LogicalVolume *logicTargetFrame = new G4LogicalVolume(solidTargetFrame, targetFrameMaterial, "logicTargetFrame");

  G4double tfholdy = (lxs->TargetChamberY - framey)/2.0 - lxs->TargetChamberThickness 
                      - lxs->TargetChamberMotorY + lxs->TargetFrameHolderGroveY;
  G4Tubs *solidTFrameHolder1 = new G4Tubs("solidTFrameHolder1", 0.0, lxs->TargetFrameHolderR, tfholdy/2.0, 0.0, 2.0*M_PI);
  G4Box *solidTFrameHolderCut = new G4Box("solidTFrameHolderCut", 2.0*lxs->TargetFrameHolderR,
                                        lxs->TargetFrameZ/2.0, lxs->TargetFrameHolderGroveY);
  G4Transform3D trfholdcut(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, tfholdy/2.0));
  G4SubtractionSolid* solidTFrameHolder = new G4SubtractionSolid("solidTFrameHolder", solidTFrameHolder1,
                                                                      solidTFrameHolderCut, trfholdcut);
  G4LogicalVolume *logicTFrameHolder = new G4LogicalVolume(solidTFrameHolder, targetFrameMaterial, "logicTFrameHolder");

  G4AssemblyVolume *targetFrameAssembly = new G4AssemblyVolume();
  G4ThreeVector traframe(lxs->TargetFrameD/2.0, 0.0, 0.0);
  targetFrameAssembly->AddPlacedVolume(logicTargetFrame, traframe, 0);
  G4ThreeVector trahold(0.0, lxs->TargetFrameHolderGroveY - (framey+tfholdy)/2.0, 0.0);
  targetFrameAssembly->AddPlacedVolume(logicTFrameHolder, trahold,
                                       new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));
  G4double motorx = 0.5*(lxs->TargetChamberMotorX - framex + lxs->TargetFrameD);
  G4ThreeVector tramotor(motorx, trahold.y() - (lxs->TargetChamberMotorY + tfholdy)/2.0, 0.0);
  targetFrameAssembly->AddPlacedVolume(logicTCMotorH, tramotor, 0);

  return targetFrameAssembly;
}



void LxTargetChamber::ConstructTargetChamberSupport(G4LogicalVolume *logicWorld, const G4String tcname, const G4double zpos)
{
// Support
  LXSetUp *lxs = LXSetUp::Instance();
  G4double ypestal = 1.5*m;
  G4double ylevel = lxs->TargetChamberY/2.0;
  G4double tblhight = -ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4double tblx = 1.5*lxs->TargetChamberX;
  G4double tblz = 2.0*lxs->TargetChamberZ;
  G4String tcsupname = tcname + "Support";
  G4AssemblyVolume *tablesupport = LxAux::BuildTable(tcsupname, tblx, tblhight, tblz, 1);
  G4ThreeVector trsupport(lxs->TargetChamberXPos, -ylevel, zpos);
  tablesupport->MakeImprint(logicWorld, trsupport, 0, 0, lxs->OverlapTest);

  G4String tcpedname = tcname + "Pedestal";
  G4LogicalVolume *pedstal = LxAux::BuildPedestal(tcsupname, tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(lxs->TargetChamberXPos, 0.5*ypestal + lxs->FloorSurfaceYpos, zpos), pedstal,
                      tcpedname, logicWorld, false, 0, lxs->OverlapTest);
}



// void LxTargetChamber::CreateMaterial()
// {
//   G4double z, a, density, fractionmass;
//   G4int ncomponents, natoms;
//   G4String symbol;
// }

