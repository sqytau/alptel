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
#include "LxBeamProfiler.hh"


void LxBeamProfiler::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();
  
  CreateMaterial();
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* bpContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerContanerMaterial);
  G4Material* bpBoxMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerContanerWallMaterial);
  
  G4Box *solidBPContainer = new G4Box("solidBeamProfilerContainer", lxs->BeamProfilerContanerX/2.0, 
                                      lxs->BeamProfilerContanerY/2.0, lxs->BeamProfilerContanerZ/2.0);
  G4LogicalVolume *logicBPContainer = new G4LogicalVolume(solidBPContainer,
                                                     bpContainerMaterial, "logicBeamProfilerContainer");

  G4Box *solidBPBox1 = new G4Box("solidBeamProfilerBox1", lxs->BeamProfilerContanerX/2.0, 
                                      lxs->BeamProfilerContanerY/2.0, lxs->BeamProfilerContanerZ/2.0);
  G4double dw = lxs->BeamProfilerContanerD;
  G4Box *solidBPBoxCut = new G4Box("solidBPBoxCut", lxs->BeamProfilerContanerX/2.0 - dw, 
                                      lxs->BeamProfilerContanerY/2.0 - dw, lxs->BeamProfilerContanerZ);
  G4SubtractionSolid* solidBPBox = new G4SubtractionSolid("solidBPBox", solidBPBox1, solidBPBoxCut);
  G4LogicalVolume *logicBPBox = new G4LogicalVolume(solidBPBox, bpBoxMaterial, "logicBeamProfilerBox");
  new G4PVPlacement (0, G4ThreeVector(), logicBPBox,
                      "BeamProfilerBox", logicBPContainer, false, 0, lxs->OverlapTest);

  G4double assbeambottom = lxs->BeamProfilerV408Y + lxs->BeamProfilerMotorGap + lxs->BeamProfilerQ545Y/2.0;
  G4double assypos = assbeambottom + lxs->BeamProfilerContanerD - lxs->BeamProfilerContanerY/2.0;
  G4double assxpos = 0.5 * (lxs->BeamProfilerV408X - lxs->BeamProfilerQ545X);
// Motors assembly
  G4AssemblyVolume *motorAssembly = ConstructMotorsAssembly();
  G4ThreeVector trmtr(assxpos, assypos, 0.0);
  motorAssembly->MakeImprint(logicBPContainer, trmtr, 0, 0, lxs->OverlapTest);
  
// PCB assembly
  G4AssemblyVolume* pcbAssembly = ConstructPCBAssembly();
  G4double gapx = lxs->BeamProfilerPCBX - lxs->BeamProfilerPCBHolderX;
  G4double pcbdx = 0.5*(lxs->BeamProfilerQ545X + lxs->BeamProfilerPCBX) + gapx - assxpos;
  G4ThreeVector trsenspcb(-pcbdx, assypos, -lxs->BeamProfilerPCBHolderZ);
  pcbAssembly->MakeImprint(logicBPContainer, trsenspcb, 0, 0, lxs->OverlapTest);
  trsenspcb.setZ(-lxs->BeamProfilerPCBHolderZ + lxs->BeamProfilerZpos[1] - lxs->BeamProfilerZpos[0]);
  pcbAssembly->MakeImprint(logicBPContainer, trsenspcb, 0, 0, lxs->OverlapTest);

//Add front and rear pannel with windows
  G4double bpfrontx = lxs->BeamProfilerContanerX - 2.0*dw;
  G4double bpfronty = lxs->BeamProfilerContanerY - 2.0*dw;
  G4Box *solidBPBoxFront1 = new G4Box("solidBPBoxFront1", bpfrontx/2.0, bpfronty/2.0, dw/2.0);
  G4Box *solidBPBoxFrontCut = new G4Box("solidBPBoxFrontCut", lxs->BeamProfilerWindowX/2.0, lxs->BeamProfilerWindowY/2.0, dw);
  G4Transform3D trfrontcut(G4RotationMatrix(), G4ThreeVector(-pcbdx, assypos, 0.0));
  G4SubtractionSolid* solidBPBoxFront = new G4SubtractionSolid("solidBPBoxFront", solidBPBoxFront1,
                                                                solidBPBoxFrontCut, trfrontcut);
  G4LogicalVolume *logicBPBoxFront = new G4LogicalVolume(solidBPBoxFront, bpBoxMaterial, "logicBPBoxFront");
  G4double frontzpos = 0.5 * (lxs->BeamProfilerContanerZ - dw) - lxs->BeamProfilerWindowZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -frontzpos), logicBPBoxFront,
                      "BeamProfilerBoxFront", logicBPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, frontzpos), logicBPBoxFront,
                      "BeamProfilerBoxFront", logicBPContainer, false, 1, lxs->OverlapTest);
//Add windows  
  G4Material* bpWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerWindowMaterial);
  G4Box *solidBPBoxWindow = new G4Box("solidBPBoxWindow", 1.1*lxs->BeamProfilerWindowX/2.0,
                                      1.1*lxs->BeamProfilerWindowY/2.0, lxs->BeamProfilerWindowZ/2.0);
  G4LogicalVolume *logicBPBoxWindow = new G4LogicalVolume(solidBPBoxWindow, bpWindowMaterial, "logicBPBoxWindow");
  G4double windowzpos = 0.5 * (lxs->BeamProfilerContanerZ - lxs->BeamProfilerWindowZ);
  new G4PVPlacement (0, G4ThreeVector(-pcbdx, assypos, -windowzpos), logicBPBoxWindow,
                      "BeamProfilerBoxWindow", logicBPContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-pcbdx, assypos, windowzpos), logicBPBoxWindow,
                      "BeamProfilerBoxWindow", logicBPContainer, false, 1, lxs->OverlapTest);
  
// PCB holder assembly
  G4AssemblyVolume* pcbHolderAssembly = ConstructPCBSupportAssembly();
  G4ThreeVector trholdpcb(-pcbdx, assypos, 0.0);
  pcbHolderAssembly->MakeImprint(logicBPContainer, trholdpcb, 0, 0, lxs->OverlapTest);

  const std::vector<G4double> &vvz = lxs->BeamProfilerZpos;
  G4double zposmean = std::accumulate(vvz.begin(), vvz.end(), 0.0)/static_cast<G4double>(vvz.size());
  new G4PVPlacement (0, G4ThreeVector(pcbdx, -assypos, zposmean), logicBPContainer,
                      "BeamProfilerContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  AddSegmentation();
  
// Support
  G4double ypestal = 1.5*m;
  G4double ylevel = assbeambottom + lxs->BeamProfilerContanerD;
  G4double tblhight = -ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4double tblx = 3.0*lxs->BeamProfilerContanerX;
  G4double tblz = 2.0*lxs->BeamProfilerContanerZ;
  G4AssemblyVolume *tablesupport = LxAux::BuildTable("ProfilerTable", tblx, tblhight, tblz, 3);
  G4ThreeVector trsupport(0.0, -ylevel, zposmean);
  tablesupport->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);

  G4LogicalVolume *pedstal = LxAux::BuildPedestal("Profiler", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.5*ypestal + lxs->FloorSurfaceYpos, zposmean), pedstal,
                      "ProfilerPedestal", fLogicWorld, false, 0, lxs->OverlapTest);
}



G4AssemblyVolume* LxBeamProfiler::ConstructPCBAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* bpPCBMaterial = G4NistManager::Instance()->FindOrBuildMaterial("FR4");
  G4Material* bpSensorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerSensorMaterial);
  G4Material* bpEnvMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerContanerMaterial);
  G4Material* bpSensMetalMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerSensorMetalization);

  G4Box *solidBPPCB1 = new G4Box("solidBPPCB1", lxs->BeamProfilerPCBX/2.0,
                                lxs->BeamProfilerPCBY/2.0, lxs->BeamProfilerPCBZ/2.0);
  G4Box *solidBPPCBCut = new G4Box("solidBPPCBCut", lxs->BeamProfilerPCBCutX/2.0,
                                lxs->BeamProfilerPCBCutY/2.0, lxs->BeamProfilerPCBZ);
  G4SubtractionSolid* solidBPPCB = new G4SubtractionSolid("solidBPPCB", solidBPPCB1, solidBPPCBCut);
  G4LogicalVolume *logicBPPCB = new G4LogicalVolume(solidBPPCB, bpPCBMaterial, "logicBPPCB");
  
  G4double bpsensorcontainerz = lxs->BeamProfilerSensorZ + lxs->BeamProfilerMetalizationZ;
  G4Box *solidBPSensorContainer = new G4Box("solidBPSensorContainer", lxs->BeamProfilerSensorX/2.0,
                                   lxs->BeamProfilerSensorY/2.0, bpsensorcontainerz/2.0);
  G4LogicalVolume *logicBPSensorContainer = new G4LogicalVolume(solidBPSensorContainer, bpEnvMaterial,
                                                                "logicBPSensorContainer");

  G4Box *solidBPSensor = new G4Box("solidBPSensor", lxs->BeamProfilerSensorX/2.0,
                                   lxs->BeamProfilerSensorY/2.0, lxs->BeamProfilerSensorZ/2.0);
  G4LogicalVolume *logicBPSensor = new G4LogicalVolume(solidBPSensor, bpSensorMaterial, "logicBPSensor");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->BeamProfilerMetalizationZ/2.0), logicBPSensor,
                      "BeamProfilerSensor", logicBPSensorContainer, false, 0, lxs->OverlapTest);

  G4Box *solidBPSensorMetal = new G4Box("solidBPSensorMetal", lxs->BeamProfilerSensorX/2.0,
                                   lxs->BeamProfilerSensorY/2.0, lxs->BeamProfilerMetalizationZ/2.0);
  G4LogicalVolume *logicBPSensorMetal = new G4LogicalVolume(solidBPSensorMetal, bpSensMetalMaterial,
                                                            "logicBPSensorMetal");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -0.5*(bpsensorcontainerz-lxs->BeamProfilerMetalizationZ)),
               logicBPSensorMetal, "BeamProfilerSensorMetal", logicBPSensorContainer, false, 0, lxs->OverlapTest);
  
  G4AssemblyVolume *bpPCBAssembly = new G4AssemblyVolume();

  G4double pcbdz = lxs->BeamProfilerPCBZ/2.0;
  G4ThreeVector trpcb(0.0, 0.0, -pcbdz);
  bpPCBAssembly->AddPlacedVolume(logicBPPCB, trpcb, 0);
  G4ThreeVector trbpsensor(0.0, 0.0, trpcb.z() - 0.5 *(lxs->BeamProfilerPCBZ + bpsensorcontainerz) );
  bpPCBAssembly->AddPlacedVolume(logicBPSensorContainer, trbpsensor, 0);
  
  return bpPCBAssembly;
}



G4AssemblyVolume* LxBeamProfiler::ConstructPCBSupportAssembly()
{ 
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* bpPCBHoldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerPCBHolderMaterial);

  G4Box *solidBPPCBHold1 = new G4Box("solidBPPCBHold1", lxs->BeamProfilerPCBHolderX/2.0,
                                lxs->BeamProfilerPCBHolderY/2.0, lxs->BeamProfilerPCBHolderZ/2.0);
  G4Tubs *solidBPPCBHoldCut1 = new G4Tubs("solidBPPCBHoldCut1", 0.0, lxs->BeamProfilerPCBHolderCutR,
                                           lxs->BeamProfilerPCBHolderZ, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidBPPCBHold = new G4SubtractionSolid("solidBPPCBHold", solidBPPCBHold1, solidBPPCBHoldCut1);
  G4LogicalVolume *logicBPPCBHold = new G4LogicalVolume(solidBPPCBHold, bpPCBHoldMaterial, "logicBPPCBHold");

  G4double hdimx = lxs->BeamProfilerQ545X + lxs->BeamProfilerPCBX - lxs->BeamProfilerPCBHolderX;
  G4double hdimy = 0.9*lxs->BeamProfilerPCBHolderY;
  G4Box *solidBPPCBHoldSide1 = new G4Box("solidBPPCBHoldSide1", hdimx/2.0, hdimy/2.0, lxs->BeamProfilerPCBHolderZ/2.0);
  hdimy -= 20.0*mm;
  G4Box *solidBPPCBHoldSideCut = new G4Box("solidBPPCBHoldSideCut", hdimx/2.0, hdimy/2.0, lxs->BeamProfilerPCBHolderZ);
  
  G4double hdxcut = 10.0*mm;
  G4Transform3D transform1(G4RotationMatrix(), G4ThreeVector(-hdxcut, 0.0, 0.0));
  G4SubtractionSolid* solidBPPCBHoldSide = new G4SubtractionSolid("solidBPPCBHoldSide",
                                                  solidBPPCBHoldSide1, solidBPPCBHoldSideCut, transform1);
  G4LogicalVolume *logicBPPCBHoldSide = new G4LogicalVolume(solidBPPCBHoldSide, bpPCBHoldMaterial, "logicBPPCBHoldSide");

  G4AssemblyVolume *bpPCBSupprtAssembly = new G4AssemblyVolume();

  G4ThreeVector trbppcbhold(0.0, 0.0, -lxs->BeamProfilerPCBHolderZ/2.0);
  bpPCBSupprtAssembly->AddPlacedVolume(logicBPPCBHold, trbppcbhold, 0);
  G4double trpcbholdsidex = 0.5*(hdimx + lxs->BeamProfilerPCBHolderX);
  G4ThreeVector trbppcbholdside(trpcbholdsidex, 0.0, -lxs->BeamProfilerPCBHolderZ/2.0);
  bpPCBSupprtAssembly->AddPlacedVolume(logicBPPCBHoldSide, trbppcbholdside, 0);
  
  return bpPCBSupprtAssembly;
}



G4AssemblyVolume* LxBeamProfiler::ConstructMotorsAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* bpV408Material = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerV408Material);
  G4Material* bpQ545Material = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerQ545Material);
  G4Material* bpaMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamProfilerAngleSupportMaterial);

  G4AssemblyVolume *motorAssembly = new G4AssemblyVolume();
  
// Motors
  G4Box *solidV408 = new G4Box("solidV408", lxs->BeamProfilerV408X/2.0, lxs->BeamProfilerV408Y/2.0, 
                               lxs->BeamProfilerV408Z/2.0);
  G4LogicalVolume *logicV408 = new G4LogicalVolume(solidV408, bpV408Material, "logicV408");
  
  G4Box *solidQ545 = new G4Box("solidV408", lxs->BeamProfilerQ545X/2.0, lxs->BeamProfilerQ545Y/2.0, 
                               lxs->BeamProfilerQ545Z/2.0);
  G4LogicalVolume *logicQ545 = new G4LogicalVolume(solidQ545, bpQ545Material, "logicQ545");
  
  G4ThreeVector trmy(0.0, 0.0, 0.5*lxs->BeamProfilerQ545Z);
  motorAssembly->AddPlacedVolume(logicQ545, trmy, 0);
  G4double mdx = 0.5*(lxs->BeamProfilerV408X - lxs->BeamProfilerQ545X);
  G4double mdy = lxs->BeamProfilerMotorGap + 0.5*(lxs->BeamProfilerV408Y+lxs->BeamProfilerQ545Y);
  G4double mdz = 0.5*lxs->BeamProfilerQ545Z;
  G4ThreeVector trmx(-mdx, -mdy, mdz);
  motorAssembly->AddPlacedVolume(logicV408, trmx, 0);
  
// Angular holder
  G4double elangledx = lxs->BeamProfilerAngleSupportX;
  G4double elangledy = lxs->BeamProfilerMotorGap + lxs->BeamProfilerQ545Y;
  G4double elanglez = 0.5*(lxs->BeamProfilerV408Z - lxs->BeamProfilerQ545Z);
  G4double eladr = lxs->BeamProfilerAngleSupportD;
  
  G4Box *solidBPAngSup1 = new G4Box("solidBPAngSup1", elangledx/2.0, elangledy/2.0, elanglez/2.0);
  G4Box *solidBPAngSupCut1 = new G4Box("solidBPAngSupCut1", elangledx, elangledy/2.0, elanglez/2.0);
  G4Transform3D trans1(G4RotationMatrix(), G4ThreeVector(0.0, eladr, eladr));
  G4SubtractionSolid* solidBPAngSup = new G4SubtractionSolid("solidBPAngSup", solidBPAngSup1, solidBPAngSupCut1, trans1);
  G4LogicalVolume *logicBPAngSup = new G4LogicalVolume(solidBPAngSup, bpaMaterial, "logicBPAngSup");

  G4double asdy = 0.5*(lxs->BeamProfilerQ545Y - elangledy);
  G4double asdz = lxs->BeamProfilerQ545Z + 0.5*elanglez;
  G4ThreeVector trmangsup(0.0, asdy, asdz);
  motorAssembly->AddPlacedVolume(logicBPAngSup, trmangsup, 0);
  
  return motorAssembly;
}
    

    
void LxBeamProfiler::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4double sensx = lxs->BeamProfilerSensorX;
  G4double sensy = lxs->BeamProfilerSensorY;

  fDetector->AddSensorSegmentation("BeamProfilerSensor", sensx, sensy, lxs->BeamProfilerNCellX, lxs->BeamProfilerNCellY);
}



G4AssemblyVolume* LxBeamProfiler::ConstructSupportAssembly()
{
//   LXSetUp *lxs = LXSetUp::Instance();
  G4AssemblyVolume *supportAssembly = new G4AssemblyVolume();
  return supportAssembly;
}



void LxBeamProfiler::CreateMaterial()
{
  G4double z, a, density, fractionmass;
  G4int ncomponents, natoms;
  G4String symbol;

  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silcone"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material *epoxy = new G4Material("BPEpoxy", density= 1.3*g/cm3, ncomponents=3);
  epoxy->AddElement(H, fractionmass=0.1310);
  epoxy->AddElement(C, fractionmass=0.5357);
  epoxy->AddElement(O, fractionmass=0.3333);
// fiber glass
  G4Material* fiberglass = new G4Material( "fiberglass_bp",density=2.61*g/cm3,ncomponents=2);
  fiberglass -> AddElement(Si, natoms=1);
  fiberglass -> AddElement(O, natoms=2);
// PCBoard material FR4
  G4Material *FR4 = new G4Material("FR4",density=1.85*g/cm3,ncomponents=2);
  FR4 ->AddMaterial( epoxy, fractionmass=0.39);
  FR4 ->AddMaterial( fiberglass, fractionmass=0.61);

}

