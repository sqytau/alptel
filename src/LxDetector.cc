//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Sphere.hh"
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


////////////////////////////////////////////////////////////////////////
///// WISDetectorTele

void WISDetectorTele::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();
  CreateMaterial();

  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* collimatorAlMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material* collimatorLeadMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

  G4double colimalx = lxs->Collimator1X;
  G4double colimaly = lxs->Collimator1Y;
  G4double colimalz = lxs->Collimator1Z;
  G4double colimalholer = lxs->Collimator1HoleR;

  G4double scontrin = lxs->SourceContainerRin;
  G4double scontrout = lxs->SourceContainerRout;
  G4double scontz = lxs->SourceContainerZ;

  G4double colimleadz = lxs->CollimatorLeadZ;
  G4double colimleadholer = lxs->CollimatorLeadHoleR;

  G4Box *solidAlCollimator0 = new G4Box("solidAlCollimator0", colimalx/2.0, colimaly/2.0, colimalz/2.0);
  G4Tubs *solidAlCollHole = new G4Tubs("solidAlCollHole", 0.0, colimalholer, colimalz, 0.0, 2.0*M_PI);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid *solidAlCollimator = new G4SubtractionSolid("solidAlCollimator",
                                 solidAlCollimator0, solidAlCollHole, transform);
  G4LogicalVolume *logicAlCollimator = new G4LogicalVolume(solidAlCollimator, collimatorAlMaterial, "logicAlCollimator");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0),
                     logicAlCollimator, "AlCollimator", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidLeadCollimator = new G4Tubs("solidLeadCollimator", colimleadholer, scontrout, colimleadz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicLeadCollimator = new G4LogicalVolume(solidLeadCollimator, collimatorLeadMaterial, "logicLeadCollimator");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -0.5*(colimalz+colimleadz)),
                     logicLeadCollimator, "LeadCollimator", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidSrContainer1 = new G4Tubs("solidSrContainer1", 0.0, scontrout, scontz/2.0, 0.0, 2.0*M_PI);
  G4Tubs *solidSrContainerCut = new G4Tubs("solidSrContainerCut", 0.0, scontrin, scontz/2.0, 0.0, 2.0*M_PI);
  G4Transform3D srcntr(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, lxs->SourceContainerTopZ));
  G4SubtractionSolid *solidSrContainer = new G4SubtractionSolid("solidSrContainer", solidSrContainer1, solidSrContainerCut, srcntr);
  G4LogicalVolume *logicSrContainer = new G4LogicalVolume(solidSrContainer, collimatorLeadMaterial, "logicSrContainer");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -0.5*(colimalz+scontz)-colimleadz),
                     logicSrContainer, "SrContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  G4LogicalVolume *tlSensor = ConstructSensor();
  G4LogicalVolume *tlPCB = ConstructPCB();

  for (size_t si = 0; si < lxs->TelescopeSensorZpos.size(); ++si) {
    G4double zpos = lxs->TelescopeSensorZpos[si] + lxs->OPPPSensorZ/2.0;
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zpos), tlSensor, "AlSensor", fLogicWorld, false, si, lxs->OverlapTest);
    zpos += (lxs->OPPPSensorZ + lxs->CarrierPCBZ)/2.0;
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zpos), tlPCB, "AlPCB", fLogicWorld, false, si, lxs->OverlapTest);
  }

  // Counting volumes in the hole of the Lead Collomator
  G4double countVolZ = 0.1 *mm;
  G4Tubs *solidCountVolume = new G4Tubs("solidCountVolume", 0.0, colimleadholer, countVolZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCountVolume = new G4LogicalVolume(solidCountVolume, environmentMaterial, "logicCountVolume");
  G4double countVolZpos = 0.5*(countVolZ - colimalz) - colimleadz;
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, countVolZpos),
                     logicCountVolume, "CounterVolume", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, countVolZpos + colimleadz - countVolZ),
                     logicCountVolume, "CounterVolume", fLogicWorld, false, 1, lxs->OverlapTest);

  ConstructShielding(fLogicWorld);
  ConstructBottomSupport(fLogicWorld);
  G4LogicalVolume *logicSr90Container = ConstructSr90Sourse();
  G4double sr90posz = -lxs->Sr90GapZ - lxs->CollimatorLeadZ -0.5*(lxs->Collimator1Z + lxs->SrSupportRingZ);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, sr90posz), logicSr90Container, "Sr90SourceAssembly", fLogicWorld, false, 0, lxs->OverlapTest);

  AddSegmentation();

}



void WISDetectorTele::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4double sensx = lxs->OPPPSensorNCellX * lxs->OPPPSensorPixelX;
  G4double sensy = lxs->OPPPSensorNCellY * lxs->OPPPSensorPixelY;

  fDetector->AddSensorSegmentation("OPPPSensitive", sensx, sensy,
                                                    lxs->OPPPSensorNCellX, lxs->OPPPSensorNCellY);
}


G4LogicalVolume* WISDetectorTele::ConstructSensor()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* sensorContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* opppSensorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPTrackerMaterial);

  G4Box *solidAlpideContainer = new G4Box("solidAlpideContainer", lxs->OPPPSensorX/2.0, lxs->OPPPSensorY/2.0, lxs->OPPPSensorZ/2.0);
  G4LogicalVolume *logicAlpideContainer = new G4LogicalVolume(solidAlpideContainer, sensorContainerMaterial,
                                                                "logicAlpideContainer");

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

  G4double lzpos = 0.0;
  G4int idet = 0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lzpos),
                    logicOPPPSensor, "OPPPSensor", logicAlpideContainer, false, idet, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector(0.0, (lxs->OPPPSensorY-sensy)/2.0,
                                       lzpos - (lxs->OPPPSensorZ-lxs->OPPPSensorPixelZ)/2.0),
                    logicOPPPSensitive, "OPPPSensitive", logicAlpideContainer, false, idet, lxs->OverlapTest);

  return logicAlpideContainer;
}



G4LogicalVolume* WISDetectorTele::ConstructPCB()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* bpPCBMaterial = G4NistManager::Instance()->FindOrBuildMaterial("FR4");

  G4Box *solidCarrierPCB1 = new G4Box("solidCarrierPCB1", lxs->CarrierPCBX/2.0,
                                lxs->CarrierPCBY/2.0, lxs->CarrierPCBZ/2.0);
  G4Box *solidCarrierPCBCut = new G4Box("solidCarrierPCBCut", lxs->CarrierPCBCutX/2.0,
                                lxs->CarrierPCBCutY/2.0, lxs->CarrierPCBZ);
  G4SubtractionSolid* solidCarrierPCB = new G4SubtractionSolid("solidCarrierPCB", solidCarrierPCB1, solidCarrierPCBCut);
  G4LogicalVolume *logicCarrierPCB = new G4LogicalVolume(solidCarrierPCB, bpPCBMaterial, "logicCarrierPCB");

  return logicCarrierPCB;
}


void WISDetectorTele::ConstructShielding(G4LogicalVolume  *logicWorld)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* shieldingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TelescopeShieldingMaterial);

  G4Box *solidShielding1 = new G4Box("solidShielding1", lxs->TeleShieldingX/2.0, lxs->TeleShieldingY/2.0, lxs->TeleShieldingZ/2.0);
  G4double scutx = lxs->TeleShieldingX - 2.0*lxs->TeleShieldingThickess;
  G4double scuty = lxs->TeleShieldingY - 2.0*lxs->TeleShieldingThickess;
  G4Box *solidShieldingCut = new G4Box("solidShieldingCut", scutx/2.0, scuty/2.0, lxs->TeleShieldingZ);
  G4Transform3D trnssens(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid* solidShielding = new G4SubtractionSolid("solidShielding", solidShielding1, solidShieldingCut, trnssens);
  G4LogicalVolume *logicShielding = new G4LogicalVolume(solidShielding, shieldingMaterial, "logicShielding");

  G4double posy = lxs->TeleShieldingY/2.0 - lxs->ShieldingGapY - lxs->TeleShieldingThickess;
  G4double posz = lxs->WoodSupportZpos - lxs->TeleShieldingZ/2.0;
  new G4PVPlacement (0, G4ThreeVector(0.0, posy, posz),
                     logicShielding, "Shielding", logicWorld, false, 0, lxs->OverlapTest);

}



void WISDetectorTele::ConstructBottomSupport(G4LogicalVolume  *logicWorld)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* plywoodMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Plywood");
  G4Material* floorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");

  G4double wsx = 2.0*lxs->TeleShieldingX;
  G4double wsy = 2.0*lxs->TeleShieldingY;
  G4Box *solidWoodSupport = new G4Box("solidWoodSupport", wsx/2.0, wsy/2.0, lxs->WoodSupportZ/2.0);
  G4LogicalVolume *logicWoodSupport = new G4LogicalVolume(solidWoodSupport, plywoodMaterial, "logicWoodSupport");

  G4Box *solidFloor = new G4Box("solidFloor", wsx/2.0, wsy/2.0, lxs->FloorZ/2.0);
  G4LogicalVolume *logicFloor = new G4LogicalVolume(solidFloor, floorMaterial, "logicFloor");

  G4double posy = wsy/4.0 - lxs->ShieldingGapY - lxs->TeleShieldingThickess;
  G4double posz = lxs->WoodSupportZpos + lxs->WoodSupportZ/2.0;
  new G4PVPlacement (0, G4ThreeVector(0.0, posy, posz), logicWoodSupport, "WoodSupport", logicWorld, false, 0, lxs->OverlapTest);
  posz += lxs->FloorGapZ + lxs->FloorZ/2.0;
  new G4PVPlacement (0, G4ThreeVector(0.0, posy, posz), logicFloor, "Floor", logicWorld, false, 0, lxs->OverlapTest);

}



G4LogicalVolume* WISDetectorTele::ConstructSr90Sourse()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* sourceMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Sr");
  G4Material* plugMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material* ringMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);

  G4Tubs *solidSr90Container = new G4Tubs("solidSr90Container", 0.0, lxs->SrSupportRingRout, lxs->SrSupportRingZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSr90Container = new G4LogicalVolume(solidSr90Container, environmentMaterial, "logicSr90Container");

  G4Tubs *solidSr90Source = new G4Tubs("solidSr90Source", 0.0, lxs->Sr90SourceR, lxs->Sr90SourceZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSr90Source = new G4LogicalVolume(solidSr90Source, sourceMaterial, "logicSr90Source");

  G4Tubs *solidSourceAlPlug = new G4Tubs("solidSourceAlPlug", 0.0, lxs->SrAlPlugR, lxs->SrAlPlugZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSourceAlPlug = new G4LogicalVolume(solidSourceAlPlug, plugMaterial, "logicSourceAlPlug");

  G4Tubs *solidSrSupportRing = new G4Tubs("solidSrSupportRing", lxs->SrSupportRingRin, lxs->SrSupportRingRout,
                                          lxs->SrSupportRingZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicSrSupportRing = new G4LogicalVolume(solidSrSupportRing, ringMaterial, "logicSrSupportRing");

  G4double posz = 0.5 * (lxs->SrSupportRingZ - lxs->Sr90SourceZ);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, posz), logicSr90Source, "Sr90Source", logicSr90Container, false, 0, lxs->OverlapTest);
  posz -= 0.5 * (lxs->SrAlPlugZ + lxs->Sr90SourceZ);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, posz), logicSourceAlPlug, "SourceAlPlug", logicSr90Container, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicSrSupportRing, "SrSupportRing", logicSr90Container, false, 0, lxs->OverlapTest);

  return logicSr90Container;
//   G4double sr90posz = -lxs->Sr90GapZ - lxs->CollimatorLeadZ -0.5*(lxs->Collimator1Z + lxs->SrSupportRingZ);
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, sr90posz), logicSr90Container, "Sr90SourceAssembly", logicWorld, false, 0, lxs->OverlapTest);

}



void WISDetectorTele::CreateMaterial()
{
  G4double z, a, density, fractionmass;
  G4int ncomponents, natoms;
  G4String symbol;

  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silcone"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
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

  G4Material *woodMaterial = new G4Material("Wood", density=0.5*g/cm3,ncomponents=4);
  woodMaterial->AddElement(C, fractionmass=0.51);
  woodMaterial->AddElement(O, fractionmass=0.42);
  woodMaterial->AddElement(H, fractionmass=0.06);
  woodMaterial->AddElement(N, fractionmass=0.01);

  G4Material *PlywoodMaterial = new G4Material("Plywood", density=0.7*g/cm3,ncomponents=2);
  PlywoodMaterial->AddMaterial(epoxy, fractionmass=13.0/28.0);
  PlywoodMaterial->AddMaterial(woodMaterial, fractionmass=15.0/28.0);

}


////////////////////////////////////////////////////////////////////////
///// WISDetectorTele

void WISDetectorTeleFrame::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();
  CreateMaterial();

  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* collimatorLeadMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");

  G4double scontrin = lxs->SourceContainerRin;
  G4double scontrout = lxs->SourceContainerRout;
  G4double scontz = lxs->SourceContainerZ;

  G4double colimleadz = lxs->CollimatorLeadZ;
  G4double colimleadholer = lxs->CollimatorLeadHoleR;

  G4Tubs *solidLeadCollimator = new G4Tubs("solidLeadCollimator", colimleadholer, scontrout, colimleadz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicLeadCollimator = new G4LogicalVolume(solidLeadCollimator, collimatorLeadMaterial, "logicLeadCollimator");

  G4double colimleadzpos = -lxs->BTargetZ - colimleadz/2.0;
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, colimleadzpos),
                     logicLeadCollimator, "LeadCollimator", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidSrContainer1 = new G4Tubs("solidSrContainer1", 0.0, scontrout, scontz/2.0, 0.0, 2.0*M_PI);
  G4Tubs *solidSrContainerCut = new G4Tubs("solidSrContainerCut", 0.0, scontrin, scontz/2.0, 0.0, 2.0*M_PI);
  G4Transform3D srcntr(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, lxs->SourceContainerTopZ));
  G4SubtractionSolid *solidSrContainer = new G4SubtractionSolid("solidSrContainer", solidSrContainer1, solidSrContainerCut, srcntr);
  G4LogicalVolume *logicSrContainer = new G4LogicalVolume(solidSrContainer, collimatorLeadMaterial, "logicSrContainer");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -lxs->BTargetZ - 0.5*(+scontz)-colimleadz),
                     logicSrContainer, "SrContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  G4LogicalVolume *tlSensor = ConstructSensor();
  G4LogicalVolume *tlPCB = ConstructPCB();
  G4LogicalVolume *tlroPCB = ConstructROPCB();
  G4LogicalVolume *tlFrame = ConstructAlFrame();
  G4double ropcbypos = 0.5*(lxs->ROPCBY + lxs->CarrierPCBY) + lxs->CarrierROGap;
  G4double frameypos = 0.5*(lxs->TeleFrameY - lxs->CarrierPCBY) - lxs->TeleFramePCBShift;

  for (size_t si = 0; si < lxs->TelescopeSensorZpos.size(); ++si) {
    G4double zpos = lxs->TelescopeSensorZpos[si] + lxs->OPPPSensorZ/2.0;
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zpos), tlSensor, "AlSensor", fLogicWorld, false, si, lxs->OverlapTest);
    zpos += (lxs->OPPPSensorZ + lxs->CarrierPCBZ)/2.0;
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zpos), tlPCB, "AlPCB", fLogicWorld, false, si, lxs->OverlapTest);
    new G4PVPlacement (0, G4ThreeVector(0.0, ropcbypos, zpos), tlroPCB, "ROPCB", fLogicWorld, false, si, lxs->OverlapTest);
    G4double framezpos = zpos + lxs->TeleFrameZ/2.0 + lxs->CarrierPCBZ/2.0;
    new G4PVPlacement (0, G4ThreeVector(0.0, frameypos, framezpos), tlFrame, "TeleFrame", fLogicWorld, false, si, lxs->OverlapTest);
  }

  // Counting volumes in the hole of the Lead Collomator
  G4double countVolZ = 0.1 *mm;
  G4Tubs *solidCountVolume = new G4Tubs("solidCountVolume", 0.0, colimleadholer, countVolZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCountVolume = new G4LogicalVolume(solidCountVolume, environmentMaterial, "logicCountVolume");
  G4double countVolZpos = 0.5*countVolZ - lxs->BTargetZ - colimleadz;
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, countVolZpos),
                     logicCountVolume, "CounterVolume", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, countVolZpos + colimleadz - countVolZ),
                     logicCountVolume, "CounterVolume", fLogicWorld, false, 1, lxs->OverlapTest);

  G4LogicalVolume *logicSr90Container = ConstructSr90Sourse();
  G4double sr90posz = colimleadzpos + colimleadz/2.0 - lxs->Sr90GapZ - 0.5*(lxs->SrSupportRingZ);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, sr90posz), logicSr90Container, "Sr90SourceAssembly", fLogicWorld, false, 0, lxs->OverlapTest);

  ConstructBottomSupport(fLogicWorld);
  G4AssemblyVolume* supportAssembly = ConstructSupportAssembly();
  G4ThreeVector trsupport(0.0, -lxs->TeleFramePCBShift - (lxs->CarrierPCBY + lxs->FrameHolderY)/2.0, 0.0);
  supportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);

  AddSegmentation();

}




void WISDetectorTeleFrame::ConstructBottomSupport(G4LogicalVolume  *logicWorld)
{
  LXSetUp *lxs = LXSetUp::Instance();
//   G4Material* plywoodMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Plywood");
  G4Material* floorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");

  G4double wsx = 2.0*lxs->TeleShieldingX;
  G4double wsy = 2.0*lxs->TeleShieldingY;
//   G4Box *solidWoodSupport = new G4Box("solidWoodSupport", wsx/2.0, wsy/2.0, lxs->WoodSupportZ/2.0);
//   G4LogicalVolume *logicWoodSupport = new G4LogicalVolume(solidWoodSupport, plywoodMaterial, "logicWoodSupport");

  G4Box *solidFloor = new G4Box("solidFloor", wsx/2.0, wsy/2.0, lxs->FloorZ/2.0);
  G4LogicalVolume *logicFloor = new G4LogicalVolume(solidFloor, floorMaterial, "logicFloor");

  G4double posy = wsy/4.0 - lxs->ShieldingGapY - lxs->TeleShieldingThickess;
//   G4double posz = lxs->WoodSupportZpos + lxs->WoodSupportZ/2.0;
//   new G4PVPlacement (0, G4ThreeVector(0.0, posy, posz), logicWoodSupport, "WoodSupport", logicWorld, false, 0, lxs->OverlapTest);
  G4double posz = lxs->FloorGapZ + lxs->FloorZ/2.0;
  new G4PVPlacement (0, G4ThreeVector(0.0, posy, posz), logicFloor, "Floor", logicWorld, false, 0, lxs->OverlapTest);

}



G4LogicalVolume* WISDetectorTeleFrame::ConstructROPCB()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* bpPCBMaterial = G4NistManager::Instance()->FindOrBuildMaterial("FR4");

  G4Box *solidROPCB = new G4Box("solidCarrierPCB1", lxs->ROPCBX/2.0,
                                lxs->ROPCBY/2.0, lxs->ROPCBZ/2.0);
  G4LogicalVolume *logicROPCB = new G4LogicalVolume(solidROPCB, bpPCBMaterial, "logicROPCB");

  return logicROPCB;
}



G4LogicalVolume* WISDetectorTeleFrame::ConstructAlFrame()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* frameMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

  G4Box *solidFrame0 = new G4Box("solidFrame0", lxs->TeleFrameX/2.0,
                                lxs->TeleFrameY/2.0, lxs->TeleFrameZ/2.0);

  G4Box *solidFrameCutAlpide = new G4Box("solidFrameCutAlpide", lxs->TeleFrameCutAlpideX/2.0,
                                lxs->TeleFrameCutAlpideY/2.0, lxs->TeleFrameZ);
  G4Box *solidFrameCutROPCB = new G4Box("solidFrameCutROPCB", lxs->ROPCBX/2.0,
                                lxs->ROPCBY/2.0, lxs->TeleFrameZ);

//   G4double cuty = 0.5*(lxs->TeleFrameCutAlpideY-lxs->TeleFrameY) + lxs->TeleFrameCutAlpideGapY;
  G4double cuty = 0.5*(lxs->CarrierPCBY-lxs->TeleFrameY) + lxs->TeleFramePCBShift;
  G4Transform3D alpcbcuttr(G4RotationMatrix(), G4ThreeVector(0.0, cuty, 0.0));
  G4SubtractionSolid *solidFrame1 = new G4SubtractionSolid("solidSrContainer", solidFrame0, solidFrameCutAlpide, alpcbcuttr);
  cuty += lxs->CarrierROGap + (lxs->CarrierPCBY+lxs->ROPCBY)/2.0;
  G4Transform3D rocbcuttr(G4RotationMatrix(), G4ThreeVector(0.0, cuty, 0.0));
  G4SubtractionSolid *solidFrame = new G4SubtractionSolid("solidSrContainer", solidFrame1, solidFrameCutROPCB, rocbcuttr);

  G4LogicalVolume *logicFrame = new G4LogicalVolume(solidFrame, frameMaterial, "logicFrame");

  return logicFrame;
}




G4AssemblyVolume* WISDetectorTeleFrame::ConstructSupportAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* supportMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material* shieldingMaterial = G4NistManager::Instance()->FindOrBuildMaterial("FR4");

  G4AssemblyVolume *teleSupportAssembly = new G4AssemblyVolume();

  G4Box *solidFrameHolder = new G4Box("solidFrameHolder", lxs->FrameHolderX/2.0,
                                lxs->FrameHolderY/2.0, lxs->FrameHolderZ/2.0);
  G4LogicalVolume *logicFrameHolder = new G4LogicalVolume(solidFrameHolder, supportMaterial, "logicFrameHolder");

  G4Box *solidSupportPole0 = new G4Box("solidSupportPole", lxs->SupportPoleX/2.0,
                                lxs->SupportPoleY/2.0, lxs->SupportPoleZ/2.0);
  G4Box *solidSupportPoleCut = new G4Box("solidSupportPoleCut", lxs->SupportPoleX/2.0 - lxs->SupportPoleThickness,
                                lxs->SupportPoleY/2.0 - lxs->SupportPoleThickness, lxs->SupportPoleZ);
  G4SubtractionSolid *solidSupportPole = new G4SubtractionSolid("solidSupportPole", solidSupportPole0, solidSupportPoleCut);
  G4LogicalVolume *logicSupportPole = new G4LogicalVolume(solidSupportPole, supportMaterial, "logicSupportPole");

  G4ThreeVector holdertr(0.0, 0.0, lxs->FrameHolderZ/2.0);
  teleSupportAssembly->AddPlacedVolume(logicFrameHolder, holdertr, 0);
  G4ThreeVector poletr(0.0, -0.5*(lxs->FrameHolderY+lxs->SupportPoleY), lxs->FloorGapZ-lxs->SupportPoleZ/2.0);
  teleSupportAssembly->AddPlacedVolume(logicSupportPole, poletr, 0);

  G4Box *solidShieldingBoxBottom0 = new G4Box("solidShieldingBoxBottom0", lxs->ShieldingBoxBottomX/2.0,
                                lxs->ShieldingBoxBottomY/2.0, lxs->ShieldingBoxBottomZ/2.0);
  G4SubtractionSolid *solidShieldingBoxBottom = new G4SubtractionSolid("solidShieldingBoxBottom",
                                                                       solidShieldingBoxBottom0, solidFrameHolder);
  G4LogicalVolume *logicShieldingBoxBottom = new G4LogicalVolume(solidShieldingBoxBottom, shieldingMaterial, "logicShieldingBoxBottom");
  G4ThreeVector shielbottomtr(0.0, 0.5*(lxs->ShieldingBoxBottomY-lxs->FrameHolderY), lxs->FrameHolderZ/2.0);
  teleSupportAssembly->AddPlacedVolume(logicShieldingBoxBottom, shielbottomtr, 0);

  return teleSupportAssembly;
}
