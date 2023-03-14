//
/// \brief Implementation of the BeamProfiler tracker class
//

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
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
#include "LxDetector.hh"
#include "LxBSMDetector.hh"
#include "LxAux.hh"



void LxBSMDetector::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  LXSetUp *lxs = LXSetUp::Instance();

  G4double bsmcalox = lxs->BSMCaloX;
  G4double bsmcaloy = lxs->BSMCaloY;
  G4double bsmcalolayerz = lxs->BSMCaloLayerZ;
  G4int    nlayers = lxs->BSMCaloNLayers;

  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* bsmcaloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BSMCaloMaterial);

  G4double bsmcaloz = nlayers * bsmcalolayerz;
  G4Box *solidBSMCaloContainer = new G4Box("solidBSMCaloContainer", bsmcalox/2.0, bsmcaloy/2.0,
                                                            bsmcaloz/2.0);
  G4LogicalVolume *logicBSMCaloContainer = new G4LogicalVolume(solidBSMCaloContainer, environmentMaterial,
                                                               "logicBSMCaloContainer");
  G4Box *solidBSMCaloLayer = new G4Box("solidBSMCaloLayer", bsmcalox/2.0, bsmcaloy/2.0,
                                                            bsmcalolayerz/2.0);
  G4LogicalVolume *logicBSMCaloLayer = new G4LogicalVolume(solidBSMCaloLayer, bsmcaloMaterial,
                                                               "logicBSMCaloLayer");

  for (G4int il = 0; il < nlayers; ++il) {
    G4double lzpos = (il + 0.5) * bsmcalolayerz - 0.5*bsmcaloz;
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lzpos), logicBSMCaloLayer, "BSMCaloLayer",
                     logicBSMCaloContainer, false, il, lxs->OverlapTest);
  }

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->BSMCaloZPos + 0.5*bsmcaloz), logicBSMCaloContainer,
                     "BSMCalo", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double ypestal =  -lxs->FloorSurfaceYpos - bsmcaloy/2.0;
  G4LogicalVolume *pedstal = LxAux::BuildPedestal("BSMCalo", 1.1*bsmcalox, ypestal, 3.0*bsmcaloz);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.5*ypestal + lxs->FloorSurfaceYpos, lxs->BSMCaloZPos + 0.5*bsmcaloz), pedstal,
                      "BSMCaloPedestal", fLogicWorld, false, 0, lxs->OverlapTest);

  AddSegmentation();
}



void LxBSMDetector::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4double bsmcalox = lxs->BSMCaloX;
  G4double bsmcaloy = lxs->BSMCaloY;
  G4int ncellx = lxs->BSMCaloNCellX;
  G4int ncelly = lxs->BSMCaloNCellY;

  fDetector->AddSensorSegmentation("BSMCaloLayer", bsmcalox, bsmcaloy, ncellx, ncelly);
}



G4AssemblyVolume* LxBSMDetector::ConstructSupportAssembly()
{
//   LXSetUp *lxs = LXSetUp::Instance();
  G4AssemblyVolume *supportAssembly = new G4AssemblyVolume();
  return supportAssembly;
}



void LxBSMDetector::CreateMaterial()
{
}

