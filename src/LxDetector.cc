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
//// LxDetectorOPPP

void LxDetectorOPPP::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();
  
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* opppDetContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPDetMaterial);
  G4Material* trackerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPTrackerMaterial);
  G4Material* caloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPCaloMaterial);

  G4double dumpMagnetZ = lxs->TypMBFieldLength;
  G4double detzpos = lxs->IPMagnetZpos + dumpMagnetZ/2.0 + lxs->OPPPDetZtoMagnet + lxs->OPPPDetZ/2.0;

  G4Box *OpppDetContainer = new G4Box("solidOpppDetContainer", lxs->OPPPDetX/2.0, lxs->OPPPDetY/2.0, lxs->OPPPDetZ/2.0);
  G4LogicalVolume *logicDetContainer = new G4LogicalVolume(OpppDetContainer, opppDetContainerMaterial, "logicDetContainer");

  G4Box *solidOpppTracker = new G4Box("solidOpppTracker", lxs->OPPPTrackerX/2.0, lxs->OPPPTrackerY/2.0, lxs->OPPPTrackerZ/2.0);
  G4LogicalVolume *logicOpppTracker = new G4LogicalVolume(solidOpppTracker, trackerMaterial, "logicOpppTracker");
  for (G4int il = 0; il < lxs->OPPPTrackerNLayers; ++il) {
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),  // Let it face backward toward beam
                      G4ThreeVector(0.0, 0.0, (-lxs->OPPPDetZ+lxs->OPPPTrackerZ)/2.0 + il*(lxs->OPPPTrackerZ+lxs->OPPPTrackerDZ)), 
                      logicOpppTracker, "OpppTracker", logicDetContainer, false, il, lxs->OverlapTest);
  }
  
  G4Box *solidOpppCalo = new G4Box("solidOpppCalo", lxs->OPPPCaloX/2.0, lxs->OPPPCaloY/2.0, lxs->OPPPCaloZ/2.0);
  G4LogicalVolume *logicOpppCalo = new G4LogicalVolume(solidOpppCalo, caloMaterial, "logicOpppCalo");
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(0.0, 0.0, (lxs->OPPPDetZ - lxs->OPPPCaloZ)/2.0), 
                    logicOpppCalo, "OpppCalo", logicDetContainer, false, 0, lxs->OverlapTest);
  
  G4double detxpos = lxs->OPPPDetX/2.0 + lxs->OPPPDetXPos;
//   G4double bpipel = lxs->OPPPDetZtoMagnet;
//   G4double thetad = atan2(detxpos, bpipel);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), 0.0), //thetad), 
                    G4ThreeVector(detxpos, 0.0, detzpos), 
                    logicDetContainer, "OpppDetContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), 0.0),  // thetad), 
                    G4ThreeVector(-detxpos, 0.0, detzpos), 
                    logicDetContainer, "OpppDetContainer", fLogicWorld, false, 1, lxs->OverlapTest);
  
  AddSegmentation();
}



void LxDetectorOPPP::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("OpppTracker", lxs->OPPPTrackerX, lxs->OPPPTrackerY, 
                                                  lxs->OPPPTrackerNCellX, lxs->OPPPTrackerNCellY);
  fDetector->AddSensorSegmentation("OpppCalo", lxs->OPPPCaloX, lxs->OPPPCaloY, 200, 20);
}


////////////////////////////////////////////////////////////////////////
/////LxDetectorCompton

void LxDetectorCompton::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();
  
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* ComptonDetContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonDetMaterial);
  G4Material* bpipeWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeWindowMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* trackerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonTrackerMaterial);
  G4Material* caloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonCaloMaterial);

  G4Box *solidComptonDetContainer = new G4Box("solidComptonDetContainer", lxs->ComptonDetX/2.0, lxs->ComptonDetY/2.0, lxs->ComptonDetZ/2.0);
  G4LogicalVolume *logicComptonDetContainer = new G4LogicalVolume(solidComptonDetContainer, 
                                                                  ComptonDetContainerMaterial, "logicComptonDetContainer");

  G4Box *solidComptonTracker = new G4Box("solidComptonTracker", lxs->ComptonTrackerX/2.0, lxs->ComptonTrackerY/2.0, lxs->ComptonTrackerZ/2.0);
  G4LogicalVolume *logicComptonTracker = new G4LogicalVolume(solidComptonTracker, trackerMaterial, "logicComptonTracker");
  for (G4int il = 0; il < lxs->ComptonTrackerNLayers; ++il) {
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),  // Let it face backward toward beam
                      G4ThreeVector(0.0, 0.0, (-lxs->ComptonDetZ+lxs->ComptonTrackerZ)/2.0 + il*(lxs->ComptonTrackerZ+lxs->ComptonTrackerDZ)), 
                      logicComptonTracker, "ComptonTracker", logicComptonDetContainer, false, il, lxs->OverlapTest);
  }
  
  G4Box *solidComptonCalo = new G4Box("solidComptonCalo", lxs->ComptonCaloX/2.0, lxs->ComptonCaloY/2.0, lxs->ComptonCaloZ/2.0);
  G4LogicalVolume *logicComptonCalo = new G4LogicalVolume(solidComptonCalo, caloMaterial, "logicComptonCalo");
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(0.0, 0.0, (lxs->ComptonDetZ - lxs->ComptonCaloZ)/2.0), 
                    logicComptonCalo, "ComptonCalo", logicComptonDetContainer, false, 0, lxs->OverlapTest);

  G4double bpipel = lxs->ComptonDetZpos - lxs->GMagnetZpos - lxs->GMagnetZ/2.0;
  G4double thetad = atan2(lxs->ComptonDetX/2.0 + 4.0*lxs->BPipeR, bpipel);
 // new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), thetad),
    new G4PVPlacement(0, 
                    G4ThreeVector((bpipel+lxs->ComptonDetZ/2.0)*sin(thetad), 0.0, 
                                  (bpipel+lxs->ComptonDetZ/2.0)*cos(thetad) + lxs->GMagnetZpos + lxs->GMagnetZ/2.0), 
                    logicComptonDetContainer, "ComptonDetContainer", fLogicWorld, false, 0, lxs->OverlapTest);
//   
//  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), thetad), 
  new G4PVPlacement(0, 
                    G4ThreeVector(-(bpipel+lxs->ComptonDetZ/2.0)*sin(thetad), 0.0, 
                                  (bpipel+lxs->ComptonDetZ/2.0)*cos(thetad) + lxs->GMagnetZpos + lxs->GMagnetZ/2.0), 
                    logicComptonDetContainer, "ComptonDetContainer", fLogicWorld, false, 1, lxs->OverlapTest);

  G4Tubs *solidBPipeComptonContainer = new G4Tubs("solidBPipeComptonContainer", 0.0, lxs->BPipeR, bpipel/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBPipeComptonContainer = new G4LogicalVolume(solidBPipeComptonContainer, vacuumMaterial, "logicBPipeComptonContainer");
  
  // Beam pipe top and bottom part 
  G4Tubs *solidBPipeComptonD = new G4Tubs("solidBPipeComptonD", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                    bpipel/2.0, lxs->DMBPipeWindow/2.0, M_PI-lxs->DMBPipeWindow);
  G4LogicalVolume *logicBPipeComptonD = new G4LogicalVolume(solidBPipeComptonD, beamPipeMaterial, "logicBPipeComptonD");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicBPipeComptonD, "BPipeComptonD", logicBPipeComptonContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0), 
                      logicBPipeComptonD, "BPipeComptonD", logicBPipeComptonContainer, false, 1, lxs->OverlapTest); 
  
  // Beam pipe windows left and right sides 
  G4Tubs *solidBPipeWindowComptonD = new G4Tubs("solidBPipeWindowComptonD", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                          bpipel/2.0, -lxs->DMBPipeWindow/2.0, lxs->DMBPipeWindow);
  G4LogicalVolume *logicBPipeWindowComptonD = new G4LogicalVolume(solidBPipeWindowComptonD, bpipeWindowMaterial, "logicBPipeWindowComptonD");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicBPipeWindowComptonD, "BPipeWindowD", logicBPipeComptonContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0), 
                      logicBPipeWindowComptonD, "BPipeWindowD", logicBPipeComptonContainer, false, 1, lxs->OverlapTest); 
  
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0,  lxs->GMagnetZpos + lxs->GMagnetZ/2.0 + bpipel/2.0), 
                     logicBPipeComptonContainer, "BPipeComptonDAssembly", fLogicWorld, false, 0, lxs->OverlapTest);
  
  AddSegmentation();
  
}



void LxDetectorCompton::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("ComptonTracker", lxs->ComptonTrackerX, lxs->ComptonTrackerY, 
                                                     lxs->ComptonTrackerNCellX, lxs->ComptonTrackerNCellY);
}


////////////////////////////////////////////////////////////////////////
/////LxDetectorCompton

void LxDetectorComptonFluka::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
//   G4Material* bpipeWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeWindowMaterial);

// Collimators
  G4Material* CollimatorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->CollimatorMaterial);
   //**************** cylindrically shaped Collimators placed between  target & Compton Magnet
  G4Tubs *solidCollimator = new G4Tubs("solidCollimator", lxs->CollimatorRin, lxs->CollimatorRout,
                                       lxs->CollimatorZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCollimator = new G4LogicalVolume(solidCollimator, CollimatorMaterial, "logicCollimator");
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->CollimatorZpos),
//                      logicCollimator, "Collimator", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->CollimatorZpos + lxs->CollimatorZ + 50.0*cm),
                     logicCollimator, "Collimator", fLogicWorld, false, 1, lxs->OverlapTest);

// // Collimators support
//   G4Material* collimatorSupportMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->CollimatorSupportMaterial);
//   G4double colsupx = 4.0*lxs->CollimatorRout;
//   G4double colsupy = -lxs->FloorSurfaceYpos - lxs->CollimatorRout/2.0;
//   G4double colsupz = 0.75*lxs->CollimatorZ;
//   G4Box *solidCollimatorSupport0 = new G4Box("solidCollimatorSupport0", colsupx/2.0, colsupy/2.0, colsupz/2.0);
//   G4Tubs *solidCollSubtr = new G4Tubs("solidCollSubtr", 0.0, lxs->CollimatorRout, lxs->CollimatorZ/2.0, 0.0, 2.0*M_PI);
//   G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(0.0, (colsupy + lxs->CollimatorRout)/2.0, 0.0));
//   G4SubtractionSolid *solidCollimatorSupport = new G4SubtractionSolid("solidCollimatorSupport",
//                                  solidCollimatorSupport0, solidCollSubtr, transform);
//   G4LogicalVolume *logicCollimatorSupport = new G4LogicalVolume(solidCollimatorSupport,
//                                         collimatorSupportMaterial, "logicCollimatorSupport");
//   new G4PVPlacement(0, G4ThreeVector(0.0, lxs->FloorSurfaceYpos+colsupy/2.0, lxs->CollimatorZpos + lxs->CollimatorZ + 50.0*cm),
//                     logicCollimatorSupport, "CollimatorSupport", fLogicWorld, false, 0, lxs->OverlapTest);

//LYSO calorimeters
 G4Material* LysoCalMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonLysoMaterial);
 G4Box *solidComptonLyso = new G4Box("solidComptonLyso", lxs->ComptonLysoX/2.0, lxs->ComptonLysoY/2.0, lxs->ComptonLysoZ/2.0);
 G4LogicalVolume *logicComptonLysoCal = new G4LogicalVolume(solidComptonLyso,
                                                                  LysoCalMaterial, "logicComptonLysoCal");
  new G4PVPlacement (0, G4ThreeVector(lxs->ComptonLysoXpos, 0.0, lxs->ComptonLysoZpos),
                     logicComptonLysoCal, "LysoCal", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-lxs->ComptonLysoXpos, 0.0, lxs->ComptonLysoZpos),
                     logicComptonLysoCal, "LysoCal", fLogicWorld, false, 1, lxs->OverlapTest);

//Beam Pipe Target - 1st collimator
  G4double pos_2ndColl = lxs->CollimatorZpos + lxs->CollimatorZ + 50.0*cm;
  G4double fHICSDumpAngle = atan2(lxs->HICSDumpFrontXPos, lxs->HICSDumpFrontZPos - lxs->IPMagnetZpos);
  G4double frontZpos = lxs->HICSDumpFrontZPos  + lxs->HICSDumpFrontXPos*tan(fHICSDumpAngle);
//  G4double lpipe_gd = lxs->CollimatorZpos - lxs->GTargetZpos - lxs->CollimatorZ/2.0 - lxs->GTargetContainerZ/2.0;
//  G4double lpipe_pos_1st = lxs->GTargetZpos + lpipe_gd/2.0 +lxs->GTargetContainerZ/2.0;
  G4double lpipe_gd = frontZpos - lxs->GTargetZpos - lxs->GTargetContainerZ/2.0;
  G4double lpipe_pos_1st = frontZpos - lpipe_gd/2.0;

 G4CutTubs *solidBeamPipeGammaT1stC = new G4CutTubs("solidBeamPipeGammaT1stCC", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG,
                         lpipe_gd/2.0, 0.0, 2.0*M_PI, G4ThreeVector(0.0, 0.0, -1.0), G4ThreeVector(-sin(fHICSDumpAngle), 0.0, cos(fHICSDumpAngle)));
//  G4Tubs *solidBeamPipeGammaT1stC = new G4Tubs("solidBeamPipeGammaT1stC", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG, lpipe_gd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGammaT1stC = new G4LogicalVolume(solidBeamPipeGammaT1stC, beamPipeMaterial, "logicBeamPipeGammaT1stC");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_1st),
                               logicBeamPipeGammaT1stC, "BeamPipeGammaT1stC", fLogicWorld, false, 0, lxs->OverlapTest);
  G4CutTubs *solidBeamPipeGammaT1stCVac = new G4CutTubs("solidBeamPipeGammaT1stCVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                         lpipe_gd/2.0, 0.0, 2.0*M_PI, G4ThreeVector(0.0, 0.0, -1.0), G4ThreeVector(-sin(fHICSDumpAngle), 0.0, cos(fHICSDumpAngle)));
//  G4Tubs *solidBeamPipeGammaT1stCVac = new G4Tubs("solidBeamPipeGammaT1stCVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, lpipe_gd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGammaT1stCVac = new G4LogicalVolume(solidBeamPipeGammaT1stCVac, vacuumMaterial, "logicBeamPipeGammaT1stCVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_1st),
                               logicBeamPipeGammaT1stCVac, "BeamPipeGammaT1stCVac", fLogicWorld, false, 0, lxs->OverlapTest);

  //Beam Pipe HICS beam dump - 2nd collimator

  G4double rearZpos = lxs->HICSDumpFrontZPos + lxs->HICSDumpZ/cos(fHICSDumpAngle) + lxs->HICSDumpFrontXPos*tan(fHICSDumpAngle);
  G4double lpipe_gd2nd = pos_2ndColl - rearZpos - lxs->CollimatorZ/2.0;
  G4double lpipe_pos_2nd = pos_2ndColl - lxs->CollimatorZ/2.0 - lpipe_gd2nd/2;
//  G4double lpipe_gd2nd = pos_2ndColl- lxs->HICSDumpFrontZPos - lxs->HICSDumpZ - lxs->CollimatorZ/2.0;
//  G4double lpipe_pos_2nd = pos_2ndColl - lxs->CollimatorZ/2.0 - lpipe_gd2nd/2;


  G4CutTubs *solidBeamPipeGamma1stC2ndC = new G4CutTubs("solidBeamPipeGamma1stC2ndC", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG,
                         lpipe_gd2nd/2.0, 0.0, 2.0*M_PI, G4ThreeVector(sin(fHICSDumpAngle), 0.0, -cos(fHICSDumpAngle)),
                        G4ThreeVector(0.0, 0.0, 1.0));

//  G4Tubs *solidBeamPipeGamma1stC2ndC = new G4Tubs("solidBeamPipeGamma1stC2ndC", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG,
//                                                  lpipe_gd2nd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma1stC2ndC = new G4LogicalVolume(solidBeamPipeGamma1stC2ndC, beamPipeMaterial, "logicBeamPipeGamma1stC2ndC");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_2nd),
                               logicBeamPipeGamma1stC2ndC, "BeamPipeGamma1stC2ndC", fLogicWorld, false, 0, lxs->OverlapTest);
  G4CutTubs *solidBeamPipeGamma1stC2ndCVac = new G4CutTubs("solidBeamPipeGamma1stC2ndCVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                         lpipe_gd2nd/2.0, 0.0, 2.0*M_PI, G4ThreeVector(sin(fHICSDumpAngle), 0.0, -cos(fHICSDumpAngle)),
                        G4ThreeVector(0.0, 0.0, 1.0));
//  G4Tubs *solidBeamPipeGamma1stC2ndCVac = new G4Tubs("solidBeamPipeGamma1stC2ndCVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, lpipe_gd2nd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma1stC2ndCVac = new G4LogicalVolume(solidBeamPipeGamma1stC2ndCVac, vacuumMaterial, "logicBeamPipeGamma1stC2ndCVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_2nd),
                               logicBeamPipeGamma1stC2ndCVac, "BeamPipeGamma1stC2ndCVac", fLogicWorld, false, 0, lxs->OverlapTest);

//Beam Pipe 1st collimator - 2nd collimator
/*  G4double lpipe_gd2nd = lxs->CollimatorZpos - lxs->GTargetZpos - lxs->CollimatorZ/2.0;
  G4double lpipe_pos_2nd = lxs->CollimatorZpos + lxs->CollimatorZ/2.0 + lpipe_gd2nd/2;
   G4Tubs *solidBeamPipeGamma1stC2ndC = new G4Tubs("solidBeamPipeGamma1stC2ndC", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG,
                                                    lpipe_gd2nd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma1stC2ndC = new G4LogicalVolume(solidBeamPipeGamma1stC2ndC, beamPipeMaterial, "logicBeamPipeGamma1stC2ndC");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_2nd),
                               logicBeamPipeGamma1stC2ndC, "BeamPipeGamma1stC2ndC", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeGamma1stC2ndCVac = new G4Tubs("solidBeamPipeGamma1stC2ndCVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, lpipe_gd2nd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma1stC2ndCVac = new G4LogicalVolume(solidBeamPipeGamma1stC2ndCVac, vacuumMaterial, "logicBeamPipeGamma1stC2ndCVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos_2nd),
                  logicBeamPipeGamma1stC2ndCVac, "BeamPipeGamma1stC2ndCVac", fLogicWorld, false, 0, lxs->OverlapTest);*/

 //round  Beam Pipe  2nd collimator - Magnet Field
//   G4double lpipe_gd3d = lxs->GMagnetZpos - lxs->GMFieldLength/2.0 - lpipe_pos_2nd - lxs->CollimatorZ - lpipe_gd2nd/2;
//   G4double lpipe_pos3d = lxs->GMagnetZpos - lxs->GMFieldLength/2.0 - lpipe_gd3d/2;
  G4double lpipe_gd3d = lxs->GMagnetZpos - lxs->FlashMFieldLength/2.0 - lpipe_pos_2nd - lxs->CollimatorZ - lpipe_gd2nd/2;
  G4double lpipe_pos3d = lxs->GMagnetZpos - lxs->FlashMFieldLength/2.0 - lpipe_gd3d/2;

  G4Tubs *solidBeamPipeGamma2ndCMagF = new G4Tubs("solidBeamPipeGamma2ndCMagF", lxs->BPipeRLG-lxs->BPipeThickness,
                                                  lxs->BPipeRLG, lpipe_gd3d/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma2ndCMagF = new G4LogicalVolume(solidBeamPipeGamma2ndCMagF, beamPipeMaterial,
                                                                    "logicBeamPipeGamma2ndCMagF");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos3d),
                               logicBeamPipeGamma2ndCMagF, "BeamPipeGamma2ndCMagF", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeGamma2ndCMagFVac = new G4Tubs("solidBeamPipeGamma2ndCMagFVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                                    lpipe_gd3d/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGamma2ndCMagFVac = new G4LogicalVolume(solidBeamPipeGamma2ndCMagFVac, vacuumMaterial,
                                                                       "logicBeamPipeGamma2ndCMagFVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos3d),
                       logicBeamPipeGamma2ndCMagFVac, "BeamPipeGamma2ndCMagFVac", fLogicWorld, false, 0, lxs->OverlapTest);

//quadratish  Beam Pipe  2nd collimator - Magnet Field
/*
////////
  G4double gmpipex = lxs->QBeamPipeContainerX;  // size of the pipe in X direction
  G4double gmpiptby = lxs->BPipeThickness; // top/bottom walls thickness
  G4double gmpipsidex = lxs->GammaBPipeWindowThickness; // side walls thickness
  G4double gfieldy = lxs->QBeamPipeContainerY;// size of the pipe in y direction

  G4Box *solidBeamPipeGamma2ndCMagF = new G4Box("solidBeamPipeGamma2ndCMagF", gmpipex/2.0, gfieldy/2.0, lpipe_gd3d/2.0);
  G4LogicalVolume *logicBeamPipeGamma2ndCMagFContainer = new G4LogicalVolume(solidBeamPipeGamma2ndCMagF,
                                                                 vacuumMaterial, "logicBeamPipeGamma2ndCMagFContainer");


  G4Box *solidGamMagPipeTopBot = new G4Box("solidGamMagPipeTopBot", gmpipex/2.0, gmpiptby/2.0, lpipe_gd3d/2.0);
  G4LogicalVolume *logicGamMagPipeTopBot = new G4LogicalVolume(solidGamMagPipeTopBot,
                                                   beamPipeMaterial, "logicGamMagPipeTopBot");

  G4Box *solidGamMagPipeSide = new G4Box("solidGamMagPipeSide", gmpipsidex/2.0, gfieldy/2.0 - gmpiptby, lpipe_gd3d/2.0);
  G4LogicalVolume *logicGamMagPipeSide = new G4LogicalVolume(solidGamMagPipeSide,
                                                   bpipeWindowMaterial, "logicGamMagPipeSide");


  new G4PVPlacement(0, G4ThreeVector(0.0, (gfieldy-gmpiptby)/2.0, 0.0), logicGamMagPipeTopBot,
                    "GMagnetFieldPipeTop", logicBeamPipeGamma2ndCMagFContainer , false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, -(gfieldy-gmpiptby)/2.0, 0.0), logicGamMagPipeTopBot,
                    "GMagnetFieldPipeBottom", logicBeamPipeGamma2ndCMagFContainer, false, 1, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector((gmpipex-gmpipsidex)/2.0, 0.0, 0.0), logicGamMagPipeSide,
                    "GMagnetFieldPipeSide", logicBeamPipeGamma2ndCMagFContainer , false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(-(gmpipex-gmpipsidex)/2.0, 0.0, 0.0), logicGamMagPipeSide,
                    "GMagnetFieldPipeSide", logicBeamPipeGamma2ndCMagFContainer, false, 1, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos3d), logicBeamPipeGamma2ndCMagFContainer,
                    "BeamPipeGamma2ndCMagF", fLogicWorld, false, 0, lxs->OverlapTest);*/
///////////
// Vacuum inside of the collimators
  G4Tubs *solidCollimatorVac = new G4Tubs("solidCollimatorVac", 0, lxs->CollimatorRin, lxs->CollimatorZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCollimatorVac = new G4LogicalVolume(solidCollimatorVac, vacuumMaterial, "logicCollimatorVac");
//  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->CollimatorZpos),
//                     logicCollimatorVac, "CollimatorVac", fLogicWorld, false, 0, lxs->OverlapTest); // we removed 1st collimator
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->CollimatorZpos + lxs->CollimatorZ + 50.0*cm),
                     logicCollimatorVac, "CollimatorVac", fLogicWorld, false, 1, lxs->OverlapTest);
///////////

/*  G4Box *solidBeamPipeGammaMagFLyso = new G4Box("solidBeamPipeGammaMagFLyso", lxs->QBeamPipeContainerX/2.0, lxs->QBeamPipeContainerY/2.0, lpipe_gd4th/2.0);
  G4LogicalVolume *logicBeamPipeGammaMagFLyso = new G4LogicalVolume(solidBeamPipeGammaMagFLyso,
                                                                  beamPipeMaterial, "logicBeamPipeGammaMagFLyso");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos4th),
                               logicBeamPipeGammaMagFLyso, "BeamPipeGammaMagFLyso", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Box *solidBeamPipeGammaMagFLysoVac = new G4Box("solidBeamPipeGammaMagFLysoVac", lxs->QBeamPipeContainerX/2.0,
                                                    lxs->QBeamPipeContainerY/2.0, lpipe_gd4th/2.0);
  G4LogicalVolume *logicBeamPipeGammaMagFLysoVac = new G4LogicalVolume(solidBeamPipeGammaMagFLysoVac,
                                                                  vacuumMaterial, "logicBeamPipeGammaMagFLysoVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos4th),
                               logicBeamPipeGammaMagFLysoVac, "BeamPipeGammaMagFLysoVac", fLogicWorld, false, 0, lxs->OverlapTest); */

//quadratish  Beam Pipe  2nd collimator - LYSO detectors
/*  G4double lpipe_gd3d = lxs->ComptonLysoZpos - lxs->ComptonLysoZ/2.0 - lpipe_pos_2nd - lxs->CollimatorZ - lpipe_gd2nd/2;
  G4double lpipe_pos3d = lxs->ComptonLysoZpos - lxs->ComptonLysoZ/2.0 - lpipe_gd3d/2;

  G4Box *QBeamPipeContainer = new G4Box("QBeamPipeContainer", lxs->QBeamPipeContainerX/2.0, lxs->QBeamPipeContainerY/2.0, lpipe_gd3d/2.0);
  G4LogicalVolume *logicQBeamPipeContainer = new G4LogicalVolume(QBeamPipeContainer,
                                                                  vacuumMaterial, "logicQBeamPipeContainer");
   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lpipe_pos3d),
                               logicQBeamPipeContainer, "QBeamPipeComptonContainer", fLogicWorld, false, 0, lxs->OverlapTest);*/


//Tracker
/*  G4Material* trackerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonTrackerMaterial);
  G4Material* ComptonDetContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptonDetMaterial);

  G4Box *solidComptonDetContainer = new G4Box("solidComptonDetContainer", lxs->ComptonDetX/2.0, lxs->ComptonDetY/2.0, lxs->ComptonDetZ/2.0);
  G4LogicalVolume *logicComptonDetContainer = new G4LogicalVolume(solidComptonDetContainer,
                                                                  ComptonDetContainerMaterial, "logicComptonDetContainer");

  G4Box *solidComptonTracker = new G4Box("solidComptonTracker", lxs->ComptonTrackerX/2.0, lxs->ComptonTrackerY/2.0, lxs->ComptonTrackerZ/2.0);
  G4LogicalVolume *logicComptonTracker = new G4LogicalVolume(solidComptonTracker, trackerMaterial, "logicComptonTracker");
  for (G4int il = 0; il < lxs->ComptonTrackerNLayers; ++il) {
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),  // Let it face backward toward beam
                      G4ThreeVector(0.0, 0.0, (-lxs->ComptonDetZ+lxs->ComptonTrackerZ)/2.0 + il*(lxs->ComptonTrackerZ+lxs->ComptonTrackerDZ)),
                      logicComptonTracker, "ComptonTracker", logicComptonDetContainer, false, il, lxs->OverlapTest);
  }
  G4double bpipel = lxs->ComptonDetZpos - lxs->GMagnetZpos - lxs->GMagnetZ/2.0;
  G4double thetad = atan2(lxs->ComptonDetX/2.0 + 4.0*lxs->BPipeR, bpipel);
 // new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, -1.0, 0.0), thetad),
    new G4PVPlacement(0,
                    G4ThreeVector((bpipel+lxs->ComptonDetZ/2.0)*sin(thetad), 0.0,
                                  (bpipel+lxs->ComptonDetZ/2.0)*cos(thetad) + lxs->GMagnetZpos + lxs->GMagnetZ/2.0),
                    logicComptonDetContainer, "ComptonDetContainer", fLogicWorld, false, 0, lxs->OverlapTest);
    new G4PVPlacement(0,
                    G4ThreeVector(-(bpipel+lxs->ComptonDetZ/2.0)*sin(thetad), 0.0,
                                  (bpipel+lxs->ComptonDetZ/2.0)*cos(thetad) + lxs->GMagnetZpos + lxs->GMagnetZ/2.0),
                    logicComptonDetContainer, "ComptonDetContainer", fLogicWorld, false, 1, lxs->OverlapTest);   */

  AddSegmentation();

  G4double supporty;
  G4AssemblyVolume* lysoSupportAssembly = ConstructSupportAssembly(supporty);
  G4double chery = lxs->ComptonElectronBeamtoStageY;
  G4double detxpos = (lxs->LYSODetTopPlateX - lxs->ComptonLysoX)/2.0 + lxs->ComptonLysoXpos;
//   G4double supypos = -(lxs->ComptonLysoY + lxs->OPPPDetTopPlateY)/2.0;
  G4double supypos = -chery + lxs->OPPPDetTopPlateY/2.0 + lxs->LYSODetSupportBallH;
  G4double detzpos = lxs->ComptonLysoZpos - (lxs->ComptonLysoZ - lxs->LYSODetTopPlateZ)/2.0;
  G4ThreeVector trsupport(detxpos, supypos, detzpos);
  lysoSupportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);
  trsupport.setX(-detxpos);
  lysoSupportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 1, lxs->OverlapTest);

//Small pedestal for LYSO detector
  G4Material* lysoPedestaMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->LYSOPedestaMaterial);
  G4double pedy = chery - lxs->ComptonLysoY/2.0 - lxs->LYSODetSupportBallH - lxs->OPPPDetTopPlateY;
  G4Box *solidLYSOPedestal = new G4Box("solidLYSOPedestal", lxs->ComptonLysoX/2.0,
                                            pedy/2.0, lxs->ComptonLysoZ/2.0);
  G4LogicalVolume *logicLYSOPedestal = new G4LogicalVolume(solidLYSOPedestal, lysoPedestaMaterial,
                                                                "logicLYSOPedestal");
  new G4PVPlacement (0, G4ThreeVector(lxs->ComptonLysoXpos, -(lxs->ComptonLysoY+pedy)/2.0, lxs->ComptonLysoZpos),
                     logicLYSOPedestal, "LysoPedestal", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-lxs->ComptonLysoXpos, -(lxs->ComptonLysoY+pedy)/2.0, lxs->ComptonLysoZpos),
                     logicLYSOPedestal, "LysoPedestal", fLogicWorld, false, 1, lxs->OverlapTest);

// Table to hold LYSO and Cherenkov detectors.
  G4double ypestal = 1.5*m;
  G4double ylevel = supporty + pedy + lxs->ComptonLysoY/2.0;
  G4double tblhight = -ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4double tblx = lxs->LYSOBasePlateX;
  G4double tblz = lxs->LYSOBasePlateZ;
  G4AssemblyVolume *tablesupport = LxAux::BuildTable("GammaSpectrTable", tblx, tblhight, tblz, 4);
  G4double tablezpos = lxs->ComptonLysoZpos - (lxs->ComptonLysoZ - lxs->LYSOBasePlateZ)/2.0;
  G4ThreeVector trtbl(0.0, -ylevel, tablezpos);
  tablesupport->MakeImprint(fLogicWorld, trtbl, 0, 0, lxs->OverlapTest);

  G4LogicalVolume *pedstal = LxAux::BuildPedestal("GammaSpectr", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.5*ypestal + lxs->FloorSurfaceYpos, tablezpos), pedstal,
                      "GammaSpectrPedestal", fLogicWorld, false, 0, lxs->OverlapTest);

}



void LxDetectorComptonFluka::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
    fDetector->AddSensorSegmentation("LysoCal", lxs->ComptonLysoX, lxs->ComptonLysoY,
                                                     lxs->ComptonLysoNCellX, lxs->ComptonLysoNCellY);
//   fDetector->AddSensorSegmentation("ComptonTracker", lxs->ComptonTrackerX, lxs->ComptonTrackerY,
//                                                      lxs->ComptonTrackerNCellX, lxs->ComptonTrackerNCellY);
}




G4AssemblyVolume* LxDetectorComptonFluka::ConstructSupportAssembly(G4double &sphight)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* lysoDetSupportMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPDetSupportMaterial);

  G4Box *solidLYSOBottomSupport = new G4Box("solidLYSOBottomSupport", lxs->LYSODetBottomSupportX/2.0,
                                            lxs->OPPPDetBottomSupportY/2.0, lxs->LYSODetBottomSupportZ/2.0);
  G4LogicalVolume *logicLYSOBottomSupport = new G4LogicalVolume(solidLYSOBottomSupport, lysoDetSupportMaterial,
                                                                "logicLYSOBottomSupport");

  G4Box *solidLYSOTopPlate = new G4Box("solidLYSOTopPlate", lxs->LYSODetTopPlateX/2.0,
                                            lxs->OPPPDetTopPlateY/2.0, lxs->LYSODetTopPlateZ/2.0);
  G4LogicalVolume *logicLYSOTopPlate = new G4LogicalVolume(solidLYSOTopPlate, lysoDetSupportMaterial,
                                                                "logicLYSOTopPlate");

  G4Sphere *solidLYSODetSupportBall = new G4Sphere("solidLYSODetSupportBall", 0.0, lxs->LYSODetSupportBallH/2.0,
                                                   0.0, 2.0*M_PI, 0.0, M_PI);
  G4LogicalVolume *logicLYSODetSupportBall = new G4LogicalVolume(solidLYSODetSupportBall, lysoDetSupportMaterial,
                                                                "logicLYSODetSupportBall");

  G4AssemblyVolume *lysoSupportAssembly = new G4AssemblyVolume();

  G4ThreeVector asstr(0.0, 0.0, 0.0);
  lysoSupportAssembly->AddPlacedVolume(logicLYSOTopPlate, asstr, 0);
  asstr.setY(-(lxs->OPPPDetBottomSupportY + lxs->OPPPDetTopPlateY)/2.0 - lxs->LYSODetSupportBallH);
  asstr.setZ((lxs->LYSODetBottomSupportZ - lxs->LYSODetTopPlateZ)/2.0);
  lysoSupportAssembly->AddPlacedVolume(logicLYSOBottomSupport, asstr, 0);

  G4double dbz = 0.5*lxs->LYSODetTopPlateZ - 1.05*lxs->LYSODetSupportBallH;
  G4double dbx = 0.5*lxs->LYSODetTopPlateX - 1.05*lxs->LYSODetSupportBallH;
  G4ThreeVector bpos(dbx, asstr.y() + (lxs->OPPPDetBottomSupportY+lxs->LYSODetSupportBallH)/2.0 , dbz);
  lysoSupportAssembly->AddPlacedVolume(logicLYSODetSupportBall, bpos, 0);
  bpos.setX(-dbx);
  lysoSupportAssembly->AddPlacedVolume(logicLYSODetSupportBall, bpos, 0);
  bpos.setX(0.0);
  bpos.setZ(-dbz);
  lysoSupportAssembly->AddPlacedVolume(logicLYSODetSupportBall, bpos, 0);

  //Hexapod
  G4double hexhight = 0.0;
  G4AssemblyVolume* hexAssembly = LxAux::BuildHexapod("LYSO", hexhight);
  G4ThreeVector hexpos(0.0, asstr.y() - (lxs->OPPPDetBottomSupportY+hexhight)/2.0, (lxs->LYSODetBottomSupportZ-lxs->LYSODetTopPlateZ)/2.0);
//  lysoSupportAssembly->AddPlacedAssembly(hexAssembly, hexpos, 0);
  LxAux::AddAssmblyVolumes(lysoSupportAssembly, hexAssembly, hexpos, 0);


  sphight = (lxs->OPPPDetTopPlateY + hexhight)/2.0 - hexpos.y();
  return lysoSupportAssembly;
}



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
