//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4EllipticalTube.hh"
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
#include "LxGammaCalo.hh"


////////////////////////////////////////////////////////////////////////
/////  LxDetectorGammaCalo

void LxDetectorGammaCalo::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* gammaDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaDumpMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* gammaCaloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaCaloMaterial);
  
  G4LogicalVolume *logicGammaDump = 0;
  G4VPhysicalVolume *edump = G4PhysicalVolumeStore::GetInstance()->GetVolume("BeamDumpAssembly", true);
  if (edump) {
    logicGammaDump = edump->GetLogicalVolume();
  } else {
    G4Tubs *solidGammaDump = new G4Tubs("solidGammaDump", 0.0, lxs->GammaDumpR, lxs->GammaDumpZ/2.0, 0.0, 2.0*M_PI);
    logicGammaDump = new G4LogicalVolume(solidGammaDump, gammaDumpMaterial, "logicGammaDump");
  }      
  if (logicGammaDump) new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->GammaDumpZpos), 
                      logicGammaDump, "GammaDump", fLogicWorld, false, 0, lxs->OverlapTest);
                     
 // cylindrically shaped 
//  G4Tubs *solidGammaCalo = new G4Tubs("solidGammaCalo", lxs->BPipeR, lxs->GammaDumpR, lxs->GammaCaloZ/2.0, 0.0, 2.0*M_PI);
//  G4LogicalVolume *logicGammaCalo = new G4LogicalVolume(solidGammaCalo, gammaCaloMaterial, "logicGammaCalo");
//  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->GammaCaloZpos), 
//                     logicGammaCalo, "GammaCalo", fLogicWorld, false, 0, lxs->OverlapTest);
                     
// Lead-glass blocks 9*9*50 cm3 //"solidLeadGlass" LeadGlass
  G4Box *solidLeadGlass = new G4Box("solidGammaCalo", lxs->LeadGlassX/2.0, lxs->LeadGlassY/2.0, lxs->LeadGlassZ/2.0);
  G4LogicalVolume *logicLeadGlass = new G4LogicalVolume(solidLeadGlass, gammaCaloMaterial, "logicLeadGlass");
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->LeadGlassY/2.0+lxs->BPipeR+lxs->BPipeThickness, lxs->GammaCaloZpos), 
                    logicLeadGlass, "GammaCalo", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness, -lxs->LeadGlassY/2.0+lxs->BPipeR+lxs->BPipeThickness, lxs->GammaCaloZpos), 
                    logicLeadGlass, "GammaCalo", fLogicWorld, false, 1, lxs->OverlapTest);     
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(-lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness, -lxs->LeadGlassY/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->GammaCaloZpos), 
                    logicLeadGlass, "GammaCalo", fLogicWorld, false, 2, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
                    G4ThreeVector(-lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->LeadGlassY/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->GammaCaloZpos), 
                    logicLeadGlass, "GammaCalo", fLogicWorld, false, 3, lxs->OverlapTest);                           
                    
                    
                    

  G4double lpipe_gd = lxs->GammaDumpZpos - lxs->ComptonDetZpos - lxs->GammaDumpZ/2.0;
  G4Tubs *solidBeamPipeGammaD = new G4Tubs("solidBeamPipeGammaD", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR, lpipe_gd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGammaD = new G4LogicalVolume(solidBeamPipeGammaD, beamPipeMaterial, "logicBeamPipeGammaD");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0), 
                               logicBeamPipeGammaD, "BeamPipeGammaD", fLogicWorld, false, 0, lxs->OverlapTest);
  
  G4Tubs *solidBeamPipeGammaDVac = new G4Tubs("solidBeamPipeGammaDVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, lpipe_gd/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeGammaDVac = new G4LogicalVolume(solidBeamPipeGammaDVac, vacuumMaterial, "logicBeamPipeGammaDVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0), 
                               logicBeamPipeGammaDVac, "BeamPipeGammaDVac", fLogicWorld, false, 0, lxs->OverlapTest);

  AddSegmentation();

}


void LxDetectorGammaCalo::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("GammaCalo", lxs->LeadGlassX, lxs->LeadGlassY, 1, 1);
}


void LxDetectorGammaCalo48::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();
  CreateMaterial();
  
  G4double bpipeRmin = 19.00 *mm;
    
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
 // G4Material* gammaDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaDumpMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* gammaCaloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaCaloMaterial);
  G4Material* LGFoilMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaCaloFoilMaterial);
  G4Material* GammaMonitorContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaMonitorDetMaterial); 
  G4Material* CalAbsorberMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->CAbsorberMaterial);
  
//Dump at the ent of the beam line
 // G4Tubs *solidGammaDumpf = new G4Tubs("solidGammaDump", 0.0, lxs->GammaDumpR, lxs->GammaDumpZ/2.0, 0.0, 2.0*M_PI);
 // G4LogicalVolume *logicGammaDump = new G4LogicalVolume(solidGammaDump, gammaDumpMaterial, "logicGammaDump");
 
  G4VPhysicalVolume *physGammaDump = G4PhysicalVolumeStore::GetInstance()->GetVolume("BeamDumpAssembly", true);
  if (physGammaDump) {
    G4LogicalVolume *logicGammaDump = physGammaDump->GetLogicalVolume();
    new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->GammaDumpZpos), 
                     logicGammaDump, "GammaDump", fLogicWorld, false, 0, lxs->OverlapTest);
  }
  
  //**************** cylindrically shaped Absorber placed between  Gamma Monitor & Dump
  G4Tubs *solidCalAbsorber = new G4Tubs("solidCalAbsorber", bpipeRmin, lxs->GammaDumpR, lxs->CAbsorberThickness/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCalAbsorber = new G4LogicalVolume(solidCalAbsorber, CalAbsorberMaterial, "logicCalAbsorber");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->CAbsorberZpos), 
                     logicCalAbsorber, "CalAbsorber", fLogicWorld, false, 0, lxs->OverlapTest);
                     
 //**************** cylindrically shaped  Gamma Monitor
//  G4Tubs *solidGammaCalo = new G4Tubs("solidGammaCalo", lxs->BPipeR, lxs->GammaDumpR, lxs->GammaCaloZ/2.0, 0.0, 2.0*M_PI);
//  G4LogicalVolume *logicGammaCalo = new G4LogicalVolume(solidGammaCalo, gammaCaloMaterial, "logicGammaCalo");
//  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->GammaCaloZpos), 
//                     logicGammaCalo, "GammaCalo", fLogicWorld, false, 0, lxs->OverlapTest);

//***************4 9*9*50 cm^3 Hermes blocks for Gamma monitor
//   G4Box *solidLeadGlass = new G4Box("solidGammaCalo", lxs->LeadGlassX/2.0, lxs->LeadGlassY/2.0, lxs->LeadGlassZ/2.0);
//   G4LogicalVolume *logicLeadGlass = new G4LogicalVolume(solidLeadGlass, gammaCaloMaterial, "logicLeadGlass");
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->LeadGlassY/2.0+lxs->BPipeR+lxs->BPipeThickness, lxs->GammaCaloZpos), 
//                     logicLeadGlass, "GammaCalo", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness, -lxs->LeadGlassY/2.0+lxs->BPipeR+lxs->BPipeThickness, lxs->GammaCaloZpos), 
//                     logicLeadGlass, "GammaCalo", fLogicWorld, false, 1, lxs->OverlapTest);     
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(-lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness, -lxs->LeadGlassY/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->GammaCaloZpos), 
//                     logicLeadGlass, "GammaCalo", fLogicWorld, false, 2, lxs->OverlapTest);
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(-lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->LeadGlassY/2.0-lxs->BPipeR-lxs->BPipeThickness, lxs->GammaCaloZpos), 
//                     logicLeadGlass, "GammaCalo", fLogicWorld, false, 3, lxs->OverlapTest);  


  

                     
// Lead-glass blocks 3.8*3.8*45 cm3 //"solidLeadGlass" LeadGlass
  		
  G4int LGidx=0;      
//  G4double shift = 12.5*mm;//17 mm
  G4double prad = bpipeRmin;
//  G4double pthkns = lxs->BPipeThickness;
  G4double lgd2 = lxs->LeadGlassX/2.0;
  
  G4double lgthkns = lxs->GammaCaloFoilThickness;
//  G4double lgdz = lxs->LeadGlassZ;
   G4double lgdz = lxs->LeadGlassZ + lgthkns;

  G4Box *solidLeadGlass = new G4Box("solidGammaCalo", lxs->LeadGlassX/2.0, lxs->LeadGlassY/2.0, lxs->LeadGlassZ/2.0);
  G4LogicalVolume *logicLeadGlass = new G4LogicalVolume(solidLeadGlass, gammaCaloMaterial, "logicLeadGlass");
//wrapping  
  G4Box *solidLGFoil = new G4Box("solidLGFoil", lgd2+lgthkns, lgd2+lgthkns, lgdz/2.0);
  G4LogicalVolume *logicLGFoil = new G4LogicalVolume(solidLGFoil, LGFoilMaterial, "logicLGFoil");
//    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), 
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -0.5*(lgdz-lxs->LeadGlassZ)), 
                      logicLeadGlass, "LeadGlass", logicLGFoil, false, 0, lxs->OverlapTest);

//4*4
//   G4double xy[4] = {(lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness+shift),(-lxs->LeadGlassX/2.0+lxs->BPipeR+lxs->BPipeThickness+shift),(lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness-shift),(-lxs->LeadGlassX/2.0-lxs->BPipeR-lxs->BPipeThickness-shift)};                   
//  for (G4int il = 0; il < 4; ++il) {
//    for (G4int jl = 0; jl < 4; ++jl) {
//      std::cout << "LG index: " << il << jl<< std::endl;
// if(!( (il==1&&jl==2)||(il==2&&jl==2)|| (il==2&&jl==1)||(il==1&&jl==1) )){ 
// LGidx=LGidx+1;
//  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                    G4ThreeVector(xy[il], xy[jl], lxs->GammaCaloZpos), 
//                      logicLGFoil, "GammaCalo", fLogicWorld, false, LGidx, lxs->OverlapTest);
//                    logicLeadGlass, "GammaCalo", fLogicWorld, false, LGidx, lxs->OverlapTest);
//                    }                               
//  			}           
//  }                 

//6*6  
//   G4double xy[6] = {(3*lgd2+prad+pthkns+shift+2*lgthkns),(lgd2+lgthkns+prad+pthkns+shift),(-lgd2-lgthkns+prad+pthkns+shift),(lgd2+lgthkns-prad-pthkns-shift),(-lgd2-lgthkns-prad-pthkns-shift),(-3*lgd2-2*lgthkns-prad-pthkns-shift)};                
//  for (G4int il = 0; il < 6; ++il) {
//    for (G4int jl = 0; jl < 6; ++jl) {
//      std::cout << "LG index: " << il << jl<< std::endl;
//  if(!( (il==2&&jl==3)||(il==2&&jl==2)|| (il==3&&jl==2)||(il==3&&jl==3) )){ 
// LGidx=LGidx+1;
//  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                    G4ThreeVector(xy[il], xy[jl], lxs->GammaCaloZpos), 
//                      logicLGFoil, "GammaCalo", fLogicWorld, false, LGidx, lxs->OverlapTest);
//                    }                               
//  			}           
//  } 
  
  //5*5  with beam pipe the same size as LG blocks
 //   G4double xy[5] = {(3*lgd2+prad+pthkns+2*lgthkns),(lgd2+lgthkns+prad+pthkns),(0.0),(-lgd2-lgthkns-prad-pthkns),(-3*lgd2-2*lgthkns-prad-pthkns)};                
//   for (G4int il = 0; il < 5; ++il) {
//     for (G4int jl = 0; jl < 5; ++jl) {
//       std::cout << "LG index: " << il << jl<< std::endl;
//   if(!( (il==2&&jl==2) )){ 
//  LGidx=LGidx+1;
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(xy[il], xy[jl], lxs->GammaCaloZpos), 
//                       logicLGFoil, "GammaCalo", fLogicWorld, false, LGidx, lxs->OverlapTest);
//                     }                               
//   			}           
//   } 
  
    //7*7  with beam pipe the same size as LG blocks
  G4int nb=7;
  G4double xsizec = nb * (lxs->LeadGlassX + 2.0*lgthkns);
  G4double ysizec = nb * (lxs->LeadGlassY + 2.0*lgthkns);
  G4Box *solidGammaMonitorC1 = new G4Box("solidGammaMonitorC1", 0.5*xsizec, 0.5*ysizec, 0.5*lgdz);
  G4Box *solidLGCut = new G4Box("solidLGCut", lgd2+lgthkns, lgd2+lgthkns, 0.6*lgdz);
  G4SubtractionSolid* solidGammaMonitorContainer = new G4SubtractionSolid("solidGammaMonitorContainer", solidGammaMonitorC1, solidLGCut);                       
  G4LogicalVolume *logicGammaMonitorContainer = new G4LogicalVolume(solidGammaMonitorContainer, 
                                                                   GammaMonitorContainerMaterial, "logicGammaMonitorContainer");

//    G4double xy[nb] = {(5*lgd2+prad+3*lgthkns),(3*lgd2+prad+2*lgthkns),(lgd2+lgthkns+prad),(0.0),(-lgd2-lgthkns-prad),(-3*lgd2-2*lgthkns-prad),(-5*lgd2-prad-3*lgthkns)};                
  for (G4int il = 0; il < nb; ++il) {
    for (G4int jl = 0; jl < nb; ++jl) {
//       std::cout << "LG index: " << il << jl<< std::endl;
      if ( il==nb/2 && jl==nb/2 ) continue; 
      G4double xpos = -0.5*xsizec + (il+0.5) * (lxs->LeadGlassX + 2.0*lgthkns);
      G4double ypos = -0.5*ysizec + (jl+0.5) * (lxs->LeadGlassY + 2.0*lgthkns);
      new G4PVPlacement(0, G4ThreeVector(xpos, ypos, 0), 
                      logicLGFoil, "GammaCalo", logicGammaMonitorContainer, false, il*nb+jl, lxs->OverlapTest);                      

//  if(!( (il==3&&jl==3) )){ 
//  LGidx=LGidx+1;
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), // Let it face backward toward beam
//                     G4ThreeVector(xy[il], xy[jl], lxs->GammaCaloZpos), 
// //                      logicLGFoil, "GammaCalo", fLogicWorld, false, LGidx, lxs->OverlapTest);
//                       logicLGFoil, "GammaCalo", logicGammaMonitorContainer, false, LGidx, lxs->OverlapTest);                      
//                     }                               
  	}           
  } 
  
   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), 
                    G4ThreeVector(0.0, 0.0, lxs->GammaCaloZpos), 
                    logicGammaMonitorContainer, "GammaMonitorContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  
//  to connect  2 beam pipes with different diameter         
/* G4Tubs *solidConnect2BeamPipeGammaD = new G4Tubs("solidConnect2BeamPipeGammaD", bpipeRmin, lxs->BPipeR, lxs->BPipeThickness, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicConnect2BeamPipeGammaD = new G4LogicalVolume(solidConnect2BeamPipeGammaD, beamPipeMaterial, "logicConnect2BeamPipeGammaD");
 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lxs->BPipeThickness), 
                               logicConnect2BeamPipeGammaD, "Connect2BeamPipeGammaD", fLogicWorld, false, 0, lxs->OverlapTest);
  
 G4double lpipe_gd = lxs->GammaDumpZpos - lxs->ComptonDetZpos - lxs->GammaDumpZ/2.0;
 G4Tubs *solidBeamPipeGammaD = new G4Tubs("solidBeamPipeGammaD", bpipeRmin-lxs->BPipeThickness, bpipeRmin, lpipe_gd/2.0, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicBeamPipeGammaD = new G4LogicalVolume(solidBeamPipeGammaD, beamPipeMaterial, "logicBeamPipeGammaD");
 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0), 
                              logicBeamPipeGammaD, "BeamPipeGammaD", fLogicWorld, false, 0, lxs->OverlapTest);
  
 G4Tubs *solidBeamPipeGammaDVac = new G4Tubs("solidBeamPipeGammaDVac", 0.0, bpipeRmin-lxs->BPipeThickness, lpipe_gd/2.0, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicBeamPipeGammaDVac = new G4LogicalVolume(solidBeamPipeGammaDVac, vacuumMaterial, "logicBeamPipeGammaDVac");

 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0), 
                              logicBeamPipeGammaDVac, "BeamPipeGammaDVac", fLogicWorld, false, 0, lxs->OverlapTest);*/

AddSegmentation();
}

void LxDetectorGammaCalo48::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("GammaCalo", lxs->LeadGlassX, lxs->LeadGlassY, 1, 1);
}



void LxDetectorGammaCalo48::CreateMaterial()
{
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;
  G4int    ncomponents;

  G4double fractionmass;
// define Elements
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);

  G4Element* Si =   new G4Element("Silicon", "Si", z=14, a=28.09 *g/mole);
  G4Element* Pb = 	new G4Element("Lead", "Pb", 82., 207.2 * g/mole);
  G4Element* K  = 	new G4Element("Potassium", "K", 19., 39.098 * g/mole);
  G4Element* As = 	new G4Element("Arsenic", "As", 33., 74.92 * g/mole);
  G4Element* Ce = 	new G4Element("Cerium", "Ce", 58., 140.12 * g/mole);

  // Materials necessary to build lead-glass
  double bigden = 1e9*g/cm3;
  G4Material* PbO = new G4Material("TF1_PbO", bigden, 2);
  PbO->AddElement(Pb, 1);
  PbO->AddElement(O, 1);

  G4Material* K2O = new G4Material("TF1_K2O", bigden, 2);
  K2O->AddElement(K, 2);
  K2O->AddElement(O, 1);

  G4Material* As2O3 = new G4Material("TF1_As2O3", bigden, 2);
  As2O3->AddElement(As, 2);
  As2O3->AddElement(O, 3);

  G4Material* SiO2 = new G4Material("F101_SiO", bigden, 2);
  SiO2->AddElement(Si, 1);
  SiO2->AddElement(O, 2);

  G4Material* Na2O = new G4Material("F101_Na2O", bigden, 2);
  Na2O->AddElement(Na, 2);
  Na2O->AddElement(O, 1);

  G4Material* TF1 = new G4Material("TF1", 3.86*g/cm3, 5);
  TF1->AddMaterial(PbO, 0.512);
  TF1->AddMaterial(SiO2, 0.413);
  TF1->AddMaterial(K2O, 0.035);
  TF1->AddMaterial(Na2O, 0.035);
  TF1->AddMaterial(As2O3, 0.005);

  //Recalculated LG composition in mass units K swapped with Na
  G4Material* LG_TF1 = new G4Material("LG_TF1", density= 3.86*g/cm3, ncomponents=6);
  LG_TF1->AddElement(Pb, fractionmass=0.475);
  LG_TF1->AddElement(O, fractionmass=0.273);
  LG_TF1->AddElement(Si, fractionmass=0.193);
  LG_TF1->AddElement(K, fractionmass=0.029);
  LG_TF1->AddElement(Na, fractionmass=0.026);
  LG_TF1->AddElement(As, fractionmass=0.004);

  G4Material* Pb3O4 = new G4Material("F101_PbO", bigden, 2);
  Pb3O4->AddElement(Pb, 3);
  Pb3O4->AddElement(O, 4);

  G4Material* Ce1 = new G4Material("F101_Ce", bigden, 1);
  Ce1->AddElement(Ce, 1);

  G4Material* F101 = new G4Material("F101", 3.86*g/cm3, 4);
  F101->AddMaterial(Pb3O4, 0.5123);
  F101->AddMaterial(SiO2, 0.4153);
  F101->AddMaterial(K2O, 0.070);
  F101->AddMaterial(Ce1, 0.002);
}



void LxDetectorGammaCalo8::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();
  CreateMaterial();

  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* gammaCaloMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaCaloMaterial);
  G4Material* LGFoilMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaCaloFoilMaterial);
  G4Material* GammaMonitorContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaMonitorDetMaterial);
  G4Material* CalSupportMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->CalSupportMaterial);
  G4Material* beamDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaBeamDumpMaterial);
  G4Material* dumpInsertMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaBeamDumpInsertMaterial);

// Construct gamma beam dump
  G4Tubs *solidBeamDump = new G4Tubs("solidBeamDump", 0.0, lxs->GammaBeamDumpR, lxs->GammaBeamDumpZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicGammaBeamDump = new G4LogicalVolume(solidBeamDump, beamDumpMaterial, "logicGammaBeamDump");

  G4Tubs *solidBeamDumpInsert = new G4Tubs("solidBeamDumpInsert", 0.0, lxs->GammaBeamDumpInsertR,
                                           lxs->GammaBeamDumpInsertZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamDumpInsert = new G4LogicalVolume(solidBeamDumpInsert, dumpInsertMaterial, "logicBeamDumpInsert");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->GammaBeamDumpInsertZ-lxs->GammaBeamDumpZ)/2.0),
                     logicBeamDumpInsert, "GammaBeamDumpInsert", logicGammaBeamDump, false, 0, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->GammaBeamDumpZpos),
                               logicGammaBeamDump, "GammaBeamDumpAssembly", fLogicWorld, false, 0, lxs->OverlapTest);

// Construct gamma beam dump support
  G4Material* GBMDSupportMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaBeamDumpSupportMaterial);
  G4double gbdsupx = 3.0*lxs->GammaBeamDumpR;
  G4double gbdsupy = -lxs->FloorSurfaceYpos - lxs->GammaBeamDumpR/2.0;
  G4double gbdsupz = 0.5*lxs->GammaBeamDumpZ;
  G4Box *solidGammaBeamDumpSupport0 = new G4Box("solidGammaBeamDumpSupport0", gbdsupx/2.0, gbdsupy/2.0, gbdsupz/2.0);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(0.0, (gbdsupy + lxs->GammaBeamDumpR)/2.0, 0.0));
  G4SubtractionSolid *solidGammaBeamDumpSupport = new G4SubtractionSolid("solidGammaBeamDumpSupport",
                                 solidGammaBeamDumpSupport0, solidBeamDump, transform);
  G4LogicalVolume *logicGammaBeamDumpSupport = new G4LogicalVolume(solidGammaBeamDumpSupport,
                                        GBMDSupportMaterial, "logicGammaBeamDumpSupport");
  new G4PVPlacement(0, G4ThreeVector(0.0, lxs->FloorSurfaceYpos+gbdsupy/2.0, lxs->GammaBeamDumpZpos),
                    logicGammaBeamDumpSupport, "GammaBeamDumpSupport", fLogicWorld, false, 0, lxs->OverlapTest);

  //12 Lead-glass blocks around the clock 3.8*3.8*45 cm3 //"solidLeadGlass" LeadGlass
  G4int NLGblocks =8;
  G4double lgthkns = lxs->GammaCaloFoilThickness;
  G4double lgdz = lxs->LeadGlassZ + lgthkns;
  G4double lgdy = lxs->LeadGlassY + 2.0*lgthkns;
  G4double lgdx = lxs->LeadGlassX + 2.0*lgthkns;

//  std::vector<double>LGphi ={0, M_PI/4.0, M_PI/2.0, 3*M_PI/4.0, M_PI, 5*M_PI/4.0, 3*M_PI/2.0, 7*M_PI/4.0};

  G4Box *solidLeadGlass = new G4Box("solidGammaCalo", lxs->LeadGlassX/2.0, lxs->LeadGlassY/2.0, lxs->LeadGlassZ/2.0);
  G4LogicalVolume *logicLeadGlass = new G4LogicalVolume(solidLeadGlass, gammaCaloMaterial, "logicLeadGlass");
//wrapping
  G4Box *solidLGFoil = new G4Box("solidLGFoil", lgdx/2.0, lgdy/2.0, lgdz/2.0);
  G4LogicalVolume *logicLGFoil = new G4LogicalVolume(solidLGFoil, LGFoilMaterial, "logicLGFoil");
//    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -0.5*(lgdz-lxs->LeadGlassZ)),
                      logicLeadGlass, "LeadGlass", logicLGFoil, false, 0, lxs->OverlapTest);


  //Gamma Monitor Container
  //round container
  G4double rcout = sqrt(pow(lxs->CalSupportR+lxs->CalSupporthickness+lgdy, 2.0) + pow(lgdx/2.0, 2.0));
//  G4Tubs *solidGammaMonitorContainer = new G4Tubs("solidGammaMonitorContainer", lxs->CalSupportR-lxs->CalSupporthickness, rcout,  lgdz/2.0,  0.0, 2.0*M_PI);

  //box container
  G4double gmcontx = lxs->CalSupportR + lgdy;
  G4double gmconty = gmcontx + 20.0*cm;
  G4Box *solidGammaMonitorContainer = new G4Box("solidGammaMonitorContainer", gmcontx, gmconty, lgdz/2.0);
  G4LogicalVolume *logicGammaMonitorContainer = new G4LogicalVolume(solidGammaMonitorContainer, 
                                                     GammaMonitorContainerMaterial, "logicGammaMonitorContainer");

   //**************** cylindrically shaped support of LG blocks
  G4Tubs *solidCalSupport = new G4Tubs("solidCalSupport", lxs->CalSupportR-lxs->CalSupporthickness,
                                       lxs->CalSupportR, lgdz/2.0,  0.0, 2.0*M_PI);
  G4LogicalVolume *logicCalSupport = new G4LogicalVolume(solidCalSupport, CalSupportMaterial, "logicCalSupport");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0),
                     logicCalSupport, "CalSupport", logicGammaMonitorContainer, false, 0, lxs->OverlapTest);

  G4double rlg=lxs->CalSupportR + lgdy/2.0;
  for (G4int i = 0; i < NLGblocks; ++i) {
    G4double LGphi = i * 2.0*M_PI/static_cast<G4double>(NLGblocks) - M_PI_2;
    new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), LGphi),
                      G4ThreeVector( rlg*cos(LGphi), rlg*sin(LGphi), 0.0),
                     logicLGFoil, "GammaCalo", logicGammaMonitorContainer, false, i, lxs->OverlapTest);
  }

// Construct gamma monitor support
  G4Material* GMSupportMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaMonitorSupportMaterial);
  G4Material* GMHolderMaterial= G4NistManager::Instance()->FindOrBuildMaterial(lxs->GammaMonitorHolderMaterial);
  G4double gmsupx = 3.0*rlg;
  G4double gmsupy = -lxs->FloorSurfaceYpos - gmconty;
  G4double gmsupz = 0.5*lxs->GammaCaloZ;
  G4Box *solidGammaMonitorSupport0 = new G4Box("solidGammaMonitorSupport0", gmsupx/2.0, gmsupy/2.0, gmsupz/2.0);
  G4LogicalVolume *logicGammaMonitorSupport = new G4LogicalVolume(solidGammaMonitorSupport0,
                                        GMSupportMaterial, "logicGammaMonitorSupport");
  new G4PVPlacement(0, G4ThreeVector(0.0, lxs->FloorSurfaceYpos+gmsupy/2.0, lxs->GammaCaloZpos),
                    logicGammaMonitorSupport, "GammaMonitorSupport", fLogicWorld, false, 0, lxs->OverlapTest);

// Construct gamma monitor metal holders
  G4double gmholddx = 1.0 *cm;
  G4double gmholdx = 1.5 * lgdy;
  G4double gmholdy = gmconty - sqrt(lxs->CalSupportR*lxs->CalSupportR - gmholdx*gmholdx);
  G4Box *solidGammaMonitorHolder0 = new G4Box("solidGammaMonitorHolder0", gmholddx/2.0, gmholdy/2.0, gmsupz/2.0);
  G4Tubs *solidCalSupportCut = new G4Tubs("solidCalSupportCut", 0.0,
                                       lxs->CalSupportR, gmsupz, 0.0, 2.0*M_PI);
  G4EllipticalTube *solidCalSupportCut1 = new G4EllipticalTube("solidCalSupportCut1", 0.3*gmholdy, 0.3*gmsupz, gmholddx);
  G4Transform3D transform0(G4RotationMatrix(), G4ThreeVector(-gmholdx, (gmconty-gmholdy/2.0), 0.0));
  G4SubtractionSolid *solidGammaMonitorHolder1 = new G4SubtractionSolid("solidGammaMonitorHolder1",
                                 solidGammaMonitorHolder0, solidCalSupportCut, transform0);
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI/2.0), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid *solidGammaMonitorHolder = new G4SubtractionSolid("solidGammaMonitorHolder",
                                 solidGammaMonitorHolder1, solidCalSupportCut1, transform1);
  G4LogicalVolume *logicGammaMonitorHolder = new G4LogicalVolume(solidGammaMonitorHolder,
                                        GMHolderMaterial, "logicGammaMonitorHolder");

//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),
//                        G4ThreeVector(-gmholdx, lxs->FloorSurfaceYpos + gmsupy + gmholdy/2.0, lxs->GammaCaloZpos),
//                        logicGammaMonitorHolder, "GammaMonitorHolder", logicGammaMonitorContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(gmholdx, lxs->FloorSurfaceYpos + gmsupy + gmholdy/2.0, lxs->GammaCaloZpos),
//                     logicGammaMonitorHolder, "GammaMonitorHolder", logicGammaMonitorContainer, false, 1, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), G4ThreeVector(-gmholdx, -(gmconty-gmholdy/2.0), 0.0),
                    logicGammaMonitorHolder, "GammaMonitorHolder", logicGammaMonitorContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(gmholdx, - (gmconty-gmholdy/2.0), 0.0),
                    logicGammaMonitorHolder, "GammaMonitorHolder", logicGammaMonitorContainer, false, 1, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->GammaCaloZpos),
                    logicGammaMonitorContainer, "GammaMonitorContainer", fLogicWorld, false, 0, lxs->OverlapTest);

//  to connect  2 beam pipes with different diameter
       
/* G4Tubs *solidConnect2BeamPipeGammaD = new G4Tubs("solidConnect2BeamPipeGammaD", lxs->BPipeRLG, lxs->BPipeR, lxs->BPipeThickness, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicConnect2BeamPipeGammaD = new G4LogicalVolume(solidConnect2BeamPipeGammaD, beamPipeMaterial, "logicConnect2BeamPipeGammaD");
 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lxs->BPipeThickness),
                               logicConnect2BeamPipeGammaD, "Connect2BeamPipeGammaD", fLogicWorld, false, 0, lxs->OverlapTest);

 G4double lpipe_gd = lxs->GammaBeamDumpZpos - lxs->ComptonDetZpos - lxs->GammaBeamDumpZ/2.0;
 G4Tubs *solidBeamPipeGammaD = new G4Tubs("solidBeamPipeGammaD", lxs->BPipeRLG-lxs->BPipeThickness, lxs->BPipeRLG, lpipe_gd/2.0, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicBeamPipeGammaD = new G4LogicalVolume(solidBeamPipeGammaD, beamPipeMaterial, "logicBeamPipeGammaD");
 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0),
                              logicBeamPipeGammaD, "BeamPipeGammaD", fLogicWorld, false, 0, lxs->OverlapTest);

 G4Tubs *solidBeamPipeGammaDVac = new G4Tubs("solidBeamPipeGammaDVac", 0.0, lxs->BPipeRLG-lxs->BPipeThickness, lpipe_gd/2.0, 0.0, 2.0*M_PI);
 G4LogicalVolume *logicBeamPipeGammaDVac = new G4LogicalVolume(solidBeamPipeGammaDVac, vacuumMaterial, "logicBeamPipeGammaDVac");
 new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->ComptonDetZpos + lpipe_gd/2.0),
                              logicBeamPipeGammaDVac, "BeamPipeGammaDVac", fLogicWorld, false, 0, lxs->OverlapTest);*/

   AddSegmentation();
}



void LxDetectorGammaCalo8::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("LeadGlass", lxs->LeadGlassX, lxs->LeadGlassY, 1, 1);
}
