//
/// \brief Implementation of the Scintillation screen and cherenkov detector
//

#include <algorithm>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4AssemblyVolume.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"
#include "LXSetUp.hh"
#include "LxAux.hh"
#include "LxDetector.hh"
#include "LxScintCherenkov.hh"


////////////////////////////////////////////////////////////////////////
//// LxSScreen


void LxSScreen::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

  // scintillator arm

  
     G4double sft = lxs->ScintFrameThickness;
     G4double scx = lxs->ScintX;
     G4double scy = lxs->ScintY;
     G4double scz1 = lxs->ScintBaseZ;
     G4double scz2 = lxs->ScintPhosphorZ;
     G4double scz3 = lxs->ScintFinishZ;
     G4double beam2stage = lxs->ComptonElectronBeamtoStageY;
     G4double sca = lxs->ScintAngle;
     G4double bsfwo = lxs->BremScintFrameWidthOffset;

     G4Material* env_mat = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
     G4Material* frame_mat = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial);
     G4Material* base_mat = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintBaseMaterial);
     G4Material* scint_mat = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintPhosphorMaterial);
     G4Material* finish_mat = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFinishMaterial);
     
     auto scintArmMotherSolid = new G4Box("scintArmMother", scx/2. + sft + lxs->ScintFrameBeamLoopX, beam2stage, beam2stage + 3*sft/2.);

     auto BremScintArmMotherSolid = new G4Box("BremScintArmMother", (scy/2. + 2*sft + bsfwo), cos(sca)*scx/2. + cos(sca)*sft + sin(sca)*sft/2., sin(sca)*scx/2. + sin(sca)*sft + cos(sca)*sft/2.);
     
  auto scintArmMotherLogical
    = new G4LogicalVolume(scintArmMotherSolid, env_mat,"scintArmMotherLogical");
  
  auto BremScintArmMotherLogical
    = new G4LogicalVolume(BremScintArmMotherSolid, env_mat,"BremScintArmMotherLogical");
  
  auto scintArmSolid 
    = new G4Box("scintArmBox", scx/2., scy/2., (scz1+scz2+scz3)/2.);
  
  auto scintArmLogical
    = new G4LogicalVolume(scintArmSolid, env_mat,"scintArmLogical");
  
  auto scintBaseSolid 
    = new G4Box("scintBaseBox", scx/2., scy/2., scz1/2.);
  
  auto scintBaseLogical
    = new G4LogicalVolume(scintBaseSolid, base_mat,"scintBaseLogical");
  
  auto scintPhosphorSolid 
    = new G4Box("scintPhosphorBox", scx/2., scy/2., scz2/2.);
  
  auto scintPhosphorLogical
    = new G4LogicalVolume(scintPhosphorSolid, scint_mat,"scintPhosphorLogical");
  
  auto scintFinishSolid 
    = new G4Box("scintFinishBox", scx/2., scy/2., scz3/2.);
  
  auto scintFinishLogical
    = new G4LogicalVolume(scintFinishSolid, finish_mat,"scintFinishLogical");
  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-(scz2+scz3)/2.),
                    scintBaseLogical,"scintBasePhysical",scintArmLogical,false,0,lxs->OverlapTest);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,(scz1 - scz3)/2.),
                    scintPhosphorLogical,"scintPhosphorPhysical",scintArmLogical,false,0,lxs->OverlapTest);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,(scz1+scz2)/2.),
                    scintFinishLogical,"scintFinishPhysical",scintArmLogical,false,0,lxs->OverlapTest);

  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                    scintArmLogical,"scintArmPhysical",scintArmMotherLogical,false,0,lxs->OverlapTest);

   auto scintFrameLegSolid 
    = new G4Box("scintFrameLegBox", sft/2., scy/2. + sft, sft/2.);

   auto scintFrameLegLogical
    = new G4LogicalVolume(scintFrameLegSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"scintFrameLegLogical");
  
   auto scintFrameBarSolid 
    = new G4Box("scintFrameBarBox", scx/2. + lxs->ScintFrameBeamLoopX/2., sft/2., sft/2.);

   auto BremScintFrameSupportLeg1Solid 
     = new G4Box("BremScintFrameSupportLeg1Box", sft/2., cos(sca)*scx/2. + cos(sca)*sft + sin(sca)*sft/2., sft/2.);   // vertical bit for brem region support

   auto BremScintFrameSupportLeg2Solid 
     = new G4Box("BremScintFrameSupportLeg2Box", sft/2., sft/2.,  sin(sca)*scx/2. + sin(sca)*sft + cos(sca)*sft/2. - sft/2.);   // vertical bit for brem region support

   auto BremScintFrameLegSolid 
    = new G4Box("scintFrameLegBox", sft/2., scy/2. + sft + bsfwo, sft/2.);

   auto BremScintFrameBarSolid 
    = new G4Box("scintFrameBarBox", scx/2., sft/2., sft/2.);

   auto scintFrameBarLogical
    = new G4LogicalVolume(scintFrameBarSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"scintFrameBarLogical");

   auto BremScintFrameLegLogical
    = new G4LogicalVolume(BremScintFrameLegSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"BremScintFrameLegLogical");

   auto BremScintFrameSupportLeg1Logical
    = new G4LogicalVolume(BremScintFrameSupportLeg1Solid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"BremScintFrameSupportLeg1Logical");

   auto BremScintFrameSupportLeg2Logical
    = new G4LogicalVolume(BremScintFrameSupportLeg2Solid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"BremScintFrameSupportLeg2Logical");

   auto BremScintFrameBarLogical
    = new G4LogicalVolume(BremScintFrameBarSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"BremScintFrameBarLogical");

     //45 degree supportive legs
       auto scintSupportLegRot = new G4RotationMatrix();
            scintSupportLegRot->rotateX(-45.*deg);
	    auto scintSupportflipY = new G4RotationMatrix();
	    scintSupportflipY->rotateY(180.*deg);

     auto scintSupportLegSolidA
       = new G4Box("scintSupportLegBoxA", sft/2. + 1.*mm, sft/2., sft/2.);

     auto scintSupportLegSolidB
       = new G4Box("scintSupportLegBoxB", sft/2., (sft/(2.*sqrt(2.)) + sqrt(2.)*beam2stage)/(2.), sft/(2.*sqrt(2.)));

     auto scintSupportLegSolidC
       = new G4Box("scintSupportLegBoxC", sft/2., sft/2., beam2stage/2.);

     auto scintSupportLegSolidD
       = new G4SubtractionSolid("scintSupportLegBoxD", scintSupportLegSolidB, scintSupportLegSolidA, scintSupportLegRot, G4ThreeVector(0., sft/(2.*sqrt(2.)) + sqrt(2.)*beam2stage/2., sft/(sqrt(2.)*2.)));
     
     // scintSupportLegRot->rotateX(90.*deg);

     auto scintSupportLegSolidE
       = new G4SubtractionSolid("scintSupportLegBoxE", scintSupportLegSolidD, scintSupportLegSolidA, scintSupportLegRot, G4ThreeVector(0., -sft/(2.*sqrt(2.)) - sqrt(2.)*beam2stage/2., sft/(sqrt(2.)*2.)));
     
     scintSupportLegRot->rotateX(270.*deg);

      auto scintSupportLegSolid
       = new G4UnionSolid("scintSupportLegBox", scintSupportLegSolidC, scintSupportLegSolidE, scintSupportLegRot, G4ThreeVector(0., beam2stage/2. - sft/2. + 1.*mm, +sft/2.));
         
     auto scintSupportLegLogical
    = new G4LogicalVolume(scintSupportLegSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintFrameMaterial),"scintSupportLegLogical");     
	 
  new G4PVPlacement(0,G4ThreeVector(-(scx + sft)/2.,0.,0.),
                    scintFrameLegLogical,"scintFrameLegPhysical",scintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(0,G4ThreeVector((scx + sft)/2.+lxs->ScintFrameBeamLoopX,0.,0.),
                    scintFrameLegLogical,"scintFrameLegPhysical",scintArmMotherLogical,false,1,lxs->OverlapTest);

  new G4PVPlacement(0,G4ThreeVector(lxs->ScintFrameBeamLoopX/2., (scy + sft)/2.,0.),
                    scintFrameBarLogical,"scintFrameBarPhysical",scintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(0,G4ThreeVector(lxs->ScintFrameBeamLoopX/2., -(scy + sft)/2.,0.),
                    scintFrameBarLogical,"scintFrameBarPhysical",scintArmMotherLogical,false,1,lxs->OverlapTest);
  
  new G4PVPlacement(0,G4ThreeVector(lxs->ScintFrameBeamLoopX/2., -(beam2stage)+sft/2.,beam2stage-sft/2.),
                    scintFrameBarLogical,"scintFrameBarPhysical",scintArmMotherLogical,false,2,lxs->OverlapTest);
  new G4PVPlacement(0,G4ThreeVector(lxs->ScintFrameBeamLoopX/2., -(beam2stage)+sft/2.,-beam2stage+sft/2.),
                    scintFrameBarLogical,"scintFrameBarPhysical",scintArmMotherLogical,false,3,lxs->OverlapTest);

  
  new G4PVPlacement(0,G4ThreeVector((scx +sft)/2. + lxs->ScintFrameBeamLoopX, -(beam2stage)+sft/2., (beam2stage)/2.),
                     scintSupportLegLogical,"scintSupportLegPhysical",scintArmMotherLogical,false,0,lxs->OverlapTest);
  
   new G4PVPlacement(0,G4ThreeVector(-(scx +sft)/2., -(beam2stage)+sft/2., (beam2stage)/2.),                                                                                              scintSupportLegLogical,"scintSupportLegPhysical",scintArmMotherLogical,false,1,lxs->OverlapTest);
  
   new G4PVPlacement(scintSupportflipY,G4ThreeVector((scx +sft)/2. + lxs->ScintFrameBeamLoopX, -(beam2stage)+sft/2., -(beam2stage)/2.),
                     scintSupportLegLogical,"scintSupportLegPhysical",scintArmMotherLogical,false,2,lxs->OverlapTest);
  
   new G4PVPlacement(scintSupportflipY,G4ThreeVector(-(scx +sft)/2., -(beam2stage)+sft/2., -(beam2stage)/2.),
                     scintSupportLegLogical,"scintSupportLegPhysical",scintArmMotherLogical,false,3,lxs->OverlapTest);
  
  auto scintRot = new G4RotationMatrix();
  auto cerenkovRot = new G4RotationMatrix();
  auto cameraRot = new G4RotationMatrix();
  auto brembasebarRot = new G4RotationMatrix();
  scintRot->rotateZ(-90.*deg);
  scintRot->rotateY(lxs->ScintAngle);
  cerenkovRot->rotateZ(-90.*deg);
  cerenkovRot->rotateY(lxs->CerenkovAngle);
  cameraRot->rotateX(lxs->ScintCameraAngle);
  brembasebarRot->rotateZ(-90.*deg);


  new G4PVPlacement(scintRot,G4ThreeVector(0.,0.,0.),
                    scintArmLogical,"BremScintArmPhysical",BremScintArmMotherLogical,false,1,lxs->OverlapTest);
  
  new G4PVPlacement(scintRot,G4ThreeVector(0., -cos(sca)*(scx+sft)/2.,  -sin(sca)*(scx+sft)/2.),
                    BremScintFrameLegLogical,"BremScintFrameLegPhysical",BremScintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(scintRot,G4ThreeVector(0., cos(sca)*(scx+sft)/2.,  +sin(sca)*(scx+sft)/2.),
                    BremScintFrameLegLogical,"BremScintFrameLegPhysical",BremScintArmMotherLogical,false,1,lxs->OverlapTest);

  new G4PVPlacement(scintRot,G4ThreeVector((scy + sft)/2.,0.,0.),
                    BremScintFrameBarLogical,"BremScintFrameBarPhysical",BremScintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(scintRot,G4ThreeVector(-(scy + sft)/2.,0.,0.),
                    BremScintFrameBarLogical,"BremScintFrameBarPhysical",BremScintArmMotherLogical,false,1,lxs->OverlapTest);
  
  new G4PVPlacement(0,G4ThreeVector(bsfwo + (scy + 3*sft)/2., 0., sin(sca)*scx/2. + sin(sca)*sft + cos(sca)*sft/2. - sft/2.),
                    BremScintFrameSupportLeg1Logical,"BremScintFrameSupportLeg1Physical",BremScintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(0,G4ThreeVector(-bsfwo - (scy + 3*sft)/2., 0., sin(sca)*scx/2. + sin(sca)*sft + cos(sca)*sft/2. - sft/2.),
                    BremScintFrameSupportLeg1Logical,"BremScintFrameSupportLeg1Physical",BremScintArmMotherLogical,false,1,lxs->OverlapTest);
  
  new G4PVPlacement(brembasebarRot,G4ThreeVector(0., -(cos(sca)*scx/2. + cos(sca)*sft + sin(sca)*sft/2. - sft/2.), sin(sca)*scx/2. + sin(sca)*sft + cos(sca)*sft/2. - sft/2.),
                    BremScintFrameLegLogical,"BremScintFrameBaseBarPhysical",BremScintArmMotherLogical,false,0,lxs->OverlapTest);
  
  new G4PVPlacement(0,G4ThreeVector(bsfwo + (scy + 3*sft)/2., -(cos(sca)*scx/2. + cos(sca)*sft + sin(sca)*sft/2. - sft/2.), -sft/2.),
                    BremScintFrameSupportLeg2Logical,"BremScintFrameSupportLeg2Physical",BremScintArmMotherLogical,false,0,lxs->OverlapTest);
  new G4PVPlacement(0,G4ThreeVector(-bsfwo - (scy + 3*sft)/2., -(cos(sca)*scx/2. + cos(sca)*sft + sin(sca)*sft/2. - sft/2.), -sft/2.),
                    BremScintFrameSupportLeg2Logical,"BremScintFrameSupportLeg2Physical",BremScintArmMotherLogical,false,1,lxs->OverlapTest);
  
  new G4PVPlacement(0,G4ThreeVector(0.,lxs->ScintXpos,lxs->ScintZpos),
                    BremScintArmMotherLogical,"BremScintArmMotherPhysical",fLogicWorld,false,0,lxs->OverlapTest);
  

//   G4double comptonelzpos = lxs->IPMagnetZpos + lxs->TypMBFieldLength/2.0  + lxs->ComptonElBackshift;
  G4double comptonelzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0  + lxs->ComptonElBackshift;

  G4double cescintz;
//   G4double dumpMagnetZ = lxs->TypMBFieldLength;
//   G4double gaptomag = 0.5 * (lxs->TypMBMagnetCoreZ - lxs->TypMBFieldLength);
  G4double dumpMagnetZ = lxs->FlashMFieldLength;
  G4double gaptomag = 0.5 * (lxs->FlashMagnetCoilZ - lxs->FlashMFieldLength);
  G4double bpipel = lxs->OPPPDetZtoMagnet - lxs->VacChambertoOPPPDetZGap - gaptomag;
  G4double scintmxpos = 0.0;
  if (abs(lxs->IPMagFieldY + 10000.0*gauss) < 1.0e-10) {
    scintmxpos = lxs->HICSScintilatorXPosMag1T;
  } else if (abs(lxs->IPMagFieldY + 20000.0*gauss) < 1.0e-10) {
    scintmxpos = lxs->HICSScintilatorXPosMag2T;
  } else {
    G4String msgstr("Field ");
    msgstr += std::to_string(lxs->IPMagFieldY/gauss) + G4String("Gauss is not supported!\n");
    G4Exception("DetectorConstruction::", "ConstructScintCerenkov()", FatalException, msgstr.c_str());
  }

  if (lxs->StickScintScreentoBeamWindow){
 
  cescintz = lxs->IPMagnetZpos + dumpMagnetZ/2.0 + bpipel/2.0 + gaptomag +
    + (lxs->OPPPDetZtoMagnet - lxs->VacChambertoOPPPDetZGap - gaptomag)/2. + lxs->ScintFrameThickness/2.;}

  else {   
    cescintz = comptonelzpos + (scz1 + scz2 + scz3)/2.0; // scy/2. + sft/sqrt(2.);}
  }

  new G4PVPlacement(0, G4ThreeVector(-(lxs->ScintX/2.0 + scintmxpos), 0.0, cescintz), scintArmMotherLogical, "HICSElectronScintillator", fLogicWorld, false, 0, lxs->OverlapTest);
  
  fDetector->AddSensorSegmentation("scintPhosphorPhysical", lxs->ScintX, lxs->ScintY, 2000, 200);

  //make camera mother volume

  auto *ScintCameraMotherBox = new G4Box("Scintillator Camera Mother Box", lxs->ScintCameraX/2. + lxs->ScintCameraShieldThickness, lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness, (lxs->ScintCameraZ + lxs->ScintCameraApertureZ + lxs->ScintCameraShieldThickness)/2.);

  auto ScintCameraMotherBoxLogical = new G4LogicalVolume(ScintCameraMotherBox,G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial),"ScintCameraMotherBoxLogical",0,0,0);

 
  //make main camera box

  auto *ScintCameraBox = new G4Box("Scintillator Camera Mother Box", lxs->ScintCameraX/2., lxs->ScintCameraY/2., lxs->ScintCameraZ/2.);

  auto ScintCameraBoxLogical = new G4LogicalVolume(ScintCameraBox,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraMaterial),"ScintCameraBoxLogical",0,0,0);
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., -lxs->ScintCameraShieldThickness/2. + lxs->ScintCameraApertureZ/2.), ScintCameraBoxLogical, "ScintCameraBoxPhysical", ScintCameraMotherBoxLogical, false,69, lxs->OverlapTest);
  
  //Make metallic aperture

    G4Tubs *ScintCameraApertureOuter = new G4Tubs("Scintillator Camera Aperture Rim", lxs->ScintCameraApertureInner/2., lxs->ScintCameraApertureOuter/2., lxs->ScintCameraApertureZ/2., 0.0, 2.0*M_PI);
    
  auto ScintCameraApertureOuterLogical = new G4LogicalVolume(ScintCameraApertureOuter,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraMaterial),"ScintCameraApertureOuterLogical",0,0,0);
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., -lxs->ScintCameraShieldThickness/2. - lxs->ScintCameraZ/2.), ScintCameraApertureOuterLogical, "ScintCameraApertureOuterPhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);

  
  // make aperture inner tube

    G4Tubs *ScintCameraApertureInner = new G4Tubs("Scintillator Camera Aperture Inner", 0., lxs->ScintCameraApertureInner/2., lxs->ScintCameraApertureZ/2., 0.0, 2.0*M_PI);
    
  auto ScintCameraApertureInnerLogical = new G4LogicalVolume(ScintCameraApertureInner,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraApertureMaterial),"ScintCameraApertureInnerLogical",0,0,0);
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., -lxs->ScintCameraShieldThickness/2. - lxs->ScintCameraZ/2.), ScintCameraApertureInnerLogical, "ScintCameraApertureInnerPhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);
  
  // make inner, outer and top shielding

  //  lxs->ScintCameraX/2. + lxs->ScintCameraShieldThickness,
  //lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness/2.,
  //(lxs->ScintCameraZ + lxs->ScintCameraApertureZ + lxs->ScintCameraShieldThickness)/2.)
  
  auto *ScintCameraShieldSide = new G4Box("Scintillator Camera Shield Side Box", lxs->ScintCameraShieldThickness/2.,   lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness, lxs->ScintCameraZ/2.);
  auto ScintCameraShieldSideLogical = new G4LogicalVolume(ScintCameraShieldSide,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraShieldMaterial),"ScintCameraShieldSideLogical",0,0,0);

  new G4PVPlacement(0, G4ThreeVector((lxs->ScintCameraX + lxs->ScintCameraShieldThickness)/2., 0., -lxs->ScintCameraShieldThickness/2. + lxs->ScintCameraApertureZ/2.), ScintCameraShieldSideLogical, "ScintCameraShieldSidePhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);    
  new G4PVPlacement(0, G4ThreeVector(-(lxs->ScintCameraX + lxs->ScintCameraShieldThickness)/2., 0., -lxs->ScintCameraShieldThickness/2. + lxs->ScintCameraApertureZ/2.), ScintCameraShieldSideLogical, "ScintCameraShieldSidePhysical", ScintCameraMotherBoxLogical, false,1, lxs->OverlapTest);
  
  auto *ScintCameraShieldTop = new G4Box("Scintillator Camera Shield Top Box", lxs->ScintCameraX/2., lxs->ScintCameraShieldThickness/2., lxs->ScintCameraZ/2.);
  auto ScintCameraShieldTopLogical = new G4LogicalVolume(ScintCameraShieldTop,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraShieldMaterial),"ScintCameraShieldTopLogical",0,0,0);
  
  new G4PVPlacement(0, G4ThreeVector(0., lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness/2., -lxs->ScintCameraShieldThickness/2. + lxs->ScintCameraApertureZ/2.), ScintCameraShieldTopLogical, "ScintCameraShieldTopPhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);



   new G4PVPlacement(0, G4ThreeVector(0., -lxs->ScintCameraY/2. - lxs->ScintCameraShieldThickness/2., -lxs->ScintCameraShieldThickness/2. + lxs->ScintCameraApertureZ/2.), ScintCameraShieldTopLogical, "ScintCameraShieldTopPhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);

  
  auto *ScintCameraShieldBack = new G4Box("Scintillator Camera Shield Back Box", lxs->ScintCameraX/2. + lxs->ScintCameraShieldThickness,   lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness, lxs->ScintCameraShieldThickness/2.);
  auto ScintCameraShieldBackLogical = new G4LogicalVolume(ScintCameraShieldBack,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraShieldMaterial),"ScintCameraShieldBackLogical",0,0,0);


  new G4PVPlacement(0, G4ThreeVector(0., 0., lxs->ScintCameraZ/2. + lxs->ScintCameraApertureZ/2.), ScintCameraShieldBackLogical, "ScintCameraShieldBackPhysical", ScintCameraMotherBoxLogical, false,0, lxs->OverlapTest);

  // place cameras
  // on top of scaffold, at 1/4 and 3/4 screen x-position 

  
  new G4PVPlacement(cameraRot, G4ThreeVector(-(scintmxpos + 3*lxs->ScintX/4.), lxs->ScintCameraYpos //+ lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness // + sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ)
					     , cescintz + 1.*m  + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. // - sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ))
					     ), ScintCameraMotherBoxLogical, "HICSScintCameraMotherBoxPhysical", fLogicWorld, false, 2, lxs->OverlapTest); 

  
  new G4PVPlacement(cameraRot, G4ThreeVector(-(scintmxpos + lxs->ScintX/4.), lxs->ScintCameraYpos// + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness // + sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ)
					     , cescintz + 1.*m + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. //- sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ))
					     ), ScintCameraMotherBoxLogical, "HICSScintCameraMotherBoxPhysical", fLogicWorld, false, 3, lxs->OverlapTest); 

  
    fDetector->AddSensorSegmentation("ScintCameraApertureInnerPhysical", lxs->ScintCameraApertureInner, lxs->ScintCameraApertureInner, 1, 1);

		    // support solid


    auto *ScintCameraPlatformBoxforsubtraction = new G4Box("Scintillator Camera Scaffolding Platform Box", lxs->ScintX/2. + 5.*cm, lxs->ScintCameraPlatformZ/2., lxs->ScintCameraSupportThickness/2.);

    auto *ScintCameraSupportLeg = new G4Tubs("Scintillator Camera Support Leg Tube", 0., lxs->ScintCameraSupportThickness/2., (lxs->CeilingSurfaceYpos - lxs->ScintCameraYpos + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness)/2.,  0.0, 2.0*M_PI);

    auto ScintCameraSupport = new G4SubtractionSolid("Scintillator Camera Support", ScintCameraSupportLeg, ScintCameraPlatformBoxforsubtraction, cameraRot, G4ThreeVector(0, 0, (lxs->CeilingSurfaceYpos - lxs->ScintCameraYpos + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness)/2. -
																					  (1- cos(lxs->ScintCameraAngle))*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.) -
																					  sin(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. -lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.)));

    auto ScintCameraSupportLogical = new G4LogicalVolume(ScintCameraSupport,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraSupportMaterial),"HICSScintCameraScaffoldLegLogical",0,0,0);
	       
  auto cameraSupportRot = new G4RotationMatrix();
    cameraSupportRot->rotateZ(1*M_PI);
    cameraSupportRot->rotateY(1*M_PI);
    cameraSupportRot->rotateX(M_PI/2.);

		    
    new G4PVPlacement(cameraSupportRot, G4ThreeVector(-(scintmxpos + lxs->ScintX/2.), lxs->CeilingSurfaceYpos - (lxs->CeilingSurfaceYpos - lxs->ScintCameraYpos + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness)/2.,
						    //				    - (lxs->CeilingSurfaceYpos-lxs->ScintCameraYpos) - lxs->ScintCameraY/2. - lxs->ScintCameraShieldThickness),
						     (cescintz + 1.*m + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. +
					       cos(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. - lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.) +
						      sin(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.))), ScintCameraSupportLogical, "HICSScintCameraSupport", fLogicWorld, false,0, lxs->OverlapTest);

		      // // make scaffold platform

   auto *ScintCameraSupportPlatform = new G4Box("Scintillator Camera Scaffolding Platform Box", lxs->ScintX/2. + 5.*cm, lxs->ScintCameraSupportThickness/2., lxs->ScintCameraPlatformZ/2.);

   auto ScintCameraSupportPlatformLogical = new G4LogicalVolume(ScintCameraSupportPlatform,G4NistManager::Instance()->FindOrBuildMaterial(lxs->ScintCameraSupportMaterial),"ScintCameraScaffoldPlatformLogical",0,0,0);

   new G4PVPlacement(cameraRot, G4ThreeVector(-(scintmxpos + lxs->ScintX/2.),
					      lxs->ScintCameraYpos -
					      cos(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.) +
					      sin(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. -lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.),
					      (cescintz + 1.*m + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. +
					       cos(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. - lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.) +
					      sin(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.))), ScintCameraSupportPlatformLogical, "HICSScintCameraSupportPlatformPhysical", fLogicWorld, false,0, lxs->OverlapTest);


   // Brem location cameras, observe normal to screen surface (coverage may be blocked by supporting table)
  // auto BremCameraRot = new G4RotationMatrix();

  // //cameraRot->rotateX( - lxs->ScintCameraAngle);

  //  BremCameraRot->rotateZ(M_PI/2.);

  //  BremCameraRot->rotateX( - lxs->ScintCameraAngle);

  //  BremCameraRot->rotateY( - lxs->ScintAngle);

  //  new G4PVPlacement(BremCameraRot, G4ThreeVector(- 60. * cm, lxs->ScintXpos - cos(lxs->ScintAngle)*lxs->ScintX/4. -  1.5 * tan(lxs->ScintAngle) * m, lxs->ScintZpos + 1.5 * m), ScintCameraMotherBoxLogical, "BremScintCameraMotherBoxPhysical", fLogicWorld, false, 0, lxs->OverlapTest); 

  //  new G4PVPlacement(BremCameraRot, G4ThreeVector(- 60. * cm, lxs->ScintXpos + cos(lxs->ScintAngle)*lxs->ScintX/4. -  1.5 * tan(lxs->ScintAngle) * m, lxs->ScintZpos + 1.5 * m), ScintCameraMotherBoxLogical, "BremScintCameraMotherBoxPhysical", fLogicWorld, false, 1, lxs->OverlapTest); 

  //  auto BremCameraSupportRot = new G4RotationMatrix();

  //  BremCameraSupportRot->rotateZ(M_PI/2.);

  //  BremCameraSupportRot->rotateX( - lxs->ScintCameraAngle);

  // // new G4PVPlacement(cameraRot, G4ThreeVector(-(scintmxpos + lxs->ScintX/4.), lxs->ScintCameraYpos// + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness // + sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ)
  // // 					     , cescintz + 1.*m + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. //- sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ))
  // // 					     ), ScintCameraMotherBoxLogical, "ScintCameraMotherBoxPhysical", fLogicWorld, false,0, lxs->OverlapTest); 

  
  //  new G4PVPlacement(BremCameraSupportRot, G4ThreeVector(- 60. * cm +
  // 					      cos(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.) - sin(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. -lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.),

  // 					      lxs->ScintXpos - 1.5 * tan(lxs->ScintAngle) * m,

  // 					      ( lxs->ScintZpos + 1.5 * m) +
  // 					       cos(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. - lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.) +
  // 					      sin(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.)), ScintCameraSupportPlatformLogical, "BremScintCameraSupportPlatformPhysical", fLogicWorld, false,0, lxs->OverlapTest);


   // Brem location cameras, observe at 1/4 and 3/4 y (height) of screen.

  auto BremCameraRot = new G4RotationMatrix();

  //cameraRot->rotateX( - lxs->ScintCameraAngle);

   BremCameraRot->rotateZ(M_PI/2.);

   BremCameraRot->rotateX( - lxs->ScintCameraAngle);

   // BremCameraRot->rotateY( - lxs->ScintAngle);

   new G4PVPlacement(BremCameraRot, G4ThreeVector(- 60. * cm, lxs->ScintXpos - cos(lxs->ScintAngle)*lxs->ScintX/4., lxs->ScintZpos + 1.5 * m), ScintCameraMotherBoxLogical, "BremScintCameraMotherBoxPhysical", fLogicWorld, false, 0, lxs->OverlapTest); 

   new G4PVPlacement(BremCameraRot, G4ThreeVector(- 60. * cm, lxs->ScintXpos + cos(lxs->ScintAngle)*lxs->ScintX/4., lxs->ScintZpos + 1.5 * m), ScintCameraMotherBoxLogical, "BremScintCameraMotherBoxPhysical", fLogicWorld, false, 1, lxs->OverlapTest); 

   auto BremCameraSupportRot = new G4RotationMatrix();

   BremCameraSupportRot->rotateZ(M_PI/2.);

   BremCameraSupportRot->rotateX( - lxs->ScintCameraAngle);

  // new G4PVPlacement(cameraRot, G4ThreeVector(-(scintmxpos + lxs->ScintX/4.), lxs->ScintCameraYpos// + lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness // + sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ)
  // 					     , cescintz + 1.*m + (lxs->ScintCameraShieldThickness - lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2. //- sin(lxs->ScintCameraAngle)*((lxs->ScintCameraZ+lxs->ScintCameraApertureZ+lxs->ScintCameraShieldThickness)/2. - lxs->ScintCameraApertureZ))
  // 					     ), ScintCameraMotherBoxLogical, "ScintCameraMotherBoxPhysical", fLogicWorld, false,0, lxs->OverlapTest); 

  
   new G4PVPlacement(BremCameraSupportRot, G4ThreeVector(- 60. * cm +
					      cos(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.) - sin(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. -lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.),

					      lxs->ScintXpos,

					      ( lxs->ScintZpos + 1.5 * m) +
					       cos(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. - lxs->ScintCameraZ - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.) +
					      sin(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness + lxs->ScintCameraSupportThickness/2.)), ScintCameraSupportPlatformLogical, "BremScintCameraSupportPlatformPhysical", fLogicWorld, false,0, lxs->OverlapTest);

   

  ConstructSupportAssembly();
  }



void LxSScreen::ConstructSupportAssembly()
{
  // Supporting Brem detector table
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

  G4double scdz = lxs->ScintBaseZ + lxs->ScintPhosphorZ + lxs->ScintFinishZ;

  G4int cc = lxs->CerenkovChannels;
  G4int ccl = lxs->CerenkovChannelLayers;
  G4double csxf = lxs->CerenkovStrawXFrequency;
  G4double cszf = lxs->CerenkovStrawZFrequency;
  G4double cslo = lxs->CerenkovStrawLayerOffset;
  G4double sir = lxs->CerenkovStrawInnerRadius;
  G4double sglt = lxs->CerenkovStrawGraphiteLayerThickness;
  G4double salt = lxs->CerenkovStrawAlLayerThickness;
  G4double sklt = lxs->CerenkovStrawKaptonLayerThickness;
  G4double splt = lxs->CerenkovStrawPolyurethaneLayerThickness;
  G4double cspb = lxs->CerenkovStrawPlateBuffer;
  G4double cbt = lxs->CerenkovBoxThickness;
  G4double cspt = lxs->CerenkovShieldingPlateThickness;

  G4double cherenkovdy = 2.0 *((cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);
  G4double cherenkovdz = 2.0 * ((ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);

  G4double scintbypos = lxs->ScintXpos - cos(lxs->ScintAngle)*lxs->ScintX/2.0 - sin(lxs->ScintAngle)*scdz/2.0;
  G4double cherbypos = lxs->CerenkovXpos - cos(lxs->CerenkovAngle)*cherenkovdy/2.0
                                        - sin(lxs->CerenkovAngle)*cherenkovdz/2.0;

  G4double tblx = 2.0*std::max(lxs->CerenkovTotalBoxHeight, lxs->ScintY);
  G4double tblz = (lxs->CerenkovZpos - lxs->ScintZpos) + sin(lxs->CerenkovAngle)*cherenkovdy/2.0 + cos(lxs->CerenkovAngle)*cherenkovdz/2.0 + sin(lxs->ScintAngle)*lxs->ScintX/2. + sin(lxs->ScintAngle)*lxs->ScintFrameThickness + cos(lxs->ScintAngle)*lxs->ScintFrameThickness/2. + 10.*cm;

  G4double baseypos = std::min(scintbypos, cherbypos);
  G4double basezpos = lxs->CerenkovZpos + sin(lxs->CerenkovAngle)*cherenkovdy/2.0 + cos(lxs->CerenkovAngle)*cherenkovdz/2.0 - tblz/2. + 5.*cm;

  // G4double basezpos = 0.5 * (lxs->ScintZpos - sin(lxs->ScintAngle)*lxs->ScintX/2.0
  //                            + lxs->CerenkovZpos - sin(lxs->CerenkovAngle)*cherenkovdz/2.0);
  
  G4double ypestal = 0.5*m;
  G4double tblhight = baseypos - ypestal - lxs->FloorSurfaceYpos;

  G4AssemblyVolume *tablesupport = LxAux::BuildTable("BremDetTable", tblx, tblhight, tblz, 4);
  G4ThreeVector trsupport(0.0, baseypos, basezpos);
  tablesupport->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);

  G4LogicalVolume *pedstal = LxAux::BuildPedestal("BremDet", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.5*ypestal + lxs->FloorSurfaceYpos, basezpos), pedstal,
                      "BremDetPedestal", fLogicWorld, false, 0, lxs->OverlapTest);

// Brem. Scintillation Cameras support
  G4double cambardy = 2.0*(lxs->ScintX/2. + 5.*cm);

  G4double cambarposx = -60. * cm + cos(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness 
                        + lxs->ScintCameraSupportThickness/2.)
                        - sin(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2. -lxs->ScintCameraZ
                                                      - lxs->ScintCameraShieldThickness + (lxs->ScintCameraShieldThickness
                                                      + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.);
  G4double cambarposy = lxs->ScintXpos; // if changing to point cameras so normal to screen, add:  - 1.5 * tan(lxs->ScintAngle) * m;

  G4double cambarposz = (lxs->ScintZpos + 1.5 * m) + cos(lxs->ScintCameraAngle)*(lxs->ScintCameraPlatformZ/2.
                            - lxs->ScintCameraZ - lxs->ScintCameraShieldThickness
                            + (lxs->ScintCameraShieldThickness + lxs->ScintCameraApertureZ + lxs->ScintCameraZ)/2.)
                            + sin(lxs->ScintCameraAngle)*(lxs->ScintCameraY/2. + lxs->ScintCameraShieldThickness
                                                          + lxs->ScintCameraSupportThickness/2.);
  tblx = 4.0 * (lxs->ScintCameraPlatformZ * sin(lxs->ScintCameraAngle)
                + lxs->ScintCameraSupportThickness * cos(lxs->ScintCameraAngle));
  tblhight = cambarposy - cambardy/2.0 - lxs->FloorSurfaceYpos;
  tblz = 3.0 * (lxs->ScintCameraPlatformZ * cos(lxs->ScintCameraAngle)
                + lxs->ScintCameraSupportThickness * sin(lxs->ScintCameraAngle));

  G4AssemblyVolume *cambartable = LxAux::BuildTable("BremScintCameraSupportTable", tblx, tblhight, tblz, 3);
  G4ThreeVector trcambarsupport(cambarposx, cambarposy-cambardy/2.0, cambarposz);
  cambartable->MakeImprint(fLogicWorld, trcambarsupport, 0, 0, lxs->OverlapTest);
}


////////////////////////////////////////////////////////////////////////
//// LxSScreen

void LxCherenkov::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

 //  // Instantiating Cerenkov detector setup, keeping commonly used quantities in local variables for brevity & speed
     
	G4double csh = lxs->CerenkovStrawHeight;

	//	G4double cdt = lxs->CerenkovDeviceOuterThickness;

  G4double sir = lxs->CerenkovStrawInnerRadius;
  G4double sglt = lxs->CerenkovStrawGraphiteLayerThickness;
  G4double salt = lxs->CerenkovStrawAlLayerThickness;
  G4double sklt = lxs->CerenkovStrawKaptonLayerThickness;
  G4double splt = lxs->CerenkovStrawPolyurethaneLayerThickness;
  G4double spt = lxs->CerenkovShieldingPlateThickness;

  auto CerenkovStrawMotherSolid = new G4Tubs("CerenkovStrawMotherSolid", 0., sir+2*(sglt+salt+sklt+splt), csh/2., 0., 2.0*M_PI); // mother volume for one complete straw - duplicate this volume for each straw

  //Creating layers of the straw substructure
  
   auto CerenkovStrawInnerSolid = new G4Tubs("CerenkovStrawInnerSolid", 0., sir, csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawGraphiteInnerLayerSolid = new G4Tubs("CerenkovStrawGraphiteInnerLayerSolid", sir, sir+sglt, csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawAlInnerLayerSolid = new G4Tubs("CerenkovStrawAlInnerLayerSolid", sir+sglt, sir+sglt+salt, csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawKaptonInnerLayerSolid = new G4Tubs("CerenkovStrawKaptonInnerLayerSolid", sir+sglt+salt, sir+sglt+salt+sklt, csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawPolyurethaneLayerSolid = new G4Tubs("CerenkovStrawPolyurethaneLayerSolid", sir+sglt+salt+sklt, sir+sglt+salt+sklt+2*splt, csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawKaptonOuterLayerSolid = new G4Tubs("CerenkovStrawKaptonOuterLayerSolid", sir+sglt+salt+sklt+2*splt, sir+sglt+salt+2*(sklt+splt), csh/2., 0., 2.0*M_PI);
   
   auto CerenkovStrawAlOuterLayerSolid = new G4Tubs("CerenkovStrawAlOuterLayerSolid", sir+sglt+salt+2*(sklt+splt), sir+sglt+2*(salt+sklt+splt), csh/2., 0., 2.0*M_PI); 

   auto CerenkovStrawGraphiteOuterLayerSolid = new G4Tubs("CerenkovStrawGraphiteOuterLayerSolid", sir+sglt+2*(salt+sklt+splt), sir+2*(sglt+salt+sklt+splt), csh/2., 0., 2.0*M_PI);

   // creating logical volumes for these straw layers including assigning materials
   
   auto CerenkovStrawMotherLogical = new G4LogicalVolume(CerenkovStrawMotherSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovMedium), "CerenkovStrawMotherLogical");

   auto CerenkovStrawInnerLogical = new G4LogicalVolume(CerenkovStrawInnerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovMedium), "CerenkovStrawInnerLogical");

   auto CerenkovStrawGraphiteInnerLayerLogical = new G4LogicalVolume(CerenkovStrawGraphiteInnerLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawGraphiteLayerMaterial), "CerenkovStrawGraphiteInnerLayerLogical");

   auto CerenkovStrawAlInnerLayerLogical = new G4LogicalVolume(CerenkovStrawAlInnerLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawAlLayerMaterial), "CerenkovStrawAlInnerLayerLogical");
   
   auto CerenkovStrawKaptonInnerLayerLogical = new G4LogicalVolume(CerenkovStrawKaptonInnerLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawKaptonLayerMaterial), "CerenkovStrawKaptonInnerLayerLogical");

   auto CerenkovStrawPolyurethaneLayerLogical = new G4LogicalVolume(CerenkovStrawPolyurethaneLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawPolyurethaneLayerMaterial), "CerenkovStrawPolyurethaneLayerLogical");

   auto CerenkovStrawKaptonOuterLayerLogical = new G4LogicalVolume(CerenkovStrawKaptonOuterLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawKaptonLayerMaterial), "CerenkovStrawKaptonOuterLayerLogical");

   auto CerenkovStrawAlOuterLayerLogical = new G4LogicalVolume(CerenkovStrawAlOuterLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawAlLayerMaterial), "CerenkovStrawAlOuterLayerLogical");
   
   auto CerenkovStrawGraphiteOuterLayerLogical = new G4LogicalVolume(CerenkovStrawGraphiteOuterLayerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovStrawGraphiteLayerMaterial), "CerenkovStrawGraphiteOuterLayerLogical");

  
   //  auto vertical =  new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0);

   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawInnerLogical, "CerenkovStrawInnerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);
   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawGraphiteInnerLayerLogical, "CerenkovStrawGraphiteInnerLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);
   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawAlInnerLayerLogical, "CerenkovStrawAlInnerLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);
   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawKaptonInnerLayerLogical, "CerenkovStrawKaptonInnerLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);

   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawPolyurethaneLayerLogical, "CerenkovStrawPolyurethaneLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);

   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawKaptonOuterLayerLogical, "CerenkovStrawKaptonOuterLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);

   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawAlOuterLayerLogical, "CerenkovStrawAlOuterLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);
   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), CerenkovStrawGraphiteOuterLayerLogical, "CerenkovStrawGraphiteOuterLayerPhysical", CerenkovStrawMotherLogical, 0, 0, lxs->OverlapTest);

 //vertically oriented

   //   new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0), G4ThreeVector(0.,0.,0.), CerenkovStrawMotherLogical, "CerenkovStrawMotherPhysical", fLogicWorld, 0,0, lxs->OverlapTest);


   // Constructing housing volumes for straws

        G4int cc = lxs->CerenkovChannels;
	G4int ccl = lxs->CerenkovChannelLayers;
	G4double cbt = lxs->CerenkovBoxThickness;
	G4double csxf = lxs->CerenkovStrawXFrequency;
	G4double cszf = lxs->CerenkovStrawZFrequency;
	G4double cslo = lxs->CerenkovStrawLayerOffset;
	G4double cspb = lxs->CerenkovStrawPlateBuffer;
	G4double cebt = lxs->CerenkovElectronicsBoardThickness;
	G4double cst = lxs->CerenkovSupportThickness;
	G4double cspt = lxs->CerenkovShieldingPlateThickness;
	
	auto CerenkovMotherSolid = new G4Box("CerenkovMotherSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt, lxs->CerenkovTotalBoxHeight/2., (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);

	//auto CerenkovMotherSolid = new G4Box("CerenkovMotherSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt, lxs->CerenkovTotalBoxHeight/2., (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt);

	auto CerenkovBoxOuterSolid = new G4Box("CerenkovBoxOuterSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt, lxs->CerenkovTotalBoxHeight/2., (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt);

	auto CerenkovBoxInnerSolid = new G4Box("CerenkovBoxInnerSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb, lxs->CerenkovTotalBoxHeight/2. - cbt,  (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb);

	auto CerenkovBoxBeamWindowSolid = new G4Box("CerenkovBoxBeamWindowSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb, lxs->CerenkovBeamWindowY/2.,  (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt - lxs->CerenkovBeamWindowThickness);

	auto CerenkovBoxSolid1 = new G4SubtractionSolid("CerenkovBoxSolid", CerenkovBoxOuterSolid, CerenkovBoxInnerSolid, 0, G4ThreeVector(0,0,0));
			
	auto CerenkovBoxSolid2 = new G4SubtractionSolid("CerenkovBoxSolid", CerenkovBoxSolid1, CerenkovBoxBeamWindowSolid, 0, G4ThreeVector(0,lxs->ComptonElectronBeamtoStageY - lxs->CerenkovTotalBoxHeight/2.,0));							    

	auto CerenkovSupportPlateSolid = new G4Box("CerenkovSupportPlateSolid",  (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb, lxs->CerenkovSupportThickness/2., (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb);	

	auto CerenkovElectronicsBoardSolid =  new G4Box("CerenkovElectronicsBoardSolid",  (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt), lxs->CerenkovElectronicsBoardThickness/2., (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt));

	auto CerenkovShieldingPlateUpperSolid = new G4Box("CerenkovShieldingPlateUpperSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt, (lxs->CerenkovTotalBoxHeight - lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2., cspt/2.);
	
	auto CerenkovShieldingPlateLowerSolid = new G4Box("CerenkovShieldingPlateLowerSolid", (cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt, (lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2., cspt/2.);
	
        auto CerenkovShieldingPlateSideSolid = new G4Box("CerenkovShieldingPlateLowerSolid", cspt/2., lxs->CerenkovTotalBoxHeight/2.,  (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt);

	
	auto CerenkovMotherLogical = new G4LogicalVolume(CerenkovMotherSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial), "CerenkovMotherLogical");
        auto CerenkovBoxLogical = new G4LogicalVolume(CerenkovBoxSolid2, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovMetal), "CerenkovBoxLogical");
	auto CerenkovBoxInnerLogical = new G4LogicalVolume(CerenkovBoxInnerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovMedium), "CerenkovBoxInnerLogical");

	auto CerenkovSupportPlateLogical = new G4LogicalVolume(CerenkovSupportPlateSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovMetal), "CerenkovSupportPlateLogical");
	auto CerenkovElectronicsBoardLogical = new G4LogicalVolume(CerenkovElectronicsBoardSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovElectronicsBoardMaterial), "CerenkovElectronicsBoardLogical");


	auto CerenkovShieldingPlateUpperLogical = new G4LogicalVolume(CerenkovShieldingPlateUpperSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovShieldingPlateMaterial), "CerenkovShieldingPlateUpperLogical");
	
	auto CerenkovShieldingPlateLowerLogical = new G4LogicalVolume(CerenkovShieldingPlateLowerSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovShieldingPlateMaterial), "CerenkovShieldingPlateLowerLogical");
	
	auto CerenkovShieldingPlateSideLogical = new G4LogicalVolume(CerenkovShieldingPlateSideSolid, G4NistManager::Instance()->FindOrBuildMaterial(lxs->CerenkovShieldingPlateMaterial), "CerenkovShieldingPlateSideLogical");

	
	
	

        G4double straw_x = -(cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2.;
	G4double straw_y = - lxs->CerenkovTotalBoxHeight/2. + csh/2. +  lxs->CerenkovSupportThickness + lxs->CerenkovSupportPosYtoStage;
	G4double straw_z = -(ccl-1) * cszf/2.;
	G4RotationMatrix* straw_vertical = new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0);

	  // iteratively building array of straws in CerenkovBoxLogical volume
	  
	for (int i=0; i<ccl; i++){
	  straw_x = -(cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. - i * cslo;
	for (int j=0; j<cc/ccl; j++){
	  new G4PVPlacement(straw_vertical, G4ThreeVector(straw_x, straw_y, straw_z), CerenkovStrawMotherLogical, "CerenkovStrawPhysical", CerenkovBoxInnerLogical, false, i*cc/ccl + j, lxs->OverlapTest);


	  // G4cout << "incepting straw number: " << i*cc/ccl + j << "at pos. (xyz): " << straw_x << ", " << straw_y << ", "<< straw_z << ", " << G4endl;
	  
	  straw_x += csxf;
//                     "ComptonCerenkovBox", fLogicWorld, false, 0, lxs->OverlapTest);
	}

	straw_z += cszf;
	}

	new G4PVPlacement(0, G4ThreeVector(0., straw_y - csh/2. - lxs->CerenkovSupportThickness/2., 0.), CerenkovSupportPlateLogical, "CerenkovSupportPlatePhysical", CerenkovBoxInnerLogical, false, 0, lxs->OverlapTest);
        
	 new G4PVPlacement(0, G4ThreeVector(0., straw_y + csh/2. + lxs->CerenkovSupportThickness/2., 0.), CerenkovSupportPlateLogical, "CerenkovSupportPlatePhysical", CerenkovBoxInnerLogical, false, 1, lxs->OverlapTest);
	 
	 new G4PVPlacement(0, G4ThreeVector(0., straw_y + csh/2. + lxs->CerenkovSupportThickness + lxs->CerenkovElectronicsBoardThickness/2., 0.), CerenkovElectronicsBoardLogical, "CerenkovElectronicsBoardPhysical", CerenkovBoxInnerLogical, false, 0, lxs->OverlapTest);  

	 // Optional Lead plates for further shielding of Cherenkov device

	 new G4PVPlacement(0, G4ThreeVector(0, lxs->CerenkovTotalBoxHeight/2. - (lxs->CerenkovTotalBoxHeight - lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2.,  (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2.), CerenkovShieldingPlateUpperLogical,
                     "CerenkovShieldingPlateUpperPhysical", CerenkovMotherLogical, false, 0, lxs->OverlapTest);

	 new G4PVPlacement(0, G4ThreeVector(0, lxs->CerenkovTotalBoxHeight/2. - (lxs->CerenkovTotalBoxHeight - lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2.,  -((ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2.)), CerenkovShieldingPlateUpperLogical,
                     "CerenkovShieldingPlateUpperPhysical", CerenkovMotherLogical, false, 1, lxs->OverlapTest);

	 new G4PVPlacement(0, G4ThreeVector(0, -lxs->CerenkovTotalBoxHeight/2. + (lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2.,  (ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2.), CerenkovShieldingPlateLowerLogical,
                     "CerenkovShieldingPlateLowerPhysical", CerenkovMotherLogical, false, 0, lxs->OverlapTest);

	 new G4PVPlacement(0, G4ThreeVector(0, -lxs->CerenkovTotalBoxHeight/2. + (lxs->ComptonElectronBeamtoStageY - lxs->CerenkovBeamWindowY/2.)/2.,  -((ccl-1)*cszf/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2.)), CerenkovShieldingPlateLowerLogical,
                     "CerenkovShieldingPlateLowerPhysical", CerenkovMotherLogical, false, 1, lxs->OverlapTest);

	 new G4PVPlacement(0, G4ThreeVector((cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2., 0, 0), CerenkovShieldingPlateSideLogical,
                     "CerenkovShieldingPlateSidePhysical", CerenkovMotherLogical, false, 0, lxs->OverlapTest);

	 new G4PVPlacement(0, G4ThreeVector(-((cc/ccl - 1)*csxf/2. + (ccl-1)*cslo/2. +  sir + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt/2.), 0, 0), CerenkovShieldingPlateSideLogical,
                     "CerenkovShieldingPlateSidePhysical", CerenkovMotherLogical, false, 1, lxs->OverlapTest);

	 
	 
	 new G4PVPlacement(0, G4ThreeVector(0,0,0), CerenkovBoxLogical,
                     "CerenkovBoxPhysical", CerenkovMotherLogical, false, 0, lxs->OverlapTest);
	 new G4PVPlacement(0, G4ThreeVector(0,0,0), CerenkovBoxInnerLogical,
                     "CerenkovInnerPhysical", CerenkovMotherLogical, false, 0, lxs->OverlapTest);
   

	 G4double cerx = ((cc/ccl - 1)*csxf/2. + (ccl-1)*cslo +  sir/2. + (sglt+salt+sklt+splt) + cspb + cbt + cspt);
         G4double cerz =  (ccl-1)*cszf/2. +  sir/2. + 2*(sglt+salt+sklt+splt) + cspb + cbt + cspt;

  // Compton electron spectrometer part
// 	 G4double comptonelzpos = lxs->IPMagnetZpos + lxs->TypMBFieldLength/2.0  + lxs->ComptonElBackshift; // lxs->OPPPDetZtoMagnet + 1.5*m
	 G4double comptonelzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0  + lxs->ComptonElBackshift;

  G4double shift_z = comptonelzpos - lxs->ComptonElectronBeamtoStageY - 3*lxs->ScintFrameThickness/2. + lxs->HICSDetBottomSupportZ; // lxs->OPPPTrackerECalZ;

  G4double cermxpos = 0.0;
  if (abs(lxs->IPMagFieldY + 10000.0*gauss) < 1.0e-10) {
    cermxpos = lxs->HICSCherenkovXPosMag1T;
  } else if (abs(lxs->IPMagFieldY + 20000.0*gauss) < 1.0e-10) {
    cermxpos = lxs->HICSCherenkovXPosMag2T;
  } else {
    G4String msgstr("Field ");
    msgstr += std::to_string(lxs->IPMagFieldY/gauss) + G4String("Gauss is not supported!\n");
    G4Exception("DetectorConstruction::", "ConstructScintCerenkov()", FatalException, msgstr.c_str());
  }
  
  // Cherenkov Detector at e-laser IP compton spectrum
  
  G4ThreeVector eccpos(-(cerx+cermxpos),  - lxs->ComptonElectronBeamtoStageY + lxs->CerenkovTotalBoxHeight/2., shift_z - cerz);
  new G4PVPlacement(0, eccpos, CerenkovMotherLogical,
                     "HICSElectronCerenkov", fLogicWorld, false, 1, lxs->OverlapTest);

   // ////// Position Cherenkov detector behind LANEX in Gamma Spectrometer
   
   G4double lysozpos = lxs->ComptonLysoZpos + lxs->ComptonLysoZ/2.0 + 100.0*mm;

   G4ThreeVector ccpos(lxs->ComptonLysoXpos - lxs->ComptonLysoX/2.0 + cerx, - lxs->ComptonElectronBeamtoStageY + lxs->CerenkovTotalBoxHeight/2., lysozpos + cerz + 15.*cm);

   new G4PVPlacement(0, ccpos, CerenkovMotherLogical, 
                      "GammaSpectrometerCerenkov", fLogicWorld, false, 2, lxs->OverlapTest);
  
   ccpos.setX(-ccpos.x());

   new G4PVPlacement(0, ccpos, CerenkovMotherLogical, 
                      "GammaSpectrometerCerenkov", fLogicWorld, false, 3, lxs->OverlapTest);
   
   // creating Brem region Cherenkov
   
   auto cerenkovRot = new G4RotationMatrix();
   cerenkovRot->rotateZ(-90.*deg);
   cerenkovRot->rotateY(lxs->CerenkovAngle);
   //   cerenkovRot->rotateZ(180.*deg);

   new G4PVPlacement(cerenkovRot, G4ThreeVector(lxs->ComptonElectronBeamtoStageY - lxs->CerenkovTotalBoxHeight/2., lxs->CerenkovXpos, lxs->CerenkovZpos), CerenkovMotherLogical,
                     "BremCerenkov", fLogicWorld, false,0, lxs->OverlapTest);
}
