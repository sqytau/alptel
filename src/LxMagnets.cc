//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>
#include <fstream>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4GenericTrap.hh"
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
#include "LxMagnets.hh"


TypMBMagnetAssembly::TypMBMagnetAssembly(const G4String mtypename) :
  fMagnetType(mtypename), fMagnetAssembly(0) 
{
}



G4LogicalVolume* TypMBMagnetAssembly::GetFieldVolume(const G4String lvname)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Box *solidTypMBFieldVol = new G4Box(lvname + G4String("solid"), lxs->TypMBMagnetHoleX/2.0, 
                                                   lxs->TypMBMagnetHoleY/2.0, lxs->TypMBFieldLength/2.0);
  G4LogicalVolume *logicTypMBFieldVol = new G4LogicalVolume(solidTypMBFieldVol, EnvironmentMaterial, lvname);
  return logicTypMBFieldVol;
}



void TypMBMagnetAssembly::ConstructMagnet()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* wireMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TypMBWireMaterial);
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);

//Container for the main central components
//   G4double magnetz = lxs->TypMBMagnetZ;
  G4double magnetz = lxs->TypMBMagnetCoreZ;
  G4Box *solidTypMBMagnetContainerOut = new G4Box("solidTypMBMagnetContainerOut", lxs->TypMBMagnetX/2.0, 
                                                  lxs->TypMBMagnetY/2.0, magnetz/2.0);
  G4Box *solidTypMBMagnetHole = new G4Box("solidTypMBMagnetHole", lxs->TypMBMagnetHoleX/2.0, 
                                                  lxs->TypMBMagnetHoleY/2.0, magnetz);
  G4SubtractionSolid* solidTypMBMagnetContainer1 = new G4SubtractionSolid("solidTypMBMagnetContainer1",
                                                  solidTypMBMagnetContainerOut, solidTypMBMagnetHole); 
  G4double hx1 = lxs->TypMBMagnetEndCapHoleX - (lxs->TypMBMagnetEndCapHoleOutX - lxs->TypMBMagnetEndCapHoleX);
  G4double hx2 = lxs->TypMBMagnetEndCapHoleOutX + (lxs->TypMBMagnetEndCapHoleOutX - lxs->TypMBMagnetEndCapHoleX);
  G4Trd *solidTypMBEndCapHole = new G4Trd("solidTypMBEndCapHole", hx1/2.0, hx2/2.0, lxs->TypMBMagnetHoleY/2.0, 
                                             lxs->TypMBMagnetHoleY/2.0, 1.5*lxs->TypMBMagnetEndCapZ);
  G4Transform3D trdhole(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0), 
                          G4ThreeVector(0.0, 0.0, (magnetz-lxs->TypMBMagnetEndCapZ)/2.0));
  G4SubtractionSolid* solidTypMBMagnetContainer = new G4SubtractionSolid("solidTypMBMagnetContainer",
                                                  solidTypMBMagnetContainer1, solidTypMBEndCapHole, trdhole); 
  G4LogicalVolume *logicTypMBMagnetContainer = new G4LogicalVolume(solidTypMBMagnetContainer, EnvironmentMaterial, 
                                                                     "logicTypMBMagnetContainer");

//Copper wireing
  G4Box *solidTypMBWireSide = new G4Box("solidTypMBWireSide", lxs->TypMBWireXW/2.0, lxs->TypMBWireXT/2.0, lxs->TypMBWireL/2.0);
  G4Box *solidTypMBWireFront = new G4Box("solidTypMBWireFront", lxs->TypMBWireH/2.0, lxs->TypMBWireXW/2.0, lxs->TypMBWireXT/2.0);
  G4double yl = lxs->TypMBWireYShift - 0.5*(lxs->TypMBWireXW + lxs->TypMBWireXT);
  if (yl < 0.0) {
    G4String msgstr("Error: ConstructDumpMagnet: wire geometry is not good!\n");
    G4Exception("DetectorConstruction::", "ConstructDumpMagnet()", FatalException, msgstr.c_str());
  }
  G4Box *solidTypMBWireV = new G4Box("solidTypMBWireV", lxs->TypMBWireXW/2.0, yl/2.0, lxs->TypMBWireXT/2.0);
  G4Tubs *solidTypMBWireLR = new G4Tubs("solidTypMBWireLR", 0.0, lxs->TypMBWireXT, lxs->TypMBWireXW/2.0, 0.0, 0.5*M_PI);
  G4Tubs *solidTypMBWireHR = new G4Tubs("solidTypMBWireHR", 0.0, lxs->TypMBWireXW, lxs->TypMBWireXT/2.0, 0.0, 0.5*M_PI);

  G4LogicalVolume *logicTypMBWireSide = new G4LogicalVolume(solidTypMBWireSide, wireMaterial, "logicTypMBWireSide");
  G4LogicalVolume *logicTypMBWireFront = new G4LogicalVolume(solidTypMBWireFront, wireMaterial, "logicTypMBWireFront");
  G4LogicalVolume *logicTypMBWireV = new G4LogicalVolume(solidTypMBWireV, wireMaterial, "logicTypMBWireV");
  G4LogicalVolume *logicTypMBWireLR = new G4LogicalVolume(solidTypMBWireLR, wireMaterial, "logicTypMBWireLR");
  G4LogicalVolume *logicTypMBWireHR = new G4LogicalVolume(solidTypMBWireHR, wireMaterial, "logicTypMBWireHR");

  G4ThreeVector vv;
  G4double xs[4] = {1.0, -1.0, -1.0, 1.0};
  G4double ys[4] = {1.0, 1.0, -1.0, -1.0};
  for (G4int ii = 0; ii < 4; ++ii) {
    new G4PVPlacement(0, G4ThreeVector(xs[ii]*(lxs->TypMBWireH + lxs->TypMBWireXW)/2.0, ys[ii]*lxs->TypMBWireXT/2.0, 0.0), 
                                logicTypMBWireSide, "TypMBWireSide", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);
    vv = G4ThreeVector(0.0, xs[ii]*(lxs->TypMBWireYShift + 0.5*lxs->TypMBWireXT), ys[ii]*(lxs->TypMBWireL + lxs->TypMBWireXT)/2.0);
    new G4PVPlacement(0, vv, logicTypMBWireFront, "TypMBWireFront", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);
    vv = G4ThreeVector(xs[ii]*(lxs->TypMBWireH + lxs->TypMBWireXW)/2.0, ys[ii]*(lxs->TypMBWireXW + yl)/2.0, 
                       (lxs->TypMBWireL + lxs->TypMBWireXT)/2.0);
    new G4PVPlacement(0, vv, logicTypMBWireV, "TypMBWireV", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);
    vv.setZ(-(lxs->TypMBWireL + lxs->TypMBWireXT)/2.0);
    new G4PVPlacement(0, vv, logicTypMBWireV, "TypMBWireV", logicTypMBMagnetContainer, false, ii+4, lxs->OverlapTest);
    
    // Joins for side copper bars 
    G4RotationMatrix *rotr = new G4RotationMatrix(0.5*M_PI, ys[ii]*0.5*M_PI, 0.0);
    vv = G4ThreeVector(xs[ii]*(lxs->TypMBWireH + lxs->TypMBWireXW)/2.0, -lxs->TypMBWireXT, ys[ii]*lxs->TypMBWireL/2.0);
    new G4PVPlacement(rotr, vv, logicTypMBWireLR, "TypMBWireLR", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);

    rotr = new G4RotationMatrix(-0.5*M_PI, ys[ii]*0.5*M_PI, 0.0);
    vv = G4ThreeVector(xs[ii]*(lxs->TypMBWireH + lxs->TypMBWireXW)/2.0, lxs->TypMBWireXT, ys[ii]*lxs->TypMBWireL/2.0);
    new G4PVPlacement(rotr, vv, logicTypMBWireLR, "TypMBWireLR", logicTypMBMagnetContainer, false, ii+4, lxs->OverlapTest);

    // Joins for transverse copper bars 
    rotr = new G4RotationMatrix((1.0-xs[ii])*0.25*M_PI, 0.0, 0.0);
    vv = G4ThreeVector(xs[ii]*lxs->TypMBWireH/2.0, lxs->TypMBWireYShift-lxs->TypMBWireXT/2.0, 
                       ys[ii]*(lxs->TypMBWireL+lxs->TypMBWireXT)/2.0);
    new G4PVPlacement(rotr, vv, logicTypMBWireHR, "TypMBWireHR", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);

    rotr = new G4RotationMatrix((xs[ii] - 3.0)*0.25*M_PI, 0.0, 0.0);
    vv = G4ThreeVector(xs[ii]*lxs->TypMBWireH/2.0, -(lxs->TypMBWireYShift-lxs->TypMBWireXT/2.0), 
                       ys[ii]*(lxs->TypMBWireL+lxs->TypMBWireXT)/2.0);
    new G4PVPlacement(rotr, vv, logicTypMBWireHR, "TypMBWireHR", logicTypMBMagnetContainer, false, ii+4, lxs->OverlapTest);
  }
  
  // Part of the magnet core
  G4double wedgez = 0.5 * lxs->TypMBCoreL * lxs->TypMBWedgeLFraction;
  G4double wedgey = lxs->TypMBWireYShift;
  G4double wedgel = sqrt(wedgey*wedgey + wedgez*wedgez);
  G4double wedgeangle = atan2(wedgey, wedgez);
  G4double magboxz = lxs->TypMBCoreL/2.0;
  G4Box *solidMagCoreBox = new G4Box("solidMagCoreBox", lxs->TypMBCoreW/2.0, lxs->TypMBCoreT/2.0, magboxz/2.0);
  G4Box *solidMagCoreCut = new G4Box("solidMagCoreCut", 0.6*lxs->TypMBCoreW, wedgey/2.0,  0.6*wedgel);
  G4Transform3D transform(G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), -wedgeangle), G4ThreeVector(0.0, 
      -0.5*(lxs->TypMBCoreT-wedgez*tan(wedgeangle)+wedgey*cos(wedgeangle)), 0.5*(magboxz-wedgez+wedgey*sin(wedgeangle))));
  G4SubtractionSolid* solidMagCoreWedge = new G4SubtractionSolid("solidMagCoreWedge", 
                                                           solidMagCoreBox, solidMagCoreCut, transform);                       
  G4LogicalVolume *logicMagCoreBox = new G4LogicalVolume(solidMagCoreWedge, magnetMaterial, "logicMagCoreBox");
  
  //small bars up and down the wire
  G4double wiretopy = lxs->TypMBMagnetY/2.0 - lxs->TypMBWireYShift 
                      - 0.5*(lxs->TypMBWireXT + lxs->TypMBWireXW) - lxs->TypMBWireGap;
  G4double wiretopz = 0.5*(lxs->TypMBMagnetCoreZ - lxs->TypMBCoreL) - lxs->TypMBMagnetEndCapZ;
  G4Box *solidTypMBWireTopBar = new G4Box("solidTypMBWireTopBar", lxs->TypMBMagnetEndCapW/2.0, wiretopy/2.0,  wiretopz/2.0);
  G4LogicalVolume *logicTypMBWireTopBar = new G4LogicalVolume(solidTypMBWireTopBar, magnetMaterial, "logicTypMBWireTopBar");
  
  for (G4int ii = 0; ii < 4; ++ii) {
    G4RotationMatrix *rotr = new G4RotationMatrix((ys[ii]-xs[ii])*M_PI/2.0, (xs[ii]<0?1.0:0.0)*M_PI, 0.0);
    // place core parts with wedge cut 
    new G4PVPlacement(rotr, G4ThreeVector(0.0, ys[ii]*(lxs->TypMBCoreT+lxs->TypMBMagnetHoleY)/2.0, xs[ii]*magboxz/2.0), 
                     logicMagCoreBox, "MagCoreBox", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);

    // place small bars up and down the wire 
    new G4PVPlacement(rotr, G4ThreeVector(0.0, ys[ii]*(lxs->TypMBMagnetY-wiretopy)/2.0, xs[ii]*(lxs->TypMBCoreL+wiretopz)/2.0), 
                     logicTypMBWireTopBar, "TypMBWireTopBar", logicTypMBMagnetContainer, false, ii, lxs->OverlapTest);
  }

  //Core endcaps  
  G4Box *solidMagEndCap1 = new G4Box("solidMagEndCap1", lxs->TypMBMagnetEndCapW/2.0, lxs->TypMBMagnetY/2.0, 
                                        lxs->TypMBMagnetEndCapZ/2.0);
//   G4Box *solidTypMBEndCapHole = new G4Box("solidTypMBEndCap1Hole", lxs->TypMBMagnetEndCapHoleX/2.0, 
//                                                   lxs->TypMBMagnetHoleY/2.0, lxs->TypMBMagnetEndCapZ);
  G4SubtractionSolid* solidTypMBCoreEndCap = new G4SubtractionSolid("solidTypMBCoreEndCap", solidMagEndCap1, 
                                                                    solidTypMBEndCapHole);                       
  G4LogicalVolume *logicTypMBCoreEndCap = new G4LogicalVolume(solidTypMBCoreEndCap, magnetMaterial, "logicTypMBCoreEndCap");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (lxs->TypMBMagnetCoreZ - lxs->TypMBMagnetEndCapZ)/2.0), 
                     logicTypMBCoreEndCap, "TypMBEndCap", logicTypMBMagnetContainer, false, 0, lxs->OverlapTest);
  
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),
                     G4ThreeVector(0.0, 0.0, -(lxs->TypMBMagnetCoreZ - lxs->TypMBMagnetEndCapZ)/2.0), 
                     logicTypMBCoreEndCap, "TypMBEndCap", logicTypMBMagnetContainer, false, 1, lxs->OverlapTest);
  
  //Side planes attached to the core structure
  G4double sidex = 0.5*(lxs->TypMBMagnetX - lxs->TypMBCoreW);
  G4Box *solidTypMBCoreSide = new G4Box("solidTypMBCoreSide", sidex/2.0, lxs->TypMBMagnetY/2.0, lxs->TypMBCoreL/2.0);
  G4LogicalVolume *logicTypMBCoreSide = new G4LogicalVolume(solidTypMBCoreSide, magnetMaterial, "logicTypMBCoreSide");
  new G4PVPlacement(0, G4ThreeVector(0.5*(lxs->TypMBMagnetX - sidex), 0.0, 0.0), 
                     logicTypMBCoreSide, "TypMBCoreSide", logicTypMBMagnetContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(-0.5*(lxs->TypMBMagnetX - sidex), 0.0, 0.0), 
                     logicTypMBCoreSide, "TypMBCoreSide", logicTypMBMagnetContainer, false, 1, lxs->OverlapTest);
  
//   //TypMB magnet endcaps  
//   G4double ecmz = 0.5*(lxs->TypMBMagnetZ - lxs->TypMBMagnetCoreZ);
//   G4Box *solidTypMBEndCap1 = new G4Box("solidTypMBEndCap1", lxs->TypMBMagnetEndCapW/2.0, lxs->TypMBMagnetY/2.0, ecmz/2.0);
//   G4Box *solidTypMBEndCap1Hole = new G4Box("solidTypMBEndCap1Hole", lxs->TypMBMagnetEndCapHoleX/2.0, 
//                                                   lxs->TypMBMagnetHoleY/2.0, lxs->TypMBMagnetEndCapZ);
//   G4SubtractionSolid* solidTypMBEndCap = new G4SubtractionSolid("solidTypMBEndCap", solidTypMBEndCap1, solidTypMBEndCap1Hole);                       
//   G4LogicalVolume *logicTypMBEndCap = new G4LogicalVolume(solidTypMBEndCap, magnetMaterial, "logicTypMBEndCap");
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.5*(lxs->TypMBMagnetZ - ecmz)), 
//                      logicTypMBEndCap, "TypMBEndCap1", logicTypMBMagnetContainer, false, 0, lxs->OverlapTest);
// 
//   // This is another not documented endcap part of funny shape, some stuff has to be estimated. 
//   G4double ecmx2 = lxs->TypMBMagnetHoleX + 2.0*lxs->TypMBWireGap;
//   G4double ecmy2 = lxs->TypMBMagnetHoleY + 2.0*lxs->TypMBWireGap;
//   G4Box *solidTypMBEndCap22 = new G4Box("solidTypMBEndCap22", ecmx2/2.0, ecmy2/2.0, ecmz/2.0);
//   G4SubtractionSolid* solidTypMBEndCap2 = new G4SubtractionSolid("solidTypMBEndCap2", solidTypMBEndCap22, solidTypMBEndCap1Hole);                       
//   G4double ecmx23 = lxs->TypMBMagnetHoleX + 2.0*lxs->TypMBMagnetHoleY;
//   G4double ecmy23 = lxs->TypMBMagnetHoleY + 2.0*lxs->TypMBMagnetHoleY;
//   G4Box *solidTypMBEndCap23 = new G4Box("solidTypMBEndCap23", ecmx23/2.0, ecmy23/2.0, lxs->TypMBWireGap/2.0);
//   G4SubtractionSolid* solidTypMBEndCap3 = new G4SubtractionSolid("solidTypMBEndCap3", solidTypMBEndCap23, solidTypMBEndCap22);                       
// 
//   G4LogicalVolume *logicTypMBEndCap2 = new G4LogicalVolume(solidTypMBEndCap2, magnetMaterial, "logicTypMBEndCap2");
//   G4LogicalVolume *logicTypMBEndCap3 = new G4LogicalVolume(solidTypMBEndCap3, magnetMaterial, "logicTypMBEndCap3");
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -0.5*(lxs->TypMBMagnetZ - ecmz)), 
//                      logicTypMBEndCap2, "TypMBEndCap2", logicTypMBMagnetContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -0.5*(lxs->TypMBMagnetZ - lxs->TypMBWireGap)), 
//                      logicTypMBEndCap3, "TypMBEndCap3", logicTypMBMagnetContainer, false, 0, lxs->OverlapTest);

  //Collect everything to the assembly
  fMagnetAssembly = new G4AssemblyVolume();
  G4ThreeVector pipetrans(0.0, 0.0, 0.0);
  fMagnetAssembly->AddPlacedVolume(logicTypMBMagnetContainer, pipetrans, 0);
}



void TypMBMagnetAssembly::CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                          const G4bool rotate)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  // Support structures
  G4double suplength = lxs->TypMBSupportH + pos.y();
  G4double magsize = lxs->TypMBMagnetY;
  G4double magwidth = lxs->TypMBMagnetX;
  if(rotate) {
    suplength = lxs->TypMBSupportV + pos.y();
    magsize = lxs->TypMBMagnetX;
    magwidth = lxs->TypMBMagnetY;
  }
  const G4int nsplane = 14;
  G4double zsf[nsplane] = {0.0, 0.02, 0.05, 0.45, 0.46, 0.47, 0.475, 0.525, 0.53, 0.54, 0.55, 0.97, 0.98, 1.0};
  G4double zsplane[nsplane];
  std::transform (zsf, zsf+nsplane, zsplane, std::bind(std::multiplies<double>(), std::placeholders::_1, suplength));
  G4double rr = lxs->TypMBSupportR;
  G4double rsouter[nsplane] = {2.5*rr, 2.5*rr, rr, rr, 0.5*rr, 0.5*rr, rr, rr, 0.5*rr, 0.5*rr, rr, rr, 1.5*rr, 1.5*rr};
  G4double rsinner[nsplane] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  G4String ssname("solidSupport");
  ssname += mname;
  G4Polycone *solidTypMBSupport = new G4Polycone(ssname, 0.0, 2.0*M_PI, nsplane, zsplane, rsinner, rsouter);
  G4String slname("logicSupport");
  slname += mname;
  G4LogicalVolume *logicTypMBSupport = new G4LogicalVolume(solidTypMBSupport, magnetMaterial, slname);

  G4String spname("Support");
  spname += mname;
  G4double dz = 0.5*lxs->TypMBCoreL - 1.05*rsouter[nsplane-1];
  G4double dx = 0.5*magwidth - 1.05*rsouter[nsplane-1];
  G4ThreeVector suppos(pos.x(), pos.y() - suplength - magsize/2.0, pos.z() + dz);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicTypMBSupport, 
                    spname, lWorld, false, 0, lxs->OverlapTest);
  suppos.setX(pos.x() - dx);
  suppos.setZ(pos.z() - dz);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicTypMBSupport, 
                    spname, lWorld, false, 1, lxs->OverlapTest);
  suppos.setX(pos.x() + dx);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicTypMBSupport, 
                    spname, lWorld, false, 2, lxs->OverlapTest);

}

/////////////////////////////////////////////////////////////////////////////////////
// TypMBMagnetSimple


TypMBMagnetSimple::TypMBMagnetSimple(const G4String mtypename) :
  TypMBMagnetAssembly(mtypename) 
{
}


void TypMBMagnetSimple::ConstructMagnet()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4double xs[4] = {1.0, -1.0, -1.0, 1.0};
  G4double ys[4] = {1.0, 1.0, -1.0, -1.0};

//Container for the main central components
//   G4double magnetz = lxs->TypMBMagnetZ;
  G4double magnetz = lxs->TypMBMagnetCoreZ;
  G4Box *solidTypMBMagSimpContainerOut = new G4Box("solidTypMBMagSimpContainerOut", lxs->TypMBMagnetX/2.0,
                                                  lxs->TypMBMagnetY/2.0, magnetz/2.0);
  G4Box *solidTypMBMagSimpHole = new G4Box("solidTypMBMagSimpHole", lxs->TypMBMagnetHoleX/2.0,
                                                  lxs->TypMBMagnetHoleY/2.0, magnetz);
  G4SubtractionSolid* solidTypMBMagSimpContainer1 = new G4SubtractionSolid("solidTypMBMagSimpContainer1",
                                                  solidTypMBMagSimpContainerOut, solidTypMBMagSimpHole);
  G4double hx1 = lxs->TypMBMagnetEndCapHoleX - (lxs->TypMBMagnetEndCapHoleOutX - lxs->TypMBMagnetEndCapHoleX);
  G4double hx2 = lxs->TypMBMagnetEndCapHoleOutX + (lxs->TypMBMagnetEndCapHoleOutX - lxs->TypMBMagnetEndCapHoleX);
  G4Trd *solidTypMBSimpEndCapHole = new G4Trd("solidTypMBSimpEndCapHole", hx1/2.0, hx2/2.0, lxs->TypMBMagnetHoleY/2.0,
                                             lxs->TypMBMagnetHoleY/2.0, 1.5*lxs->TypMBMagnetEndCapZ);
  G4Transform3D trdhole(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
                          G4ThreeVector(0.0, 0.0, (magnetz-lxs->TypMBMagnetEndCapZ)/2.0));
  G4SubtractionSolid* solidTypMBMagSimpContainer = new G4SubtractionSolid("solidTypMBMagSimpContainer",
                                                  solidTypMBMagSimpContainer1, solidTypMBSimpEndCapHole, trdhole);
  G4LogicalVolume *logicTypMBMagSimpContainer = new G4LogicalVolume(solidTypMBMagSimpContainer, EnvironmentMaterial,
                                                                     "logicTypMBMagSimpContainer");

  // Part of the magnet core
  G4Box *solidMagSimpCoreBox1 = new G4Box("solidMagSimpCoreBox1", lxs->TypMBMagnetX/2.0, lxs->TypMBMagnetY/2.0,
                                         lxs->TypMBCoreL/2.0);
  G4SubtractionSolid* solidMagSimpCoreBox = new G4SubtractionSolid("solidMagSimpCoreBox",
                                                  solidMagSimpCoreBox1, solidTypMBMagSimpHole);
  G4LogicalVolume *logicMagSimpCoreBox = new G4LogicalVolume(solidMagSimpCoreBox, magnetMaterial, "logicMagSimpCoreBox");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                     logicMagSimpCoreBox, "MagSimpCoreBox", logicTypMBMagSimpContainer, false, 0, lxs->OverlapTest);

  //small bars up and down the wire
  G4double wiretopy = lxs->TypMBMagnetY/2.0 - lxs->TypMBWireYShift
                      - 0.5*(lxs->TypMBWireXT + lxs->TypMBWireXW) - lxs->TypMBWireGap;
  G4double wiretopz = 0.5*(lxs->TypMBMagnetCoreZ - lxs->TypMBCoreL) - lxs->TypMBMagnetEndCapZ;
  G4Box *solidTypMBSimpWireTopBar = new G4Box("solidTypMBSimpWireTopBar", lxs->TypMBMagnetEndCapW/2.0,
                                               wiretopy/2.0,  wiretopz/2.0);
  G4LogicalVolume *logicTypMBSimpWireTopBar = new G4LogicalVolume(solidTypMBSimpWireTopBar, magnetMaterial,
                                                                  "logicTypMBSimpWireTopBar");
  for (G4int ii = 0; ii < 4; ++ii) {
    G4RotationMatrix *rotr = new G4RotationMatrix((ys[ii]-xs[ii])*M_PI/2.0, (xs[ii]<0?1.0:0.0)*M_PI, 0.0);
    new G4PVPlacement(rotr, G4ThreeVector(0.0, ys[ii]*(lxs->TypMBMagnetY-wiretopy)/2.0, xs[ii]*(lxs->TypMBCoreL+wiretopz)/2.0),
                     logicTypMBSimpWireTopBar, "TypMBSimpWireTopBar", logicTypMBMagSimpContainer, false, ii, lxs->OverlapTest);
  }

  //Core endcaps
  G4Box *solidMagSimpEndCap1 = new G4Box("solidMagSimpEndCap1", lxs->TypMBMagnetEndCapW/2.0, lxs->TypMBMagnetY/2.0,
                                        lxs->TypMBMagnetEndCapZ/2.0);
  G4SubtractionSolid* solidTypMBSimpCoreEndCap = new G4SubtractionSolid("solidTypMBSimpCoreEndCap", solidMagSimpEndCap1,
                                                                        solidTypMBSimpEndCapHole);
  G4LogicalVolume *logicTypMBSimpCoreEndCap = new G4LogicalVolume(solidTypMBSimpCoreEndCap, magnetMaterial,
                                                                  "logicTypMBSimpCoreEndCap");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (lxs->TypMBMagnetCoreZ - lxs->TypMBMagnetEndCapZ)/2.0),
                     logicTypMBSimpCoreEndCap, "TypMBSimpEndCap", logicTypMBMagSimpContainer, false, 0, lxs->OverlapTest);

  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI),
                     G4ThreeVector(0.0, 0.0, -(lxs->TypMBMagnetCoreZ - lxs->TypMBMagnetEndCapZ)/2.0),
                     logicTypMBSimpCoreEndCap, "TypMBSimpEndCap", logicTypMBMagSimpContainer, false, 1, lxs->OverlapTest);

  // Support structures
  const G4int nsplane = 6;
  G4double zsf[nsplane] = {0.0, 0.01, 0.03, 0.97, 0.98, 1.0};
  G4double zsplane[nsplane];
  std::transform (zsf, zsf+nsplane, zsplane, std::bind(std::multiplies<double>(), std::placeholders::_1, lxs->TypMBSupportH));
  G4double rr = lxs->TypMBSupportR;
  G4double rsouter[nsplane] = {2.5*rr, 2.5*rr, rr, rr, 0.7*rr, 0.7*rr};
  G4double rsinner[nsplane] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4Polycone *solidTypMBSimpleSupport = new G4Polycone("solidTypMBSimpleSupport", 0.0, 2.0*M_PI, nsplane,
                                                       zsplane, rsinner, rsouter);
  G4LogicalVolume *logicTypMBSimpleSupport = new G4LogicalVolume(solidTypMBSimpleSupport, magnetMaterial,
                                                                 "logicTypMBSimpleSupport");

  //Collect everything to the assembly
  fMagnetAssembly = new G4AssemblyVolume();
  G4ThreeVector pipetrans(0.0, 0.0, 0.0);
  fMagnetAssembly->AddPlacedVolume(logicTypMBMagSimpContainer, pipetrans, 0);

  G4RotationMatrix *rot1 = new G4RotationMatrix(0.0, 0.5*M_PI, 0.0);
  pipetrans = G4ThreeVector(lxs->TypMBSupportPosDx/2.0, -lxs->TypMBSupportH - lxs->TypMBMagnetY/2.0,
                                                               lxs->TypMBSupportPosDz/2.0);
  fMagnetAssembly->AddPlacedVolume(logicTypMBSimpleSupport, pipetrans, rot1);
  pipetrans = G4ThreeVector(lxs->TypMBSupportPosDx/2.0, -lxs->TypMBSupportH - lxs->TypMBMagnetY/2.0,
                                                              -lxs->TypMBSupportPosDz/2.0);
  fMagnetAssembly->AddPlacedVolume(logicTypMBSimpleSupport, pipetrans, rot1);
  pipetrans = G4ThreeVector(-lxs->TypMBSupportPosDx/2.0, -lxs->TypMBSupportH - lxs->TypMBMagnetY/2.0, 0.0);
  fMagnetAssembly->AddPlacedVolume(logicTypMBSimpleSupport, pipetrans, rot1);

}



/////////////////////////////////////////////////////////////////////////////////////
// PDSMagnetAssembly


PDSMagnetAssembly::PDSMagnetAssembly(const G4String mtypename) :
  fMagnetType(mtypename) 
{
}


void PDSMagnetAssembly::ConstructMagnet()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* wireMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TypMBWireMaterial);
  
  G4double magnetz = lxs->GMagnetZ - 2.0*lxs->GMagCoilThickness;
  G4double dxpos = lxs->GMagnetCutX - lxs->GMagnetX;
  G4Box *solidMagSimpContainerOut = new G4Box("solidMagSimpContainerOut", lxs->GMagnetX/2.0, 
                                                  lxs->GMagnetY/2.0,   magnetz/2.0);
  G4Box *solidMagSimpHole = new G4Box("solidMagSimpHole", lxs->GMagnetX/2.0, 
                                                  lxs->GMagnetCutY/2.0, magnetz);
  G4Transform3D trcut(G4RotationMatrix(), G4ThreeVector(dxpos, 0.0, 0.0));
  G4SubtractionSolid* solidGMagnetYoke = new G4SubtractionSolid("solidGMagnetYoke", solidMagSimpContainerOut, solidMagSimpHole, trcut); 
  G4LogicalVolume *logicGMagnetYoke = new G4LogicalVolume(solidGMagnetYoke, magnetMaterial, "logicGMagnetYoke");

//Core 
  G4double corex = lxs->GMagnetCoreX;
  G4double corey = 0.5*(lxs->GMagnetCutY - lxs->GMagneteffY);
  G4Box *solidGMagnetCore = new G4Box("solidGMagnetCore", corex/2.0, corey/2.0,  magnetz/2.0);
  G4LogicalVolume *logicGMagnetCore = new G4LogicalVolume(solidGMagnetCore, magnetMaterial,
                                                                 "logicGMagnetCore");
//Coil
  G4double coilz = lxs->GMagnetZ; 
  G4double coily = lxs->GMagCoilH;  
  G4double coilx = corex + 2.0*lxs->GMagCoilThickness;                                                               
  G4Box *solidGMagnetCoilOut = new G4Box("solidGMagnetCoilOut",  coilx/2.0, coily/2.0,   coilz/2.0);
  G4Box *solidGMagnetCoilHole = new G4Box("solidGMagnetCoilHole", corex/2.0, corey/2.0,  magnetz/2);
//  G4Transform3D trcut(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid* solidGMagnetCoil = new G4SubtractionSolid("solidGMagnetCoil", solidGMagnetCoilOut, solidGMagnetCoilHole); 
  
  G4LogicalVolume *logicGMagnetCoil = new G4LogicalVolume(solidGMagnetCoil, wireMaterial,
                                                                 "logicGMagnetCoil");
                                                                 
    //Collect everything to the assembly
  fMagnetAssembly = new G4AssemblyVolume();
  
  G4ThreeVector vtr(0.0, 0.0, 0.0);
  fMagnetAssembly->AddPlacedVolume(logicGMagnetYoke, vtr, 0);
  G4double coreXpos = 0.5*(corex - lxs->GMagnetX);
  vtr.setX(coreXpos);
  vtr.setY((lxs->GMagnetCutY-corey)/2.0);
  fMagnetAssembly->AddPlacedVolume(logicGMagnetCore, vtr, 0);
  vtr.setY(-vtr.y());
  fMagnetAssembly->AddPlacedVolume(logicGMagnetCore, vtr, 0);
  vtr.setY((lxs->GMagnetCutY-coily)/2.0);
  fMagnetAssembly->AddPlacedVolume(logicGMagnetCoil, vtr, 0);
  vtr.setY(-vtr.y());
  fMagnetAssembly->AddPlacedVolume(logicGMagnetCoil, vtr, 0);

}



G4LogicalVolume* PDSMagnetAssembly::GetFieldVolume(const G4String lvname)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Box *solidGMFieldVol = new G4Box(lvname + G4String("solid"), lxs->GMFieldX/2.0, 
                                                   lxs->GMagneteffY/2.0, lxs->GMFieldLength/2.0);
  G4LogicalVolume *logicGMFieldVol = new G4LogicalVolume(solidGMFieldVol, EnvironmentMaterial, lvname);
  return logicGMFieldVol;
}



void PDSMagnetAssembly::CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname, 
                                 G4bool rotate) 
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  // Support structures
  G4double magnetz = lxs->GMagnetZ - 2.0*lxs->GMagCoilThickness;
  G4double suplength = -lxs->FloorSurfaceYpos - lxs->GMagnetY/2.0 + pos.y();
  if(rotate) {
    G4String msgstr("Error: PDSMagnetAssembly::CostructSupport: rotation of this magnet is not supported!\n");
    G4Exception("PDSMagnetAssembly::", "CostructSupport()", FatalException, msgstr.c_str());
  }
  const G4int nsplane = 14;
  G4double zsf[nsplane] = {0.0, 0.02, 0.05, 0.45, 0.46, 0.47, 0.475, 0.525, 0.53, 0.54, 0.55, 0.97, 0.98, 1.0};
  G4double zsplane[nsplane];
  std::transform (zsf, zsf+nsplane, zsplane, std::bind(std::multiplies<double>(), std::placeholders::_1, suplength));
  G4double rr = lxs->TypMBSupportR;
  G4double rsouter[nsplane] = {2.5*rr, 2.5*rr, rr, rr, 0.5*rr, 0.5*rr, rr, rr, 0.5*rr, 0.5*rr, rr, rr, 1.5*rr, 1.5*rr};
  G4double rsinner[nsplane] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  G4String ssname("solidSupport");
  ssname += mname;
  G4Polycone *solidPDSMagSupport = new G4Polycone(ssname, 0.0, 2.0*M_PI, nsplane, zsplane, rsinner, rsouter);
  G4String slname("logicSupport");
  slname += mname;
  G4LogicalVolume *logicPDSMagSupport = new G4LogicalVolume(solidPDSMagSupport, magnetMaterial, slname);

  G4String spname("Support");
  spname += mname;
  G4double dz = 0.5*magnetz - 1.05*rsouter[nsplane-1];
  G4double dx = 0.5*lxs->GMagnetX - 1.05*rsouter[nsplane-1];
  G4ThreeVector suppos(pos.x(), pos.y() - suplength - lxs->GMagnetY/2.0, pos.z() + dz);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicPDSMagSupport, 
                    spname, lWorld, false, 0, lxs->OverlapTest);
  suppos.setX(pos.x() - dx);
  suppos.setZ(pos.z() - dz);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicPDSMagSupport, 
                    spname, lWorld, false, 1, lxs->OverlapTest);
  suppos.setX(pos.x() + dx);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), 0.5*M_PI), suppos, logicPDSMagSupport, 
                    spname, lWorld, false, 2, lxs->OverlapTest);
}



/////////////////////////////////////////////////////////////////////////////////////
// FlashMagnetAssembly


FlashMagnetAssembly::FlashMagnetAssembly(const G4String mtypename) :
  fMagnetType(mtypename)
{}



G4LogicalVolume* FlashMagnetAssembly::GetFieldVolume(const G4String lvname)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* EnvironmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Box *solidFlashMFieldVol = new G4Box(lvname + G4String("solid"), lxs->FlashMFieldX/2.0,
                                                  lxs->FlashMagneteffY/2.0, lxs->FlashMFieldLength/2.0);
  G4LogicalVolume *logicFlashMFieldVol = new G4LogicalVolume(solidFlashMFieldVol, EnvironmentMaterial, lvname);
  return logicFlashMFieldVol;
}



void FlashMagnetAssembly::ConstructMagnet()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* magMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* wireMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TypMBWireMaterial);

  G4double hmagy = lxs->FlashMagnetY/2.0;
  G4Box *solidFlashMagnetHalf1 = new G4Box("solidFlashMagnetHalf1", lxs->FlashMagnetX/2.0, hmagy/2.0,
                                            lxs->FlashMagnetZ/2.0);
  G4Box *solidFlashMMiddleCut = new G4Box("solidFlashMMiddleCut", lxs->FlashMMiddleCutX/2.0,
                                                  lxs->FlashMMiddleCutY/2.0, lxs->FlashMagnetZ);
  G4double acutx = 2.0 * lxs->FlashMagnetCutD * sqrt(2.0);
  G4double acuty = 2.0 * lxs->FlashMagnetCutD/sqrt(2.0);
  G4Box *solidFlashMACut = new G4Box("solidFlashMACut", acutx/2.0, acuty/2.0, lxs->FlashMagnetZ);

  G4Transform3D trmcut(G4RotationMatrix(), G4ThreeVector(0.0, -hmagy/2.0, 0.0));
  G4SubtractionSolid* solidFlashMagnetHalf2 = new G4SubtractionSolid("solidFlashMagnetHalf2",
                                                  solidFlashMagnetHalf1, solidFlashMMiddleCut, trmcut);

  G4Transform3D tracut1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), -M_PI/4.0),
                         G4ThreeVector(lxs->FlashMagnetX/2.0, hmagy/2.0, 0.0));
  G4SubtractionSolid* solidFlashMagnetHalf3 = new G4SubtractionSolid("solidFlashMagnetHalf3",
                                                  solidFlashMagnetHalf2, solidFlashMACut, tracut1);
  G4Transform3D tracut2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI/4.0),
                         G4ThreeVector(-lxs->FlashMagnetX/2.0, hmagy/2.0, 0.0));
  G4SubtractionSolid* solidFlashMagnetHalf = new G4SubtractionSolid("solidFlashMagnetHalf",
                                                  solidFlashMagnetHalf3, solidFlashMACut, tracut2);
  G4LogicalVolume *logicFlashMagnetHalf = new G4LogicalVolume(solidFlashMagnetHalf, magMaterial, "logicFlashMagnetHalf");
// Core
  G4double corey = 0.5*(lxs->FlashMMiddleCutY - lxs->FlashMFieldGapY);
  G4Box *solidFlashMagnetCore1 = new G4Box("solidFlashMagnetCore1", lxs->FlashMagnetCoreX/2.0,
                                            corey/2.0, lxs->FlashMagnetCoreZ/2.0);
  G4double corecutd = 20.0*mm;
  G4double acorecutx = 2.0 * corecutd * sqrt(2.0);
  G4double acorecuty = 2.0 * corecutd/sqrt(2.0);
  G4Box *solidFlashMagnetCoreCut = new G4Box("solidFlashMagnetCoreCut", acorecutx/2.0,
                                                  acorecuty/2.0, lxs->FlashMagnetCoreZ);
  G4Transform3D tracorecut1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI/4.0),
                         G4ThreeVector(lxs->FlashMagnetCoreX/2.0, -corey/2.0, 0.0));
  G4SubtractionSolid* solidFlashMagnetCore2 = new G4SubtractionSolid("solidFlashMagnetCore2",
                                     solidFlashMagnetCore1, solidFlashMagnetCoreCut, tracorecut1);
  G4Transform3D tracorecut2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), -M_PI/4.0),
                         G4ThreeVector(-lxs->FlashMagnetCoreX/2.0, -corey/2.0, 0.0));
  G4SubtractionSolid* solidFlashMagnetCore = new G4SubtractionSolid("solidFlashMagnetCore",
                                     solidFlashMagnetCore2, solidFlashMagnetCoreCut, tracorecut2);
  G4LogicalVolume *logicFlashMagnetCore = new G4LogicalVolume(solidFlashMagnetCore, magMaterial, "logicFlashMagnetCore");

// Coil
  G4double lcoilz = lxs->FlashMagnetCoreZ + 2.0*lxs->FlashMagnetCoilGapH;
  G4Box *solidFlashMagnetCoilZ = new G4Box("solidFlashMagnetCoilZ", lxs->FlashMagnetCoilH/2.0,
                                            lxs->FlashMagnetCoilD/2.0, lcoilz/2.0);
  G4LogicalVolume *logicFlashMagnetCoilZ = new G4LogicalVolume(solidFlashMagnetCoilZ,
                                                               wireMaterial, "logicFlashMagnetCoilZ");
  G4double hcoilx = lxs->FlashMagnetCoreX + 2.0*lxs->FlashMagnetCoilGapH;
  G4Box *solidFlashMagnetCoilX = new G4Box("solidFlashMagnetCoilX", hcoilx/2.0,
                                            lxs->FlashMagnetCoilD/2.0, lxs->FlashMagnetCoilH/2.0);
  G4LogicalVolume *logicFlashMagnetCoilX = new G4LogicalVolume(solidFlashMagnetCoilX,
                                                               wireMaterial, "logicFlashMagnetCoilX");
  G4Tubs *solidFlashMCoilR = new G4Tubs("solidFlashMCoilR", 0.0, lxs->FlashMagnetCoilH,
                                         lxs->FlashMagnetCoilD/2.0, 0.0, 0.5*M_PI);
  G4LogicalVolume *logicFlashMCoilR = new G4LogicalVolume(solidFlashMCoilR,
                                                               wireMaterial, "logicFlashMCoilR");
//Collect everything to the assemblies
  G4AssemblyVolume  *magnetHalfAssembly = new G4AssemblyVolume();
  G4ThreeVector trmh(0.0, hmagy/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetHalf, trmh, 0);
  G4ThreeVector trmcoreh(0.0, (lxs->FlashMMiddleCutY - corey)/2.0, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetCore, trmcoreh, 0);

  G4double coillx = 0.5*(lxs->FlashMagnetCoilH + lxs->FlashMagnetCoreX) + lxs->FlashMagnetCoilGapH;
  G4double coilly = 0.5*(lxs->FlashMMiddleCutY - lxs->FlashMagnetCoilD);
  G4ThreeVector trmcoill(coillx, coilly, 0.0);
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetCoilZ, trmcoill, 0);
  trmcoill.setX(-coillx);
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetCoilZ, trmcoill, 0);
  G4double coilhz = 0.5*(lxs->FlashMagnetZ+lxs->FlashMagnetCoilH)+lxs->FlashMagnetCoilGapH;
  G4ThreeVector trmcoilh(0.0, coilly, coilhz);
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetCoilX, trmcoilh, 0);
  trmcoilh.setZ(-trmcoilh.z());
  magnetHalfAssembly->AddPlacedVolume(logicFlashMagnetCoilX, trmcoilh, 0);

  G4double coilrx = coillx - lxs->FlashMagnetCoilH/2.0;
  G4double coilrz = coilhz - lxs->FlashMagnetCoilH/2.0;
  for (int ii = 0; ii < 4; ++ii) {
    G4ThreeVector trmcoilr(pow(-1.0, ii&1)*coilrx, coilly, pow(-1.0, (ii>>1))*coilrz);
    G4RotationMatrix *crrot = new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0);
    crrot->rotateY(pow(-1.0,ii)*((ii&1)+((ii>>1)&1))*M_PI/2.0);
    magnetHalfAssembly->AddPlacedVolume(logicFlashMCoilR, trmcoilr, crrot);
  }

  fMagnetAssembly = new G4AssemblyVolume();
  G4ThreeVector trmhassmbly(0.0, 0.0, 0.0);
//   fMagnetAssembly->AddPlacedAssembly(magnetHalfAssembly, trmhassmbly, 0);
//   fMagnetAssembly->AddPlacedAssembly(magnetHalfAssembly, trmhassmbly,
//                                      new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));
  LxAux::AddAssmblyVolumes(fMagnetAssembly, magnetHalfAssembly, trmhassmbly, 0);
  LxAux::AddAssmblyVolumes(fMagnetAssembly, magnetHalfAssembly, trmhassmbly, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI));

}



void FlashMagnetAssembly::CostructSupport(G4LogicalVolume* lWorld, const G4ThreeVector& pos, const G4String& mname,
                                 G4bool rotate)
{
// Support
  LXSetUp *lxs = LXSetUp::Instance();
  G4double ypestal = 1.5*m;
  G4double ylevel = lxs->FlashMagnetY/2.0 - pos.y();
  G4double tblx = lxs->FlashMagnetX - 2.0*lxs->FlashMagnetCutD;
  G4double tblz = lxs->FlashMagnetZ - 2.0*lxs->FlashMagnetCutD;
  if (rotate) {
    ylevel = lxs->FlashMagnetX/2.0 - pos.y();
    tblx = lxs->FlashMagnetY - 2.0*lxs->FlashMagnetCutD;
  }
  G4double tblhight = -ylevel - ypestal - lxs->FloorSurfaceYpos;
  G4String supname = mname + "Support";
  G4AssemblyVolume *tablesupport = LxAux::BuildTable(supname, tblx, tblhight, tblz, 3, 0.0);
  G4ThreeVector trsupport(pos.x(), -ylevel, pos.z());
  tablesupport->MakeImprint(lWorld, trsupport, 0, 0, lxs->OverlapTest);

  G4String pedname = mname + "Support";
  G4LogicalVolume *pedstal = LxAux::BuildPedestal("Profiler", tblx, ypestal, tblz);
  new G4PVPlacement (0, G4ThreeVector(pos.x(), 0.5*ypestal + lxs->FloorSurfaceYpos, pos.z()), pedstal,
                      pedname, lWorld, false, 0, lxs->OverlapTest);
}



PDSMagnetField::PDSMagnetField(const G4ThreeVector pos): G4MagneticField(), fPosition(pos),
  fMin(G4ThreeVector(0.0, 0.0, 0.0)),
  fMax(G4ThreeVector(0.0, 0.0, 0.0)),
  fBField(G4ThreeVector(0.0, 0.0, 0.0))
{
  LXSetUp *lxs = LXSetUp::Instance();
//   fMin = G4ThreeVector(-lxs->GMFieldX/2.0, -lxs->GMagneteffY/2.0, -lxs->GMFieldLength/2.0);
//   fMax = G4ThreeVector(lxs->GMFieldX/2.0, lxs->GMagneteffY/2.0, lxs->GMFieldLength/2.0);

  G4double magnetz = lxs->GMagnetZ - 2.0*lxs->GMagCoilThickness;
  G4double magnetx = lxs->GMagnetCoreX;
  fMin = G4ThreeVector(-magnetx/2.0, -lxs->GMagneteffY/2.0, -magnetz/2.0);
  fMax = G4ThreeVector(magnetx/2.0, lxs->GMagneteffY/2.0, magnetz/2.0);
//   fBField = G4ThreeVector(0.0, lxs->GMagFieldY, 0.0);
  fBField = G4ThreeVector(0.0, 1.24491*tesla, 0.0);

// Parameter for field approximation using Fermi Dirac like function
  std::vector<double> paramx{1.24491*tesla, -100.0*mm, 28.6578*mm, 100.0*mm, 28.9133*mm};
  fInvExpParams[0] = paramx;
  std::vector<double> paramz{1.24491*tesla, -615.807*mm, 28.6578*mm, 621.643*mm, 28.9133*mm};
  fInvExpParams[2] = paramz;
  
  fMin.setX(paramx[1]);
  fMin.setZ(paramz[1]);
  fMax.setX(paramx[3]);
  fMax.setZ(paramz[3]);
  
  LoadFieldData("/home/sqy/work/luxe/magnet/field_measurements/Feldkarten_TDC1_i300a.txt");
}



void PDSMagnetField::GetFieldValue(const G4double ppos[4], G4double *MagField) const 
{
  MagField[0] = MagField[1] = MagField[2] = 0.0;
//   if(ppos[0] != ppos[0]) return;

//   bool invol = true;
//   for (int ii = 2; ii >= 0; --ii) {
//     G4double cmpos = ppos[ii] - fPosition[ii];
//     invol = invol && (cmpos >= fMin[ii]) && (cmpos <= fMax[ii]);
//     if (!invol) return;
//   }

  MagField[0] = fBField.x();
//   MagField[1] = fBField.y();
  MagField[1] = fInvExpParams.at(0).at(0) * FField(ppos[2]-fPosition[2], 2) * FField(ppos[0]-fPosition[0], 0);
//   MagField[1] = DataField(0.0, ppos[2]-fPosition[2]) * FField(ppos[0]-fPosition[0], 0);
  MagField[2] = fBField.z();
  
//   std::cout << "Magnet pos: " << ppos[0] << "  " << ppos[1] << "  " << ppos[2] << std::endl;
  
//   MagField[1] = DataField(ppos[0]-fPosition[0], ppos[2]-fPosition[2]);
}



// G4double PDSMagnetField::FField(const G4double x, const G4int ii) const
// {
//   G4double tt = 20.0;
//   G4double ff1 = 1.0 + exp( (fMin[ii] - x)/tt );
//   G4double ff2 = 1.0 + exp( (x - fMax[ii])/tt );
//   return 1.0 / (ff1 * ff2);
// }


G4double PDSMagnetField::FField(const G4double x, const G4int ii) const
{
  const std::vector<double> &pp = fInvExpParams.at(ii);
  G4double ff1 = 1.0 + exp( (pp[1] - x)/pp[2] );
  G4double ff2 = 1.0 + exp( (x - pp[3])/pp[4] );
  return 1.0 / (ff1 * ff2);
}



void PDSMagnetField::LoadFieldData(const G4String fname)
{
  std::fstream  fdataf;

  fdataf.open(fname, std::ios::in);
  if (!fdataf.is_open()) {
    G4String msgstr("Failed to load magnetic field data from the file ");
    msgstr += fname;
    G4Exception("PDSMagnetField::", "LoadFieldData()", FatalException, msgstr.c_str());
  }

  double x, z, bf;
  while (!fdataf.eof()) {
      fdataf >> x >> z >> bf;
      if (fdataf.eof()) break;
      fFieldData.push_back(G4ThreeVector(x*mm, z*mm, bf*tesla));
//       std::cout << x << " " << z << "  " << bf << std::endl;
  }
  fdataf.close();
  G4cout << fFieldData.size() << " Field data points were loaded from the file " << fname << G4endl;
}



G4double PDSMagnetField::DataField(const G4double xx, const G4double yy) const
{
  std::vector<G4double> drv;
  std::for_each(fFieldData.begin(), fFieldData.end(), 
                [&](const G4ThreeVector &vv){drv.push_back( pow(vv.x()-xx, 2.0) + pow(vv.y()-yy, 2.0 ));} );
  std::vector<size_t> vind(drv.size());
  size_t ind = 0;
  std::generate(vind.begin(), vind.end(), [&](){return ind++;});
  point_cmp<G4double> pcmp(&drv);
  const size_t nsort = 10;
  std::partial_sort (vind.begin(), vind.begin()+nsort, vind.end(), pcmp);
   
  size_t ind3 = 2;
  while (CheckCollinear(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]]) 
//          && (!CheckInside(fFieldData[vind[0]], fFieldData[vind[1]], fFieldData[vind[ind3]], xx, yy)) 
        ) {
    ind3++;
    if (ind3 > nsort) { 
      G4String msgstr("Number of elements for partial sorting is not sufficient");
      G4Exception("PDSMagnetField::", "DataField()", FatalException, msgstr.c_str());
    }
  }

//   std::cout << "DataField: " << xx << "  " << yy << "  " << drv[vind[0]] 
//             << "  " << drv[vind[1]] << "  " << drv[vind[ind3]] << std::endl;
  return FindZ(fFieldData.at(vind[0]), fFieldData.at(vind[1]), fFieldData.at(vind[ind3]), xx, yy);
}



G4double PDSMagnetField::FindZ(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3, 
                               const double xx, const double yy) const
{
// Makes linear interpolation finding Z coordinate for a given (x,y) 
// as a point on the surface defined by three 3D points r1, r2, r3.
  G4ThreeVector dr1 = r2 -r1;
  G4ThreeVector dr2 = r3 -r1;

  G4ThreeVector normv = dr1.cross(dr2).unit();

  return r1.z() - (normv.x()*(xx-r1.x()) + normv.y()*(yy-r1.y())) / normv.z();
}



G4bool PDSMagnetField::CheckCollinear(const G4ThreeVector &r1, const G4ThreeVector &r2, const G4ThreeVector &r3) const
{
// Check if 2D projection to xy-plane of three 3D points are on one line.
  const double eps = 1.0e-15;
  G4ThreeVector dr1 = r2 -r1;
  G4ThreeVector dr2 = r3 -r1;
  dr1.setZ(0.0);  
  dr2.setZ(0.0);
  return ((dr1.cross(dr2).mag2()) < eps);
}

