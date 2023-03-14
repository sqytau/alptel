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
#include "LxAux.hh"
#include "LXSetUp.hh"
#include "LxInteractionChambers.hh"


//////////////////////////////////////////////////////////////////////////////////////////////
// CircularChamber

CircularChamber::CircularChamber(const G4String mtypename) :
  InteractionChamberAssembly(mtypename)
{
  ConstructInteractionChamber();  
}



void CircularChamber::ConstructInteractionChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* IPChamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->IPChamberMaterial);

  fIPChamberAssembly = new G4AssemblyVolume();
  
  G4Tubs *solidIPContainer = new G4Tubs("solidIPContainer", 0.0, lxs->IPRmax, lxs->IPHight/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPContainer = new G4LogicalVolume(solidIPContainer, vacuumMaterial, "logicIPContainer");
  
  const std::vector<G4double> &rhmin = lxs->IPChamberHolesRmin;
  const std::vector<G4double> &cutangles = lxs->IPChamberHolesPhi;
  
// Construct side surface with holes for in/out pipes
  for (G4int ii = 0; ii < cutangles.size(); ++ii) {

    G4double phi0 = (ii==0) ? (cutangles.back() - 2.0*M_PI) : cutangles[ii-1];
    G4double phi1 = (ii==cutangles.size()-1) ? (cutangles[0] + 2.0*M_PI) : cutangles[ii+1];

    G4double phi = cutangles[ii];
    G4double phis = 0.5*(phi0-phi);
    G4double dphi = 0.5*(phi1 - phi0);
    G4String secname("IPSide"); 
    secname += G4UIcommand::ConvertToString(ii);

    G4double hholecut = (lxs->IPRmax - sqrt(lxs->IPRmin*lxs->IPRmin - rhmin[ii]*rhmin[ii]));
    G4double dx = 0.5*hholecut + sqrt(lxs->IPRmin*lxs->IPRmin - rhmin[ii]*rhmin[ii]);
    G4Tubs *solidIPHCut = new G4Tubs(G4String("solidCut")+secname, 0.0, rhmin[ii], hholecut, 0.0, 2.0*M_PI);
    
    G4Tubs *solidIPSide = new G4Tubs(G4String("solid")+secname, lxs->IPRmin, lxs->IPRmax, lxs->IPHight/2.0, phis, dphi);
    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), 0.5*M_PI), G4ThreeVector(dx, 0.0, 0.0));
    G4SubtractionSolid* solidIPSideWithHole = new G4SubtractionSolid(G4String("solid")+secname+G4String("WithHole"),
                                                                       solidIPSide, solidIPHCut, transform);                       
    G4LogicalVolume *logicIPSideWithHole = new G4LogicalVolume(solidIPSideWithHole, IPChamberMaterial,
                                                                 G4String("logic")+secname+G4String("WithHole"));
    new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, -1.0), phi), G4ThreeVector(0.0, 0.0, 0.0), 
                         logicIPSideWithHole, secname+G4String("WithHole"), logicIPContainer, false, 0, lxs->OverlapTest);
  }

// Construct in/out pipes
  G4double cutha = asin(lxs->IPHight/(2.0*lxs->IPRmax));
  G4Tubs *solidIPChamberPipeCut = new G4Tubs("solidIPChamberPipeCut", lxs->IPRmin*cos(cutha), lxs->IPRmax,
                                                lxs->IPHight/2.0, -cutha, 2.0*cutha); 
  
  for (G4int ii = 0; ii < cutangles.size(); ++ii) {
    G4double pipermax = rhmin[ii] + lxs->BPipeThickness;
    G4double hholecut1 = (lxs->IPRmax - sqrt(lxs->IPRmax*lxs->IPRmax - pipermax*pipermax));
    G4double pipelength = lxs->IPChamberPipeLength + hholecut1;
    G4double phi = cutangles[ii];

    G4String pname("IPChambePipe"); 
    pname += G4UIcommand::ConvertToString(ii);

    G4Tubs *solidIPChamberPipe1 = new G4Tubs(G4String("solid")+pname, rhmin[ii], pipermax, 
                                                                    0.5*pipelength, 0.0, 2.0*M_PI);
    G4Tubs *solidIPChamberPipeVac1 = new G4Tubs(G4String("solidVac")+pname, 0.0, rhmin[ii], 
                                                                    0.5*pipelength, 0.0, 2.0*M_PI);
    G4double dx = 0.5*pipelength - hholecut1 + lxs->IPRmax;

    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), 0.5*M_PI), G4ThreeVector(0.0, 0.0, dx));
    G4SubtractionSolid* solidIPChamberPipe = new G4SubtractionSolid(G4String("solid")+pname+G4String("WithCut"),
                                                            solidIPChamberPipe1, solidIPChamberPipeCut, transform);                       

    G4SubtractionSolid* solidIPChamberPipeVac = new G4SubtractionSolid(G4String("solidVac")+pname+G4String("WithCut"),
                                                            solidIPChamberPipeVac1, solidIPChamberPipeCut, transform);                       
    
    G4LogicalVolume *logicIPChamberPipe = new G4LogicalVolume(solidIPChamberPipe, IPChamberMaterial,
                                                                G4String("logic")+pname);
    G4LogicalVolume *logicIPChamberPipeVac = new G4LogicalVolume(solidIPChamberPipeVac, vacuumMaterial,
                                                                G4String("logicVac")+pname);

    G4RotationMatrix *prot = new G4RotationMatrix(0.5*M_PI, 0.5*M_PI, -phi-0.5*M_PI);
    G4ThreeVector pipetrans(dx*cos(phi), dx*sin(phi), 0.0);
    fIPChamberAssembly->AddPlacedVolume(logicIPChamberPipe, pipetrans, prot);
    fIPChamberAssembly->AddPlacedVolume(logicIPChamberPipeVac, pipetrans, prot);
    
    // Flange
    G4Tubs *solidIPChamberPipeFlange = new G4Tubs(G4String("solidFlange")+pname, pipermax, 
                            pipermax+lxs->IPChamberPipeFlangedR, 0.5*lxs->IPChamberPipeFlangeThickness, 0.0, 2.0*M_PI);
    G4LogicalVolume *logicIPChamberPipeFlange = new G4LogicalVolume(solidIPChamberPipeFlange, IPChamberMaterial,
                                                                G4String("logicFlange")+pname);

    G4double dxfl = lxs->IPRmax + lxs->IPChamberPipeLength - 0.5*lxs->IPChamberPipeFlangeThickness;
    G4ThreeVector flangerans(dxfl*cos(phi), dxfl*sin(phi), 0.0);
    fIPChamberAssembly->AddPlacedVolume(logicIPChamberPipeFlange, flangerans, prot);
  }
  
  //Bottom
  G4Tubs *solidIPChamberBottom = new G4Tubs("solidIPChamberBottom", 0.0, lxs->IPRmax, 
                                                 0.5*lxs->IPChamberBottomThickness, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPChamberBottom = new G4LogicalVolume(solidIPChamberBottom, IPChamberMaterial,
                                                                "logicIPChamberBottom");
  G4ThreeVector trbtm(0.0, 0.0, -0.5*(lxs->IPChamberBottomThickness+lxs->IPHight));
  fIPChamberAssembly->AddPlacedVolume(logicIPChamberBottom, trbtm, 0);
  
  //Top flange and cover
  G4Tubs *solidIPChamberTopFlange = new G4Tubs("solidIPChamberTopFlange", lxs->IPRmax, lxs->IPChamberTopR, 
                                                 0.5*lxs->IPChamberTopThickness, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPChamberTopFlange = new G4LogicalVolume(solidIPChamberTopFlange, IPChamberMaterial,
                                                                "logicIPChamberTopFlange");
  G4ThreeVector trtopflg(0.0, 0.0, 0.5*(lxs->IPHight-lxs->IPChamberTopThickness));
  fIPChamberAssembly->AddPlacedVolume(logicIPChamberTopFlange, trtopflg, 0);

  G4Tubs *solidIPChamberTop = new G4Tubs("solidIPChamberTop", 0.0, lxs->IPChamberTopR, 
                                                 0.5*lxs->IPChamberTopThickness, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPChamberTop = new G4LogicalVolume(solidIPChamberTop, IPChamberMaterial,
                                                                "logicIPChamberTop");
  G4ThreeVector trtop(0.0, 0.0, 0.5*(lxs->IPChamberTopThickness+lxs->IPHight));
  fIPChamberAssembly->AddPlacedVolume(logicIPChamberTop, trtop, 0);
  
  //Add circular container to the assembly
  G4ThreeVector tr0(0.0, 0.0, 0.0);
  fIPChamberAssembly->AddPlacedVolume(logicIPContainer, tr0, 0);
  
  // Support structures
  const G4int nsplane = 14;
  G4double zsf[nsplane] = {0.0, 0.02, 0.05, 0.45, 0.46, 0.47, 0.475, 0.525, 0.53, 0.54, 0.55, 0.97, 0.98, 1.0};
  G4double zsplane[nsplane];
  std::transform (zsf, zsf+nsplane, zsplane, std::bind(std::multiplies<double>(), std::placeholders::_1, lxs->IPChamberSupportH));
  G4double rr = lxs->IPChamberSupportR;
  G4double rsouter[nsplane] = {2.5*rr, 2.5*rr, 1.5*rr, rr, 0.5*rr, 0.5*rr, rr, rr, 0.5*rr, 0.5*rr, rr, rr, 1.5*rr, 1.5*rr};
  G4double rsinner[nsplane] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4Polycone *solidIPCSupport = new G4Polycone("solidIPChamberSupport", 0.0, 2.0*M_PI, nsplane, zsplane, rsinner, rsouter);
  G4LogicalVolume *logicIPChamberSupport = new G4LogicalVolume(solidIPCSupport, IPChamberMaterial, "logicIPChamberSupport");

  for (G4int ii = 0; ii < 3; ++ii) {
    G4double phi = static_cast<G4double>(ii) * 2.0*M_PI/3.0 - M_PI/6.0;
    G4ThreeVector vv(lxs->IPChamberSupportRpos * cos(phi), lxs->IPChamberSupportRpos * sin(phi), lxs->FloorSurfaceYpos);
    fIPChamberAssembly->AddPlacedVolume(logicIPChamberSupport, vv, 0);
  }
  
}


//////////////////////////////////////////////////////////////////////////////////////////////
// LxInteractionChamber

LxInteractionChamber::LxInteractionChamber(const G4String mtypename) :
  InteractionChamberAssembly(mtypename) 
{
  ConstructInteractionChamber();  
}



void LxInteractionChamber::ConstructInteractionChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* tauIChamberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TAUIChamberMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);

  G4Box *solidTAUICContainer = new G4Box("solidTAUICContainer", lxs->TAUIChamberX/2.0, lxs->TAUIChamberY/2.0, lxs->TAUIChamberZ/2.0);
  G4LogicalVolume *logicTAUICContainer = new G4LogicalVolume(solidTAUICContainer, vacuumMaterial, "logicTAUICContainer");

  //side panels
  G4Box *solidTAUICSide = new G4Box("solidTAUICSide", lxs->TAUIChamberSideThickness/2.0, lxs->TAUIChamberY/2.0, 
                                    lxs->TAUIChamberZ/2.0);
  G4LogicalVolume *logicTAUICSide = new G4LogicalVolume(solidTAUICSide, tauIChamberMaterial, "logicTAUICSide");
  G4ThreeVector sidepos((lxs->TAUIChamberX-lxs->TAUIChamberSideThickness)/2.0, 0.0, 0.0);
  new G4PVPlacement (0, sidepos, logicTAUICSide, "TAUICSide", logicTAUICContainer, false, 0, lxs->OverlapTest); 
  sidepos.setX(-sidepos.getX());
  new G4PVPlacement (0, sidepos, logicTAUICSide, "TAUICSide", logicTAUICContainer, false, 1, lxs->OverlapTest); 

  //front and rear panels
  G4double beamtocy = lxs->TAUIChamberYpos;
  G4double dxf = lxs->TAUIChamberX - 2.0* lxs->TAUIChamberSideThickness;
  G4Box *solidTAUICFront1 = new G4Box("solidTAUICFront1", dxf/2.0, lxs->TAUIChamberY/2.0, lxs->TAUIChamberFrontThickness/2.0);
  G4Tubs *solidTAUICFPCut = new G4Tubs("solidTAUICFPCut", 0.0, lxs->TAUIChamberElPipeRIn,
                                       lxs->TAUIChamberFrontThickness, 0.0, 2.0*M_PI);
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0), G4ThreeVector(-lxs->TAUIChamberXpos, -beamtocy, 0.0));
  G4SubtractionSolid* solidTAUICFront = new G4SubtractionSolid("solidTAUICFront", solidTAUICFront1, solidTAUICFPCut, transform1);
  G4LogicalVolume *logicTAUICFront = new G4LogicalVolume(solidTAUICFront, tauIChamberMaterial, "logicTAUICFront");
  //rear panel
  G4double ipmbphx = lxs->FlashMFieldX/2.0 - lxs->IPMAgnetBeamPipeXGap;
  G4Box *solidTAUICFPRearCut = new G4Box("solidTAUICFPRearCut", ipmbphx-lxs->BPipeThickness,
                               lxs->FlashMagneteffY/2.0 - lxs->BPipeThickness, lxs->TAUIChamberFrontThickness);
  G4SubtractionSolid* solidTAUICFPRear = new G4SubtractionSolid("solidTAUICFPRear", solidTAUICFront1,
                                                                  solidTAUICFPRearCut, transform1);
  G4LogicalVolume *logicTAUICFPRear = new G4LogicalVolume(solidTAUICFPRear, tauIChamberMaterial, "logicTAUICFPRear");
  // Place front and rear panels
  G4ThreeVector frontpos(0.0, 0.0, (lxs->TAUIChamberZ - lxs->TAUIChamberFrontThickness)/2.0);
  new G4PVPlacement (0, frontpos, logicTAUICFPRear, "TAUICRear", logicTAUICContainer, false, 0, lxs->OverlapTest);
  frontpos.setZ(-frontpos.getZ());
  new G4PVPlacement (0, frontpos, logicTAUICFront, "TAUICFront", logicTAUICContainer, false, 0, lxs->OverlapTest);

  // IP volume for tracks interception
  G4Tubs *solidIPVolume = new G4Tubs("solidIPVolume", 0.0, lxs->IPVolumeR, lxs->IPVolumeZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicIPVolume = new G4LogicalVolume(solidIPVolume, vacuumMaterial, "logicIPVolume");
  new G4PVPlacement (0, G4ThreeVector(-lxs->TAUIChamberXpos, -lxs->TAUIChamberYpos, -lxs->TAUIChamberZpos), logicIPVolume, "IPVolume",
                     logicTAUICContainer, false, 0, lxs->OverlapTest);

  ConstructMirrors(logicTAUICContainer);

  //top and bottom panels
  G4Box *solidTAUICBottom = new G4Box("solidTAUICBottom", lxs->TAUIChamberBottomX/2.0, lxs->TAUIChamberBottomY/2.0,
                                       lxs->TAUIChamberBottomZ/2.0);
  G4LogicalVolume *logicTAUICBottom = new G4LogicalVolume(solidTAUICBottom, tauIChamberMaterial, "logicTAUICBottom");
  G4Box *solidTAUICTop = new G4Box("solidTAUICTop", lxs->TAUIChamberBottomX/2.0, lxs->TAUIChamberBottomY/2.0,
                                       lxs->TAUIChamberBottomZ/2.0);
  G4LogicalVolume *logicTAUICTop = new G4LogicalVolume(solidTAUICTop, tauIChamberMaterial, "logicTAUICTop");

  // In/out pipes
  G4Tubs *solidTAUIChamberBPipe = new G4Tubs("solidTAUIChamberBPipe", lxs->TAUIChamberElPipeRIn, lxs->TAUIChamberElPipeROut,
                                                                    0.5*lxs->TAUIChamberPipeLength, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicTAUIChamberBPipe = new G4LogicalVolume(solidTAUIChamberBPipe, beamPipeMaterial, "logicTAUIChamberBPipe");
  G4Tubs *solidTAUIChamberBPipeVac = new G4Tubs("solidTAUIChamberBPipeVac", 0.0, lxs->TAUIChamberElPipeRIn,
                                                                    0.5*lxs->TAUIChamberPipeLength, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicTAUIChamberBPipeVac = new G4LogicalVolume(solidTAUIChamberBPipeVac, vacuumMaterial, "logicTAUIChamberBPipeVac");
  // Flange
  G4Tubs *solidTAUIChamberBPipeFlange = new G4Tubs("solidTAUIChamberBPipeFlange", lxs->TAUIChamberElPipeROut,
                            lxs->TAUIChamberElPipeROut+lxs->TAUIChamberPipeFlangedH, 0.5*lxs->TAUIChamberPipeFlangedL, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicTAUIChamberBPipeFlange = new G4LogicalVolume(solidTAUIChamberBPipeFlange, beamPipeMaterial,
                                                                "logicTAUIChamberBPipeFlange");
  // Out pipe
  G4Box *solidTAUIChamberBPipeOut1 = new G4Box("solidTAUIChamberBPipeOut1", ipmbphx, lxs->FlashMagneteffY/2.0,
                                                0.5*lxs->TAUIChamberPipeLength);
  G4Box *solidTAUIChamberBPipeOutCut = new G4Box("solidTAUIChamberBPipeOutCut", ipmbphx-lxs->BPipeThickness,
                             lxs->FlashMagneteffY/2.0 - lxs->BPipeThickness, lxs->TAUIChamberPipeLength);
  G4SubtractionSolid* solidTAUIChamberBPipeOut = new G4SubtractionSolid("solidTAUIChamberBPipeOut",
                                               solidTAUIChamberBPipeOut1, solidTAUIChamberBPipeOutCut);
  G4LogicalVolume *logicTAUIChamberBPipeOut = new G4LogicalVolume(solidTAUIChamberBPipeOut, beamPipeMaterial,
                                                                  "logicTAUIChamberBPipeOut");
  G4Box *solidTAUIChamberBPipeOutVac = new G4Box("solidTAUIChamberBPipeOutVac", ipmbphx-lxs->BPipeThickness,
                             lxs->FlashMagneteffY/2.0 - lxs->BPipeThickness, 0.5*lxs->TAUIChamberPipeLength);
  G4LogicalVolume *logicTAUIChamberBPipeOutVac = new G4LogicalVolume(solidTAUIChamberBPipeOutVac, vacuumMaterial,
                                                                    "logicTAUIChamberBPipeOutVac");
  //support
  const G4int nsplane = 4;
  G4double zsf[nsplane] = {0.0, 0.97, 0.98, 1.0};
  G4double zsplane[nsplane];
  std::transform (zsf, zsf+nsplane, zsplane, std::bind(std::multiplies<double>(), std::placeholders::_1, lxs->TAUIChamberSupportH));
  G4double rr = lxs->TAUIChamberSupportR;
  G4double rsouter[nsplane] = {rr, rr, 0.5*rr, 0.5*rr};
  G4double rsinner[nsplane] = {0.0, 0.0, 0.0, 0.0};
  G4Polycone *solidTAUICSupport = new G4Polycone("solidTAUICSupport", 0.0, 2.0*M_PI, nsplane, zsplane, rsinner, rsouter);
  G4LogicalVolume *logicTAUICSupport = new G4LogicalVolume(solidTAUICSupport, tauIChamberMaterial, "logicTAUICSupport");

  //Add everything to the assembly  
  fIPChamberAssembly = new G4AssemblyVolume();
  G4RotationMatrix *srot = new G4RotationMatrix(0.0, 0.5*M_PI, 0.0);
  for (G4int ii = 0; ii < 4; ++ii) {
    G4double sdx = pow(-1.0, ii&1) * (lxs->TAUIChamberX/2.0 - lxs->TAUIChamberSupportR);
    G4double sdz = pow(-1.0, (ii>>1)&1) * (lxs->TAUIChamberZ/2.0 - lxs->TAUIChamberSupportR);
//     G4ThreeVector spos(sdx, -lxs->TAUIChamberSupportH - lxs->TAUIChamberBeamBottomY, sdz);
//     G4ThreeVector spos(0.0, -lxs->TAUIChamberSupportH - lxs->TAUIChamberBeamBottomY, 0.0);
    G4ThreeVector spos(sdx, -lxs->TAUIChamberSupportH - lxs->TAUIChamberY/2.0 - lxs->TAUIChamberBottomY, sdz);
    fIPChamberAssembly->AddPlacedVolume(logicTAUICSupport, spos, srot);
  }

  beamtocy = 0.0;
  G4ThreeVector ipcpos(0.0, beamtocy, 0.0);
  fIPChamberAssembly->AddPlacedVolume(logicTAUICContainer, ipcpos, 0);
  ipcpos.setY((lxs->TAUIChamberY + lxs->TAUIChamberBottomY)/2.0 + beamtocy);
  fIPChamberAssembly->AddPlacedVolume(logicTAUICTop, ipcpos, 0);
  ipcpos.setY(-(lxs->TAUIChamberY + lxs->TAUIChamberBottomY)/2.0 + beamtocy);
  fIPChamberAssembly->AddPlacedVolume(logicTAUICBottom, ipcpos, 0);

  G4ThreeVector ipcbppos(-lxs->TAUIChamberXpos, -lxs->TAUIChamberYpos, (lxs->TAUIChamberZ + lxs->TAUIChamberPipeLength)/2.0);
//   fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipe, ipcbppos, 0);
//   fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeVac, ipcbppos, 0);
  fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeOut, ipcbppos, 0);
  fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeOutVac, ipcbppos, 0);
  ipcbppos.setZ(-ipcbppos.getZ());
  fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipe, ipcbppos, 0);
  fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeVac, ipcbppos, 0);
  ipcbppos.setZ((lxs->TAUIChamberZ - lxs->TAUIChamberPipeFlangedL)/2.0 + lxs->TAUIChamberPipeLength);
//   fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeFlange, ipcbppos, 0);
  ipcbppos.setZ(-ipcbppos.getZ());
  fIPChamberAssembly->AddPlacedVolume(logicTAUIChamberBPipeFlange, ipcbppos, 0);

// Support for IC internal breadboard
  G4double tblhight = -lxs->TAUIChamberBeamBottomY - lxs->ICBBoardTableGapY - lxs->FloorSurfaceYpos;
  G4AssemblyVolume *icbboardtable = LxAux::BuildTable("ICBBoardTable", lxs->ICBBoardTableX, tblhight, lxs->ICBBoardTableZ, 4,
                                                      lxs->ICBBoardTableY, lxs->ICBBoardTableLegR);
  G4double ylevel = lxs->TAUIChamberY/2.0 + lxs->TAUIChamberBottomY + lxs->ICBBoardTableGapY;
  G4ThreeVector trbbtable(0.0, -ylevel, 0.0);
//   fIPChamberAssembly->AddPlacedAssembly(icbboardtable, trbbtable, 0);
  LxAux::AddAssmblyVolumes(fIPChamberAssembly, icbboardtable, trbbtable, 0);

  G4AssemblyVolume *icbboardSupport = LxAux::BuildTable("ICBBoardSupport", lxs->ICBBoardTableX,
                                                        lxs->ICBBoardTableGapY, lxs->ICBBoardTableZ, 4, 0.0);
  G4double yslevel = lxs->TAUIChamberY/2.0 + lxs->TAUIChamberBottomY;
  G4ThreeVector trbbsupport(0.0, -yslevel, 0.0);
//   fIPChamberAssembly->AddPlacedAssembly(icbboardSupport, trbbsupport, 0);
  LxAux::AddAssmblyVolumes(fIPChamberAssembly, icbboardSupport, trbbsupport, 0);

}



void LxInteractionChamber::ConstructMirrors(G4LogicalVolume *lgTAUICContainer)
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* tauICBBoardMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TAUIChamberBBoardMaterial);

  G4Box *solidTAUICBBoard = new G4Box("solidTAUICBBoard", lxs->TAUIChamberBBoardX/2.0, lxs->TAUIChamberBBoardY/2.0, lxs->TAUIChamberBBoardZ/2.0);
  G4LogicalVolume *logicTAUICBBoard = new G4LogicalVolume(solidTAUICBBoard, tauICBBoardMaterial, "logicTAUICBBoard");

  G4double bbposy = 0.5*(lxs->TAUIChamberY - lxs->TAUIChamberBBoardY) - lxs->TAUIChamberBBoardGapY;
  new G4PVPlacement (0, G4ThreeVector(0.0, -bbposy, 0.0), logicTAUICBBoard, "ICBBoard", lgTAUICContainer, false, 0, lxs->OverlapTest);

// Construct mirrors
  G4double mh = -lxs->TAUIChamberYpos + lxs->TAUIChamberY/2.0 - lxs->TAUIChamberBBoardY - lxs->TAUIChamberBBoardGapY;
  G4AssemblyVolume *mirrorAssembly;

  for (size_t mi = 0; mi < lxs->TAUIChamberMirrorR.size(); ++mi) {
    mirrorAssembly = ConstructMirrorAssembly(lxs->TAUIChamberMirrorR[mi], lxs->TAUIChamberMirrorD[mi], mh, lxs->TAUIChamberMirrorName[mi]);
    G4ThreeVector mpos(lxs->TAUIChamberMirrorXpos[mi]-lxs->TAUIChamberXpos, -lxs->TAUIChamberYpos,
                       lxs->TAUIChamberMirrorZpos[mi]-lxs->TAUIChamberZpos);
    mirrorAssembly->MakeImprint(lgTAUICContainer, mpos,
                    new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), lxs->TAUIChamberMirrorTheta[mi]), 0, lxs->OverlapTest);
  }
}



G4AssemblyVolume* LxInteractionChamber::ConstructMirrorAssembly(const G4double r, const G4double d, const G4double h, const G4String mname)
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* tauICMirrorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TAUIChamberMirrorMaterial);
  G4Material* tauICMirrorRingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TAUIChamberMirrorRingMaterial);
  G4Material* tauICMirrorHoldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->TAUIChamberMirrorHolderMaterial);

  G4String mrr("ICMirror");
  G4String mring("ICMirrorRing");
  G4String hold("ICMirrorHold");
  G4String sld("solid");
  G4String lgc("logic");

  G4double rh = lxs->TAUIChamberMirrorRingH;

  G4Tubs *solidICMirror = new G4Tubs(sld + mrr + mname, 0.0, r-rh, d/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicICMirror = new G4LogicalVolume(solidICMirror, tauICMirrorMaterial, lgc + mrr + mname);

//   G4Tubs *solidICMirrorRing = new G4Tubs(sld + mring + mname, r-rh, r, d/2.0, 0.0, 2.0*M_PI);
  G4Tubs *solidICMirrorRing = new G4Tubs(sld + mring + mname, r-rh, r, d/2.0, 0.0, 4.0*M_PI/3.0);
  G4LogicalVolume *logicICMirrorRing = new G4LogicalVolume(solidICMirrorRing, tauICMirrorRingMaterial, lgc + mring + mname);

  G4Tubs *solidICMirrorHold = new G4Tubs(sld + hold + mname, 0.0, d/2.0, (h-r)/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicICMirrorHold = new G4LogicalVolume(solidICMirrorHold, tauICMirrorHoldMaterial, lgc + hold + mname);

  G4AssemblyVolume *mirrorAssembly = new G4AssemblyVolume();

  G4ThreeVector mpos(0.0, 0.0, 0.0);
  mirrorAssembly->AddPlacedVolume(logicICMirror, mpos, 0);
  mirrorAssembly->AddPlacedVolume(logicICMirrorRing, mpos, new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 4.0*M_PI/3.0));
  G4ThreeVector hpos(0.0, -(h+r)/2.0, 0.0);
  mirrorAssembly->AddPlacedVolume(logicICMirrorHold, hpos, new G4RotationMatrix(G4ThreeVector(1.0, 0.0, 0.0), M_PI/2.0));

  return mirrorAssembly;
}



//////////////////////////////////////////////////////////////////////////////////////////////
// SimpleIPBox

SimpleIPBox::SimpleIPBox(const G4String mtypename) :
  InteractionChamberAssembly(mtypename) 
{
  ConstructInteractionChamber();  
}



void SimpleIPBox::ConstructInteractionChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();
  
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  
  G4Box *ipContainer = new G4Box("solidIPContainer", lxs->IPBoxX/2.0, lxs->IPBoxY/2.0, lxs->IPContainerZ/2.0);
  G4LogicalVolume *ipLogicalContainer = new G4LogicalVolume(ipContainer, worldMaterial, "logicIPContainer");
  
  G4Box *ipContainerVacuum = new G4Box("solidIPVacuumContainer", lxs->IPBoxX/2.0, lxs->IPBoxY/2.0, lxs->IPBoxZ/2.0+lxs->IPBoxThickness);
  G4LogicalVolume *ipLogicalContainerVacuum = new G4LogicalVolume(ipContainerVacuum, vacuumMaterial, "logicIPContainerVacuum");

  G4Box *solidBoxOuter = new G4Box("solidIPOuterBox", lxs->IPBoxX/2.0, lxs->IPBoxY/2.0, lxs->IPBoxZ/2.0);
  G4Box *solidBoxInner = new G4Box("solidIPInnerBox", (lxs->IPBoxX/2.0 - lxs->IPBoxThickness), (lxs->IPBoxY/2.0 - lxs->IPBoxThickness), 
                                                                                               lxs->IPBoxZ/2.0+ lxs->IPBoxThickness);
  G4SubtractionSolid* ipbox_xy = new G4SubtractionSolid("solidIPBoxXY", solidBoxOuter, solidBoxInner);                       
  G4LogicalVolume *logicIPBoxXY = new G4LogicalVolume(ipbox_xy, beamPipeMaterial, "logicIPBoxXY");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicIPBoxXY, "IPBoxXY", ipLogicalContainerVacuum, false, 0, lxs->OverlapTest); 

  // Hole for electron beam pipe 
  G4Box *solidIPFront = new G4Box("solidIPFront0", lxs->IPBoxX/2.0, lxs->IPBoxY/2.0, lxs->IPBoxThickness/2.0);
  G4Tubs *solidIPFrontBH = new G4Tubs("solidIPFrontBH", 0.0, lxs->BPipeR-lxs->BPipeThickness, lxs->IPBoxThickness, 0.0, 2.0*M_PI);
  G4SubtractionSolid* ipfront = new G4SubtractionSolid("solidIPFront", solidIPFront, solidIPFrontBH);                       

  // Hole for laser beam pipe 
  G4double dx = (lxs->IPBoxZ + lxs->IPBoxThickness)/2.0 * tan(lxs->BeamCrossingAngle);
  G4double dl = (lxs->BPipeR - lxs->BPipeThickness)*tan(lxs->BeamCrossingAngle) + 0.5*(lxs->IPBoxThickness/cos(lxs->BeamCrossingAngle));
  G4Tubs *solidIPFrontLH = new G4Tubs("solidIPFrontLH", 0.0, lxs->BPipeR-lxs->BPipeThickness, 1.3*dl, 0.0, 2.0*M_PI);
  G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), -lxs->BeamCrossingAngle), G4ThreeVector(dx, 0.0, 0.0));
  G4SubtractionSolid* ipfrontbhlh = new G4SubtractionSolid("solidIPFront", ipfront, solidIPFrontLH, transform);                       
  
  // Place front and back cover of the box 
  G4LogicalVolume *logicIPFront = new G4LogicalVolume(ipfrontbhlh, beamPipeMaterial, "logicIPFront");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, -(lxs->IPBoxZ + lxs->IPBoxThickness)/2.0), 
                                 logicIPFront, "IPFrontB", ipLogicalContainerVacuum, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), G4ThreeVector(0.0, 0.0, (lxs->IPBoxZ + lxs->IPBoxThickness)/2.0), 
                                 logicIPFront, "IPFrontF", ipLogicalContainerVacuum, false, 1, lxs->OverlapTest); 

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), 
                                 ipLogicalContainerVacuum, "IPBoxVac", ipLogicalContainer, false, 0, lxs->OverlapTest); 
  
  // Electron beam pipe pieces 
  G4Tubs *solidElectronPipe = new G4Tubs("solidElectronPipe", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR, 
                                         0.25*(lxs->IPContainerZ-lxs->IPBoxZ-2.0*lxs->IPBoxThickness), 0.0, 2.0*M_PI);
  G4LogicalVolume *logicElectronPipe = new G4LogicalVolume(solidElectronPipe, beamPipeMaterial, "logicElectronPipe");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0), 
                                 logicElectronPipe, "ElectronPipe", ipLogicalContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), 
                                 G4ThreeVector(0.0, 0.0, -(lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0), 
                                 logicElectronPipe, "ElectronPipe", ipLogicalContainer, false, 1, lxs->OverlapTest); 

  G4Tubs *solidElectronPipeVac = new G4Tubs("solidElectronPipeVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, 
                                         0.25*(lxs->IPContainerZ-lxs->IPBoxZ-2.0*lxs->IPBoxThickness), 0.0, 2.0*M_PI);
  G4LogicalVolume *logicElectronPipeVac = new G4LogicalVolume(solidElectronPipeVac, vacuumMaterial, "logicElectronPipeVac");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0), 
                                 logicElectronPipeVac, "ElectronPipeVac", ipLogicalContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI), 
                                 G4ThreeVector(0.0, 0.0, -(lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0), 
                                 logicElectronPipeVac, "ElectronPipeVac", ipLogicalContainer, false, 1, lxs->OverlapTest); 
  
  // Laser beam pipe pieces 
  G4CutTubs *solidLaserPipe = new G4CutTubs("solidLaserPipe", lxs->LaserPipeR-lxs->LaserPipeThickness, lxs->LaserPipeR, 
                              0.25*(lxs->IPContainerZ-lxs->IPBoxZ-2.0*lxs->IPBoxThickness)/cos(lxs->BeamCrossingAngle), 
                              0.0, 2.0*M_PI, G4ThreeVector(-sin(lxs->BeamCrossingAngle), 0, -cos(lxs->BeamCrossingAngle)),
                              G4ThreeVector(sin(lxs->BeamCrossingAngle), 0, cos(lxs->BeamCrossingAngle)));
  G4LogicalVolume *logicLaserPipe = new G4LogicalVolume(solidLaserPipe, beamPipeMaterial, "logicLaserPipe");
  G4double dxl = (lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0 * tan(lxs->BeamCrossingAngle);
  G4double dzl = (lxs->IPBoxZ + lxs->IPContainerZ + 2.0*lxs->IPBoxThickness)/4.0;
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), lxs->BeamCrossingAngle), G4ThreeVector(-dxl, 0.0, dzl),
                     logicLaserPipe, "LaserPipe", ipLogicalContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI+lxs->BeamCrossingAngle), G4ThreeVector(dxl, 0.0, -dzl),
                     logicLaserPipe, "LaserPipe", ipLogicalContainer, false, 1, lxs->OverlapTest); 
      
  G4CutTubs *solidLaserPipeVac = new G4CutTubs("solidLaserPipeVac", 0.0, lxs->LaserPipeR-lxs->LaserPipeThickness, 
                              0.25*(lxs->IPContainerZ-lxs->IPBoxZ-2.0*lxs->IPBoxThickness)/cos(lxs->BeamCrossingAngle), 
                              0.0, 2.0*M_PI, G4ThreeVector(-sin(lxs->BeamCrossingAngle), 0, -cos(lxs->BeamCrossingAngle)),
                              G4ThreeVector(sin(lxs->BeamCrossingAngle), 0, cos(lxs->BeamCrossingAngle)));
  G4LogicalVolume *logicLaserPipeVac = new G4LogicalVolume(solidLaserPipeVac, vacuumMaterial, "logicLaserPipeVac");
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), lxs->BeamCrossingAngle), G4ThreeVector(-dxl, 0.0, dzl),
                     logicLaserPipeVac, "LaserPipeVac", ipLogicalContainer, false, 0, lxs->OverlapTest); 
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), M_PI+lxs->BeamCrossingAngle), G4ThreeVector(dxl, 0.0, -dzl),
                     logicLaserPipeVac, "LaserPipeVac", ipLogicalContainer, false, 1, lxs->OverlapTest); 
  
  fIPChamberAssembly = new G4AssemblyVolume();
  G4ThreeVector iptr(0.0, 0.0, 0.0);
  fIPChamberAssembly->AddPlacedVolume(ipLogicalContainer, iptr, 0);

}




