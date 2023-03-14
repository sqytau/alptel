//
/// \brief Implementation of the LxAux class
//

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "LXSetUp.hh"
#include "LxAux.hh"


////////////////////////////////////////////////////////////////////////
//// LxAux

G4AssemblyVolume *LxAux::BuildTable(const G4String tname, const G4double xs, const G4double ys, 
                            const G4double zs, const int nleg, const G4double ytop, const G4double rleg)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* tableTopPlateMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPBasePlateMaterial);

  G4double syplate = ytop;
  if (syplate < 0.0)  syplate = lxs->OPPPBasePlateY;
  
  G4AssemblyVolume *tableAssembly = new G4AssemblyVolume();

  G4String legname = tname + "Leg";
  G4String sld = "solid";
  G4String lgc = "logic";

  if (syplate > 1.0e-15) {
    G4String topname = tname + "Top";
    G4Box *solidTopPlate = new G4Box(sld+topname, xs/2.0, syplate/2.0, zs/2.0);
    G4LogicalVolume *logicTopPlate = new G4LogicalVolume(solidTopPlate, tableTopPlateMaterial, lgc+topname);
    G4ThreeVector trplate(0.0, -syplate/2.0, 0.0);
    tableAssembly->AddPlacedVolume(logicTopPlate, trplate, 0);
  }

  G4double rout = rleg;
  if (rout < 0.0) rout = lxs->OPPPBasePlateR;
  G4double blegz = ys - syplate;
  G4Tubs *solidPlateLeg = new G4Tubs(sld+legname, 0.0, rout, blegz/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicPlateLeg = new G4LogicalVolume(solidPlateLeg, tableTopPlateMaterial, lgc+legname);

  int corrnleg = nleg;
  if (nleg > 4) corrnleg = 4;
  if (xs < 6.0 * rout && corrnleg > 2 ) corrnleg = 2;
  if (zs < 6.0 * rout ) corrnleg = 1;

  std::vector<G4ThreeVector> lmove;
  G4double dx = xs/2.0 - 1.5*rout;
  G4double dy = syplate + blegz/2.0;
  G4double dz = zs/2.0 - 1.5*rout;
  switch (corrnleg) {
    case 1 : lmove.push_back(G4ThreeVector(0.0, -dy, 0.0)); break;
    case 2 : lmove.push_back(G4ThreeVector(0.0, -dy, -dz));
             lmove.push_back(G4ThreeVector(0.0, -dy, dz)); break;
    case 3 : lmove.push_back(G4ThreeVector(-dx, -dy, -dz));
             lmove.push_back(G4ThreeVector(dx,  -dy, -dz));
             lmove.push_back(G4ThreeVector(0.0, -dy, dz));  break;
    case 4 : lmove.push_back(G4ThreeVector(-dx, -dy, -dz));
             lmove.push_back(G4ThreeVector(-dx, -dy, dz));
             lmove.push_back(G4ThreeVector(dx, -dy, -dz));
             lmove.push_back(G4ThreeVector(dx, -dy, dz));  break;
  }

  G4RotationMatrix *lrot = new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0);
  for (auto &trv : lmove) {
    tableAssembly->AddPlacedVolume(logicPlateLeg, trv, lrot);
  }

  return tableAssembly;
}



G4LogicalVolume *LxAux::BuildPedestal(const G4String pname, const G4double xs, const G4double ys, 
                                         const G4double zs)
{
//   LXSetUp *lxs = LXSetUp::Instance();
  G4Material* pedestalMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");

  G4String volname = pname + "Pedestal";
  G4String sld = "solid";
  G4String lgc = "logic";

  G4Box *solidPedestal = new G4Box(sld+volname, xs/2.0, ys/2.0, zs/2.0);
  G4LogicalVolume *logicPedestal = new G4LogicalVolume(solidPedestal, pedestalMaterial, lgc+volname);

  return logicPedestal;
}



G4AssemblyVolume *LxAux::BuildHexapod(const G4String pname, G4double &hexhight)
{
  //Hexapod
  // hexhight is first subtracted from the original hight of the hexapod (lxs->OPPPHexapodY)
  // and then updated to its hight as in the assmbly returned. Assmbly id centered in Y.
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppHexapodMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPHexapodMaterial);

  G4String topvolname = pname + "HexTop";
  G4String botvolname = pname + "HexBottom";
  G4String midvolname = pname + "HexMiddle";
  G4String sld = "solid";
  G4String lgc = "logic";

  G4Tubs *solidOPPPHexTop = new G4Tubs(sld + topvolname, 0.0, lxs->OPPPHexapodUpR, lxs->OPPPHexapodUpH/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPHexTop = new G4LogicalVolume(solidOPPPHexTop, opppHexapodMaterial, lgc + topvolname);
  G4Tubs *solidOPPPHexBottom = new G4Tubs(sld + botvolname, 0.0, lxs->OPPPHexapodDownR, lxs->OPPPHexapodDownH/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPHexBottom = new G4LogicalVolume(solidOPPPHexBottom, opppHexapodMaterial, lgc + botvolname);
  G4double hexmh = lxs->OPPPHexapodY - lxs->OPPPHexapodUpH - lxs->OPPPHexapodDownH - hexhight;
  G4double hexmrup = 0.9 * lxs->OPPPHexapodUpR;
  G4double hexmrdwn =  0.9 * lxs->OPPPHexapodDownR;
  hexhight = lxs->OPPPHexapodUpH + lxs->OPPPHexapodDownH + hexmh;

//   G4Cons *solidOPPPHexMiddle = new G4Cons("solidOPPPHexMiddle", hexmrdwn-lxs->OPPPHexapodThickness, hexmrdwn,
//                                           hexmrup-lxs->OPPPHexapodThickness, hexmrup, hexmh/2.0, 0.0, 2.0*M_PI);
//   G4LogicalVolume *logicOPPPHexMiddle = new G4LogicalVolume(solidOPPPHexMiddle, opppHexapodMaterial, "logicOPPPHexMiddle");

  G4double hsphi = M_PI/6.0;
  G4ThreeVector lupposl(hexmrup*cos(-hsphi/2.0),  hexhight/2.0 - lxs->OPPPHexapodUpH, hexmrup*sin(-hsphi/2.0));
  G4ThreeVector ldwposl(hexmrdwn*cos(hsphi/2.0), -hexhight/2.0 + lxs->OPPPHexapodDownH, hexmrdwn*sin(hsphi/2.0));
  G4ThreeVector lposdifl = lupposl - ldwposl;
  G4double hexsidel = lposdifl.mag();
  std::vector<G4double> rohexleg{0.0, 7.5,  7.5,  25.0, 25.0, 7.5,  7.5,  0.0};
  std::vector<G4double> rihexleg{0.0, 5.0,  5.0,  22.0, 22.0, 5.0,  5.0,  0.0};
  std::vector<G4double>  lhexleg{0.0, 0.04, 0.14, 0.16, 0.84, 0.86, 0.96, 1.0};
  std::for_each(lhexleg.begin(), lhexleg.end(), [=](G4double &x) { x *= hexsidel;});

  G4Polycone *solidOPPPHexMiddle = new G4Polycone(sld + midvolname, 0.0, 2.0*M_PI, rohexleg.size(), &lhexleg[0], &rihexleg[0], &rohexleg[0]);
  G4LogicalVolume *logicOPPPHexMiddle = new G4LogicalVolume(solidOPPPHexMiddle, opppHexapodMaterial, lgc + midvolname);

  G4AssemblyVolume *oppHexapodAssembly = new G4AssemblyVolume();

  G4ThreeVector hexpos(0.0, (hexhight - lxs->OPPPHexapodUpH)/2.0, 0.0);
  oppHexapodAssembly->AddPlacedVolume(logicOPPPHexTop, hexpos, new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));
  hexpos.setY(-(hexhight - lxs->OPPPHexapodDownH)/2.0);
  oppHexapodAssembly->AddPlacedVolume(logicOPPPHexBottom, hexpos, new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));

  for (int nl = 0; nl < 3; ++nl) {
    G4double phinl = nl * 2.0*M_PI/3.0;
    G4ThreeVector luppos(hexmrup*cos(phinl-hsphi/2.0),  hexhight/2.0 - lxs->OPPPHexapodUpH, hexmrup*sin(phinl-hsphi/2.0));
    G4ThreeVector ldwpos(hexmrdwn*cos(phinl+hsphi/2.0), -hexhight/2.0 + lxs->OPPPHexapodDownH, hexmrdwn*sin(phinl+hsphi/2.0));
    G4ThreeVector lposdif = luppos - ldwpos;
    G4RotationMatrix *rotlhex = new G4RotationMatrix();
    rotlhex->rotateY(lposdif.theta());
    rotlhex->rotateZ(lposdif.phi());
    oppHexapodAssembly->AddPlacedVolume(logicOPPPHexMiddle, ldwpos, rotlhex);
  }

  for (int nl = 0; nl < 3; ++nl) {
    G4double phinl = nl * 2.0*M_PI/3.0 - 2.0*hsphi;
    G4ThreeVector luppos(hexmrup*cos(phinl+hsphi/2.0),  hexhight/2.0 - lxs->OPPPHexapodUpH, hexmrup*sin(phinl+hsphi/2.0));
    G4ThreeVector ldwpos(hexmrdwn*cos(phinl-hsphi/2.0), -hexhight/2.0 + lxs->OPPPHexapodDownH, hexmrdwn*sin(phinl-hsphi/2.0));
    G4ThreeVector lposdif = luppos - ldwpos;
    G4RotationMatrix *rotlhex = new G4RotationMatrix();
    rotlhex->rotateY(lposdif.theta());
    rotlhex->rotateZ(lposdif.phi());
    oppHexapodAssembly->AddPlacedVolume(logicOPPPHexMiddle, ldwpos, rotlhex);
  }

  return oppHexapodAssembly;
}



G4AssemblyVolume* LxAux::BuildConeHexapod(const G4String pname, G4double &hexhight)
{
  //Hexapod with cone
  // hexhight is first subtracted from the original hight of the hexapod (lxs->OPPPHexapodY)
  // and then updated to its hight as in the assmbly returned. Assmbly id centered in Y.
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppHexapodMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPHexapodMaterial);

  G4String topvolname = pname + "HexTop";
  G4String botvolname = pname + "HexBottom";
  G4String midvolname = pname + "HexMiddle";
  G4String sld = "solid";
  G4String lgc = "logic";

  G4Tubs *solidOPPPHexTop = new G4Tubs(sld + topvolname, 0.0, lxs->OPPPHexapodUpR, lxs->OPPPHexapodUpH/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPHexTop = new G4LogicalVolume(solidOPPPHexTop, opppHexapodMaterial, lgc + topvolname);
  G4Tubs *solidOPPPHexBottom = new G4Tubs(sld + botvolname, 0.0, lxs->OPPPHexapodDownR, lxs->OPPPHexapodDownH/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPHexBottom = new G4LogicalVolume(solidOPPPHexBottom, opppHexapodMaterial, lgc + botvolname);
  G4double hexmh = lxs->OPPPHexapodY - lxs->OPPPHexapodUpH - lxs->OPPPHexapodDownH - hexhight;
  G4double hexmrup = 0.8 * lxs->OPPPHexapodUpR;
  G4double hexmrdwn =  0.8 * lxs->OPPPHexapodDownR;
  G4Cons *solidOPPPHexMiddle = new G4Cons(sld + midvolname, hexmrdwn-lxs->OPPPHexapodThickness, hexmrdwn,
                                          hexmrup-lxs->OPPPHexapodThickness, hexmrup, hexmh/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicOPPPHexMiddle = new G4LogicalVolume(solidOPPPHexMiddle, opppHexapodMaterial, lgc + midvolname);

  G4AssemblyVolume *oppHexapodAssembly = new G4AssemblyVolume();

  hexhight = lxs->OPPPHexapodUpH + lxs->OPPPHexapodDownH + hexmh;
  G4ThreeVector hexpos(0.0, (hexhight - lxs->OPPPHexapodUpH)/2.0, 0.0);
  oppHexapodAssembly->AddPlacedVolume(logicOPPPHexTop, hexpos, new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));
  hexpos.setY(hexpos.y() - (lxs->OPPPHexapodUpH + hexmh)/2.0);
  oppHexapodAssembly->AddPlacedVolume(logicOPPPHexMiddle, hexpos, new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));
  hexpos.setY(-(hexhight - lxs->OPPPHexapodDownH)/2.0);
  oppHexapodAssembly->AddPlacedVolume(logicOPPPHexBottom, hexpos, new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), M_PI/2.0));

  return oppHexapodAssembly;
}



void LxAux::AddAssmblyVolumes(G4AssemblyVolume* vasm, G4AssemblyVolume* cpyasm, G4ThreeVector &translation, G4RotationMatrix *rotation)
{
  G4RotationMatrix vrot;
  if (rotation) vrot = *rotation;
  G4Transform3D transformation(vrot, translation);

  std::vector<G4AssemblyTriplet>::iterator avitr = cpyasm->GetTripletsIterator();
  for( std::size_t i = 0; i < cpyasm->TotalTriplets(); ++i, ++avitr)
  {
    G4AssemblyTriplet triplet = *avitr;
    G4Transform3D Ta(*triplet.GetRotation(), triplet.GetTranslation());
    if ( triplet.IsReflection() )  { Ta = Ta * G4ReflectZ3D(); }

    G4Transform3D Tfinal = transformation * Ta;

    if ( triplet.GetVolume() )
    {
      vasm->AddPlacedVolume(triplet.GetVolume(), Tfinal);
    }
    else if ( triplet.GetAssembly() )
    {
      // Place volumes in this assembly with composed transformation
      G4Exception("LxAux::AddAssmblyVolumes(..)", "AssmblyVolume", FatalException, "Triplet with assembly is not supported");
//       G4ThreeVector atranslate(Tfinal.getTranslation());
//       AddAssmblyVolumes(vasm, triplet.GetAssembly(), atranslate, new G4RotationMatrix(Tfinal.getRotation()));
    }
    else
    {
      G4Exception("LxAux::AddAssmblyVolumes(..)", "GeomVol0003", FatalException, "Triplet has no volume and no assembly");
    }
  }
}

