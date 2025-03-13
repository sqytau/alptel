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

