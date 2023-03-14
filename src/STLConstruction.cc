//
/// \brief Implementation of the DetectorConstruction class
//

#include <vector>
#include <string>
#include <cstdlib>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "G4Exception.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include "LXSetUp.hh"
#include "STLConstruction.hh"
#include "LxDetector.hh"

#include "CADMesh.hh"

////////////////////////////////////////////////////////////////////////
//// STLConstruction


STLConstruction::STLConstruction(DetectorConstruction *detc): LxDetector(detc)
{
  std::vector<std::string> files{ 
      "BackWall_Concrete.stl",
      "BeamWallA_Concrete.stl",
      "BeamWallB_Concrete.stl",
      "BeamWallC_Concrete.stl",
      "CeillingA_Concrete.stl",
      "CeillingB_Concrete.stl",
      "ConcreteBlockShaft_Concrete.stl",
      "Crane_Concrete.stl",
      "Floor_Concrete.stl",
      "OuterWallA_Concrete.stl",
      "OuterWallB_Concrete.stl",
      "RadiationDoor_Concrete.stl",
      "RampA_Aluminium.stl",
      "RampB_Aluminium.stl",
      "StairCaseA_Concrete.stl",
      "StairCaseB_Concrete.stl"
  };
  
  fFileList = files;
}

void STLConstruction::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();
  
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* concreteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");
  G4Material* AluminiumMaterial = G4NistManager::Instance()->FindOrBuildMaterial("Aluminium");
  
  char *envpath = getenv("LUXE_DIR");
  if (!envpath) {
    G4String msgstr("Envitonment variable LUXE_DIR is not defined! Using ./");
    G4Exception("DetectorConstruction::", "ConstructWalls()", JustWarning, msgstr.c_str());
  }
  std::string path=std::string(envpath ? envpath : ".")+std::string("/stl/");
  G4cout<<"Debug::Build wall in::"+path<<G4endl;

  for (auto filesi : fFileList){	
	G4Material* mat;
	std::string str="";
	
	if(filesi.find("_Concrete")==std::string::npos){
		mat=AluminiumMaterial;
		str=filesi.substr (0,filesi.find("_Aluminum"));
		G4cout<<"Debug::Build wall in Aluminum::"+str<<G4endl;
	}
	else{
		mat=concreteMaterial;
		str=filesi.substr (0,filesi.find("_Concrete"));
		G4cout<<"Debug::Build wall in Concrete::"+str<<G4endl;
	}
	G4cout<<path+filesi<<G4endl;
    auto element = CADMesh::TessellatedMesh::FromSTL(path+filesi);
    //element->SetScale(1.0);

    auto element_logical = new G4LogicalVolume( element->GetSolid() , mat, "logical_"+str, 0, 0, 0);
    new G4PVPlacement( 0, G4ThreeVector(0, 0, 0), element_logical, str, fLogicWorld, false, 0 ,lxs->OverlapTest);
  }
}


