//
/// \brief Implementation of the LxUG03 class
//

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "LXSetUp.hh"
#include "LxDetector.hh"
#include "LxUG03.hh"


////////////////////////////////////////////////////////////////////////
//// LxUG03

void LxUG03::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  LXSetUp *lxs = LXSetUp::Instance();
  G4double floory = lxs->FloorSurfaceYpos;
  G4double shieldfez = 300.0*mm;
  G4double shieldpez = 150.0*mm;
  G4double shieldy = 2500.0*mm;

//Primary solids
  std::map<std::string, std::tuple<std::string, std::vector<G4double> > > ugb;
  ugb["FloorA"]       = std::make_tuple("Box", std::vector<G4double>({10000.0*mm, 1500.0*mm, 40000.0*mm}));
  ugb["StairCaseB"]   = std::make_tuple("Box", std::vector<G4double>({2000.0*mm, 6350.0*mm, 5000.0*mm}));
  ugb["StairCaseA"]   = std::make_tuple("Box", std::vector<G4double>({2850.0*mm, 6350.0*mm, 1350.0*mm}));
  ugb["OuterWallA"]   = std::make_tuple("Box", std::vector<G4double>({1000.0*mm, 6350.0*mm, 25000.0*mm}));
  ugb["OuterWallB"]   = std::make_tuple("Box", std::vector<G4double>({1000.0*mm, 4311.0*mm, 11500.0*mm}));
  ugb["BackWall"]     = std::make_tuple("Box", std::vector<G4double>({5550.0*mm, 4311.0*mm, 1000.0*mm}));
  ugb["BeamWallC"]    = std::make_tuple("Box", std::vector<G4double>({1500.0*mm, 4311.0*mm, 12500.0*mm}));
  ugb["BeamWallB"]    = std::make_tuple("Box", std::vector<G4double>({1350.0*mm, 6350.0*mm, 9012.0*mm}));
  ugb["BeamWallA1"]   = std::make_tuple("Box", std::vector<G4double>({1350.0*mm, 6350.0*mm, 7012.0*mm}));
  ugb["BeamWallACut"] = std::make_tuple("Box", std::vector<G4double>({2700.0*mm, 8000.0*mm, 1000.0*mm}));
  ugb["FloorB"]       = std::make_tuple("Box", std::vector<G4double>({1350.0*mm, 1088.0*mm, 11512.0*mm}));
  ugb["CeillingA"]    = std::make_tuple("Box", std::vector<G4double>({9500.0*mm, 2000.0*mm, 29000.0*mm}));
  ugb["CeillingB"]    = std::make_tuple("Box", std::vector<G4double>({8050.0*mm, 2000.0*mm, 11500.0*mm}));
  ugb["CeillingC"]    = std::make_tuple("Box", std::vector<G4double>({7050.0*mm, 2000.0*mm, 1000.0*mm}));
  ugb["VehicleRack1"]  = std::make_tuple("Box", std::vector<G4double>({800.0*mm, 2500.0*mm, 2000.0*mm}));
  ugb["VehicleRackCut"]  = std::make_tuple("Box", std::vector<G4double>({794.0*mm, 2494.0*mm, 1994.0*mm}));
  ugb["VehicleRackIn"]  = std::make_tuple("Box", std::vector<G4double>({794.0*mm, 2494.0*mm, 1994.0*mm}));
  ugb["XFELDumpPLC1"]  = std::make_tuple("Box", std::vector<G4double>({3500.0*mm, 2000.0*mm, 1000.0*mm}));
  ugb["XFELDumpPLCCut"]  = std::make_tuple("Box", std::vector<G4double>({3494.0*mm, 1994.0*mm, 994.0*mm}));
  ugb["XFELDumpPLCIn"]  = std::make_tuple("Box", std::vector<G4double>({3494.0*mm, 1994.0*mm, 994.0*mm}));
  ugb["XFELDumpPLCShieldingFe"]  = std::make_tuple("Box", std::vector<G4double>({5000.0*mm, shieldy, shieldfez}));
  ugb["XFELDumpPLCShieldingPE1"]  = std::make_tuple("Box", std::vector<G4double>({5000.0*mm, shieldy, shieldpez}));
  ugb["XFELDumpPLCShieldingPE2"]  = std::make_tuple("Box", std::vector<G4double>({5000.0*mm, shieldy, shieldpez}));
  
//Composit solids
  std::map<std::string, std::tuple<G4int, std::string, std::string> > compositsolid;
  compositsolid["BeamWallA"]    = std::make_tuple(0, "BeamWallA1", "BeamWallACut");
  compositsolid["VehicleRack"]    = std::make_tuple(0, "VehicleRack1", "VehicleRackCut");
  compositsolid["XFELDumpPLC"]    = std::make_tuple(0, "XFELDumpPLC1", "XFELDumpPLCCut");

//Logic volume materials
  std::map<std::string, std::string > material;
  material["FloorA"]       = "ShieldingConcrete";
  material["StairCaseB"]   = "ShieldingConcrete";
  material["StairCaseA"]   = "ShieldingConcrete";
  material["OuterWallA"]   = "ShieldingConcrete";
  material["OuterWallB"]   = "ShieldingConcrete";
  material["BackWall"]     = "ShieldingConcrete";
  material["BeamWallC"]    = "ShieldingConcrete";
  material["BeamWallB"]    = "ShieldingConcrete";
  material["BeamWallA"]    = "ShieldingConcrete";
  material["FloorB"]       = "ShieldingConcrete";
  material["CeillingA"]    = "ShieldingConcrete";
  material["CeillingB"]    = "ShieldingConcrete";
  material["CeillingC"]    = "ShieldingConcrete";
  material["VehicleRack"]  = "Iron";
  material["XFELDumpPLC"]  = "Iron";
  material["VehicleRackIn"]  = "Air20";
  material["XFELDumpPLCIn"]  = "Air20";
  material["XFELDumpPLCShieldingFe"]  = "Iron";
  material["XFELDumpPLCShieldingPE1"]  = "G4_POLYETHYLENE";
  material["XFELDumpPLCShieldingPE2"]  = "G4_POLYETHYLENE";

//Transformations for placements and for composit solids
  G4double yrot = -6.610000000000011*deg;
  G4double shieldzpos = -10570.0*mm;
  G4double shieldxpos = -3670.0*mm;
  G4double pedz = 0.5*(shieldfez+shieldpez) * cos(yrot);
  G4double pedx = 0.5*(shieldfez+shieldpez) * sin(fabs(yrot));

  std::map<std::string, std::tuple<G4double, G4double, G4double, G4double, G4double, G4double> > trns;
  trns["FloorA"]     = std::make_tuple(0.0, yrot, 0.0, -1070.52*mm, (-3275.0+47.0)*mm, -1040.21*mm);
  trns["StairCaseB"] = std::make_tuple(0.0, yrot, 0.0, (-3292.7085503990165 + 100.0)*mm,
                                       707.333*mm, (-19609.546097159982 + 1200.0)*mm);
  trns["StairCaseA"] = std::make_tuple(0.0, yrot, 0.0, (-5499.451656730336 + 100.0)*mm,
                                       701.502*mm, (-17521.06469811116 + 1200.0)*mm);
  trns["OuterWallA"] = std::make_tuple(0.0, yrot, 0.0, -5997.300876172876*mm, 709.0*mm,
                                       -4510.529795282763*mm);
  trns["OuterWallB"] = std::make_tuple(0.0, yrot, 0.0, -3889.771965292013*mm, -301.47*mm,
                                       13631.394944194511*mm);
  trns["BackWall"]   = std::make_tuple(0.0, yrot, 0.0, -6.684735700617555*mm, -307.287*mm,
                                       18445.934715733307*mm);
  trns["BeamWallC"]  = std::make_tuple(0.0, yrot, 0.0, 2844.0701965404655*mm, -310.857*mm,
                                       12350.514755128403*mm);
  trns["BeamWallB"]  = std::make_tuple(0.0, yrot, 0.0, 1834.0301742113813*mm, 697.0*mm,
                                       1631.5840052750361*mm);
  trns["BeamWallA"]  = std::make_tuple(0.0, -0.23098034806726217*deg, 0.0, 1303.320168838236*mm, 697.0*mm,
                                       -6275.462965825784*mm);
  trns["FloorB"]    = std::make_tuple(0.0, -0.23098034806726217*deg, 0.0, 1260.3765850945992*mm,
                                       (-1934.2499992931325+47.0)*mm, -15549.791176070099*mm);
  trns["CeillingA"] = std::make_tuple(0.0, yrot, 0.0, -1994.8760499745401*mm, 4892.02*mm,
                                      -6933.98412989145*mm);
  trns["CeillingB"] = std::make_tuple(0.0, yrot, 0.0, -385.1920100823297*mm, 2865.16*mm,
                                      13230.532752802359*mm);
  trns["CeillingC"] = std::make_tuple(0.0, yrot, 0.0, -490.34424265166126*mm, 2862.0*mm,
                                      6950.248013723539*mm);
  trns["VehicleRack"]   = std::make_tuple(0.0, yrot, 0.0, -5950.0*mm, (floory + 2500.0/2.0)*mm, -12000.0*mm);
  trns["VehicleRackIn"] = std::make_tuple(0.0, yrot, 0.0, -5950.0*mm, (floory + 2500.0/2.0)*mm, -12000.0*mm);
  trns["XFELDumpPLC"]   = std::make_tuple(0.0, yrot, 0.0, -3400.0*mm, (floory + 2000.0/2.0)*mm, -13500.0*mm);
  trns["XFELDumpPLCIn"] = std::make_tuple(0.0, yrot, 0.0, -3400.0*mm, (floory + 2000.0/2.0)*mm, -13500.0*mm);

  trns["XFELDumpPLCShieldingFe"] = std::make_tuple(0.0, yrot, 0.0, shieldxpos, (floory + shieldy/2.0)*mm, shieldzpos);
  trns["XFELDumpPLCShieldingPE1"] = std::make_tuple(0.0, yrot, 0.0, shieldxpos-pedx, (floory + shieldy/2.0)*mm, shieldzpos-pedz);
  trns["XFELDumpPLCShieldingPE2"] = std::make_tuple(0.0, yrot, 0.0, shieldxpos+pedx, (floory + shieldy/2.0)*mm, shieldzpos+pedz);

  trns["Cut_BeamWallA1_BeamWallACut"] = std::make_tuple(0.0, 6.37901965193275*deg, 0.0, 0.0, 0.0, 3930.53688727682*mm);
  trns["Cut_VehicleRack1_VehicleRackCut"] = std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  trns["Cut_XFELDumpPLC1_XFELDumpPLCCut"] = std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  std::map<std::string, G4VSolid*>  volsolid;

  for (const auto &ssolid : ugb) {
    std::string vname = std::string("solid")+ssolid.first;
    const auto &spar = ssolid.second;
    const std::string &solidtype = std::get<0>(spar);
    const auto &parv = std::get<1>(spar);
    if (solidtype == "Box") {
      volsolid[ssolid.first] = new G4Box(vname.c_str(), parv.at(0)/2.0, parv.at(1)/2.0, parv.at(2)/2.0);
    } else {
      G4String msgstr("Solid type ");
      msgstr += solidtype + G4String(" is not supported! Ignore it\n");
      G4Exception("LxUG03::", "Construct()", JustWarning, msgstr.c_str());
    }
  }

  for (const auto &csolid : compositsolid) {
    const G4int &soptype = std::get<0>(csolid.second);
    const std::string &solv1 = std::get<1>(csolid.second);
    const std::string &solv2 = std::get<2>(csolid.second);
    if (soptype == 0) { // subtraction
      std::string trnsname = std::string("Cut_") + solv1 + std::string("_") + solv2;
      G4double rx, ry, rz, tx, ty, tz;
      std::tie (rx, ry, rz, tx, ty, tz) = trns[trnsname];

      G4RotationMatrix rotm  = G4RotationMatrix();
      rotm.rotateX(rx);
      rotm.rotateY(ry);
      rotm.rotateZ(rz);
      G4Transform3D transform1(rotm, G4ThreeVector(tx, ty, tz));
      G4String vsname = G4String("solid") + csolid.first;
      volsolid[csolid.first] = new G4SubtractionSolid(vsname, volsolid[solv1], volsolid[solv2], transform1);
    } else {
      G4String msgstr("Composit solid type ");
      msgstr += std::to_string(soptype) + G4String(" is not supported! Ignore it\n");
      G4Exception("LxUG03::", "Construct()", JustWarning, msgstr.c_str());
    }
  }

  for (const auto &volm : material) {
    G4Material* vmt = G4NistManager::Instance()->FindOrBuildMaterial(volm.second);
    G4String lvname = G4String("logic") + volm.first;
    G4LogicalVolume *lv = new G4LogicalVolume(volsolid[volm.first], vmt, lvname);

    G4double rx, ry, rz, tx, ty, tz;
    std::tie (rx, ry, rz, tx, ty, tz) = trns[volm.first];
    G4RotationMatrix *rotm  = new G4RotationMatrix();
    rotm->rotateX(rx);
    rotm->rotateY(ry);
    rotm->rotateZ(rz);
    new G4PVPlacement(rotm, G4ThreeVector(tx, ty, tz), lv, volm.first, fLogicWorld, false, 0, lxs->OverlapTest);
  }

}



void LxUG03::CreateMaterial()
{

}

