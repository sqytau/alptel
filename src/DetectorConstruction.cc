//
/// \brief Implementation of the DetectorConstruction class
//

#include <algorithm>
#include <functional>
#include <utility>
#include <typeinfo>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4GenericTrap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
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
#include "LxDetector.hh"

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fAbsorberMaterial(0),fWorldMaterial(0),fDefaultWorld(true),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fSolidAbsorber(0),fLogicAbsorber(0),fPhysiAbsorber(0),
 fDetectorMessenger(0)
{
  // default parameter values of the calorimeter
  LXSetUp *lxs = LXSetUp::Instance();

  fAbsorberThickness = lxs->BTargetZ;
  fAbsorberSizeXY    = lxs->BTargetX;
  fZposAbs           = lxs->BTargetZpos;
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   (LXSetUp::Instance()->EnvironmentMaterial);
  SetAbsorberMaterial("Tungsten");
//   SetAbsorberMaterial("Air20");
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);
}



DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
  for (auto itr : fDetList) {
    if (itr.second) {delete itr.second;} 
  } 
}


void DetectorConstruction::AddSensorSegmentation(const G4String sname, const G4double xsize, 
                                                 const G4double ysize, const G4int nx, const G4int ny)
{
  EventAction *evatmp = new EventAction();
  int res = evatmp->TestSegmentationEncoding(nx, ny);
  delete evatmp;
  if (!res) {
    fSensors[sname] = std::make_tuple (xsize, ysize, nx, ny);
  } else {
    G4String msgstr("Failed to add sensitive detector ");
    msgstr += sname;
    G4Exception("DetectorConstruction::", "AddSensorSegmentation()", FatalException, msgstr.c_str());
  }
}



G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}



void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials

  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //

  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
  G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
  G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
  G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

  //
  // define simple materials
  //

  new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
  new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
  new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

  G4Material* lAr = new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
  lAr->AddElement(Ar, natoms=1);

  new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
  new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
  new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
  new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //

  G4Material* H2O = new G4Material("Water",density= 1.000*g/cm3,ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78*eV);

  G4Material* CH = new G4Material("Plastic",density= 1.04*g/cm3,ncomponents=2);
  CH->AddElement(C, natoms=1);
  CH->AddElement(H, natoms=1);

  G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* Air20 = 
    new G4Material("Air20", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20->AddElement(N, fractionmass=0.7);
  Air20->AddElement(O, fractionmass=0.3);

  G4Material* Air20Refractive =
    new G4Material("Air20Refractive", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
  Air20Refractive->AddElement(N, fractionmass=0.7);
  Air20Refractive->AddElement(O, fractionmass=0.3);

  //Graphite
  //
  G4Material* Graphite = new G4Material("Graphite", density= 1.7*g/cm3, ncomponents=1);
  Graphite->AddElement(C, fractionmass=1.);

  //Havar
  //
  G4Element* Cr = new G4Element("Chrome", "Cr", z=25, a=  51.996*g/mole);
  G4Element* Fe = new G4Element("Iron"  , "Fe", z=26, a=  55.845*g/mole);
  G4Element* Co = new G4Element("Cobalt", "Co", z=27, a=  58.933*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28, a=  58.693*g/mole);
  G4Element* W  = new G4Element("Tungsten","W", z=74, a= 183.850*g/mole);

  G4Material* Havar = new G4Material("Havar", density= 8.3*g/cm3, ncomponents=5);
  Havar->AddElement(Cr, fractionmass=0.1785);
  Havar->AddElement(Fe, fractionmass=0.1822);
  Havar->AddElement(Co, fractionmass=0.4452);
  Havar->AddElement(Ni, fractionmass=0.1310);
  Havar->AddElement(W , fractionmass=0.0631);

  //
  // examples of gas
  //
  G4Material* ArgonGas = new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);

  G4Material* HeliumGas = new G4Material("HeliumGas", z=2, a=4.003*g/mole, density = 0.166322 *mg/cm3, kStateGas, 293.15*kelvin, 1*atmosphere);

  new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, 293.15*kelvin, 1*atmosphere);

  G4Material* CO2 =
    new G4Material("CarbonicGas", density= 1.977*mg/cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* ArCO2 =
    new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
  ArCO2->AddElement (Ar,  fractionmass=0.7844);
  ArCO2->AddMaterial(CO2, fractionmass=0.2156);

  //another way to define mixture of gas per volume
  G4Material* NewArCO2 =
    new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
  NewArCO2->AddElement (Ar, natoms=8);
  NewArCO2->AddElement (C,  natoms=2);
  NewArCO2->AddElement (O,  natoms=4);

  G4Material* ArCH4 = 
    new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
  ArCH4->AddElement (Ar, natoms=93);
  ArCH4->AddElement (C,  natoms=7);
  ArCH4->AddElement (H,  natoms=28);

  G4Material* XeCH = 
    new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                   kStateGas, 293.15*kelvin, 1*atmosphere);
  XeCH->AddElement (Xe, natoms=875);
  XeCH->AddElement (C,  natoms=225);
  XeCH->AddElement (H,  natoms=700);

  G4Material* steam = new G4Material("WaterSteam", density= 1.0*mg/cm3, ncomponents=1);
  steam->AddMaterial(H2O, fractionmass=1.);
  steam->GetIonisation()->SetMeanExcitationEnergy(71.6*eV);

  //
  // example of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1, a=1.01*g/mole,density, kStateGas,temperature,pressure);

  // Define Stainless steal

  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.09 *g/mole);
  G4Element* Mn  = new G4Element("Manganese","Mn", z=25, a=54.938044 *g/mole);

  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
  StainlessSteel->AddElement(C, fractionmass=0.001);
  StainlessSteel->AddElement(Si, fractionmass=0.007);
  StainlessSteel->AddElement(Cr, fractionmass=0.18);
  StainlessSteel->AddElement(Mn, fractionmass=0.01);
  StainlessSteel->AddElement(Fe, fractionmass=0.712);
  StainlessSteel->AddElement(Ni, fractionmass=0.09);  

  G4Element* Al = 	new G4Element("Aluminum", "Al", 13., 26.9815386 * g/mole);
  G4Element* Ce = 	new G4Element("Cerium", "Ce", 58., 140.12 * g/mole);
  G4Element* Ca = 	new G4Element("Calcium", "Ca", 20.0, 40.078 * g/mole);

     // Concrete, must  check recipe for concrete

  G4double crdensity = 2.5*g/cm3;
  G4Material* ShieldingConcrete = new G4Material("ShieldingConcrete", crdensity, 6);
  ShieldingConcrete->AddElement(O,  0.52);
  ShieldingConcrete->AddElement(Si, 0.325);
  ShieldingConcrete->AddElement(Ca, 0.06);
  ShieldingConcrete->AddElement(Na, 0.015);
  ShieldingConcrete->AddElement(Fe, 0.04);
  ShieldingConcrete->AddElement(Al, 0.04);

  G4Material* fPstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  fPstyrene->AddElement(C, 8);
  fPstyrene->AddElement(H, 8);
  G4Material* lanex;
  lanex = G4NistManager::Instance()->FindOrBuildMaterial("G4_GADOLINIUM_OXYSULFIDE");

  LXSetUp *lxs = LXSetUp::Instance();

  G4Element* elLu = new G4Element("Lutetium", "Lu", 71, 176.944  *g/mole);
  G4Element* elY  = new G4Element( "Yttrium",  "Y", 39,  88.90585*g/mole);

  G4Material* lysoMaterial = new G4Material("LYSO", 7.1*g/cm3, 5, kStateSolid, 300.0 * kelvin);
  lysoMaterial->AddElement(elLu, 71.43/100.0);
  lysoMaterial->AddElement(elY,   4.03/100.0);
  lysoMaterial->AddElement(Si,  6.37/100.0);
  lysoMaterial->AddElement(O,  18.14/100.0);
  lysoMaterial->AddElement(Ce,  0.02/100.0);

  G4Material* lanexGammaDetMaterial = new G4Material("LANEX", density= 7.5*g/cm3, ncomponents=2);
  G4Material* Gos_base = G4NistManager::Instance()->FindOrBuildMaterial("G4_GADOLINIUM_OXYSULFIDE");
  G4Element* Tb  = G4NistManager::Instance()->FindOrBuildElement(65);
  lanexGammaDetMaterial->AddElement(Tb, 0.05*perCent);
  lanexGammaDetMaterial->AddMaterial(Gos_base, 99.95*perCent);

  density     = 1.250 *mg/cm3;
  temperature = 300.*kelvin;
  pressure    = 1.0*atmosphere;
  G4Material* CO = new G4Material("CarbonMonoxide", density, ncomponents=2, kStateGas, temperature, pressure);
  CO->AddElement(C, natoms=1);
  CO->AddElement(O, natoms=1);

  density     = 8.376E-05 *g/cm3;
  pressure    = 1.0*atmosphere;
  G4Material* H2Gass = new G4Material("HydrogenGass", density, ncomponents=1, kStateGas, temperature, pressure);
  H2Gass->AddElement(H, natoms=2);

  density     = 2.89e-15 *mg/cm3;
  pressure    = 1.0e-6 *pascal;
  G4Material* xfelVacuum = new G4Material("XFELVacuum", density, ncomponents=2, kStateGas, temperature, pressure);
  xfelVacuum->AddMaterial(CO, 0.78);
  xfelVacuum->AddMaterial(H2Gass, 0.22);
}



void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fZstartAbs = fZposAbs-0.5*fAbsorberThickness;
  fZendAbs   = fZposAbs+0.5*fAbsorberThickness;

  if (fDefaultWorld) {
     fWorldSizeZ = 1.5*m;
     fWorldSizeXY= 1.5*m;
  }
}


  
G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{ 
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // World
  //
  fSolidWorld = new G4Box("World", fWorldSizeXY/2.0, fWorldSizeXY/2.0, fWorldSizeZ/2.0); 
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fWorldMaterial, "World");
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, true);
  
  ConstructLuxeDetectors();

//   AssignRegions();
// ////////////////////
//   DumpBFieldModel();
  
//   G4double maxStep = 10.0*nm;
//   G4double maxStep = 1.0*um;
//   fStepLimit = new G4UserLimits(maxStep);
//   fStepLimit->SetMaxAllowedStep(maxStep);
//   fLogicAbsorber->SetUserLimits(fStepLimit);
  
  PrintCalorParameters();         
  
  return fPhysiWorld;
}



void DetectorConstruction::ConstructLuxeDetectors()
{

//   fDetList["WISDetectorTele"] = new WISDetectorTele(this);
  fDetList["WISDetectorTeleFrame"] = new WISDetectorTeleFrame(this);

  for (auto &nd : fDetList) { nd.second->Construct(); }

  LXSetUp *lxs = LXSetUp::Instance();
  fAbsorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial("StainlessSteel");

//   G4Box* solidCollimator2 = new G4Box("solidCollimator2", lxs->Collimator2X/2.0, lxs->Collimator2Y/2.0, lxs->BTargetZ/2.0);

  G4Box* solidCollimator2a = new G4Box("solidCollimator2a", lxs->Collimator2X/2.0, lxs->Collimator2Y/2.0, lxs->BTargetZ/2.0);
  G4Tubs *solidColl2Hole = new G4Tubs("solidColl2Hole", 0.0, lxs->Collimator1HoleR, lxs->BTargetZ, 0.0, 2.0*M_PI);
  G4Transform3D transform(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid *solidCollimator2 = new G4SubtractionSolid("solidCollimator2",
                                 solidCollimator2a, solidColl2Hole, transform);

  G4LogicalVolume* logicCollimator2 = new G4LogicalVolume(solidCollimator2, fAbsorberMaterial, "logicCollimator2");
  G4double collimator2Zpos = -lxs->BTargetZ/2.0;
  fPhysiAbsorber = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, collimator2Zpos), logicCollimator2, "Collimator2",
                                                                      fLogicWorld, false, 0, lxs->OverlapTest);

}



void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << fWorldMaterial    << G4endl;
  G4cout << "\n" << fAbsorberMaterial << G4endl;
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMaterial->GetName();
  G4cout << ". The transverse size (XY) of the world is " 
         << G4BestUnit(fWorldSizeXY,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(fAbsorberThickness,"Length")
         << " of " << fAbsorberMaterial->GetName();
  G4cout << ". The transverse size (XY) is " 
         << G4BestUnit(fAbsorberSizeXY,"Length") << G4endl;
  G4cout << " Z position of the middle of the absorber "
         << G4BestUnit(fZposAbs,"Length");
  G4cout << G4endl;
}



void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fAbsorberMaterial != pttoMaterial) {
    fAbsorberMaterial = pttoMaterial;                  
    if(fLogicAbsorber) fLogicAbsorber->SetMaterial(fAbsorberMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}



void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && fWorldMaterial != pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
//     if(fLogicMagnet) fLogicMagnet->SetMaterial(fWorldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}



void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  LXSetUp *lxs = LXSetUp::Instance();
  fAbsorberThickness = val;
  lxs->BTargetZ = fAbsorberThickness;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetAbsorberSizeXY(G4double val)
{
  fAbsorberSizeXY = val;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldSizeZ = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetWorldSizeXY(G4double val)
{
  fWorldSizeXY = val;
  fDefaultWorld = false;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fZposAbs  = val;
  ComputeCalorParameters();
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
}



void DetectorConstruction::AssignRegions()
{ 
  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;

  LXSetUp *lxs = LXSetUp::Instance();

  regName = "BDumpRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  if (!region) {
    region = new G4Region(regName);
  }
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(lxs->BeamDumpProductionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(lxs->BeamDumpProductionCut, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(lxs->BeamDumpProductionCut, G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(lxs->BeamDumpProductionCut); // same cuts for gamma, proton, e- and e+
  region->SetProductionCuts(cuts);
  G4LogicalVolume *bdump = G4LogicalVolumeStore::GetInstance()->GetVolume("logicBeamDump");
  if (bdump) {
    bdump->SetRegion(region);
    region->AddRootLogicalVolume(bdump);
  }

}


