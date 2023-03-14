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
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4IntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4DormandPrince745.hh"
#include "G4CashKarpRKF45.hh"

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
#include "LxMagnets.hh"
#include "LxVacChambers.hh"
#include "LxInteractionChambers.hh"
#include "LxDetector.hh"
#include "LxGammaCalo.hh"
#include "LxTrackerOPPP.hh"
#include "LxECal.hh"
#include "LxHICSBeam.hh"
#include "LxBeamProfiler.hh"
#include "LxScintCherenkov.hh"
#include "LxBSMDetector.hh"
#include "LxTargetChamber.hh"

#include "STLConstruction.hh"
#include "LxUG03.hh"
#include "LxBFields.hh"


DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fAbsorberMaterial(0),fWorldMaterial(0),fDefaultWorld(true),
 fSolidWorld(0),fLogicWorld(0),fPhysiWorld(0),
 fSolidAbsorber(0),fLogicAbsorber(0),fPhysiAbsorber(0),
 fDetectorMessenger(0), fMagnetZPos(0), fMagnetSizeZ(0), fTargetChamber(0)
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

  if (lxs->ScintCerenkovPhysics){

    // give argon, Helium, air optical properties
     G4double rIndexArgon[] = {1.000282, 1.000294, 1.000318, 1.0003675};
  G4double rEnergiesArgon[] = {2.184*eV, 3.982*eV, 5.9379*eV, 8.093*eV};

  G4double rIndexAir[] = {1.00027362, 1.00027670, 1.00028208, 1.00029024, 1.00030205};
  G4double rEnergiesAir[] = {1.*eV, 2.*eV, 3.*eV, 4.*eV, 5.*eV};

     G4double rIndexHelium[] = {1.00003469, 1.00003473, 1.000034886, 1.000035032};
  G4double rEnergiesHelium[] = {0.60216*eV, 1.0015*eV, 1.9927*eV, 2.5825*eV};

  G4MaterialPropertiesTable* fMPTArgon = new G4MaterialPropertiesTable();
  const G4int vectorsizeArgon = sizeof(rIndexArgon)/sizeof(G4double);
  fMPTArgon->AddProperty("RINDEX", rEnergiesArgon, rIndexArgon, vectorsizeArgon);

  ArgonGas->SetMaterialPropertiesTable(fMPTArgon);

  G4MaterialPropertiesTable* fMPTAir = new G4MaterialPropertiesTable();
  const G4int vectorsizeAir = sizeof(rIndexAir)/sizeof(G4double);
  fMPTAir->AddProperty("RINDEX", rEnergiesAir, rIndexAir, vectorsizeAir);

  Air20Refractive->SetMaterialPropertiesTable(fMPTAir);

  G4MaterialPropertiesTable* fMPTHelium = new G4MaterialPropertiesTable();
  const G4int vectorsizeHelium = sizeof(rIndexHelium)/sizeof(G4double);
  fMPTHelium->AddProperty("RINDEX", rEnergiesHelium, rIndexHelium, vectorsizeHelium);

  HeliumGas->SetMaterialPropertiesTable(fMPTHelium);

  // give Polystyrene optical properties

  G4double wls_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
  const G4int wlsnum = sizeof(wls_Energy)/sizeof(G4double);

  G4double rIndexPstyrene[]={ 1.5, 1.5, 1.5, 1.5};
  assert(sizeof(rIndexPstyrene) == sizeof(wls_Energy));
  G4double absorption1[]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  assert(sizeof(absorption1) == sizeof(wls_Energy));
  G4double scintilFast[]={0.00, 0.00, 1.00, 1.00};
  assert(sizeof(scintilFast) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* fMPTPStyrene = new G4MaterialPropertiesTable();
  fMPTPStyrene->AddProperty("RINDEX",wls_Energy,rIndexPstyrene,wlsnum);
  fMPTPStyrene->AddProperty("ABSLENGTH",wls_Energy,absorption1,wlsnum);
  fMPTPStyrene->AddProperty("FASTCOMPONENT",wls_Energy, scintilFast,wlsnum);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  fMPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);

  // Set the Birks Constant for the Polystyrene scintillator
  fPstyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

// Give Gadox optical properties
    G4MaterialPropertiesTable* MPT_lanex = new G4MaterialPropertiesTable();
    G4double lanex_energies[] = {1*eV, 3*eV, 5*eV};
    G4int vectorsizegadox = sizeof(lanex_energies)/sizeof(G4double);
    G4double lanex_refractiveIndex[] = {2.3, 2.3, 2.3};
    G4double lanex_absorption[] = {420*cm, 420*cm, 1*mm};
    G4double lanex_scints[] = {0.5, 5, 0};
    MPT_lanex->AddProperty("RINDEX", lanex_energies, lanex_refractiveIndex, vectorsizegadox);
    MPT_lanex->AddProperty("ABSLENGTH", lanex_energies, lanex_absorption, vectorsizegadox);
    MPT_lanex->AddProperty("SLOWCOMPONENT", lanex_energies, lanex_scints, vectorsizegadox);
    MPT_lanex->AddConstProperty("SCINTILLATIONYIELD",89000./MeV);
    MPT_lanex->AddConstProperty("SLOWTIMECONSTANT",380.0*microsecond);
    MPT_lanex->AddConstProperty("RESOLUTIONSCALE",1.0);
    G4double lanex_mie[] = {60.3e-3*mm, 60.3e-3*mm, 60.3e-3*mm};
    MPT_lanex->AddProperty("MIEHG", lanex_energies, lanex_mie, vectorsizegadox);
    MPT_lanex->AddConstProperty("MIEHG_FORWARD", 0.91);
    MPT_lanex->AddConstProperty("MIEHG_BACKWARD", 0.91);
    MPT_lanex->AddConstProperty("MIEHG_FORWARD_RATIO", 1.0);
    lanex->SetMaterialPropertiesTable(MPT_lanex);

  }

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
  
// // ConstructInfrastructure();
//   ConstructIPChamber();
//   ConstructBeamDump();
//   ConstructShielding();
//   ConstructBeamPipeToDump();
//   ConstructBremsTarget();
//   ConstructDumpMagnet();
//   ConstructVacuumChamber();
//   ConstructBeamPipeToIP();
//   ConstructIPMagnet();
  ConstructLuxeDetectors();
//   ConstructGammaMagnet();
//   ConstructGammaTarget();
//   ConstructComptShielding();
// //   ConstructCleaningMagnet();
//   AssignRegions();
// ////////////////////
  DumpBFieldModel();
  
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
//   fDetList["ComptonDet"] = new LxDetectorComptonFluka(this);
//   fDetList["GammaDet"] = new LxDetectorGammaCalo8(this);
//   fDetList["OpppTracker"] = new LxTrackerOPPP(this);
//   fDetList["LxECal"] = new LxECal(this);
//   fDetList["LxHICSBeamDump"] = new LxHICSBeam(this);
//   fDetList["BeamProfiler"] = new LxBeamProfiler(this);
//   fDetList["ScintillatorScreen"] = new LxSScreen(this);
//   fDetList["Cherenkov"] = new LxCherenkov(this);
//   fDetList["UG03"] = new LxUG03(this);
//   fDetList["BSMDetector"] = new LxBSMDetector(this);

//   fDetList["WISDetectorTele"] = new WISDetectorTele(this);
  fDetList["WISDetectorTeleFrame"] = new WISDetectorTeleFrame(this);

//   fDetList["STLComponents"] = new STLConstruction(this);

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



void DetectorConstruction::ConstructInfrastructure()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* concreteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");
  
  G4Box *solidFloor = new G4Box("solidFloor", fWorldSizeXY/2.0, lxs->FloorY/2.0, fWorldSizeZ/2.0);
  G4LogicalVolume *logicFloor = new G4LogicalVolume(solidFloor, concreteMaterial, "logicFloor");
  G4double ypos = lxs->FloorSurfaceYpos - lxs->FloorY/2.0;   // lxs->FloorSurfaceYpos < 0
  new G4PVPlacement(0, G4ThreeVector(0.0, ypos, 0), logicFloor, "Floor", fLogicWorld, false, 0, lxs->OverlapTest);
}



void DetectorConstruction::ConstructIPChamber()
{
  LXSetUp *lxs = LXSetUp::Instance();
  InteractionChamberAssembly *ipc = 0;
//   G4ThreeVector tr0(0.0, 0.0, 0.0);
  G4ThreeVector tr0(lxs->TAUIChamberXpos, lxs->TAUIChamberYpos, lxs->TAUIChamberZpos);
//   ipc = new CircularChamber("CircularChamber");
//   ipc->GetAssembly()->MakeImprint(fLogicWorld, tr0, new G4RotationMatrix(-0.5*M_PI, 0.5*M_PI, 0.0), 0, lxs->OverlapTest);
//   ipc = new SimpleIPBox("SimpleBox");
//   ipc->GetAssembly()->MakeImprint(fLogicWorld, tr0, new G4RotationMatrix(0.0, 0.0, 0.0), 0, lxs->OverlapTest);
  ipc = new LxInteractionChamber("LxInteractionChamber");
  ipc->GetAssembly()->MakeImprint(fLogicWorld, tr0, new G4RotationMatrix(0.0, 0.0, 0.0), 0, lxs->OverlapTest);
}



void DetectorConstruction::ConstructDumpMagnet(const G4String magType)
{
  LxMagnetAssembly *mag = 0; 
  LXSetUp *lxs = LXSetUp::Instance();

  G4String mtype(magType);
//   G4String mtype("TypMBMagnetSimple");
  auto magnetitr = fMagentsTypeCollection.find(mtype);
  if (magnetitr == fMagentsTypeCollection.end()) {
    mag = new FlashMagnetAssembly(mtype);
//     mag = new TypMBMagnetAssembly(mtype);
//     mag = new TypMBMagnetSimple(mtype);
    fMagentsTypeCollection[mtype] = mag;
  } else {
    mag = magnetitr->second;
  }

  G4LogicalVolume *fieldvol = mag->GetFieldVolume("logicDumpMagnetField");
//   G4ThreeVector b_field(0.0, lxs->DumpMagFieldY, 0.0);
//   ConstructMagnet(fieldvol, b_field.rotateZ(-0.5*M_PI));

  const G4Box  *mfbox = dynamic_cast<G4Box*>(fieldvol->GetSolid());
  G4double fhlength = mfbox->GetZHalfLength();
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
//   // Beam pipe cylindrical
//   G4Tubs *solidDMBPipe = new G4Tubs("solidDMBPipe", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                        fhlength, 0.0, 2.0*M_PI);
//   G4LogicalVolume *logicDMBPipe = new G4LogicalVolume(solidDMBPipe, beamPipeMaterial, "logicDMBPipe");
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicDMBPipe, "DMBPipe", fieldvol,
//                                          false, 0, lxs->OverlapTest);
//   G4Tubs *solidDMBPipeVac = new G4Tubs("solidDMBPipeVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
//                                           fhlength, 0.0, 2.0*M_PI);
//   G4LogicalVolume *logicDMBPipeVac = new G4LogicalVolume(solidDMBPipeVac, vacuumMaterial, "logicDMBPipeVac");
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicDMBPipeVac, "DMBPipeVac", fieldvol,
//                                          false, 0, lxs->OverlapTest);
  // Beam pipe rectangular
  G4Box *solidDMBPipeContainer = new G4Box("solidDMBPipeContainer", 2.0*(lxs->BPipeR+lxs->BPipeThickness),
                                           lxs->BPipeR, fhlength);
  G4LogicalVolume *logicDMBPipeContainer = new G4LogicalVolume(solidDMBPipeContainer,
                                                               vacuumMaterial, "logicDMBPipeContainer");
  G4Box *solidDMBPipeOuter = new G4Box("solidDMBPipeOuter", 2.0*(lxs->BPipeR+lxs->BPipeThickness),
                                       lxs->BPipeR, fhlength);
  G4Box *solidDMBPipeInner = new G4Box("solidDMBPipeInner", 2.0*lxs->BPipeR+lxs->BPipeThickness, 
                                         (lxs->BPipeR - lxs->BPipeThickness), fhlength-lxs->BPipeThickness);
  G4SubtractionSolid* ipbox_xy = new G4SubtractionSolid("solidDMBPipeXY", solidDMBPipeOuter, solidDMBPipeInner);

// Round hole for the income pipe
  G4Tubs *hin_b = new G4Tubs("solidSHoleIn", 0.0, lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeThickness, 0.0, 2.0*M_PI);
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0), 
                           G4ThreeVector(lxs->BPipeR+2.0*lxs->BPipeThickness, 0.0, -fhlength+lxs->BPipeThickness/2.0));
  G4SubtractionSolid* ipbox_xyh1 = new G4SubtractionSolid("solidDMBPipe1", ipbox_xy, hin_b, transform1);

// Rectangular hole for the outgoing pipe
  G4Box *out_cut = new G4Box("solidoutcut", 2.0*lxs->BPipeR+lxs->BPipeThickness, 
                                         (lxs->BPipeR - lxs->BPipeThickness), lxs->BPipeThickness);
  G4Transform3D transform2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0), 
                           G4ThreeVector(0.0, 0.0, fhlength-lxs->BPipeThickness/2.0));
  G4SubtractionSolid* solidDMBPipe = new G4SubtractionSolid("solidDMBPipe", ipbox_xyh1, out_cut, transform2);
  G4LogicalVolume *logicDMBPipe = new G4LogicalVolume(solidDMBPipe, beamPipeMaterial, "logicDMBPipe");
  
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicDMBPipe, "DMBPipe", 
                     logicDMBPipeContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(-(lxs->BPipeR+2.0*lxs->BPipeThickness)-lxs->DumpMagnetYpos, 0.0, 0.0), 
                     logicDMBPipeContainer, "DMBPipeContainer", fieldvol, false, 0, lxs->OverlapTest);

  fMagnetZPos = 2.0*fhlength;
  fMagnetSizeZ = lxs->DumpMagnetZpos;
  G4ThreeVector trm(lxs->DumpMagnetXpos, lxs->DumpMagnetYpos, lxs->DumpMagnetZpos);
//   mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 0, 0, lxs->OverlapTest);
  mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 
                              new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.5*M_PI), 0, lxs->OverlapTest);
  mag->CostructSupport(fLogicWorld, trm, "DumpMagnet", true);
  
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, -1.0), 0.5*M_PI), trm, fieldvol, 
                    "DumpMagnetField", fLogicWorld, false, 0, lxs->OverlapTest);

// Add local magnetic field
  G4String magid = "Brems";
  LxBField *ipfield = ComposeFieldObject(magid, trm);
  G4FieldManager* fieldMgr = new G4FieldManager(ipfield);
  ConfigureFieldManager(fieldMgr, ipfield);
 
  fieldvol->SetFieldManager(fieldMgr, true);
  G4AutoDelete::Register(ipfield);
  G4AutoDelete::Register(fieldMgr);

  // Add field to neighboring volumes
  if (ipfield->IsZLimited()) {
    std::vector<G4String> fieldvolname{"logicBeamPipeTMVac", "logicBeamSplitContainer"};
    AddFieldToLogVolumes(fieldMgr, fieldvolname);
  } else {G4cout << "Brems magnetic field is not limited in z, do not add neighboring volumes.\n";}

}



void DetectorConstruction::ConstructBeamDump()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* beamDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamDumpMaterial);
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* dumpInsertMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamDumpInsertMaterial);
  
  G4double dx_dump_pipe = lxs->BeamDumpFrontXpos - lxs->BeamDumpR * cos(lxs->BeamDumpAngle) - lxs->BPipeR;
  G4LogicalVolume *logicBeamDump = 0;
  if (dx_dump_pipe < 0.0) {
  // needs to cut beam dump
    G4Tubs *solidBeamDump0 = new G4Tubs("solidBeamDump0", 0.0, lxs->BeamDumpR, lxs->BeamDumpZ/2.0, 0.0, 2.0*M_PI);
    G4double lbpipecut = lxs->BeamDumpZ/cos(lxs->BeamDumpAngle) + 2.0*lxs->BPipeR*tan(lxs->BeamDumpAngle);
    G4Tubs *solidBeamPipeHole = new G4Tubs("solidBeamPipeHole", 0.0, lxs->BPipeR, 0.5*lbpipecut, 0.0, 2.0*M_PI);

//     G4Transform3D transform(G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), lxs->BeamDumpAngle),
//                             G4ThreeVector(lxs->BeamDumpXpos/cos(lxs->BeamDumpAngle), 0.0, 0.0));
    G4Transform3D transform(G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), lxs->BeamDumpAngle),
                            G4ThreeVector(0.0, lxs->BeamDumpXpos/cos(lxs->BeamDumpAngle), 0.0));

    G4SubtractionSolid* solidBeamDump = new G4SubtractionSolid("solidBeamDump", solidBeamDump0,
                                                               solidBeamPipeHole, transform);                       
    logicBeamDump = new G4LogicalVolume(solidBeamDump, beamDumpMaterial, "logicBeamDump");
//      G4String msgstr("Error: beam dump overlaps with beam pipe!\n");
//      G4Exception("DetectorConstruction::", "ConstructBeamDump()", FatalException, msgstr.c_str());
  } else {
    G4Tubs *solidBeamDump = new G4Tubs("solidBeamDump", 0.0, lxs->BeamDumpR, lxs->BeamDumpZ/2.0, 0.0, 2.0*M_PI);
    logicBeamDump = new G4LogicalVolume(solidBeamDump, beamDumpMaterial, "logicBeamDump");
  }

  G4Tubs *solidBeamDumpInsert = new G4Tubs("solidBeamDumpInsert", 0.0, lxs->BeamDumpInsertR,
                                           lxs->BeamDumpInsertZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamDumpInsert = new G4LogicalVolume(solidBeamDumpInsert, dumpInsertMaterial,
                                                             "logicBeamDumpInsert");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, (lxs->BeamDumpInsertZ-lxs->BeamDumpZ)/2.0),
                   logicBeamDumpInsert, "BeamDumpInsert", logicBeamDump, false, 0, lxs->OverlapTest);

//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), lxs->BeamDumpAngle),
//                     G4ThreeVector(-lxs->BeamDumpXpos, 0.0, lxs->BeamDumpZpos),
//                                logicBeamDump, "BeamDumpAssembly", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), lxs->BeamDumpAngle),
                    G4ThreeVector(0.0, -lxs->BeamDumpXpos, lxs->BeamDumpZpos),
                               logicBeamDump, "BeamDumpAssembly", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double lbpipe = (2.0*lxs->BeamDumpR * tan(lxs->BeamDumpAngle) + lxs->BeamDumpZ) * cos(lxs->BeamDumpAngle);
  G4Tubs *solidBeamPipeNextToDump = new G4Tubs("solidBeamPipeNextToDump", lxs->BPipeR-lxs->BPipeThickness,
                                               lxs->BPipeR, 0.5*lbpipe, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeNextToDump = new G4LogicalVolume(solidBeamPipeNextToDump, beamPipeMaterial,
                                                                 "logicBeamPipeNextToDump");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->BeamDumpZpos),
                      logicBeamPipeNextToDump, "BeamPipeNextToDump", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeNextToDumpVac = new G4Tubs("solidBeamPipeNextToDumpVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                                  0.5*lbpipe, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeNextToDumpVac = new G4LogicalVolume(solidBeamPipeNextToDumpVac,
                                                                    vacuumMaterial, "logicBeamPipeNextToDumpVac");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, lxs->BeamDumpZpos),
                      logicBeamPipeNextToDumpVac, "BeamPipeNextToDumpVac", fLogicWorld, false, 0, lxs->OverlapTest);

}



void DetectorConstruction::ConstructShielding()
{
  LXSetUp *lxs = LXSetUp::Instance();
    
  G4Material* shildingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ShieldingMaterial);
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* shildingAbsMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ShieldingAbsorberMaterial);
  G4Material* environmentMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);

  G4double ydeep = 2.0 * (lxs->BeamDumpR + lxs->ShieldingDeepMargine);
  G4double bdumpxproj = 2.0*lxs->BeamDumpR * cos(lxs->BeamDumpAngle) + lxs->BeamDumpZ * sin(lxs->BeamDumpAngle); 
  G4double xdeep, xdeeppos;
//   if (lxs->BeamDumpFrontXpos + lxs->BPipeR < 0.5*bdumpxproj) {
  if (lxs->BeamDumpFrontXpos + lxs->BPipeR < lxs->BeamDumpR * cos(lxs->BeamDumpAngle)) {
    xdeep = bdumpxproj + 2.0*lxs->ShieldingDeepMargine;
    xdeeppos = lxs->BeamDumpXpos;
  } else {
    xdeep = lxs->BeamDumpXpos + lxs->BPipeR + 0.5*bdumpxproj + 2.0*lxs->ShieldingDeepMargine;
    xdeeppos = 0.5 * xdeep - lxs->BPipeR - lxs->ShieldingDeepMargine;
  }
  G4double ddd = 0.5*(2.0*lxs->BeamDumpR * sin(lxs->BeamDumpAngle) + lxs->BeamDumpZ * cos(lxs->BeamDumpAngle)) 
                 + lxs->ShieldingDeepMargine; 
  G4double zdeeppos =  0.5 * (lxs->ShieldingDeepZ - lxs->ShieldingZ);

  G4double zshieldpos = lxs->BeamDumpZpos - lxs->ShieldingDeepZ + 0.5*lxs->ShieldingZ + ddd;
  lxs->ShieldingZpos = zshieldpos;
  
  G4double walldx = lxs->ShieldingX/2.0 - lxs->ShieldingDXBeamWall;
  G4double floordy = -(lxs->ShieldingY/2.0 + lxs->FloorSurfaceYpos);
  
  G4Box *solidShielding0 = new G4Box("solidShielding0", lxs->ShieldingX/2.0, lxs->ShieldingY/2.0, lxs->ShieldingZ/2.0);
  G4Tubs *solidShieldingH = new G4Tubs("solidShieldingH", 0.0, lxs->BPipeR, 0.55*lxs->ShieldingZ, 0.0, 2.0*M_PI);
//   G4Box *solidShieldingDeep = new G4Box("solidShieldingDeep", xdeep/2.0, ydeep/2.0, 
//                                         lxs->ShieldingDeepZ/2.0 + lxs->ShieldingDeepMargine);
  G4Box *solidShieldingDeep = new G4Box("solidShieldingDeep", ydeep/2.0, xdeep/2.0, 
                                        lxs->ShieldingDeepZ/2.0 + lxs->ShieldingDeepMargine);
  
  G4Transform3D transformh(G4RotationMatrix(), G4ThreeVector(walldx, floordy, 0.0));
  G4SubtractionSolid* solidShielding1 = new G4SubtractionSolid("solidShielding1", solidShielding0, 
                                                                solidShieldingH, transformh);
//   G4Transform3D transformd(G4RotationMatrix(), G4ThreeVector(-xdeeppos + walldx, 0.0, zdeeppos));
  G4Transform3D transformd(G4RotationMatrix(), G4ThreeVector(walldx, -xdeeppos+floordy, zdeeppos));
  G4SubtractionSolid* solidShielding = new G4SubtractionSolid("solidShielding", solidShielding1, 
                                                               solidShieldingDeep, transformd);
  
  G4LogicalVolume *logicShielding = new G4LogicalVolume(solidShielding, shildingMaterial, "logicShielding");
  new G4PVPlacement(0, G4ThreeVector(-walldx, -floordy, zshieldpos), 
                    logicShielding, "Shielding", fLogicWorld, false, 0, lxs->OverlapTest);

// Add front lead absorber bar
  G4double shabsorby = (floordy-xdeeppos) + lxs->ShieldingY/2.0 - xdeep/2.0;
//   shabsorby += lxs->ShieldingDeepMargine + lxs->BeamDumpR/2.0;
  G4Box *solidShieldAbsorber = new G4Box("solidShieldAbsorber", lxs->ShieldingAbsorberX/2.0,
                                         shabsorby/2.0, lxs->ShieldingAbsorberZ/2.0);
  G4LogicalVolume *logicShieldAbsorber = new G4LogicalVolume(solidShieldAbsorber, shildingAbsMaterial,
                                                             "logicShieldAbsorber");
  new G4PVPlacement(0, G4ThreeVector(walldx, 0.5*(shabsorby-lxs->ShieldingY),
                                     0.5*(lxs->ShieldingAbsorberZ-lxs->ShieldingZ)),
                    logicShieldAbsorber, "ShieldingAbsorber", logicShielding, false, 0, lxs->OverlapTest);
// Add air slit to the front absorber bar
  G4Box *solidShieldAbsSlit = new G4Box("solidShieldAbsSlit", lxs->ShieldingAbsorberSlitX/2.0,
                                         shabsorby/2.0, lxs->ShieldingAbsorberSlitZ/2.0);
  G4LogicalVolume *logicShieldAbsSlit = new G4LogicalVolume(solidShieldAbsSlit, environmentMaterial,
                                                             "logicShieldAbsSlit");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.5*(lxs->ShieldingAbsorberSlitZ - lxs->ShieldingAbsorberZ)),
                    logicShieldAbsSlit, "ShieldingAbsorberSlit", logicShieldAbsorber, false, 0, lxs->OverlapTest);

// Add lead absorber bar to cover the gap between shielding and dump
  G4double shabsorbtopy = lxs->ShieldingDeepMargine + lxs->BeamDumpR/2.0;
  G4Box *solidShieldAbsorberTop1 = new G4Box("solidShieldAbsorberTop1", lxs->ShieldingAbsorberX/2.0,
                                         shabsorbtopy/2.0, lxs->ShieldingAbsorberTopZ/2.0);
  G4Box *solidShieldAbsorberTopCut = new G4Box("solidShieldAbsorberTopCut", lxs->ShieldingAbsorberSlitX/2.0,
                                         shabsorbtopy, lxs->ShieldingAbsorberSlitZ/2.0);
  G4Transform3D trtopslitcut(G4RotationMatrix(),
                  G4ThreeVector(0.0, 0.0, 0.5*(lxs->ShieldingAbsorberSlitZ - lxs->ShieldingAbsorberTopZ)));
  G4SubtractionSolid* solidShieldAbsorberTop = new G4SubtractionSolid("solidShieldAbsorberTop",
                                 solidShieldAbsorberTop1, solidShieldAbsorberTopCut, trtopslitcut);
  G4LogicalVolume *logicShieldAbsorberTop = new G4LogicalVolume(solidShieldAbsorberTop, shildingAbsMaterial,
                                                             "logicShieldAbsorberTop");
  G4double shabstopposy = shabsorby + shabsorbtopy/2.0 + lxs->FloorSurfaceYpos;
  G4double shabstopposz = zshieldpos + 0.5*(lxs->ShieldingAbsorberTopZ - lxs->ShieldingZ);
  new G4PVPlacement(0, G4ThreeVector(0.0, shabstopposy, shabstopposz),
                    logicShieldAbsorberTop, "ShieldAbsorberTop", fLogicWorld, false, 0, lxs->OverlapTest);

// Add beam pipe through the shielding
  G4double lshildpipe = lxs->ShieldingZ - lxs->ShieldingDeepZ + lxs->ShieldingDeepMargine;
  G4double shildpipezpos = zshieldpos + 0.5 * (lxs->ShieldingZ - lshildpipe);
  
  G4Tubs *solidShieldingPipe = new G4Tubs("solidShieldingPipe", lxs->BPipeR - lxs->BPipeThickness, 
                                          lxs->BPipeR, lshildpipe/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicShieldingPipe = new G4LogicalVolume(solidShieldingPipe, beamPipeMaterial, "logicShieldingPipe");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, shildpipezpos), logicShieldingPipe, "ShieldingPipe", fLogicWorld, 
                    false, 0, lxs->OverlapTest);

  G4Tubs *solidShieldingPipeVac = new G4Tubs("solidShieldingPipeVac", 0.0, lxs->BPipeR - lxs->BPipeThickness, 
                                             lshildpipe/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicShieldingPipeVac = new G4LogicalVolume(solidShieldingPipeVac, vacuumMaterial,
                                                               "logicShieldingPipeVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, shildpipezpos), logicShieldingPipeVac, "ShieldingPipeVac", fLogicWorld, 
                    false, 0, lxs->OverlapTest);
    
}



void DetectorConstruction::ConstructBeamPipeToDump()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* beamPipeDumpMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamDumpPipeMaterial);

//   G4double magnetz = lxs->DumpMagnetZ;
//   G4double magnetz = lxs->TypMBFieldLength;
  G4double magnetz = lxs->FlashMFieldLength;
  
  G4Box *solidBeamSplitContainer = new G4Box("solidBeamSplitContainer", 2.0*(lxs->BPipeR+lxs->BPipeThickness),
                                             lxs->BPipeR, lxs->BeamPipeSplitterZ/2.0);
  //This is to have vacuum insede of the splitter
  G4LogicalVolume *logicBeamSplitContainer = new G4LogicalVolume(solidBeamSplitContainer, vacuumMaterial,
                                                                 "logicBeamSplitContainer");
  G4Box *solidBeamSplitOuter = new G4Box("solidBeamSplitOuter", 2.0*(lxs->BPipeR+lxs->BPipeThickness), lxs->BPipeR, 
                                                                lxs->BeamPipeSplitterZ/2.0);
  G4Box *solidBeamSplitInner = new G4Box("solidBeamSplitInner", 2.0*lxs->BPipeR+lxs->BPipeThickness, 
                                         (lxs->BPipeR - lxs->BPipeThickness), lxs->BeamPipeSplitterZ/2.0-lxs->BPipeThickness);
  G4SubtractionSolid* ipbox_xy = new G4SubtractionSolid("solidBeamSplitXY", solidBeamSplitOuter, solidBeamSplitInner);
  G4Tubs *hin_b = new G4Tubs("solidSHoleIn", 0.0, lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeThickness, 0.0, 2.0*M_PI);

// // This is round hole for the income pipe
//   G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
//          G4ThreeVector(lxs->BPipeR+2.0*lxs->BPipeThickness, 0.0, -(lxs->BeamPipeSplitterZ-lxs->BPipeThickness)/2.0));
//   G4SubtractionSolid* ipbox_xyh1 = new G4SubtractionSolid("solidBeamSplit1", ipbox_xy, hin_b, transform1);

// This is rectangular hole for the income pipe
  G4Box *in_cut = new G4Box("solidincut", 2.0*lxs->BPipeR+lxs->BPipeThickness,
                                         (lxs->BPipeR - lxs->BPipeThickness), lxs->BPipeThickness);
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
                           G4ThreeVector(0.0, 0.0, -(lxs->BeamPipeSplitterZ-lxs->BPipeThickness)/2.0));
  G4SubtractionSolid* ipbox_xyh1 = new G4SubtractionSolid("solidBeamSplit1", ipbox_xy, in_cut, transform1);

  G4Transform3D transform2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
                                            G4ThreeVector(lxs->BPipeR+2.0*lxs->BPipeThickness, 0.0,
                                                  (lxs->BeamPipeSplitterZ-lxs->BPipeThickness)/2.0));
  G4SubtractionSolid* ipbox_xyh2 = new G4SubtractionSolid("solidBeamSplit2", ipbox_xyh1, hin_b, transform2);
  G4Transform3D transform3(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
          G4ThreeVector(-(lxs->BPipeR+2.0*lxs->BPipeThickness), 0.0, (lxs->BeamPipeSplitterZ-lxs->BPipeThickness)/2.0));
  G4SubtractionSolid* solidBeamSplit = new G4SubtractionSolid("solidBeamSplit", ipbox_xyh2, hin_b, transform3);
  G4LogicalVolume *logicBeamSplit = new G4LogicalVolume(solidBeamSplit, beamPipeMaterial, "logicBeamSplit");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                               logicBeamSplit, "BeamSplit", logicBeamSplitContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(-(lxs->BPipeR+2.0*lxs->BPipeThickness), 0.0,
//                                      lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ/2.0),
//                     logicBeamSplitContainer, "BeamSplitContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.5*M_PI),
                    G4ThreeVector(0.0, -(lxs->BPipeR+2.0*lxs->BPipeThickness),
                                     lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ/2.0),
                          logicBeamSplitContainer, "BeamSplitContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double dxdump = lxs->BeamDumpXpos - 0.5*lxs->BeamDumpZ * tan(lxs->BeamDumpAngle);
  G4double dxsplit = 2.0*(lxs->BPipeR+2.0*lxs->BPipeThickness);
  G4double lpipe = lxs->BeamDumpZpos - lxs->DumpMagnetZpos - 0.5*magnetz 
        - 0.5*lxs->BeamDumpZ*cos(lxs->BeamDumpAngle) - lxs->BeamDumpR*sin(lxs->BeamDumpAngle) - lxs->BeamPipeSplitterZ;

  G4Tubs *solidStrightPipe = new G4Tubs("solidStrightPipe", lxs->BPipeR-lxs->BPipeThickness,
                                        lxs->BPipeR, lpipe/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeMB = new G4LogicalVolume(solidStrightPipe, beamPipeMaterial, "logicBeamPipeMB");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + lpipe/2.0),
                               logicBeamPipeMB, "BeamPipeMB", fLogicWorld, false, 0, lxs->OverlapTest);
  
  G4Tubs *solidStrightPipeVac = new G4Tubs("solidStrightPipeVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                           lpipe/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeMBVac = new G4LogicalVolume(solidStrightPipeVac, vacuumMaterial, "logicBeamPipeMBVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + lpipe/2.0),
                               logicBeamPipeMBVac, "BeamPipeMBVac", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double zpipedump =  lxs->BeamDumpZpos - lxs->DumpMagnetZpos - 0.5*magnetz
                      - 0.5*lxs->BeamDumpZ*cos(lxs->BeamDumpAngle) - lxs->BeamPipeSplitterZ;
  G4double lpipedump = sqrt(pow(zpipedump, 2.0) + pow(dxdump-dxsplit, 2.0));
  G4double dumpipeangle = atan2(dxdump-dxsplit, zpipedump);

//   G4CutTubs *solidPipeToDump = new G4CutTubs("solidPipeToDump", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                          lpipedump/2.0, 0.0, 2.0*M_PI,
//                                          G4ThreeVector(-(dxdump-dxsplit)/lpipedump, 0.0, -zpipedump/lpipedump),
//          G4ThreeVector(sin(dumpipeangle - lxs->BeamDumpAngle), 0.0, cos(dumpipeangle - lxs->BeamDumpAngle)));
  G4CutTubs *solidPipeToDump = new G4CutTubs("solidPipeToDump", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                         lpipedump/2.0, 0.0, 2.0*M_PI,
                                         G4ThreeVector(0.0, -(dxdump-dxsplit)/lpipedump, -zpipedump/lpipedump),
          G4ThreeVector(0.0, sin(dumpipeangle - lxs->BeamDumpAngle), cos(dumpipeangle - lxs->BeamDumpAngle)));

  G4LogicalVolume *logicBeamPipeMD = new G4LogicalVolume(solidPipeToDump, beamPipeDumpMaterial, "logicBeamPipeMD");
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), atan((dxdump-dxsplit)/zpipedump)), 
//                     G4ThreeVector(-(dxdump+dxsplit)/2.0, 0.0, 
//                     lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + zpipedump/2.0), 
//                                logicBeamPipeMD, "BeamPipeMD", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), atan((dxdump-dxsplit)/zpipedump)), 
                    G4ThreeVector(0.0, -(dxdump+dxsplit)/2.0, 
                    lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + zpipedump/2.0), 
                               logicBeamPipeMD, "BeamPipeMD", fLogicWorld, false, 0, lxs->OverlapTest);

//   G4CutTubs *solidPipeToDumpVac = new G4CutTubs("solidPipeToDumpVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, 
//                                          lpipedump/2.0, 0.0, 2.0*M_PI,
//                                          G4ThreeVector(-(dxdump-dxsplit)/lpipedump, 0.0, -zpipedump/lpipedump),
//                                          G4ThreeVector(sin(dumpipeangle - lxs->BeamDumpAngle), 0.0, cos(dumpipeangle - lxs->BeamDumpAngle)));
  G4CutTubs *solidPipeToDumpVac = new G4CutTubs("solidPipeToDumpVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, 
                                         lpipedump/2.0, 0.0, 2.0*M_PI,
                                         G4ThreeVector(0.0, -(dxdump-dxsplit)/lpipedump, -zpipedump/lpipedump),
                        G4ThreeVector(0.0, sin(dumpipeangle - lxs->BeamDumpAngle), cos(dumpipeangle - lxs->BeamDumpAngle)));
  
  G4LogicalVolume *logicBeamPipeMDVac = new G4LogicalVolume(solidPipeToDumpVac, vacuumMaterial, "logicBeamPipeMDVac");
//   new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 1.0, 0.0), atan((dxdump-dxsplit)/zpipedump)), 
//                     G4ThreeVector(-(dxdump+dxsplit)/2.0, 0.0, 
//                     lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + zpipedump/2.0), 
//                                logicBeamPipeMDVac, "BeamPipeMDVac", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(-1.0, 0.0, 0.0), atan((dxdump-dxsplit)/zpipedump)), 
                    G4ThreeVector(0.0, -(dxdump+dxsplit)/2.0, 
                    lxs->DumpMagnetZpos + 0.5*magnetz + lxs->BeamPipeSplitterZ + zpipedump/2.0), 
                               logicBeamPipeMDVac, "BeamPipeMDVac", fLogicWorld, false, 0, lxs->OverlapTest);
}



void DetectorConstruction::ConstructBeamPipeToIP()
{
  LXSetUp *lxs = LXSetUp::Instance();

//   G4double magnetz = lxs->DumpMagnetZ;
//   G4double magnetz = lxs->TypMBFieldLength;
  G4double magnetz = lxs->FlashMFieldLength;

  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  G4double lpipe_sip = -lxs->ShieldingZpos - lxs->IPContainerZ/2.0 - lxs->ShieldingZ/2.0 + lxs->TAUIChamberZpos;
//   G4Tubs *solidBeamPipeSIP = new G4Tubs("solidBeamPipeSIP", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                         lpipe_sip/2.0, 0.0, 2.0*M_PI);
  G4Cons *solidBeamPipeSIP = new G4Cons("solidBeamPipeSIP", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR, lxs->TAUIChamberElPipeRIn,
                                        lxs->TAUIChamberElPipeROut, lpipe_sip/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeSIP = new G4LogicalVolume(solidBeamPipeSIP, beamPipeMaterial, "logicBeamPipeSIP");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(lpipe_sip + lxs->IPContainerZ)/2.0 + lxs->TAUIChamberZpos),
                               logicBeamPipeSIP, "BeamPipeSIP", fLogicWorld, false, 0, lxs->OverlapTest);

//   G4Tubs *solidBeamPipeSIPVac = new G4Tubs("solidBeamPipeSIPVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
//                                            lpipe_sip/2.0, 0.0, 2.0*M_PI);
  G4Cons *solidBeamPipeSIPVac = new G4Cons("solidBeamPipeSIPVac", 0.0, lxs->BPipeR-lxs->BPipeThickness, 0.0,
                                           lxs->TAUIChamberElPipeRIn, lpipe_sip/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeSIPVac = new G4LogicalVolume(solidBeamPipeSIPVac, vacuumMaterial, "logicBeamPipeSIPVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(lpipe_sip + lxs->IPContainerZ)/2.0 + lxs->TAUIChamberZpos),
                               logicBeamPipeSIPVac, "BeamPipeSIPVac", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double lpipe_ipm = lxs->IPMagnetZpos - magnetz/2.0 - lxs->IPContainerZ/2.0 - lxs->TAUIChamberZpos;
//   G4Tubs *solidBeamPipeIPM = new G4Tubs("solidBeamPipeIPM", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                         lpipe_ipm/2.0, 0.0, 2.0*M_PI);
  G4double ipmbphx = lxs->FlashMFieldX/2.0 - lxs->IPMAgnetBeamPipeXGap;
  G4Box *solidBeamPipeIPM1 = new G4Box("solidBeamPipeIPM1", ipmbphx, lxs->FlashMagneteffY/2.0,
                                                lpipe_ipm/2.0);
  G4Box *solidBeamPipeIPMCut = new G4Box("solidBeamPipeIPMCut", ipmbphx-lxs->BPipeThickness,
                             lxs->FlashMagneteffY/2.0 - lxs->BPipeThickness, lpipe_ipm);
  G4SubtractionSolid* solidBeamPipeIPM = new G4SubtractionSolid("solidBeamPipeIPM",
                                               solidBeamPipeIPM1, solidBeamPipeIPMCut);
  G4LogicalVolume *logicBeamPipeIPM = new G4LogicalVolume(solidBeamPipeIPM, beamPipeMaterial, "logicBeamPipeIPM");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (lpipe_ipm + lxs->IPContainerZ)/2.0 + lxs->TAUIChamberZpos),
                               logicBeamPipeIPM, "BeamPipeIPM", fLogicWorld, false, 0, lxs->OverlapTest);

//   G4Tubs *solidBeamPipeIPMVac = new G4Tubs("solidBeamPipeIPMVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
//                                             lpipe_ipm/2.0, 0.0, 2.0*M_PI);
  G4Box *solidBeamPipeIPMVac = new G4Box("solidBeamPipeIPMVac", ipmbphx-lxs->BPipeThickness,
                             lxs->FlashMagneteffY/2.0 - lxs->BPipeThickness, lpipe_ipm/2.0);
  G4LogicalVolume *logicBeamPipeIPMVac = new G4LogicalVolume(solidBeamPipeIPMVac, vacuumMaterial, "logicBeamPipeIPMVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (lpipe_ipm + lxs->IPContainerZ)/2.0 + lxs->TAUIChamberZpos),
                               logicBeamPipeIPMVac, "BeamPipeIPMVac", fLogicWorld, false, 0, lxs->OverlapTest);

}



void DetectorConstruction::ConstructVacuumChamber()
{
  LxVacChamberAssembly *vc = 0;
  LXSetUp *lxs = LXSetUp::Instance();
  
//   G4double dumpMagnetZ = lxs->TypMBFieldLength;
//   G4double gaptomag = 0.5 * (lxs->TypMBMagnetCoreZ - lxs->TypMBFieldLength);
  G4double dumpMagnetZ = lxs->FlashMFieldLength;
  G4double gaptomag = 0.5 * (lxs->FlashMagnetCoilZ - lxs->FlashMagnetZ);
  G4double bpipel = lxs->OPPPDetZtoMagnet - lxs->VacChambertoOPPPDetZGap - gaptomag;
  
//   G4String vctype("PipeVacuumChamber");
//   vc = new PipeChamberAssembly(vctype);
  G4String vctype("DetectorVacuumChamber");
  vc = new DetChamberAssembly(vctype);
   
  G4ThreeVector trm(0.0, 0.0, lxs->IPMagnetZpos + dumpMagnetZ/2.0 + bpipel/2.0 + gaptomag);
  vc->GetAssembly()->MakeImprint(fLogicWorld, trm, 0, 0, lxs->OverlapTest);
}



void DetectorConstruction::ConstructBremsTarget()
{
  LXSetUp *lxs = LXSetUp::Instance();

  if (!fTargetChamber) fTargetChamber = new LxTargetChamber();
  G4LogicalVolume *logicBremsTargetChamber = fTargetChamber->GetTargetChamber();

//   G4double dumpMagnetZ = lxs->DumpMagnetZ;
//   G4double dumpMagnetZ = lxs->TypMBMagnetZ;
  G4double dumpMagnetZ = lxs->FlashMagnetCoilZ;

  fZposAbs = lxs->DumpMagnetZpos - 0.5*(dumpMagnetZ + lxs->TargetChamberZ) - lxs->BTargetChamberMagnetGapZ;
  lxs->BTargetZpos = fZposAbs;
  fZstartAbs = fZposAbs-0.5*fAbsorberThickness;
  fZendAbs   = fZposAbs+0.5*fAbsorberThickness;

  new G4PVPlacement(0, G4ThreeVector(lxs->TargetChamberXPos, 0.0, lxs->BTargetZpos),
                       logicBremsTargetChamber, "BremsTargetChamber", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  fAbsorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BTargetMaterial);
  G4LogicalVolume *logicTargetVolume = fTargetChamber->GetTargetVolume("BTarget");
  const G4Box  *tvolbox = dynamic_cast<G4Box*>(logicTargetVolume->GetSolid());
  G4double tvolhx = tvolbox->GetXHalfLength();
  G4double tvolhy = tvolbox->GetYHalfLength();

  G4Box* solidBremsTarget = new G4Box("solidBremsTarget", tvolhx, tvolhy, lxs->BTargetZ/2.0);
  G4LogicalVolume* logicBremsTarget = new G4LogicalVolume(solidBremsTarget, fAbsorberMaterial, "logicBremsTarget");
  fPhysiAbsorber = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicBremsTarget, "BremsTarget",
                                                                      logicTargetVolume, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->BTargetZpos),
                       logicTargetVolume, "BremsTargetContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  fTargetChamber->ConstructTargetChamberSupport(fLogicWorld, "BremsTargetChamber", lxs->BTargetZpos);

  fSolidAbsorber = solidBremsTarget;
  fLogicAbsorber = logicBremsTarget;

// beam pipes
//   G4double dumpMagnetfZ = lxs->TypMBFieldLength;
  G4double dumpMagnetfZ = lxs->FlashMFieldLength;

  G4double lpipe_tm = lxs->DumpMagnetZpos - lxs->BTargetZpos -  0.5*(lxs->TargetChamberZ + dumpMagnetfZ);
  G4Tubs *solidBeamPipeTM = new G4Tubs("solidBeamPipeTM", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                       lpipe_tm/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeTM = new G4LogicalVolume(solidBeamPipeTM, beamPipeMaterial, "logicBeamPipeTM");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->BTargetZpos + (lxs->TargetChamberZ + lpipe_tm)/2.0),
                               logicBeamPipeTM, "BeamPipeTM", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeTMVac = new G4Tubs("solidBeamPipeTMVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                          lpipe_tm/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeTMVac = new G4LogicalVolume(solidBeamPipeTMVac, vacuumMaterial, "logicBeamPipeTMVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->BTargetZpos + (lxs->TargetChamberZ + lpipe_tm)/2.0),
                               logicBeamPipeTMVac, "BeamPipeTMVac", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double lpipe_inc = lxs->BTargetZpos + fWorldSizeZ/2.0 - lxs->TargetChamberZ/2.0;
  G4Tubs *solidBeamPipeInc = new G4Tubs("solidBeamPipeInc", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                        lpipe_inc/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeInc = new G4LogicalVolume(solidBeamPipeInc, beamPipeMaterial, "logicBeamPipeInc");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->BTargetZpos - (lxs->TargetChamberZ + lpipe_inc)/2.0),
                               logicBeamPipeInc, "BeamPipeInc", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeIncVac = new G4Tubs("solidBeamPipeIncVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                           lpipe_inc/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeIncVac = new G4LogicalVolume(solidBeamPipeIncVac, vacuumMaterial, "logicBeamPipeIncVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->BTargetZpos - (lxs->TargetChamberZ + lpipe_inc)/2.0),
                               logicBeamPipeIncVac, "BeamPipeIncVac", fLogicWorld, false, 0, lxs->OverlapTest);

}



void DetectorConstruction::ConstructIPMagnet(const G4String magType)
{
  LxMagnetAssembly *mag = 0;
  LXSetUp *lxs = LXSetUp::Instance();

  G4String mtype(magType);
//   G4String mtype("TypMBMagnetSimple");
  auto magnetitr = fMagentsTypeCollection.find(mtype);
  if (magnetitr == fMagentsTypeCollection.end()) {
    mag = new FlashMagnetAssembly(mtype);
//     mag = new TypMBMagnetAssembly(mtype);
//     mag = new TypMBMagnetSimple(mtype);
    fMagentsTypeCollection[mtype] = mag;
  } else {
    mag = magnetitr->second;
  }

  G4LogicalVolume *fieldvol = mag->GetFieldVolume("logicIPMagnetField");

  const G4Box  *mfbox = dynamic_cast<G4Box*>(fieldvol->GetSolid());
  G4double fhlength = mfbox->GetZHalfLength();
  G4double fhwidth = mfbox->GetXHalfLength();
  G4double fhhight = mfbox->GetYHalfLength();
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

// Beam pipe rectangular
//   G4double ipmbphx = lxs->TypMBMagnetHoleX/2.0 - lxs->IPMAgnetBeamPipeXGap;
  G4double ipmbphx = fhwidth - lxs->IPMAgnetBeamPipeXGap;

// Container
  G4Box *solidIPMBPipeContainer = new G4Box("solidIPMBPipeContainer", ipmbphx, fhhight, fhlength);
  G4LogicalVolume *logicIPMBPipeContainer = new G4LogicalVolume(solidIPMBPipeContainer,
                                                                 vacuumMaterial, "logicIPMBPipeContainer");

  G4Box *soliIPDMBPipeOuter = new G4Box("soliIPDMBPipeOuter", ipmbphx, fhhight, fhlength);
  G4Box *solidIPMBPipeInner = new G4Box("solidIPMBPipeInner", ipmbphx-lxs->BPipeThickness,
                                         fhhight - lxs->BPipeThickness, fhlength-lxs->BPipeThickness);
  G4SubtractionSolid* ipbox_xy = new G4SubtractionSolid("solidIPMBPipeXY", soliIPDMBPipeOuter, solidIPMBPipeInner);

// Round hole for the income pipe
//   G4Tubs *hin_b = new G4Tubs("solidIPMSHoleIn", 0.0, lxs->BPipeR-lxs->BPipeThickness,
//                              lxs->BPipeThickness, 0.0, 2.0*M_PI);
// Rectangular hole for the income pipe
  G4Box *hin_b = new G4Box("solidIPMSHoleIn", ipmbphx-lxs->BPipeThickness,
                                         fhhight - lxs->BPipeThickness, lxs->BPipeThickness);
  G4Transform3D transform1(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
                           G4ThreeVector(0.0, 0.0, -fhlength+lxs->BPipeThickness/2.0));
  G4SubtractionSolid* ipbox_xyh1 = new G4SubtractionSolid("solidIPMBPipe1", ipbox_xy, hin_b, transform1);

// Rectangular hole for the outgoing pipe
  G4Box *out_cut = new G4Box("solidipmoutcut", ipmbphx-lxs->BPipeThickness,
                                         fhhight - lxs->BPipeThickness, lxs->BPipeThickness);
  G4Transform3D transform2(G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.0),
                           G4ThreeVector(0.0, 0.0, fhlength-lxs->BPipeThickness/2.0));
  G4SubtractionSolid* solidIPMBPipe = new G4SubtractionSolid("solidIPMBPipe", ipbox_xyh1, out_cut, transform2);
  G4LogicalVolume *logicIPMBPipe = new G4LogicalVolume(solidIPMBPipe, beamPipeMaterial, "logicIPMBPipe");

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicIPMBPipe, "IPMBPipe",
                     logicIPMBPipeContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0),
                     logicIPMBPipeContainer, "IPMBPipeContainer", fieldvol, false, 0, lxs->OverlapTest);

  G4ThreeVector trm(lxs->IPMagnetXpos, lxs->IPMagnetYpos, lxs->IPMagnetZpos);
  mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 0, 0, lxs->OverlapTest);
  mag->CostructSupport(fLogicWorld, trm, "IPMagnet", false);
  new G4PVPlacement(0, trm, fieldvol, "IPMagnetField", fLogicWorld, false, 0, lxs->OverlapTest);

// Add local magnetic field
  G4String magid = "IP";
  LxBField *ipfield = ComposeFieldObject(magid, trm);
  G4FieldManager* fieldMgr = new G4FieldManager(ipfield);
  ConfigureFieldManager(fieldMgr, ipfield);
 
  fieldvol->SetFieldManager(fieldMgr, true);
  G4AutoDelete::Register(ipfield);
  G4AutoDelete::Register(fieldMgr);

  // Add field to neighboring volumes
  if (ipfield->IsZLimited()) {
    std::vector<G4String> fieldvolname{"logicVCMagFieldJoin", "logicVCMagFieldJoin1",
                                       "logicVCContainer", "logicBeamPipeIPMVac"};
    AddFieldToLogVolumes(fieldMgr, fieldvolname);
  } else {G4cout << "IP magnetic field is not limited in z, do not add neighboring volumes.\n";}
}



void DetectorConstruction::ConstructCleaningMagnet()
{
  LXSetUp *lxs = LXSetUp::Instance();

  G4Material* magnetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->MagnetMaterial);
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* bpipeWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeWindowMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

  G4Box *solidCleaningMContainer = new G4Box("solidCleaningMContainer", (lxs->DumpCleanMagX/2.0 + lxs->DumpCleanMagThickness), 
                                              (lxs->DumpCleanMagY/2.0 + lxs->DumpCleanMagThickness), lxs->DumpCleanMagZ/2.0);
  G4LogicalVolume *logicCleaningMContainer = new G4LogicalVolume(solidCleaningMContainer,
                                                                 fWorldMaterial, "logicCleaningMContainer");

  G4Box *solidCleaningMBoxOuter = new G4Box("solidCleaningMBoxOuter", (lxs->DumpCleanMagX/2.0 + lxs->DumpCleanMagThickness), 
                                               (lxs->DumpCleanMagY/2.0 + lxs->DumpCleanMagThickness), lxs->DumpCleanMagZ/2.0);
  G4Box *solidCleaningMBoxInner = new G4Box("solidCleaningMBoxInner", lxs->DumpCleanMagX/2.0, lxs->DumpCleanMagY/2.0, 1.1*lxs->DumpCleanMagZ/2.0);
  G4SubtractionSolid* ipbox_xy = new G4SubtractionSolid("solidCleaningMBoxXY", solidCleaningMBoxOuter, solidCleaningMBoxInner);                       
  G4LogicalVolume *logicCleaningMBoxXY = new G4LogicalVolume(ipbox_xy, magnetMaterial, "logicCleaningMBoxXY");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicCleaningMBoxXY, "CleaningMBoxXY", logicCleaningMContainer, false, 0, lxs->OverlapTest);

//   ///////// Beam pipe windos are in horizontal plane  
//   // Beam pipe top and bottom part 
//   G4Tubs *solidCleaningMBPipe = new G4Tubs("solidCleaningMBPipe", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                     lxs->DumpCleanMagZ/2.0, lxs->DMBPipeWindow/2.0, M_PI-lxs->DMBPipeWindow);
//   G4LogicalVolume *logicCleaningMBPipe = new G4LogicalVolume(solidCleaningMBPipe, beamPipeMaterial, "logicCleaningMBPipe");
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicCleaningMBPipe, "CleaningMBPipe", logicCleaningMContainer, false, 0, lxs->OverlapTest); 
//   new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0), 
//                       logicCleaningMBPipe, "CleaningMBPipe", logicCleaningMContainer, false, 1, lxs->OverlapTest); 
//   
//   // Beam pipe windows left and right sides 
//   G4Tubs *solidCleaningMBPipeWindow = new G4Tubs("solidCleaningMBPipeWindow", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
//                                           lxs->DumpCleanMagZ/2.0, -lxs->DMBPipeWindow/2.0, lxs->DMBPipeWindow);
//   G4LogicalVolume *logicCleaningMBPipeWindow = new G4LogicalVolume(solidCleaningMBPipeWindow, bpipeWindowMaterial, "logicCleaningMBPipeWindow");
//   new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicCleaningMBPipeWindow, "CleaningMBPipeWindow", logicCleaningMContainer, 
//                                                false, 0, lxs->OverlapTest); 
//   new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI), G4ThreeVector(0.0, 0.0, 0.0), 
//                       logicCleaningMBPipeWindow, "CleaningMBPipeWindow", logicCleaningMContainer, false, 1, lxs->OverlapTest); 

  ///////// Beam pipe windos are in vertical plane
  // Beam pipe top and bottom part 
  G4Tubs *solidCleaningMBPipe = new G4Tubs("solidCleaningMBPipe", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                    lxs->DumpCleanMagZ/2.0, lxs->DMBPipeWindow/2.0, M_PI-lxs->DMBPipeWindow);
  G4LogicalVolume *logicCleaningMBPipe = new G4LogicalVolume(solidCleaningMBPipe, beamPipeMaterial, "logicCleaningMBPipe");
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.5*M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicCleaningMBPipe, "CleaningMBPipe", logicCleaningMContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0),  1.5*M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicCleaningMBPipe, "CleaningMBPipe", logicCleaningMContainer, false, 1, lxs->OverlapTest);

  // Beam pipe windows left and right sides
  G4Tubs *solidCleaningMBPipeWindow = new G4Tubs("solidCleaningMBPipeWindow", lxs->BPipeR-lxs->BPipeThickness, lxs->BPipeR,
                                          lxs->DumpCleanMagZ/2.0, -lxs->DMBPipeWindow/2.0, lxs->DMBPipeWindow);
  G4LogicalVolume *logicCleaningMBPipeWindow = new G4LogicalVolume(solidCleaningMBPipeWindow,
                                                                   bpipeWindowMaterial, "logicCleaningMBPipeWindow");
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 0.5*M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicCleaningMBPipeWindow, "CleaningMBPipeWindow", logicCleaningMContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement (new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), 1.5*M_PI), G4ThreeVector(0.0, 0.0, 0.0),
                      logicCleaningMBPipeWindow, "CleaningMBPipeWindow", logicCleaningMContainer, false, 1, lxs->OverlapTest);

  G4Tubs *solidCleaningMBPipeVac = new G4Tubs("solidCleaningMBPipeVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                              lxs->DumpCleanMagZ/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicCleaningMBPipeVac = new G4LogicalVolume(solidCleaningMBPipeVac,
                                                                vacuumMaterial, "logicCleaningMBPipeVac");
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, 0.0), logicCleaningMBPipeVac, "CleaningMBPipeVac",
                     logicCleaningMContainer, false, 0, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->DumpCleanMagZpoz),
                               logicCleaningMContainer, "DumpCleaningMagnet", fLogicWorld, false, 0, lxs->OverlapTest);

  ConstructMagnet(logicCleaningMContainer, G4ThreeVector(lxs->DumpCleanMagField, 0.0, 0.0));
}




void DetectorConstruction::ConstructGammaTarget()
{
  LXSetUp *lxs = LXSetUp::Instance();

  if (!fTargetChamber) fTargetChamber = new LxTargetChamber();
  G4LogicalVolume *logicGammaTargetChamber = fTargetChamber->GetTargetChamber();
  new G4PVPlacement(0, G4ThreeVector(lxs->TargetChamberXPos, 0.0, lxs->GTargetZpos),
                       logicGammaTargetChamber, "GammaTargetChamber", fLogicWorld, false, 1, lxs->OverlapTest);

  G4Material *gammaAbsorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->GTargetMaterial);

  G4LogicalVolume *logicTargetVolume = fTargetChamber->GetTargetVolume("GammaTarget");
  const G4Box  *tvolbox = dynamic_cast<G4Box*>(logicTargetVolume->GetSolid());
  G4double tvolhx = tvolbox->GetXHalfLength();
  G4double tvolhy = tvolbox->GetYHalfLength();

  if (lxs->GTargetType == LXSetUp::tTargetType::twire) {
    fGSolidAbsorber = new G4Tubs("GsolidAbsorber", 0, lxs->GTargetZ/2.0, tvolhy, 0, 2.0*M_PI);

    fGLogicAbsorber = new G4LogicalVolume(fGSolidAbsorber, gammaAbsorberMaterial,  "GlogicAbsorber");
    fGPhysiAbsorber = new G4PVPlacement(new G4RotationMatrix(0.0, M_PI/2.0, 0.0), G4ThreeVector(0.0, 0.0,  0.0),
                   fGLogicAbsorber,  "GAbsorber", logicTargetVolume, false, 0, lxs->OverlapTest);
  } else {
    fGSolidAbsorber = new G4Box("GsolidAbsorber",  tvolhx, tvolhy, lxs->GTargetZ/2.0);
    fGLogicAbsorber = new G4LogicalVolume(fGSolidAbsorber,  gammaAbsorberMaterial,  "GlogicAbsorber");

    fGPhysiAbsorber = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), fGLogicAbsorber,
                                "GAbsorber", logicTargetVolume, false, 0, lxs->OverlapTest);
  }

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->GTargetZpos),
                       logicTargetVolume, "GammaTargetContainer", fLogicWorld, false, 0, lxs->OverlapTest);
  fTargetChamber->ConstructTargetChamberSupport(fLogicWorld, "GammaTargetChamber", lxs->GTargetZpos);
// Beam pipes
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);

//   G4double dumpMagnetZ = lxs->TypMBFieldLength;
  G4double dumpMagnetZ = lxs->FlashMFieldLength;
  G4double opppdetzpos = lxs->IPMagnetZpos + dumpMagnetZ/2.0 + lxs->OPPPDetZtoMagnet;
  G4double lpipe_dt = lxs->GTargetZpos - opppdetzpos -  0.5*lxs->TargetChamberZ;
  G4Tubs *solidBeamPipeOPPPDGT = new G4Tubs("solidBeamPipeOPPPDGT", lxs->BPipeR-lxs->BPipeThickness,
                                            lxs->BPipeR, lpipe_dt/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeOPPPDGT = new G4LogicalVolume(solidBeamPipeOPPPDGT, beamPipeMaterial,
                                                              "logicBeamPipeOPPPDGT");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->GTargetZpos - 0.5 * (lxs->TargetChamberZ + lpipe_dt)),
                    logicBeamPipeOPPPDGT, "BeamPipeOPPPDGT", fLogicWorld, false, 0, lxs->OverlapTest);

  G4Tubs *solidBeamPipeOPPPDGTVac = new G4Tubs("solidBeamPipeOPPPDGTVac", 0.0, lxs->BPipeR-lxs->BPipeThickness,
                                               lpipe_dt/2.0, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeOPPPDGTVac = new G4LogicalVolume(solidBeamPipeOPPPDGTVac, vacuumMaterial,
                                                                 "logicBeamPipeOPPPDGTVac");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lxs->GTargetZpos - 0.5 * (lxs->TargetChamberZ + lpipe_dt)),
                    logicBeamPipeOPPPDGTVac, "BeamPipeOPPPDGTVac", fLogicWorld, false, 0, lxs->OverlapTest);
}



void DetectorConstruction::ConstructGammaMagnet()
{
  LxMagnetAssembly *mag = 0;
  LXSetUp *lxs = LXSetUp::Instance();

  G4String mtype("FlashMagnet");
  auto magnetitr = fMagentsTypeCollection.find(mtype);
  if (magnetitr == fMagentsTypeCollection.end()) {
    mag = new FlashMagnetAssembly(mtype);
//     mag = new PDSMagnetAssembly(mtype);
    fMagentsTypeCollection[mtype] = mag;
  } else {
    mag = magnetitr->second;
  }

   G4LogicalVolume *fieldvol = mag->GetFieldVolume("logicGMagnetField");
//   ConstructMagnet(fieldvol, G4ThreeVector(0.0, lxs->GMagFieldY, 0.0));

// Beam pipe through the field
  G4Material* beamPipeMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeMaterial);
  G4Material* vacuumMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeVacuumMaterial);
  G4Material* bpipeWindowMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeWindowMaterial);

  G4double gmpipex = lxs->QBeamPipeContainerX;  // size of the pipe in X direction
  G4double gmpiptby = lxs->GChamberTBWallThickness; // top/bottom walls thickness
  G4double gmpipsidex = lxs->GammaBPipeWindowThickness; // side walls thickness, beam pipe window thickness

  const G4Box  *mfbox = dynamic_cast<G4Box*>(fieldvol->GetSolid());
  G4double gfieldz = mfbox->GetZHalfLength();
  G4double gfieldy = mfbox->GetYHalfLength();
  G4double gfieldx = mfbox->GetXHalfLength();

  G4Box *solidGammaMagnetPipeContainer = new G4Box("solidGammaMagnetPipeContainer", gmpipex/2.0,
                                                   gfieldy, gfieldz);
  G4LogicalVolume *logicGammaMagnetPipeContainer = new G4LogicalVolume(solidGammaMagnetPipeContainer,
                                                   vacuumMaterial, "logicGammaMagnetPipeContainer");
//Top, bottom
  G4Box *solidGamMagPipeTopBot = new G4Box("solidGamMagPipeTopBot", gmpipex/2.0, gmpiptby/2.0, gfieldz);
  G4LogicalVolume *logicGamMagPipeTopBot = new G4LogicalVolume(solidGamMagPipeTopBot,
                                                   beamPipeMaterial, "logicGamMagPipeTopBot");
 //Sides
  G4double swy = gfieldy - gmpiptby;  // hight of the pipe or chamber
  G4Box *solidGamMagPipeSide = new G4Box("solidGamMagPipeSide", gmpipsidex/2.0, swy, gfieldz);
  G4LogicalVolume *logicGamMagPipeSide = new G4LogicalVolume(solidGamMagPipeSide,
                                                   beamPipeMaterial, "logicGamMagPipeSide");
//Front window
  G4Box *solidGamMagPipeFrontWind = new G4Box("solidGamMagPipeFrontWind", gmpipex/2.0 , swy, gmpipsidex/2.0);
  G4LogicalVolume *logicGamMagPipeFrontWind = new G4LogicalVolume(solidGamMagPipeFrontWind,
                                                   bpipeWindowMaterial, "logicGamMagPipeFrontWind");
//Rear part
  G4Box *solidGamMagPipeRearBox = new G4Box("solidGamMagPipeRearBox", gmpipex/2.0 , swy, gmpiptby/2.0);
  G4double incutr = lxs->BPipeRLG-lxs->BPipeThickness;
  if (incutr >= swy) incutr = 0.96*swy;
  G4Tubs *solidGamMagPipeRearBeamPipeCut = new G4Tubs("solidGamMagPipeRearBeamPipeCut", 0.0,
                                    incutr, gmpiptby, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidGamMagPipeRear = new G4SubtractionSolid("solidGamMagPipeRear",
                                                 solidGamMagPipeRearBox, solidGamMagPipeRearBeamPipeCut);
  G4LogicalVolume *logicGamMagPipeRear = new G4LogicalVolume(solidGamMagPipeRear,
                                                   beamPipeMaterial, "logicGamMagPipeRear");

  new G4PVPlacement(0, G4ThreeVector(0.0, gfieldy-gmpiptby/2.0, 0.0), logicGamMagPipeTopBot,
                    "GMagnetFieldPipeTop", logicGammaMagnetPipeContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, -(gfieldy-gmpiptby/2.0), 0.0), logicGamMagPipeTopBot,
                    "GMagnetFieldPipeBottom", logicGammaMagnetPipeContainer, false, 1, lxs->OverlapTest);

//   new G4PVPlacement(0, G4ThreeVector((gmpipex-gmpipsidex)/2.0, 0.0, 0.0), logicGamMagPipeSide,
//                     "GMagnetFieldPipeSide", logicGammaMagnetPipeContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(-(gmpipex-gmpipsidex)/2.0, 0.0, 0.0), logicGamMagPipeSide,
//                     "GMagnetFieldPipeSide", logicGammaMagnetPipeContainer, false, 1, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, gfieldz-gmpipsidex/2.0), logicGamMagPipeFrontWind,
//                     "GamMagPipeFrontWind", logicGammaMagnetPipeContainer, false, 1, lxs->OverlapTest);
// Added chamber in front
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(gfieldz-gmpiptby/2.0)), logicGamMagPipeRear,
                    "GamMagPipeRear", logicGammaMagnetPipeContainer, false, 1, lxs->OverlapTest);                                   

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logicGammaMagnetPipeContainer,
                    "GMagnetFieldBeamPipe", fieldvol, false, 0, lxs->OverlapTest);

//To  make wide beam pipe inside the magnet
//   WideBeamPipeContainerX
  G4double wpipecontainerx =  gfieldx - gmpipex/2.0 - lxs->GammaMagnetBeamPipeXGap;
  G4Box *solidGammaMagnetWidePipeContainer = new G4Box("solidGammaMagnetWidePipeContainer", wpipecontainerx/2.0,
                                                   gfieldy, gfieldz);
  G4LogicalVolume *logicGammaMagnetWidePipeContainer = new G4LogicalVolume(solidGammaMagnetWidePipeContainer,
                                                   vacuumMaterial, "logicGammaMagnetWidePipeContainer");
  G4Box *solidGamMagWidePipeTopBot = new G4Box("solidGamMagWidePipeTopBot", wpipecontainerx/2.0,
                                               gmpiptby/2.0, gfieldz);
  G4LogicalVolume *logicGamMagWidePipeTopBot = new G4LogicalVolume(solidGamMagWidePipeTopBot,
                                                   beamPipeMaterial, "logicGamMagWidePipeTopBot");
  G4Box *solidGamMagWidePipeSide = new G4Box("solidGamMagWidePipeSide", gmpipsidex/2.0, gfieldy - gmpiptby, gfieldz);
  G4LogicalVolume *logicGamMagWidePipeSide = new G4LogicalVolume(solidGamMagWidePipeSide,
                                                   beamPipeMaterial, "logicGamMagWidePipeSide");
  G4Box *solidGamMagWidePipeFront = new G4Box("solidGamMagWidePipeFront", wpipecontainerx/2.0 - gmpipsidex/2.0,
                                              gfieldy - gmpiptby, gmpipsidex/2.0);
  G4LogicalVolume *logicGamMagWidePipeFront = new G4LogicalVolume(solidGamMagWidePipeFront,
                                                   bpipeWindowMaterial, "logicGamMagWidePipeFront");
  G4Box *solidGamMagWidePipeRear = new G4Box("solidGamMagWidePipeRear", wpipecontainerx/2.0 - gmpipsidex/2.0,
                                             gfieldy - gmpiptby, gmpiptby/2.0);
  G4LogicalVolume *logicGamMagWidePipeRear = new G4LogicalVolume(solidGamMagWidePipeRear,
                                                   beamPipeMaterial, "logicGamMagWidePipeRear");

//Placements
//   new G4PVPlacement(0, G4ThreeVector((gmpipex-gmpipsidex)/2.0, 0.0, lxs->GMagnetZpos), logicGamMagPipeSide,
//                     "GMagnetWidePipeSide", logicGammaMagnetWidePipeContainer, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(-(lxs->GMagnetCutX)/2.0, 0.0, lxs->GMagnetZpos), logicGamMagPipeSide,
//                     "GMagnetWidePipeSide", logicGammaMagnetWidePipeContainer, false, 1, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, gfieldy-gmpiptby/2.0, 0.0), logicGamMagWidePipeTopBot,
                    "GMagnetFieldWidePipeTop", logicGammaMagnetWidePipeContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, -(gfieldy-gmpiptby/2.0), 0.0), logicGamMagWidePipeTopBot,
                    "GMagnetFieldWidePipeBottom", logicGammaMagnetWidePipeContainer, false, 1, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector((wpipecontainerx-gmpipsidex)/2.0, 0.0, 0.0), logicGamMagWidePipeSide,
                     "GMagnetFieldWidePipeSide", logicGammaMagnetWidePipeContainer, false, 0, lxs->OverlapTest);
//  new G4PVPlacement(0, G4ThreeVector(-gmpipsidex/2.0, 0.0, gfieldz-gmpipsidex/2.0), logicGamMagWidePipeFront,
//                    "GamMagWidePipeFront", logicGammaMagnetWidePipeContainer, false, 0, lxs->OverlapTest);
// Added chamber in front
  new G4PVPlacement(0, G4ThreeVector(-gmpipsidex/2.0, 0.0, -(gfieldz-gmpiptby/2.0)), logicGamMagWidePipeRear,
                    "GamMagWidePipeRear", logicGammaMagnetWidePipeContainer, false, 0, lxs->OverlapTest);

  new G4PVPlacement(0, G4ThreeVector((gmpipex + wpipecontainerx)/2.0, 0.0, 0.0),
               logicGammaMagnetWidePipeContainer, "GMagnetWideBeamPipeL", fieldvol, false, 0, lxs->OverlapTest);

  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI),
                    G4ThreeVector(-(gmpipex + wpipecontainerx)/2.0, 0.0, 0.0),
               logicGammaMagnetWidePipeContainer, "GMagnetWideBeamPipeR", fieldvol, false, 1, lxs->OverlapTest);

//To make vacuum chamber
  G4double vclength = lxs->ComptonLysoZpos - lxs->GMagnetZpos - lxs->ComptonLysoZ/2.0
                      - gfieldz - lxs->GamVacChamberGap;
  G4double magnetx = wpipecontainerx + gmpipex/2.0;
  G4double detx =  lxs->ComptonLysoX + gmpipsidex;

  G4Trd *solidGVCContainer = new G4Trd("solidGVCContainer", magnetx, detx, gfieldy, gfieldy, vclength/2.0);

  G4LogicalVolume *logicGVCContainer = new G4LogicalVolume(solidGVCContainer, vacuumMaterial, "logicGVCContainer");

 //Top and bottom wals of the chamber
  G4Trd *solidGamChamTopBot = new G4Trd("solidGamChamTopBot", magnetx, detx, gmpiptby/2.0, gmpiptby/2.0, vclength/2.0);
  G4LogicalVolume *logicGamChamTopBot = new G4LogicalVolume(solidGamChamTopBot,
                                                   beamPipeMaterial, "logicGamChamTopBot");
  new G4PVPlacement(0, G4ThreeVector(0.0, gfieldy-gmpiptby/2.0, 0.0), logicGamChamTopBot,
                    "GamChamTop", logicGVCContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, -(gfieldy-gmpiptby/2.0), 0.0), logicGamChamTopBot,
                    "GamChamBottom", logicGVCContainer, false, 1, lxs->OverlapTest);

 //Window of the chamber
  G4double wtheta = atan2((detx-magnetx), vclength);
  G4double vcwndz = gmpipsidex;
  G4double vcwallx = detx - vcwndz * tan(wtheta);
  G4double BPipeRProf = lxs->ComptonLysoXpos - lxs->ComptonLysoX/2.0; //1.0 *cm;
  G4Box *solidVCWindow1 = new G4Box("solidVCWindow1", vcwallx, swy, vcwndz/2.0);
  G4Tubs *solidCVBeamPipeHole = new G4Tubs("solidCVBeamPipeHole", 0.0, BPipeRProf-lxs->BPipeThickness, vcwndz, 0.0, 2.0*M_PI);
  G4SubtractionSolid* solidVCWindow = new G4SubtractionSolid("solidVCWindow", solidVCWindow1, solidCVBeamPipeHole);
  G4LogicalVolume *logicGVCWindow = new G4LogicalVolume(solidVCWindow, bpipeWindowMaterial, "logicGVCWindow");
  G4ThreeVector wndtrans(0.0, 0.0, (vclength-vcwndz)/2.0);
  new G4PVPlacement(0, wndtrans, logicGVCWindow, "GVCWindow", logicGVCContainer, false, 0, lxs->OverlapTest);

 //Side wall of the chamber
  G4double   swdx = 0.5*gmpipsidex / cos(wtheta);
  G4Para *solidGamChamSide = new G4Para("solidGamChamSide", swdx/2.0, swy, (vclength-vcwndz)/2.0, 0.0, wtheta, 0.0);
  G4LogicalVolume *logicGamChamSideWall = new G4LogicalVolume(solidGamChamSide,
                                                   beamPipeMaterial, "logicGamChamSideWall");

  new G4PVPlacement(0, G4ThreeVector(0.5*(vclength-vcwndz)*tan(wtheta)-swdx/2.0 + magnetx, 0.0, -vcwndz/2.0),
                    logicGamChamSideWall, "GamChamSide", logicGVCContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(new G4RotationMatrix(G4ThreeVector(0.0, 0.0, 1.0), M_PI),
                        G4ThreeVector(-0.5*(vclength-vcwndz)*tan(wtheta)+swdx/2.0 - magnetx, 0.0, -vcwndz/2.0),
                        logicGamChamSideWall, "GamChamSide", logicGVCContainer, false, 1, lxs->OverlapTest);

//Short piece of beampipe to join vacuum chamber with the pipe toward photon detectors
  G4double lbpipe = lxs->BeamPipe2ProfilerZ;
  G4Tubs *solidBeamPipeVCG = new G4Tubs("solidBeamPipeVCG", BPipeRProf-lxs->BPipeThickness, BPipeRProf,
                                           0.5*lbpipe, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeVCG = new G4LogicalVolume(solidBeamPipeVCG, beamPipeMaterial, "logicBeamPipeVCG");
  G4Tubs *solidBeamPipeVCGVac = new G4Tubs("solidBeamPipeVCGVac", 0.0, BPipeRProf-lxs->BPipeThickness,
                                           0.5*lbpipe, 0.0, 2.0*M_PI);
//Lid at the end of the beam pipe
  G4Material* bpipeLidMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->BeamPipeLidMaterial);
  G4LogicalVolume *logicBeamPipeVCGVac = new G4LogicalVolume(solidBeamPipeVCGVac,
                                                             vacuumMaterial, "logicBeamPipeVCGVac");
  G4Tubs *solidBeamPipeVCGLid = new G4Tubs("solidBeamPipeVCGLid", 0.0, BPipeRProf,
                                           0.5*lxs->BeamPipeLidThickness, 0.0, 2.0*M_PI);
  G4LogicalVolume *logicBeamPipeVCGLid = new G4LogicalVolume(solidBeamPipeVCGLid,
                                                             bpipeLidMaterial, "logicBeamPipeVCGLid");

  G4double chpos = lxs->GMagnetZpos + vclength/2.0 + gfieldz;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, chpos), logicGVCContainer,
                    "GMagnetWideChamber", fLogicWorld, false, 0, lxs->OverlapTest);

  G4double  bppos =chpos + (vclength+lbpipe)/2.0;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bppos), logicBeamPipeVCG,
                    "BeamPipeVCGProf", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bppos), logicBeamPipeVCGVac,
                    "BeamPipeVCGProfVac", fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bppos + (lbpipe + lxs->BeamPipeLidThickness)/2.0), logicBeamPipeVCGLid,
                    "BeamPipeVCGProfLid", fLogicWorld, false, 0, lxs->OverlapTest);

// Add magnet and field volume to the world
//   G4ThreeVector trm(lxs->GMagnetXpos+0.5*(lxs->GMagnetX - mfbox->GetXHalfLength()),
//                     lxs->GMagnetYpos, lxs->GMagnetZpos);
  G4ThreeVector trm(lxs->GMagnetXpos, lxs->GMagnetYpos, lxs->GMagnetZpos);
  mag->GetAssembly()->MakeImprint(fLogicWorld, trm, 0, 0, lxs->OverlapTest);
  mag->CostructSupport(fLogicWorld, trm, "GMagnet");

  new G4PVPlacement(0, trm, fieldvol, "GMagnetField", fLogicWorld, false, 0, lxs->OverlapTest);

// Add local magnetic field
  G4String magid = "Gamma";
  LxBField *gammamagf = ComposeFieldObject(magid, trm);
  G4FieldManager* fieldMgr = new G4FieldManager(gammamagf);
  ConfigureFieldManager(fieldMgr, gammamagf);

  fieldvol->SetFieldManager(fieldMgr, true);
  G4AutoDelete::Register(gammamagf);
  G4AutoDelete::Register(fieldMgr);

  // Add field to neighboring volumes in case it uses not constant model
  if (gammamagf->IsZLimited()) {
    std::vector<G4String> fieldvolname{"logicGVCContainer", "logicBeamPipeGamma2ndCMagFVac"};
    AddFieldToLogVolumes(fieldMgr, fieldvolname);
  } else {G4cout << "Gamma magnetic field is not limited in z, do not add neighboring volumes.\n";}

}



void DetectorConstruction::AddFieldToLogVolumes(G4FieldManager* fieldMgr, const std::vector<G4String> vname)
{
  G4LogicalVolumeStore *lvstor = G4LogicalVolumeStore::GetInstance();
  for (auto &vitr : vname) {
    G4LogicalVolume  *lv = lvstor->GetVolume (vitr);
    if (lv) {
      lv->SetFieldManager(fieldMgr, true);
      G4cout << "Add magnetic field to " << vitr << G4endl;
    } else G4cout << "!!!!!!!!!!!!!!!! Field was not added to volume " << vitr << " !!!!!!!!!!!!!!!!\n";
  }
}



void DetectorConstruction::ConstructComptShielding()
{
  LXSetUp *lxs = LXSetUp::Instance();
  
  G4Material* concreteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("ShieldingConcrete");
  G4Material* shildingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptShieldingMaterial);
  G4Material* shildingMaterial2 = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptShieldingMaterial2);
  G4Material* shildingPlateMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ComptShieldingPlateMaterial);

  G4Box *solidShielding0 = new G4Box("solidComptShielding0", lxs->ComptShieldingX/2.0, lxs->ComptShieldingY/2.0, 
                                     lxs->ComptShieldingZ/2.0);
  G4Tubs *solidShieldingH = new G4Tubs("solidComptShieldingH", 0.0, lxs->BPipeRLG, 0.55*lxs->ComptShieldingZ, 0.0, 2.0*M_PI);

// Concrete support for the shielding (not finished)
//   G4double  floordysup = -(lxs->ComptShieldingSuppY/2.0 + lxs->FloorSurfaceYpos);
//   - lxs->ComptonLysoY/2.0 - lxs->LYSODetSupportBallH - lxs->OPPPDetTopPlateY
//   G4Box *solidShieldingSupp = new G4Box("solidShieldingSupp", lxs->ComptShieldingX/2.0, lxs->ComptShieldingSuppY/2.0, 
//                                      3*lxs->ComptShieldingZ/2.0);  
//   G4LogicalVolume *logicShieldingSupp = new G4LogicalVolume(solidShieldingSupp, concreteMaterial, "logicComptShieldingSupp");   
//   G4double tableypos = lxs->FloorSurfaceYpos + lxs->ComptShieldingSuppY/2.0;                                                  
//   new G4PVPlacement(0, G4ThreeVector(0.0, -floordysup, lxs->ComptShieldingZpos), logicShieldingSupp,
//                     "ComptShieldingSupport", fLogicWorld, false, 0, lxs->OverlapTest);

// 3-layer shielding
  G4double floordy = -(lxs->ComptShieldingY/2.0 + lxs->FloorSurfaceYpos);
  G4Transform3D transformh(G4RotationMatrix(), G4ThreeVector(0.0, floordy, 0.0));
  G4SubtractionSolid* solidShielding1 = new G4SubtractionSolid("solidComptShielding", solidShielding0,
                                                               solidShieldingH, transformh);

//   G4LogicalVolume *logicShielding = new G4LogicalVolume(solidShielding1, shildingMaterial, "logicComptShielding");
//   G4LogicalVolume *logicShielding2 = new G4LogicalVolume(solidShielding1, shildingMaterial2, "logicComptShielding2");
//   G4double zs = lxs->ComptShieldingZ;
//   new G4PVPlacement(0, G4ThreeVector(0.0, -floordy, lxs->ComptShieldingZpos), logicShielding,
//                     "ComptShieldingFe", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, -floordy, lxs->ComptShieldingZpos - zs), logicShielding2,
//                     "ComptShielding1Al", fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(0.0, -floordy, lxs->ComptShieldingZpos + zs), logicShielding2,
//                     "ComptShielding2Al", fLogicWorld, false, 0, lxs->OverlapTest);


  //Lead plates imbedded in front of the shielding
  G4double ShieldingPlateX = lxs->ComptShieldingX;//100.0* cm;
  G4double ShieldingPlateY = 10.0* cm;
  G4double ShieldingPlateZ = 30.0* cm;
  G4double ShieldingPlatePosZ =  lxs->ComptShieldingZpos-lxs->ComptShieldingZ/2.0+ShieldingPlateZ/2.0;
  G4Box *solidShieldingPlate0 = new G4Box("solidShieldingPlate0", ShieldingPlateX/2.0, ShieldingPlateY/2.0,
                                     ShieldingPlateZ/2.0);
  G4Tubs *solidShieldingT = new G4Tubs("solidComptShieldingT", 0.0, lxs->BPipeRLG, 0.55*ShieldingPlateZ, 0.0, 2.0*M_PI);
  G4Transform3D transformt(G4RotationMatrix(), G4ThreeVector(0.0, 0.0, 0.0));
  G4SubtractionSolid* solidShieldingPlate = new G4SubtractionSolid("solidShieldingPlate", solidShieldingPlate0,
                                                               solidShieldingT, transformt);
  G4LogicalVolume *logicShieldingPlate = new G4LogicalVolume(solidShieldingPlate, shildingPlateMaterial, "logicComptShieldingPlate");

  // 1-layer concrete shielding
    G4Transform3D transformp(G4RotationMatrix(), G4ThreeVector(0.0, floordy, -lxs->ComptShieldingZ/2.0+ShieldingPlateZ/2.0));
    G4SubtractionSolid* solidShieldingP = new G4SubtractionSolid("solidShieldingP", solidShielding1,
                                                               solidShieldingPlate0, transformp); // cut the place for the plate
    G4LogicalVolume *logicShielding = new G4LogicalVolume(solidShieldingP, concreteMaterial, "logicComptShielding");
    new G4PVPlacement(0, G4ThreeVector(0.0, -floordy, lxs->ComptShieldingZpos), logicShielding,
                     "ComptShieldingConcrete", fLogicWorld, false, 0, lxs->OverlapTest);
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, ShieldingPlatePosZ), logicShieldingPlate,
                     "ComptShieldingPlate", fLogicWorld, false, 0, lxs->OverlapTest);   //insert plate placement

//Pipe through the shielding
//   G4double lshildpipe = lxs->ShieldingZ - lxs->ShieldingDeepZ + lxs->ShieldingDeepMargine;
//   G4double shildpipezpos = lxs->ComptShieldingZpos + 0.5 * (lxs->ShieldingZ - lshildpipe);
//
//   G4Tubs *solidShieldingPipe = new G4Tubs("solidShieldingPipe", lxs->BPipeR - lxs->BPipeThickness, lxs->BPipeR, lshildpipe/2.0, 0.0, 2.0*M_PI);
//   G4LogicalVolume *logicShieldingPipe = new G4LogicalVolume(solidShieldingPipe, beamPipeMaterial, "logicShieldingPipe");
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, shildpipezpos), logicShieldingPipe, "ShieldingPipe", fLogicWorld, false, 0, lxs->OverlapTest);
//
//   G4Tubs *solidShieldingPipeVac = new G4Tubs("solidShieldingPipeVac", 0.0, lxs->BPipeR - lxs->BPipeThickness, lshildpipe/2.0, 0.0, 2.0*M_PI);
//   G4LogicalVolume *logicShieldingPipeVac = new G4LogicalVolume(solidShieldingPipeVac, vacuumMaterial, "logicShieldingPipeVac");
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, shildpipezpos), logicShieldingPipeVac, "ShieldingPipeVac", fLogicWorld, false, 0, lxs->OverlapTest);
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



void DetectorConstruction::SetGammaBeamDumpLength(G4double val)
{
  LXSetUp::Instance()->SetGammaBeamDumpZ(val);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::SetGammaBeamDumpMaterial(G4String val)
{
  LXSetUp *lxs = LXSetUp::Instance();
  lxs->GammaBeamDumpMaterial = val;
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


void DetectorConstruction::ConstructMagnet(G4LogicalVolume *mfContainer, const G4ThreeVector mf)
{
//   G4double magnetFieldValue = 14000.0*gauss;
//   G4double magnetFieldValue = 0.0*gauss;
  
//   G4ThreeVector  fieldVector( 0.0, -magnetFieldValue, 0.0);  
  G4ThreeVector  fieldVector(mf);  
  G4MagneticField *magField = new G4UniformMagField( fieldVector );    
  G4FieldManager  *localFieldMgr = new G4FieldManager (magField); 
  localFieldMgr->CreateChordFinder(magField);

  mfContainer->SetFieldManager(localFieldMgr, true);   
}   



void DetectorConstruction::ConfigureFieldManager(G4FieldManager *fieldMgr, G4MagneticField *mfield)
{
  G4double minStep = 0.010*mm;
  auto pEquation = new G4Mag_UsualEqRhs(mfield);
//   G4int nvar = pEquation->GetNumberOfVariables();
  G4int nvar = 8;
  auto pStepper = new G4DormandPrince745( pEquation, nvar );
  auto pIntgrationDriver = new G4IntegrationDriver<G4DormandPrince745>(minStep, pStepper, nvar);
//   G4int nvar = 6; //it is default value in constructor;
//   auto pStepper = new G4CashKarpRKF45( pEquation, nvar );
//   auto pIntgrationDriver = new G4IntegrationDriver<G4CashKarpRKF45>(minStep, pStepper, nvar);
//   G4int nvar = 6; //it is default value in constructor;
//   auto pStepper = new ClassicalRK4( pEquation, nvar );
//   auto pStepper = new G4NystromRK4( pEquation );
  G4ChordFinder *pChordFinder = new G4ChordFinder(pIntgrationDriver);
  fieldMgr->SetChordFinder( pChordFinder );

  fieldMgr->SetMinimumEpsilonStep( 1.0e-5 );
  fieldMgr->SetMaximumEpsilonStep( 1.0e-4 );
  fieldMgr->SetDeltaOneStep( 0.5e-3 * mm );
  fieldMgr->SetDeltaIntersection(0.5e-3 * mm);

//   // magnetic field
//   G4FieldManager fieldMgr = new G4FieldManager();
//   fieldMgr->SetDetectorField(mfield);
//   fieldMgr->CreateChordFinder(mfield);
//   G4bool forceToAllDaughters = true;
//   fieldvol->SetFieldManager(fieldMgr, forceToAllDaughters);

}



LxBField* DetectorConstruction::ComposeFieldObject(const G4String magid, const G4ThreeVector magpos)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4double bxval, byval, bzval;
  bxval = byval = bzval = 0.0;

  if (magid == "IP") byval = lxs->IPMagFieldY;
  if (magid == "Brems") bxval = lxs->DumpMagFieldY;
  if (magid == "Gamma") byval = lxs->GMagFieldY;

  typedef VFieldComponentDistrib<FieldDistribution, FieldDistribution, FieldDistribution>  LxTFieldComponent;

  std::vector<FieldDistribution*> bcj(9,0);
  const auto msit = fBFieldModelsInfo.find(magid);
  if (msit != fBFieldModelsInfo.end()) {
    const auto &msetv = msit->second;
    for (auto svitr = msetv.cbegin(); svitr != msetv.cend(); ++svitr) {
      const G4String &fcomponent = std::get<0>(*svitr);
      const G4String &coordinate = std::get<1>(*svitr);
      const G4String &dmodel = std::get<2>(*svitr);
      const G4String &params = std::get<3>(*svitr);

      if (fcomponent == "Bx") {
        if (coordinate == "x") {
          bcj[0] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[1] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[2] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          bxval = G4UIcommand::ConvertToDimensionedDouble(params.data());
        }
      }

      if (fcomponent == "By") {
        if (coordinate == "x") {
          bcj[3] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[4] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[5] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          byval = G4UIcommand::ConvertToDimensionedDouble(params.data());
        }
      }

      if (fcomponent == "Bz") {
        if (coordinate == "x") {
          bcj[6] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "y") {
          bcj[7] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "z") {
          bcj[8] = CreateFieldDistribution(dmodel, params);
        }
        if (coordinate == "bvalue") {
          bzval = G4UIcommand::ConvertToDimensionedDouble(params.data());
        }
      }

    }
  }

  std::for_each(bcj.begin(), bcj.end(), [](FieldDistribution* &bd) {if (!bd) bd = new FieldDistribution(); });

  LxTFieldComponent *bx = new LxTFieldComponent(bcj[0], bcj[1], bcj[2], bxval);
  LxTFieldComponent *by = new LxTFieldComponent(bcj[3], bcj[4], bcj[5], byval);
  LxTFieldComponent *bz = new LxTFieldComponent(bcj[6], bcj[7], bcj[8], bzval);

  LxBField *bfield = new LxBField(bx, by, bz, magpos);

  G4cout << "=========  " << magid << " magnet settings" << "  =========" << G4endl;
  G4cout << "Bx: x: " << typeid(*(bcj[0])).name() 
          << "  y: " << typeid(*(bcj[1])).name() 
          << "  z: " << typeid(*(bcj[2])).name() << "  Value: " << bxval/tesla << "T" << G4endl;
  G4cout << "By: x: " << typeid(*(bcj[3])).name() 
          << "  y: " << typeid(*(bcj[4])).name() 
          << "  z: " << typeid(*(bcj[5])).name() << "  Value: " << byval/tesla << "T"  << G4endl;
  G4cout << "Bz: x: " << typeid(*(bcj[6])).name() 
          << "  y: " << typeid(*(bcj[7])).name() 
          << "  z: " << typeid(*(bcj[8])).name() << "  Value: " << bzval/tesla << "T"  << G4endl;
  G4cout << "====================================" << G4endl;

 return bfield;
}



FieldDistribution* DetectorConstruction::CreateFieldDistribution(const G4String &fmodel, const G4String &params)
{
  FieldDistribution *bfd = 0;
  std::stringstream paramstr(params.data());
  G4String vunits; 
  if (fmodel == "const") {
    G4double cmin, cmax, vu;
    paramstr >> cmin >> cmax >> vunits;
    vu = G4UIcommand::ValueOf(vunits);
    bfd = new FieldConst(cmin*vu, cmax*vu);
  }
  if (fmodel == "f_fd") {
    G4double cmin, cmax, tmin, tmax, vu;
    paramstr >> cmin >> cmax >> tmin >> tmax >> vunits;
    vu = G4UIcommand::ValueOf(vunits);
    bfd = new FieldFD(cmin*vu, cmax*vu, tmin*vu, tmax*vu);
  }
  if (fmodel == "f_err") {
    G4double cmin, cmax, tmin, tmax, vu;
    paramstr >> cmin >> cmax >> tmin >> tmax >> vunits;
    vu = G4UIcommand::ValueOf(vunits);
    bfd = new FieldErrF(cmin*vu, cmax*vu, tmin*vu, tmax*vu);
  }
  return bfd;
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
  
  regName = "BShieldingRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  if (!region) {
    region = new G4Region(regName);
  }
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(lxs->ShieldingProductionCut, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(lxs->ShieldingProductionCut, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(lxs->ShieldingProductionCut, G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(lxs->ShieldingProductionCut); // same cuts for gamma, proton, e- and e+
  region->SetProductionCuts(cuts);
  G4LogicalVolume *bshilding = G4LogicalVolumeStore::GetInstance()->GetVolume("logicShielding");
  if (bshilding) {
    bshilding->SetRegion(region);
    region->AddRootLogicalVolume(bshilding);
  }
  
  regName = "GammaDumpRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  if (!region) {
    region = new G4Region(regName);
  }
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(lxs->GammaDumpProductionCut, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(lxs->GammaDumpProductionCut, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(lxs->GammaDumpProductionCut, G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(lxs->GammaDumpProductionCut); // same cuts for gamma, proton, e- and e+
  region->SetProductionCuts(cuts);
  G4LogicalVolume *gammadump = G4LogicalVolumeStore::GetInstance()->GetVolume("logicGammaBeamDump");
  if (gammadump) {
    gammadump->SetRegion(region);
    region->AddRootLogicalVolume(gammadump);
  }
    
  regName = "HICSDumpRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  if (!region) {
    region = new G4Region(regName);
  }
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(lxs->HICSDumpProductionCut, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(lxs->HICSDumpProductionCut, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(lxs->HICSDumpProductionCut, G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(lxs->HICSDumpProductionCut); // same cuts for gamma, proton, e- and e+
  region->SetProductionCuts(cuts);
  G4LogicalVolume *hicsdump = G4LogicalVolumeStore::GetInstance()->GetVolume("logicHICSDump");
  if (hicsdump) {
    hicsdump->SetRegion(region);
    region->AddRootLogicalVolume(hicsdump);
  }

}


void DetectorConstruction::SetAbsorberType(G4String val)
{
  if (val == "foil") {
    LXSetUp::Instance()->GTargetType = LXSetUp::tTargetType::tfoil;
  } else if (val == "wire") {
    LXSetUp::Instance()->GTargetType = LXSetUp::tTargetType::twire;
  } else {
    G4cout << "DetectorConstruction::SetAbsorberType: <" << val << ">"
           << " is not supported!"
           << G4endl;
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}



void DetectorConstruction::DumpBFieldModel()
{
  G4cout << "======= Following setting are used to create magnetic fields =======\n"; 
  for (const auto &mitr : fBFieldModelsInfo) {
    const auto &tinfo = mitr.second;
    G4cout << "MagnetID: " << mitr.first << G4endl;
    for (const auto &vitr : tinfo) {
      G4cout << "     " << std::get<0>(vitr) << "  " << std::get<1>(vitr)
             << "  " << std::get<2>(vitr) << "  " << std::get<3>(vitr) << G4endl;
    }
  }
  G4cout << "====================================================================\n"; 
}

