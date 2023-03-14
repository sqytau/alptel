//
/// \brief Implementation of the LxECal class
//

#include <algorithm>
#include <functional>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
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
#include "LxAux.hh"
#include "LxDetector.hh"
#include "LxECal.hh"


////////////////////////////////////////////////////////////////////////
//// LxECal

void LxECal::Construct()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  CreateMaterial();

  LXSetUp *lxs = LXSetUp::Instance();
//   G4String TBeam_senrio = "T AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS";
  G4String TBeam_senrio = "AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS AS";

  G4Material* ecalContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* Silicon = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  G4Material* Aluminium = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
//   G4Material* Tungsten = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

  G4Material* PLASTIC_SC = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* Wabsorber_PL = G4NistManager::Instance()->FindOrBuildMaterial("Wabsorber_PL");
  G4Material* Wabsorber_MGS = G4NistManager::Instance()->FindOrBuildMaterial("Wabsorber_MGS");
  G4Material* FanoutMatB = G4NistManager::Instance()->FindOrBuildMaterial("BackFanoutMaterial");
  G4Material* FanoutMatF = G4NistManager::Instance()->FindOrBuildMaterial("FrontFanoutMaterial");
  G4Material* C_fiber = G4NistManager::Instance()->FindOrBuildMaterial("ECalCarbonFiber");

  G4double base_airx = lxs->ECalX;
  G4double base_airy = lxs->ECalY;

  G4double Lcal_absorber_pitch = 1  *mm;
  G4double hAbsorberPitchDZ = Lcal_absorber_pitch / 2.0;

  G4double Lcal_tungsten_thickness = 3.5 *mm;
  G4double Lcal_tungsten_hdz = Lcal_tungsten_thickness / 2.0;
  G4double hTungstenDZ = Lcal_tungsten_hdz;      // half thickness absorber

  G4double Lcal_silicon_thickness  = 0.32 *mm;
  G4double Lcal_silicon_hdz  = Lcal_silicon_thickness / 2.0;
  G4double hSiliconDZ = Lcal_silicon_hdz;       // half thickness of the silicon

  G4double Lcal_epoxy_heightF  = 0.065*mm;
  G4double Lcal_kapton_heightF = 0.050*mm;
  G4double Lcal_copper_heightF = 0.035*mm;
  G4double Lcal_fanoutF_thickness = Lcal_epoxy_heightF + Lcal_kapton_heightF + Lcal_copper_heightF;
  G4double Lcal_fanoutF_hdz = Lcal_fanoutF_thickness / 2.0;
  G4double hFanoutFrontDZ = Lcal_fanoutF_hdz;

  G4double Lcal_epoxy_heightB  = 0.060 *mm;
  G4double Lcal_kapton_heightB = 0.065 *mm;
  G4double Lcal_copper_heightB = 0.025 *mm;
  G4double Lcal_fanoutB_thickness = Lcal_epoxy_heightB + Lcal_kapton_heightB + Lcal_copper_heightB;
  G4double Lcal_fanoutB_hdz = Lcal_fanoutB_thickness / 2.0;
  G4double hFanoutBackDZ  = Lcal_fanoutB_hdz;       // half thickness fanout back

  G4double Lcal_pad_metal_thickness = 0.02 *mm;
  G4double hMetalDZ = Lcal_pad_metal_thickness/2.0;
  G4double zpos_PSC = 1130.0 *mm; // for 2016 from 4th Telescope plane to Plastic scintilator
  G4double zpos_TR_PSC = 20.0 *mm; // for 2016 from 4th Telescope plane to TRiger Plastic scintilator

  G4int n_layers = lxs->ECalNLayers;

  G4double DUTextrahz = 0.002 *mm;

  G4double airhz =  hAbsorberPitchDZ + hTungstenDZ;  // should be 2.25 mm half a slut
  G4double airhz1mm = hAbsorberPitchDZ; // should be 0.5 mm
  G4double airhz_PSC =  10.0*mm ;
  G4double airhz_TR_PSC =  3.75*mm ;

  fECalLayerZ = airhz * 2.0 + DUTextrahz;
  G4double base_airz = n_layers * (airhz * 2.0 + DUTextrahz);
//   G4double base_airz = n_layers * airhz * 2.0;

  G4Box *solidECalContainer = new G4Box("solidECalContainer", base_airx/2.0, base_airy/2.0, base_airz/2.0);
  G4LogicalVolume *logicECalContainer = new G4LogicalVolume(solidECalContainer,
                                                            ecalContainerMaterial, "logicECalContainer");

  //Test Plane "Telescope"
  G4Box *solidBaseUnit = new G4Box ( "solidBaseUnit", base_airx/2.0, base_airy/2.0, airhz );
  G4LogicalVolume *logicBaseUnit = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnit");
  G4LogicalVolume *logicBaseUnit_T = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnit_T");
  G4LogicalVolume *logicBaseUnit_S = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnit_S");
  G4LogicalVolume *logicBaseUnit_PL = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnit_PL");
  G4LogicalVolume *logicBaseUnitONLYAB = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnitONLYAB");
  G4LogicalVolume *logicBaseUnitONLYAB_PL = new G4LogicalVolume (solidBaseUnit, Air, "logicBaseUnitONLYAB_PL");

  G4Box *solidBaseUnitFor1mm = new G4Box ( "solidBaseUnitFor1mm", base_airx/2.0, base_airy/2.0, airhz1mm );
  G4LogicalVolume *logicBaseUnitFor1mm = new G4LogicalVolume (solidBaseUnitFor1mm, Air, "logicBaseUnitFor1mm");
//   G4LogicalVolume *logicBaseUnitFor1mm_NoS = new G4LogicalVolume (solidBaseUnitFor1mm, Air, "logicBaseUnitFor1mm_NoS");

  ///for Plastic SC
  G4Box *solidBase_PSC = new G4Box ( "solidBase_PSC", base_airx/2.0, base_airy/2.0, airhz_PSC );
  G4LogicalVolume *logicBaseUnit_PSC = new G4LogicalVolume (solidBase_PSC, PLASTIC_SC , "logicBaseUnit_PSC");

  //for trigger SC
  G4Box *solidBase_TR_PSC = new G4Box ( "solidBase_TR_PSC", base_airx/2.0, base_airy/2.0, airhz_TR_PSC );
  G4LogicalVolume *logicBaseUnit_TR_PSC = new G4LogicalVolume (solidBase_TR_PSC, PLASTIC_SC , "logicBaseUnit_TR_PSC");

  G4double zposAbs0 = -airhz + hTungstenDZ;
  G4double zposAbs0_MSG = -airhz + hTungstenDZ*1.02;

  G4Box *solidAbs0 = new G4Box("solidAbs0", base_airx/2.0, base_airy/2.0, hTungstenDZ );
  G4LogicalVolume *logicAbs0_PL = new G4LogicalVolume(solidAbs0, Wabsorber_PL , "logicAbs0_PL");
//   G4LogicalVolume *logicAbs0 = new G4LogicalVolume(solidAbs0, Tungsten, "logicAbs0");

  // nominal thiknes of MSG plate is ~3.57
  G4Box *solidAbs0_MGS = new G4Box("solidAbs0_MGS", base_airx/2.0, base_airy/2.0, 1.02*hTungstenDZ); 
  G4LogicalVolume *logicAbs0_MGS = new G4LogicalVolume(solidAbs0_MGS, Wabsorber_MGS, "logicAbs0_MGS");

 new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zposAbs0_MSG), logicAbs0_MGS, "Absorber0_MGS"  ,logicBaseUnit,
                    false, 0, lxs->OverlapTest);
 new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zposAbs0)    , logicAbs0_PL,  "Absorber0_PL"   ,logicBaseUnit_PL,
                    false, 0, lxs->OverlapTest);
 new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zposAbs0_MSG), logicAbs0_MGS, "Absorber0_MGS"  ,logicBaseUnitONLYAB,
                    false, 0, lxs->OverlapTest);
 new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, zposAbs0)    , logicAbs0_PL,  "Absorber0_PL"   ,logicBaseUnitONLYAB_PL,
                    false, 0, lxs->OverlapTest);

//----------------------------------------------
//	carbon fiber suppurt
//---------------------------------------------

  G4double CFhz = 0.395 *mm;

  G4Box *solidCF = new G4Box ( "solidCF", base_airx/2.0, base_airy/2.0, CFhz );
  G4LogicalVolume *logicCF = new G4LogicalVolume (solidCF, C_fiber, "logicCF");

//---------------------------------------------
//		sensor stracture from front to back  :
//			1. 0.150 kapton front.
//			1.2. 0.010 epoxy glue
//			2.1. 0.020 AL on sensor
//			2.2. 0.320 Si sensor
//			2.3. 0.020 AL on sensor
//			3.1. 0.040 epoxy (condactive - not take into acount)
//			3.2. 0.025 cupper (on kapton)
//			3.3. 0.065 kapton back
// 			3.4. 0.020 epoxy glue
//			total 0.670
//---------------------------------------------
  G4double sensor_x = lxs->ECalSensorPixelX * lxs->ECalSensorNCellX;
  G4double sensor_y = lxs->ECalSensorPixelY * lxs->ECalSensorNCellY;

  G4Box *solidECalSensor = new G4Box("solidECalSensor", sensor_x/2.0, sensor_y/2.0, hSiliconDZ);
  G4Box *solidMetal = new G4Box("solidMetal", sensor_x/2.0, sensor_y/2.0, hMetalDZ);
  G4Box *solidFanoutFrnt  = new G4Box("solidFanoutFrnt", sensor_x/2.0, sensor_y/2.0, hFanoutFrontDZ);
  G4Box *solidFanoutBack  = new G4Box("solidFanoutBack", sensor_x/2.0 , sensor_y/2.0, hFanoutBackDZ);

  G4LogicalVolume *logicECalSensor = new G4LogicalVolume(solidECalSensor, Silicon, "logicECalSensor");
  G4LogicalVolume *logicMetalV = new G4LogicalVolume(solidMetal, Aluminium, "logicMetalV");
  G4LogicalVolume *logicFanoutFrnt = new G4LogicalVolume(solidFanoutFrnt, FanoutMatF, "logicFanoutFront");
  G4LogicalVolume *logicFanoutBack = new G4LogicalVolume(solidFanoutBack, FanoutMatB, "logicFanoutFront");

//-------------------------------------------------------------
//		placment of sensor part in CF
//-------------------------------------------------------------
  G4double ypos = 0.0;

  G4double SensorAtCF = -CFhz + hFanoutFrontDZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, ypos, SensorAtCF), logicFanoutFrnt , "FunOut0", logicCF,
                      false, 0, lxs->OverlapTest);

  SensorAtCF = SensorAtCF + hFanoutFrontDZ +  hMetalDZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, ypos, SensorAtCF), logicMetalV , "PadMetal0", logicCF,
                      false, 0, lxs->OverlapTest);

  SensorAtCF = SensorAtCF + hMetalDZ + hSiliconDZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, ypos, SensorAtCF), logicECalSensor, "ECalSensor", logicCF,
                     false, 0, lxs->OverlapTest);

  SensorAtCF = SensorAtCF + hSiliconDZ + hMetalDZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, ypos, SensorAtCF), logicMetalV, "PadMetal1", logicCF, 
                     false, 0, lxs->OverlapTest);

  SensorAtCF =  SensorAtCF + hMetalDZ + hFanoutBackDZ;
  new G4PVPlacement (0, G4ThreeVector(0.0, ypos, SensorAtCF), logicFanoutBack, "FunOut1", logicCF, 
                      false, 0, lxs->OverlapTest);

  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, airhz - CFhz), logicCF, "CF0",logicBaseUnit, 
                     false, 0, lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, airhz - CFhz), logicCF, "CF0",logicBaseUnit_PL, 
                     false, 0,lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, airhz1mm - CFhz), logicCF, "CF0", logicBaseUnitFor1mm, 
                     false, 0,lxs->OverlapTest);
  new G4PVPlacement (0, G4ThreeVector(0.0, 0.0, airhz - CFhz), logicCF, "CF0", logicBaseUnit_S, 
                     false, 0,lxs->OverlapTest);

  int Is_SC = 0;
  int Is_TSC = 0;
  int Is_stag = 0;
  G4double zyposLC = 0.0 *mm;
  G4double zposLC = -base_airz/2.0;

  G4double ypos_stag[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  std::istringstream iss((std::string)TBeam_senrio);
  G4int iplacelayer = 0; // for sensors
  G4int iplaceElements = 0; // for all
  std::string delimiter = ":";
  while ((iss)&& (iplacelayer < n_layers)) {
//     std::cout<<"!!!!!!!!!!!!!!!Put layer!!!!!!!!!!!!!!\n";
    std::stringstream placement_name("");
    std::string Osub;
    iss >> Osub;
//     std::cout << "Substring:" << Osub << std::endl;
    std::string tmpsub = Osub.substr(0, Osub.find(delimiter));
    std::string sub;
    int i_air = 0;

    if ( tmpsub.length() < Osub.length()) { // check if the ":" for air gup size for tracker is there
      Osub.erase(0, 1 + tmpsub.length());
      i_air = std::atoi(Osub.c_str()) *mm;
      sub = tmpsub;
    } else {
      sub = tmpsub;
    }

    if (sub == "S") {
      placement_name << "DUTAS" << iplacelayer;
      new  G4PVPlacement (0, G4ThreeVector(0.0, zyposLC+ypos_stag[iplacelayer], zposLC + airhz),logicBaseUnit_S,
                          placement_name.str().c_str(), logicECalContainer, 0, iplacelayer, lxs->OverlapTest);
//       G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<< " at z = " 
//       << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC += 2.0*airhz + DUTextrahz;
      ++iplacelayer;
    }  else if (sub == "T") {
      placement_name << "DUTAS" << iplacelayer;
      new  G4PVPlacement (0, G4ThreeVector(0.0, zyposLC+ypos_stag[iplacelayer], zposLC + airhz),logicBaseUnit_T,
                          placement_name.str().c_str(), logicECalContainer, 0, iplacelayer, lxs->OverlapTest);
//       std::cout << " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<< " at z = " << zposLC
//              <<" layer with name "<< placement_name.str().c_str() << std::endl;
      zposLC += 2.0*airhz + DUTextrahz;
      ++iplacelayer;
    }	else if (sub == "AS") {
      placement_name << "DUTAS" << iplacelayer;
      new  G4PVPlacement (0, G4ThreeVector(0.0, zyposLC+ypos_stag[iplacelayer], zposLC + airhz),logicBaseUnit, 
                           placement_name.str().c_str(), logicECalContainer, 0, iplacelayer, lxs->OverlapTest);
//       std::cout << " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<< " at z = " 
//              << zposLC<<" layer with name "<< placement_name.str().c_str() << std::endl;
      zposLC += 2.0*airhz + DUTextrahz;
      ++iplacelayer;
    } else if(sub == "AS_PL") {
      placement_name << "DUTASP" << iplacelayer;
      new  G4PVPlacement(0, G4ThreeVector(0.0, zyposLC+ypos_stag[iplacelayer], zposLC + airhz  ),logicBaseUnit_PL,
                           placement_name.str().c_str(), logicECalContainer, 0, iplacelayer,lxs->OverlapTest);
//       G4cout<< " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<<" at z = " << zposLC
//             <<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC += 2.0*airhz + DUTextrahz;
      ++iplacelayer;
    } else if (sub == "A_PL") {
      placement_name << "ABP" << iplaceElements;
      new  G4PVPlacement (0, G4ThreeVector(0.0, zyposLC, zposLC +airhz), logicBaseUnitONLYAB_PL, 
                             placement_name.str().c_str(), logicECalContainer, 0, iplaceElements, lxs->OverlapTest);
//       G4cout<< " placed "<<iplaceElements <<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC += 2.0*airhz + DUTextrahz;
      //iplaceA++;
    } else if (sub == "A") {
      placement_name << "AB" << iplaceElements;
      new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zposLC + airhz ),logicBaseUnitONLYAB, 
                           placement_name.str().c_str(), logicECalContainer, 0, iplaceElements,lxs->OverlapTest);
//       G4cout<< " placed "<<iplaceElements <<" at z = " << zposLC<<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC = zposLC + 2.0*airhz + DUTextrahz;
//    iplaceA++;
    } else if (sub == "JSM") {
      placement_name << "DUTJSM" << iplacelayer;
      new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC+ypos_stag[iplacelayer], zposLC+ airhz1mm ),logicBaseUnitFor1mm,
                           placement_name.str().c_str(), logicECalContainer, 0, iplacelayer, lxs->OverlapTest);
//       G4cout << " placed "<<iplaceElements << " as sensor number  : " << iplacelayer<<" at z = " << zposLC
//              <<" layer with name "<< placement_name.str().c_str() << G4endl;
      zposLC= zposLC+ 2.0*airhz1mm + DUTextrahz + i_air; //add the air gup after
      ++iplacelayer;
    } else if (sub == "SC") {
      if(Is_SC == 0 ) {
        placement_name << "SC" ;
        new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zpos_PSC ),logicBaseUnit_PSC, placement_name.str().c_str(), 
                             logicECalContainer, 0, 1, lxs->OverlapTest);
//         G4cout<< " placed SC" <<" at z = " << zpos_PSC <<" layer with name "<< placement_name.str().c_str() << G4endl;
        Is_SC =1; // can by placed only 1 time
      }
    } else if(sub == "TSC") {
      if(Is_TSC == 0 ){
        placement_name << "TR_SC" ;
        new  G4PVPlacement ( 0, G4ThreeVector( 0., zyposLC, zpos_TR_PSC ),logicBaseUnit_TR_PSC, 
                             placement_name.str().c_str(), logicECalContainer, 0, 1, lxs->OverlapTest);
//         G4cout<< " placed SC" <<" at z = " << zpos_TR_PSC <<" layer with name "<< placement_name.str().c_str() << G4endl;
        Is_TSC =1; // can by placed only 1 time
      }
    } else if(sub == "stag") {
      if(Is_stag == 0 ) {
//         G4cout<< "Stag should applay before all sensor layers " << G4endl;
        double TB_stag[8] = {0.0, 0.0, 0.2*mm, -0.7 *mm, 1.5 *mm, -1.0 *mm, 0.0 ,0.0};// asked by marina
        // proper stagging
        //double TB_stag[8] = {-0.11 * mm, -1.26* mm, 0.46*mm, -0.275644 *mm, 1.87705 *mm, -1.2183 *mm, 0.53323 *mm,0 *mm};
        //double TB_stag[8] = {0.11 * mm, 1.26* mm, -0.46*mm, 0.275644 *mm, -1.87705 *mm, 1.2183 *mm, -0.53323 *mm,0 *mm};
        for(int istag = 0; istag < 8; istag++) {
          ypos_stag[istag] = TB_stag[istag];
//           G4cout << "Stugging: layer " << istag <<" movment :  " << ypos_stag[istag] << G4endl;
        }
        Is_stag =1; // can by placed only 1 time
      }
    } else {
//       G4cout << "Substring: " << sub <<" in layer " << iplaceElements <<" is not recognized we will brack"<< G4endl;
      break;
    }
    ++iplaceElements;
  }

  G4double shift_x = base_airx/2.0 + lxs->ECalXpos;
  G4double trackerzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet + lxs->OPPPDetTopPlateZ/2.0;
  G4double shift_z = trackerzpos + lxs->OPPPDetTopPlateZ/2.0 + lxs->OPPPTrackerECalZ + base_airz/2.0;
  new G4PVPlacement(0, G4ThreeVector(shift_x, 0.0, shift_z), logicECalContainer, "OPPPECal", 
                    fLogicWorld, false, 0, lxs->OverlapTest);
//   new G4PVPlacement(0, G4ThreeVector(-shift_x, 0.0, shift_z), logicECalContainer, "OPPPECal", 
//                     fLogicWorld, false, 1, lxs->OverlapTest);

  G4LogicalVolume *logicECalPCBContainer = ConstructPCB(base_airx, base_airz);
  new G4PVPlacement(0, G4ThreeVector(shift_x, 0.5*(base_airy + lxs->OPPPECalPCBY), shift_z), logicECalPCBContainer,
                                "OPPPECalPCBContainer", fLogicWorld, false, 0, lxs->OverlapTest);

  G4AssemblyVolume *casingAssembly = ConstructCasingAssembly(base_airz);
  G4double casingx = shift_x;
  G4double casingz = shift_z + lxs->ECalCasingGapZ/2.0;
  G4ThreeVector trcasing(casingx, 0.0, casingz);
  casingAssembly->MakeImprint(fLogicWorld, trcasing, 0, 0, lxs->OverlapTest);

  G4AssemblyVolume *supportAssembly = ConstructSupportAssembly();

  G4double detxpos = lxs->OPPPDetTopPlateX/2.0 + lxs->OPPPDetXPos;
  G4double supypos = -(lxs->ECalY + lxs->OPPPDetTopPlateY)/2.0;
  G4double detzpos = shift_z - base_airz/2.0 + lxs->ECalTopPlateZ/2.0;
  G4ThreeVector trsupport(detxpos, supypos, detzpos);
  supportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 0, lxs->OverlapTest);
//   trsupport.setX(-detxpos);
//   supportAssembly->MakeImprint(fLogicWorld, trsupport, 0, 1, lxs->OverlapTest);

  ConstructShielding();
  ConstructDumpShielding();
  AddSegmentation();
}



void LxECal::AddSegmentation()
{
  LXSetUp *lxs = LXSetUp::Instance();
  fDetector->AddSensorSegmentation("ECalSensor", lxs->ECalX, lxs->ECalY, 
                                                  lxs->ECalSensorNCellX, lxs->ECalSensorNCellY);
}



G4AssemblyVolume* LxECal::ConstructSupportAssembly()
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* opppDetSupportMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPDetSupportMaterial);

  G4Box *solidECalTopPlate = new G4Box("solidECalTopPlate", lxs->OPPPDetTopPlateX/2.0, 
                                            lxs->OPPPDetTopPlateY/2.0, lxs->ECalTopPlateZ/2.0);
  G4LogicalVolume *logicECalTopPlate = new G4LogicalVolume(solidECalTopPlate, opppDetSupportMaterial,
                                                                "logicECalTopPlate");
  
  G4double trackerh = lxs->OPPPColdPlateY/2.0 + lxs->OPPPStaveLiftY + lxs->OPPPDetTopPlateY + lxs->OPPPDetSupportBallH;
  G4double ECalSupportBallH = trackerh - lxs->ECalY/2.0 - lxs->OPPPDetTopPlateY;
  G4Sphere *solidECalSupportBall = new G4Sphere("solidECalSupportBall", 0.0, ECalSupportBallH/2.0,
                                                   0.0, 2.0*M_PI, 0.0, M_PI);
  G4LogicalVolume *logicECalSupportBall = new G4LogicalVolume(solidECalSupportBall, opppDetSupportMaterial,
                                                                "logicECalSupportBall");

  G4AssemblyVolume *ecalSupportAssembly = new G4AssemblyVolume();
  G4ThreeVector asstr(0.0, 0.0, 0.0);
  ecalSupportAssembly->AddPlacedVolume(logicECalTopPlate, asstr, 0);

  G4double dbz = 0.5*lxs->ECalTopPlateZ - 1.05*ECalSupportBallH;
  G4double dbx = 0.5*lxs->OPPPDetTopPlateX - 1.05*ECalSupportBallH;
  G4ThreeVector bpos(dbx, asstr.y() - (lxs->OPPPDetTopPlateY+ECalSupportBallH)/2.0 , dbz);
  ecalSupportAssembly->AddPlacedVolume(logicECalSupportBall, bpos, 0);
  bpos.setX(-dbx);
  ecalSupportAssembly->AddPlacedVolume(logicECalSupportBall, bpos, 0);
  bpos.setX(0.0);
  bpos.setZ(-dbz);
  ecalSupportAssembly->AddPlacedVolume(logicECalSupportBall, bpos, 0);

  return ecalSupportAssembly;
}



G4AssemblyVolume* LxECal::ConstructCasingAssembly(const G4double ecalz)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* ecalCasingMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ECalCasingMaterial);

  double casingz = ecalz + lxs->ECalCasingGapZ;
  G4Box *solidECalCasingTop = new G4Box("solidECalCasingTop", lxs->ECalX/2.0,
                                            lxs->ECalCasingTopY/2.0, casingz/2.0);
  G4LogicalVolume *logicECalCasingTop = new G4LogicalVolume(solidECalCasingTop, ecalCasingMaterial,
                                                                "logicECalCasingTop");

  G4double casingy = lxs->ECalY + lxs->ECalCasingTopY + lxs->OPPPECalPCBY;
  G4Box *solidECalCasingSide = new G4Box("solidECalCasingSide", lxs->ECalCasingSideX/2.0,
                            casingy/2.0, (casingz+lxs->ECalCasingBackZ)/2.0);
  G4LogicalVolume *logicECalCasingSide = new G4LogicalVolume(solidECalCasingSide, ecalCasingMaterial,
                                                                "logicECalCasingSide");

  G4Box *solidECalCasingBack = new G4Box("solidECalCasingBack", lxs->ECalX/2.0,
                                         casingy/2.0, lxs->ECalCasingBackZ/2.0);
  G4LogicalVolume *logicECalCasingBack = new G4LogicalVolume(solidECalCasingBack, ecalCasingMaterial,
                                                                "logicECalCasingBack");
  G4AssemblyVolume *ecalCasingAssembly = new G4AssemblyVolume();
  G4double sideypos = 0.5 * (casingy - lxs->ECalY);
  G4ThreeVector toppos(0.0, sideypos + (casingy - lxs->ECalCasingTopY)/2.0, 0.0);
  ecalCasingAssembly->AddPlacedVolume(logicECalCasingTop, toppos, 0);
  G4ThreeVector sidepos(0.5*(lxs->ECalX+lxs->ECalCasingSideX), sideypos, lxs->ECalCasingBackZ/2.0);
  ecalCasingAssembly->AddPlacedVolume(logicECalCasingSide, sidepos, 0);
  sidepos.setX(-sidepos.x());
  ecalCasingAssembly->AddPlacedVolume(logicECalCasingSide, sidepos, 0);
  G4ThreeVector backpos(0.0, sideypos, 0.5*(casingz+lxs->ECalCasingBackZ));
  ecalCasingAssembly->AddPlacedVolume(logicECalCasingBack, backpos, 0);

  return ecalCasingAssembly;
}



G4LogicalVolume* LxECal::ConstructPCB(const G4double ecalx, const G4double ecalz)
{
  LXSetUp *lxs = LXSetUp::Instance();
  G4Material* envMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* ecalPCBMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->OPPPECalPCBMaterial);

  G4double pcbcontx = ecalx;
  G4double pcbcontz = ecalz;
  G4Box *solidECalPCBContainer = new G4Box("solidECalPCBContainer", pcbcontx/2.0, lxs->OPPPECalPCBY/2.0, pcbcontz/2.0);
  G4LogicalVolume *logicECalPCBContainer = new G4LogicalVolume(solidECalPCBContainer,
                                                            envMaterial, "logicECalPCBContainer");
  G4Box *solidECalPCB = new G4Box("solidECalPCB", pcbcontx/2.0, lxs->OPPPECalPCBY/2.0, lxs->OPPPECalPCBZ/2.0);
  G4LogicalVolume *logicECalPCB = new G4LogicalVolume(solidECalPCB, ecalPCBMaterial, "logicECalPCB");

  G4int n_layers = lxs->ECalNLayers-1;
  G4double offsetz = fECalLayerZ - lxs->OPPPECalPCBZ;
  G4double layerz = fECalLayerZ;
  for (G4int il = 0; il < n_layers; ++il) {
    G4double lzpos = il*layerz + offsetz + (lxs->OPPPECalPCBZ - pcbcontz)/2.0;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, lzpos), logicECalPCB,
                                "OPPPECalPCB", logicECalPCBContainer, false, il, lxs->OverlapTest);
  }

  return logicECalPCBContainer;
}



void LxECal::ConstructShielding()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

// Shielding along the beam pipe
  G4Material* ecalPipeShieldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ECalPipeShieldMaterial);
  G4Box *solidECalPipeShield = new G4Box("solidECalPipeShield", lxs->ECalPipeShieldX/2.0,
                                            lxs->ECalPipeShieldY/2.0, lxs->ECalPipeShieldZ/2.0);
  G4LogicalVolume *logicECalPipeShield = new G4LogicalVolume(solidECalPipeShield, ecalPipeShieldMaterial,
                                                                "logicECalPipeShield");
  G4double bpshielxpos = lxs->ECalXpos - lxs->ECalCasingSideX - lxs->ECalPipeShieldX/2.0;
  G4double bpshielzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet
                         - lxs->ECalPipeShieldShiftZ + lxs->ECalPipeShieldZ/2.0;
  new G4PVPlacement(0, G4ThreeVector(bpshielxpos, 0.0, bpshielzpos), logicECalPipeShield, "ECalPipeShield",
                    fLogicWorld, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(-bpshielxpos, 0.0, bpshielzpos), logicECalPipeShield, "ECalPipeShield",
                    fLogicWorld, false, 1, lxs->OverlapTest);
}



void LxECal::ConstructDumpShielding()
{
  const G4VPhysicalVolume *physicalWorld = fDetector->GetphysiWorld();
  G4LogicalVolume   *fLogicWorld = physicalWorld->GetLogicalVolume();

  LXSetUp *lxs = LXSetUp::Instance();

// Shielding which protect from the beam dump
  G4Material* ecalDumpShieldContainerMaterial = G4NistManager::Instance()->FindOrBuildMaterial(lxs->EnvironmentMaterial);
  G4Material* ecalDumpShieldMaterialInner = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ECalDumpShieldMaterialInner);
  G4Material* ecalDumpShieldMaterialOuter = G4NistManager::Instance()->FindOrBuildMaterial(lxs->ECalDumpShieldMaterialOuter);

  G4double shield_z = lxs->ECalDumpShieldInnerZ + 2.0 * lxs->ECalDumpShieldOuterZ;

  G4Box *solidECalDumpShieldContainer = new G4Box("solidECalDumpShieldContainer", lxs->ECalDumpShieldX/2.0,
                                            lxs->ECalDumpShieldY/2.0, shield_z/2.0);
  G4LogicalVolume *logicECalDumpShieldContainer = new G4LogicalVolume(solidECalDumpShieldContainer, ecalDumpShieldContainerMaterial,
                                                                "logicECalDumpShieldContainer");

  G4Box *solidECalDumpShieldInner = new G4Box("solidECalDumpShieldInner", lxs->ECalDumpShieldX/2.0,
                                            lxs->ECalDumpShieldY/2.0, lxs->ECalDumpShieldInnerZ/2.0);
  G4LogicalVolume *logicECalDumpShieldInner = new G4LogicalVolume(solidECalDumpShieldInner, ecalDumpShieldMaterialInner,
                                                                "logicECalDumpShieldInner");

  G4Box *solidECalDumpShieldOuter = new G4Box("solidECalDumpShieldOuter", lxs->ECalDumpShieldX/2.0,
                                            lxs->ECalDumpShieldY/2.0, lxs->ECalDumpShieldOuterZ/2.0);
  G4LogicalVolume *logicECalDumpShieldOuter = new G4LogicalVolume(solidECalDumpShieldOuter, ecalDumpShieldMaterialOuter,
                                                                "logicECalDumpShieldOuter");

  new G4PVPlacement(0, G4ThreeVector(), logicECalDumpShieldInner, "ECalDumpShieldInner",
                    logicECalDumpShieldContainer, false, 0, lxs->OverlapTest);

  G4double dz = 0.5*(lxs->ECalDumpShieldInnerZ + lxs->ECalDumpShieldOuterZ);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -dz), logicECalDumpShieldOuter, "ECalDumpShieldOuter",
                    logicECalDumpShieldContainer, false, 0, lxs->OverlapTest);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, dz), logicECalDumpShieldOuter, "ECalDumpShieldOuter",
                    logicECalDumpShieldContainer, false, 1, lxs->OverlapTest);

  G4double ecaldsXpos = lxs->ECalDumpShieldXpos + lxs->ECalDumpShieldX/2.0;
  G4double ecaldsYpos = lxs->ECalDumpShieldYpos;
  G4double ecalz = lxs->ECalNLayers * fECalLayerZ;
  G4double trackerzpos = lxs->IPMagnetZpos + lxs->FlashMFieldLength/2.0 + lxs->OPPPDetZtoMagnet + lxs->OPPPDetTopPlateZ/2.0;
  G4double ecaldsZpos = ecalz + trackerzpos + lxs->OPPPDetTopPlateZ/2.0 + lxs->OPPPTrackerECalZ;
  ecaldsZpos += lxs->ECalDumpShieldZpos + shield_z/2.0;

  new G4PVPlacement(0, G4ThreeVector(ecaldsXpos, ecaldsYpos, ecaldsZpos), logicECalDumpShieldContainer, "ECalDumpShiel",
                    fLogicWorld, false, 0, lxs->OverlapTest);

  G4double supx = 0.9 * lxs->ECalDumpShieldX;
  G4double supy = -lxs->FloorSurfaceYpos - lxs->ECalDumpShieldY/2.0;
  G4double supz = 0.9 * shield_z;
  G4LogicalVolume *pedstal = LxAux::BuildPedestal("ECalDumpShiel", supx, supy, supz);

  new G4PVPlacement(0, G4ThreeVector(ecaldsXpos, -(supy+lxs->ECalDumpShieldY)/2.0, ecaldsZpos), pedstal, "ECalDumpShielPedestal",
                    fLogicWorld, false, 0, lxs->OverlapTest);
}



void LxECal::CreateMaterial()
{
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  
// 
  G4int ncomponents;
  G4double fractionmass;

  G4NistManager *materials = G4NistManager::Instance();
  
  G4Material* Copper = materials->FindOrBuildMaterial("G4_Cu");
  G4Material* Kapton = materials->FindOrBuildMaterial("G4_KAPTON");
  G4Material* Tungsten = materials->FindOrBuildMaterial("G4_W");
  G4Material* Ni = materials->FindOrBuildMaterial("G4_Ni");
  G4Material* Carbon = materials->FindOrBuildMaterial("G4_C");
  
      // EPOXY - build up from elements
  G4Element* H  =  new G4Element("Hydrogen", symbol="H", z= 1.0, a= 1.01*g/mole);
  G4Element* C  =  new G4Element("Carbon"  , symbol="C", z= 6.0, a= 12.01*g/mole);
  G4Element* O  =  new G4Element("Oxygen"  , symbol="O", z= 8.0, a= 16.00*g/mole);
    
  G4Material* Epoxy = materials->FindOrBuildMaterial("Epoxy");
  if (!Epoxy) {
    Epoxy = new G4Material("Epoxy", density= 1.3*g/cm3, ncomponents=3);
    Epoxy->AddElement(H, fractionmass=0.1310);
    Epoxy->AddElement(C, fractionmass=0.5357);
    Epoxy->AddElement(O, fractionmass=0.3333);
  }
    
  G4double epoxydens =  Epoxy->GetDensity();
  G4double kaptondens = Kapton->GetDensity();
  G4double copperdens = Copper->GetDensity();

  // FrontFanoutMaterial
  G4double Lcal_epoxy_propF    = 1.0;
  G4double Lcal_copper_propF   = 0.5;
  G4double Lcal_epoxy_heightF  = 0.065*mm;
  G4double Lcal_kapton_heightF = 0.050*mm;
  G4double Lcal_copper_heightF = 0.035*mm;
  G4double Lcal_fanoutF_thickness = Lcal_epoxy_heightF + Lcal_kapton_heightF + Lcal_copper_heightF;
  
  G4double epoxydensF  = epoxydens  * Lcal_epoxy_propF;
  G4double copperdensF = copperdens * Lcal_copper_propF;
  G4double epoxyfracF  = Lcal_epoxy_heightF  / Lcal_fanoutF_thickness;
  G4double kaptonfracF = Lcal_kapton_heightF / Lcal_fanoutF_thickness;
  G4double copperfracF = Lcal_copper_heightF / Lcal_fanoutF_thickness;
  G4double frontDensity = (epoxydensF  * epoxyfracF + kaptondens  * kaptonfracF + copperdensF * copperfracF);

  G4Material* FanoutMatF = new G4Material("FrontFanoutMaterial", frontDensity, ncomponents=3);
  FanoutMatF->AddMaterial(Copper, fractionmass=epoxyfracF);
  FanoutMatF->AddMaterial(Kapton, fractionmass=kaptonfracF);
  FanoutMatF->AddMaterial(Epoxy,  fractionmass=copperfracF);
  
  // BackFanoutMaterial
  G4double Lcal_epoxy_propB    = 1.0;
  G4double Lcal_copper_propB   = 1.0;
  G4double Lcal_epoxy_heightB  = 0.060 *mm;
  G4double Lcal_kapton_heightB = 0.065 *mm;
  G4double Lcal_copper_heightB = 0.025 *mm;
  G4double Lcal_fanoutB_thickness = Lcal_epoxy_heightB + Lcal_kapton_heightB + Lcal_copper_heightB;
                                    
  G4double epoxydensB = epoxydens * Lcal_epoxy_propB;
  G4double copperdensB = copperdens * Lcal_copper_propB;
  G4double epoxyfracB  = Lcal_epoxy_heightB  / Lcal_fanoutB_thickness;
  G4double kaptonfracB = Lcal_kapton_heightB / Lcal_fanoutB_thickness;
  G4double copperfracB = Lcal_copper_heightB / Lcal_fanoutB_thickness;
  G4double backDensity = (2.0* epoxydensB  * epoxyfracB + kaptondens  * kaptonfracB + copperdensB * copperfracB);
  G4Material* FanoutMatB = new G4Material("BackFanoutMaterial", backDensity, ncomponents=3);
  FanoutMatB->AddMaterial(Copper, fractionmass=epoxyfracB);
  FanoutMatB->AddMaterial(Kapton, fractionmass=kaptonfracB);
  FanoutMatB->AddMaterial(Epoxy,  fractionmass=copperfracB);
  
  G4Material* C_fiber = new G4Material("ECalCarbonFiber", 1.6*g/cm3, ncomponents=2);
  C_fiber->AddMaterial(Carbon, fractionmass=50.0*perCent);
  C_fiber->AddMaterial(Epoxy, fractionmass=50.0*perCent);
  
  //Absorber MGS
  G4Material* Wabsorber_MGS = new G4Material("Wabsorber_MGS", 17.7*g/cm3, ncomponents=3);
  Wabsorber_MGS->AddMaterial(Tungsten, fractionmass=93.0*perCent);
  Wabsorber_MGS->AddMaterial(Ni, fractionmass=5.25*perCent);
  Wabsorber_MGS->AddMaterial(Copper, fractionmass=1.75*perCent);

  // for the 2014 TB plate
  G4Material* Wabsorber_PL = new G4Material("Wabsorber_PL", 18.0*g/cm3 ,ncomponents=3);
  Wabsorber_PL->AddMaterial(Tungsten, fractionmass=95.0*perCent);
  Wabsorber_PL->AddMaterial(Ni, fractionmass=2.5*perCent);
  Wabsorber_PL->AddMaterial(Copper, fractionmass=2.5*perCent);

  // fiber glass
  G4int natoms;
  G4Element* Si = new G4Element("Silcone"  ,symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Material* fiberglass = new G4Material( "fiberglass",density=2.61*g/cm3,ncomponents=2);
  fiberglass->AddElement(Si, natoms=1);
  fiberglass->AddElement(O, natoms=2);
  // PCBoard material FR4
  G4Material *FR4 = new G4Material("ECAL_FR4",density=1.93*g/cm3,ncomponents=3);
  FR4->AddMaterial( Epoxy, fractionmass=0.34);
  FR4->AddMaterial( fiberglass, fractionmass=0.56);
  FR4->AddMaterial( Copper, fractionmass=0.1);

}
