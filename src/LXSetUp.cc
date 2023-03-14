
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "LXSetUp.hh"


LXSetUp *LXSetUp::flxSetup = 0;


LXSetUp::LXSetUp()
{
  OverlapTest = true;

  EnvironmentMaterial = "Air20";
//   EnvironmentMaterial = "Galactic";
  FloorSurfaceYpos = -2478.0 *mm; //-1700.0 *mm;
  FloorY = 500.0 *mm;

  //  CeilingSurfaceYpos = 3866.0 *mm; // approximately

  CeilingSurfaceYpos = 1816.0 *mm; // approx. for Scint screen camera (IP) ceiling mount. downstream (lower) part of chamber 

  IPBoxX = 45.0 *cm;
  IPBoxY = 25.0 *cm;
//   IPContainerZ = 100.0 *cm;
  IPBoxZ = 60.0 *cm;
  IPBoxThickness = 5.0 *mm;

  BPipeR = 24.13 *mm;          // DN40;
  BPipeThickness = 1.651 *mm;  // seems to be the thinnest for DN40  //2.0 *mm;

  BeamPipeMaterial = "Aluminium";  // "StainlessSteel";  //"Iron";
  BeamPipeWindowMaterial = "Aluminium";

//   BeamPipeWindowMaterial = "Aluminium";
//   BeamPipeWindowMaterial = "Beryllium";
//   BeamPipeWindowMaterial = "Galactic";

  BeamCrossingAngle = 0.3 *rad;
  BeamPipeVacuumMaterial = "XFELVacuum"; //"Galactic";
  LaserPipeR = 25.0 *mm;
  LaserPipeThickness = 2.0 *mm;

  MagnetMaterial = "Iron";
  DumpMagnetXpos = 0.0 *mm; // -50.0 *mm;
  DumpMagnetYpos = -100.0 *mm;
//   DumpMagnetZpos = -800.0 *cm;
  DumpMagnetZpos = -650.0 *cm;
  DumpMagnetX = 40.0*cm;
  DumpMagnetY = 20.0 *cm;
  DumpMagnetZ = 100.0 *cm;
//   DumpMagnetZ = 1520.0 *mm;
  DumpMagnetThickness = 50.0 *cm;

//   DumpMagFieldY = 18000.0*gauss; // for 16.5 GeV // 20000.0*gauss; // for 17.5 GeV

  DumpMagFieldY = 0.0*gauss;
  IPMagFieldY = -10000.0*gauss; //-16000.0*gauss;
//   IPMagFieldY = -20000.0*gauss; //-16000.0*gauss;

  DMBPipeWindow = M_PI/4.0 *rad;
  CerenkovchannelY = 5. * mm;
  CerenkovBoxThickness = 1. * cm;

//since Magnet / Brem detectors are rotated 90 degrees, X dimensions here (Scint/Cerenkov) map to global Y dimensions & vice versa  
  ScintCerenkovPhysics = false;
  ScintPhysicsYield = 0.01;
  ScintAngle = 15. *deg;
  ScintXpos = (-500. - 110.)*mm;
  ScintZpos = DumpMagnetZpos + (1520./2. + 275.) *mm; //  defining in relation to position of the beam-dump magnet and its length   
  CerenkovAngle = 0. *deg;
  CerenkovXpos = -605. *mm;
  CerenkovZpos = DumpMagnetZpos + (1520./2. + 950.) *mm; 
  ScintX = 1000. *mm;
  ScintY = 100. *mm;
  ScintBaseZ = 0.19 *mm;
  ScintPhosphorZ = 0.3 *mm;
  ScintFinishZ = 0.01 *mm;
  ScintFrameThickness = 10. *mm;
  ScintFrameMaterial = "Aluminium";
  ScintFrameBeamLoopX = 15. * mm; // for HICS region
  BremScintFrameWidthOffset = 15.*cm; // increase in width on either side of Brem Scint frame for camera visibility 
  
  CerenkovChannels = 240;         
  CerenkovStrawHeight = 200. *mm;
  CerenkovChannelLayers = 4;
  CerenkovStrawInnerRadius = 4.02/2. * mm;
  CerenkovStrawGraphiteLayerThickness = 6. * 0.001 * mm;
  CerenkovStrawAlLayerThickness = 0.2 * 0.001 * mm;
  CerenkovStrawKaptonLayerThickness = 25 * 0.001 * mm;
  CerenkovStrawPolyurethaneLayerThickness = 5 * 0.001 * mm;
  CerenkovElectronicsBoardThickness = 2. * mm;
  CerenkovShieldingPlateThickness = 0.0001 * mm; // 0.5 * cm; // same material / thickness as EcalPipeShield 
  CerenkovBoxThickness = 0.9 * cm;
  CerenkovBeamWindowThickness = 0.15 * mm;
    
  CerenkovStrawXFrequency = 16. * mm;
  CerenkovStrawZFrequency = 16. * mm;
  CerenkovStrawLayerOffset = 4. * mm; // = CerenkovStrawXFrequency/ CerenkovChannelLayers and approx total straw diameter
  CerenkovStrawPlateBuffer = 1. * cm;
  CerenkovSupportThickness = 1. * cm;
  CerenkovSupportPosYtoStage = CerenkovBoxThickness + 5. * mm;
  CerenkovBeamWindowY = 15. * mm;
  CerenkovUpperSupporttoElectronicsBoard = 5. * mm;
  CerenkovTotalBoxHeight = CerenkovStrawHeight + 2*CerenkovBoxThickness + 2* CerenkovSupportThickness + CerenkovElectronicsBoardThickness + CerenkovUpperSupporttoElectronicsBoard + 5. * mm;
 CerenkovTotalBoxLength =  2*((CerenkovChannels/CerenkovChannelLayers - 1)*CerenkovStrawXFrequency/2. + (CerenkovChannelLayers-1)*CerenkovStrawLayerOffset/2. +  CerenkovStrawInnerRadius + 2*(CerenkovStrawGraphiteLayerThickness+CerenkovStrawAlLayerThickness+CerenkovStrawKaptonLayerThickness+CerenkovStrawPolyurethaneLayerThickness) + CerenkovStrawPlateBuffer + CerenkovBoxThickness + CerenkovShieldingPlateThickness);


  
  ScintBaseMaterial = "G4_POLYETHYLENE";
  ScintPhosphorMaterial = "G4_GADOLINIUM_OXYSULFIDE";
  ScintFinishMaterial = "G4_POLYETHYLENE";
  CerenkovMetal = "Aluminium";
  CerenkovMedium = "Air20Refractive"; // "ArgonGas";
  CerenkovStrawGraphiteLayerMaterial = "G4_GRAPHITE";
  CerenkovStrawAlLayerMaterial = "Aluminium";
  CerenkovStrawKaptonLayerMaterial = "G4_KAPTON";
  CerenkovStrawPolyurethaneLayerMaterial = "G4_POLYETHYLENE";
  CerenkovElectronicsBoardMaterial = "G4_KAPTON";
  CerenkovShieldingPlateMaterial = "Air20"; // "Lead"; // same material / thickness as EcalPipeShield 

  ScintCameraX = 29. *mm;
  ScintCameraY = 29. *mm;
  ScintCameraZ = 42. *mm;
  ScintCameraAngle = 22. *deg;
  ScintCameraApertureZ = 12. *mm;
  ScintCameraApertureOuter = 29. *mm;
  ScintCameraApertureInner = 20. *mm;
  ScintCameraShieldThickness = 2. *cm;
  ScintCameraSupportThickness = 5. *cm;
  ScintCameraPlatformZ = 15. *cm;
  StickScintScreentoBeamWindow = false;

 
  ScintCameraMaterial = "Aluminium";
  ScintCameraShieldMaterial = "Lead";
  ScintCameraApertureMaterial = "Air20";
  ScintCameraSupportMaterial = "Aluminium";
    
  BeamDumpR = 13.0 *cm;
  BeamDumpZ = 50.0 *cm;
  BeamDumpFrontZpos = 350 *cm + DumpMagnetZpos + DumpMagnetZ/2.0;
  BeamDumpFrontXpos = 13.5 *cm;
  BeamDumpYpos = 0.0 *mm;
  BeamDumpAngle = atan2(BeamDumpFrontXpos, fabs(BeamDumpFrontZpos - DumpMagnetZpos));
  BeamDumpZpos = BeamDumpFrontZpos + 0.5*BeamDumpZ*cos(BeamDumpAngle);
  BeamDumpXpos = BeamDumpFrontXpos + 0.5*BeamDumpZ*sin(BeamDumpAngle);
  DumpBeamHoleRin = BPipeR - BPipeThickness;
  DumpBeamHoleRout = BPipeR - BPipeThickness;
//   BeamDumpMaterial = "Iron";
  BeamDumpMaterial = "Copper";

  BeamDumpInsertR = 6.5 *cm;
  BeamDumpInsertZ = 20.0 *cm;
  BeamDumpInsertMaterial = "Aluminium";

  BeamDumpProductionCut = 5.0 *cm;

  BeamDumpPipeMaterial = "Aluminium";
  BeamDumpPipeR = 0.5*63.0 *mm;
  BeamDumpPipeThickness = 1.0 *mm;
  
//   DumpBeamHoleRin = 1.0 *mm;
//   DumpBeamHoleRout = 1.0 *mm;
  
//   DumpCleanMagZpoz = -220.0 *cm;
  DumpCleanMagZpoz = -290.0 *cm;
  DumpCleanMagZ = 100.0 *cm;
  DumpCleanMagX = 20.0 *cm;
  DumpCleanMagY = 20.0 *cm;
  DumpCleanMagThickness = 20.0 *cm;
  DumpCleanMagField = 14000.0*gauss;
  
  ShieldingMaterial = "Iron"; //"ShieldingConcrete";
//   ShieldingZpos = -200.0 *cm;
//   ShieldingZpos = -140.0 *cm;
  ShieldingZpos = 0.0 *cm;  // this is not used, it is calculated in DetectorConstruction::ConstructShielding() 
  ShieldingX = 300.0 *cm; //122.0 *cm; //200.0 *cm;
  ShieldingY = 400.0 *cm; // 350 *cm;
  ShieldingDXBeamWall = 61.5 *cm;
  // ShieldingY = 2.0* std::fabs(FloorSurfaceYpos);
//   ShieldingZ = 80.0 *cm;
//   ShieldingZ = 60.0 *cm;
  ShieldingZ = 150.0 *cm;
  ShieldingDeepZ = 100.0 *cm;
  ShieldingDeepMargine = 2.0 *cm;
  ShieldingProductionCut = 5.0 *cm;
  ShieldingAbsorberX = 20.0 *cm;
  ShieldingAbsorberZ = 20.0 *cm;
  ShieldingAbsorberMaterial = "Lead";
  ShieldingAbsorberSlitX = 2.0 *cm;
  ShieldingAbsorberSlitZ = 5.0 *cm;
  ShieldingAbsorberTopZ = 30.0 *cm;

//   BeamPipeSplitterZ = 200.0 *cm;  // 1.4 Tesla
  BeamPipeSplitterZ = 100.0 *cm;     
  
  IPMagnetXpos = 0.0 *mm;
  IPMagnetYpos = 0.0 *mm;
//   IPMagnetZpos = 150.0 *cm;
  IPMagnetZpos = (150.0 + 55.0) *cm;

//   OPPPDetZpos = 300.0 *cm;
  OPPPDetZtoMagnet = 130.0 *cm;  // considere the distance from the end of field area of the magnet
  OPPPDetXPos = BPipeR + 5.0*mm; //57.0 *mm;  // actually the distance from the beam axis to the detector
  OPPPDetYLowPos = 150.0 *mm;
  OPPPDetX = 100.0 *cm;
  OPPPDetY = 10.0 *cm;
  OPPPDetZ = 70.0 *cm;
//   OPPPDetMaterial = "Galactic";
  OPPPDetMaterial = "Air20";
  
  OPPPTrackerX = OPPPDetX;
  OPPPTrackerY = OPPPDetY;
  OPPPTrackerZ = 0.05 *mm;
  OPPPTrackerDZ = 10.0 *cm;
  OPPPTrackerNLayers = 4;
  OPPPTrackerMaterial = "Silicon";
  OPPPTrackerNCellX = 34200; //29.24 um //20000; //50um
  OPPPTrackerNCellY = 3720;  //26.88 um //2000; //50um

  OPPPCaloX = OPPPDetX;
  OPPPCaloY = OPPPDetY;
  OPPPCaloZ = OPPPDetZ - OPPPTrackerNLayers * (OPPPTrackerZ + OPPPTrackerDZ);
  OPPPCaloMaterial = "Tungsten";

  BTargetX = 10.0 *cm;
  BTargetY = 10.0 *cm;
  BTargetZ = 15.0*mm; //0.035 *mm;
  BTargetChamberMagnetGapZ = 20.0 *cm; // distance from the target chamber to the magnet
  BTargetMaterial = "Tungsten";

  TargetChamberX = 300.0 *mm;
  TargetChamberY = 300.0 *mm;
  TargetChamberZ = 10.0 *cm;
  TargetChamberXPos = 8.0 *cm;
  TargetChamberThickness = 5.0 *mm;
  TargetChamberMaterial = "Aluminium";
  TargetChamberMotorMaterial = "StainlessSteel";
  TargetChamberMotorX = 150.0 *mm;
  TargetChamberMotorY = 15.0 *mm;
  TargetChamberMotorZ = 45.0 *mm;
  TargetFrameD = 15.0 *mm;
  TargetFrameZ = 10.0 *mm;
  TargetFrameMaterial = "Aluminium";
  TargetFrameHolderR = 1.5 * TargetFrameZ;
  TargetFrameHolderGroveY = 5.0 *mm;
  ChamberTargetVolX = BTargetX;
  ChamberTargetVolY = BTargetY;
  ChamberTargetVolZ = TargetFrameZ;

  BTargetZpos = DumpMagnetZpos - 0.5*(DumpMagnetZ + TargetChamberZ) - BTargetChamberMagnetGapZ;

      // Photon detector system

   // Pipe
  BPipeRLG = BPipeR; //2.5 *cm; //19.00 *mm;          // Special size of the pipe in FDS ;

  GTargetX = 5.0 *cm;
  GTargetY = 5.0 *cm;
  GTargetZ = 0.01 *mm; //0.2 *mm;
  GTargetZpos = 6.5 *m;// 3.86 from IP to OPPPdet + 2m 
  GTargetMaterial = "Tungsten"; // "G4_KAPTON"; //
//   GTargetType = twire;
  GTargetType = tfoil;

  GTargetContainerX = 2.0 * GTargetX;
  GTargetContainerY = 2.0 * GTargetX;
  GTargetContainerZ = (GTargetZ < 1.0*cm) ? 10.0 *cm : 2.0 * GTargetZ;
  


  CollimatorMaterial = "Lead";
  CollimatorSupportMaterial = "ShieldingConcrete";
  CollimatorZ = 50.0 *cm; 
  CollimatorZpos = GTargetZpos + CollimatorZ/2.0 + 50.0*cm;
  CollimatorRout = 5.0 *cm;
  CollimatorRin = 0.4 *cm;
  
  GMagnetXpos = 0.0 *cm;
  GMagnetYpos = 0.0 *cm;
  GMagnetY = 97.0 *cm;
  GMagnetX = 73.75 *cm;
  GMagnetZ = 128.0 *cm;
  GMagnetZpos = CollimatorZpos + 1.5*CollimatorZ + GMagnetZ/2.0 + 2.0*50.0*cm; //10.0 *m;
  GMagFieldY = 1.0*tesla; //14000.0 *gauss;//10000.0 *gauss;
  GMagCoilThickness = 150.0* mm;
  GMagCoilH = 160.0* mm;
  GMagnetCutX= 36.0*cm;
  GMagnetCutY= 47.0*cm;
  GMagneteffY = 5.0*cm;
  GMagnetCoreX = 18.0*cm;
  GMFieldX = 30.0*cm;
  GMFieldLength = 980.0 *mm;

  QBeamPipeContainerX = GMagnetCoreX;  //5.0 *cm;
  GChamberTBWallThickness = 5.0 *mm;
  GammaBPipeWindowThickness = 0.5 *mm;
  WideBeamPipeContainerX = GMagnetCutX - GMFieldX; //orig GMagnetCutX - GMagnetCoreX;
  GammaMagnetBeamPipeXGap = 1.0 *mm;

  ComptonDetX = 100.0 *cm;
  ComptonDetY = 10.0 *cm;
  ComptonDetZ = 50.0 *cm;
  ComptonDetXpos = 10.0 *cm + ComptonDetX/2.0;
  ComptonDetYpos = 0.0 *cm;
  ComptonDetZpos = GMagnetZpos + 0.75 *m + 0.5*(ComptonDetZ + GMFieldLength);
//   ComptonDetMaterial = "Galactic";
  ComptonDetMaterial = "Air20";

  ComptonLysoX = 100.0*cm; //8.0 *cm;
  ComptonLysoY = 5.0 *cm;
  ComptonLysoZ = 0.5 *mm;
//  ComptonLysoDZ = 2.0 *cm;
  ComptonLysoNCellX = 2000; //80; //2 x 1 mm (finer segmentation in x; the deflection direction) giving 25 x 80 bins
  ComptonLysoNCellY = 100;
  ComptonLysoMaterial = "G4_GADOLINIUM_OXYSULFIDE"; // "LANEX"; //"LYSO";
  ComptonLysoZpos = GMagnetZpos + 0.75 *m + 0.5*(ComptonLysoZ + GMFieldLength);
  ComptonLysoXpos = 1.0 *cm + ComptonLysoX/2.0;
  LYSODetTopPlateX = 520.0 *mm;
  LYSODetTopPlateZ = 70.0 *mm;
  LYSODetBottomSupportX = 520.0 *mm;
  LYSODetBottomSupportZ = 350.0 *mm;
  LYSODetSupportBallH = 20.0 *mm;
  LYSOPedestaMaterial = "Aluminium";
  LYSOBasePlateX = 2.0*(LYSODetTopPlateX + 1.0 *cm);
  LYSOBasePlateY = 100.0 *mm;
  LYSOBasePlateZ = 500.0 *mm;

  ComptonTrackerX = ComptonDetX;
  ComptonTrackerY = ComptonDetY;
  ComptonTrackerZ = 0.05 *mm;
  ComptonTrackerDZ = 5.0 *cm;
  ComptonTrackerNLayers = 4;
  ComptonTrackerMaterial = "Silicon";
  ComptonTrackerNCellX = 34200; //29.24 um 20000; //50um
  ComptonTrackerNCellY = 3720;  //26.88 um 2000; //50um

  ComptonCaloX = ComptonDetX;
  ComptonCaloY = ComptonDetY;
  ComptonCaloZ = ComptonDetZ - ComptonTrackerNLayers * (ComptonTrackerZ + ComptonTrackerDZ);
  ComptonCaloMaterial = "Tungsten";

 // Gamma Beam Dump
//   GammaDumpZpos = 12.0 *m;
//   GammaDumpZpos =  ComptShieldingZpos + 0.5*(ComptShieldingZ) + 2.5 *m;
//   GammaDumpZ = 50.0 *cm; //100.0 *cm;
//   GammaDumpR = 30.0 *cm;
//   GammaDumpMaterial = "Iron";
  GammaDumpProductionCut = 1.0 *mm;// 5.0 *cm;
  GammaBeamDumpZpos = ComptonDetZpos + 3.0 *m;
  GammaBeamDumpZ = 100.0 *cm; //50.0 *cm;
  GammaBeamDumpR = 30.0 *cm;
  GammaBeamDumpMaterial = "Copper";
  GammaBeamDumpInsertR = BPipeRLG; // 13.0  6.5 
  GammaBeamDumpInsertZ = 15.0 *cm; //40.0 *cm; //20
  GammaBeamDumpInsertMaterial = "Air20"; //"Galactic"; // "Aluminium";
  GammaBeamDumpSupportMaterial = "ShieldingConcrete";

// Gamma Monitor
  GammaMonitorDetMaterial = "Air20"; //"Galactic";
  GammaMonitorSupportMaterial = "ShieldingConcrete";
  GammaMonitorHolderMaterial = "StainlessSteel";
  GammaCaloZ = 50.0 *cm; //0.5 *cm;//10.0 * 0.3504 *cm;
  GammaCaloZpos = GammaBeamDumpZpos - 0.5*(GammaBeamDumpZ + GammaCaloZ) - 100.0*mm; //- 20.0*mm;
  GammaCaloMaterial = "LG_TF1"; //"TF1"; //"F101"; //"Silicon";"Tungsten"
  GammaCaloFoilMaterial = "Aluminium";
  GammaCaloFoilThickness = 0.016 *mm; // standard household foil
  LeadGlassX = 3.8 *cm;
  LeadGlassY = 3.8 *cm;
  LeadGlassZ = 45.0 *cm;

  CalSupportMaterial = "StainlessSteel";
  CalSupportR = 12.0 *cm;
  CalSupporthickness = 1.0 *cm;
  CAbsorberMaterial = "Lead"; //"Aluminium";
  CAbsorberThickness = 40.0 *mm;
  CAbsorberZpos = GammaBeamDumpZpos - 0.5*(GammaBeamDumpZ + CAbsorberThickness) - 20.0*mm;

  ComptShieldingMaterial = "Iron";
  ComptShieldingMaterial2 = "Aluminium";
  ComptShieldingZpos = ComptonDetZpos + 1.3 *m;
  ComptShieldingX = 360.0 *cm;
  // ComptShieldingY = 2.0* std::fabs(FloorSurfaceYpos);
  ComptShieldingY = 300 *cm;
  ComptShieldingZ = 80.0 *cm;
  ComptShieldingProductionCut = 5.0 *cm;
  ComptShieldingSuppY = 1.0 *m;
  ComptShieldingPlateMaterial = "Lead";

  IPChamberMaterial = "StainlessSteel";
  IPRmin = 452.44 *mm;
  IPRmax = 457.20 *mm;
  IPHight = 269.88 *mm;
  IPRHoleCut = 76.505 *mm;
  IPChamberHolesRmin = {BPipeR-BPipeThickness, 76.505 *mm, 76.505 *mm, 76.505 *mm, 
                        BPipeR-BPipeThickness, 76.505 *mm, 76.505 *mm, 76.505 *mm};
  IPChamberHolesPhi = {0.0, M_PI/4.0, 2.0*M_PI/4.0, 3.0*M_PI/4.0, 4.0*M_PI/4.0, 5.0*M_PI/4.0, 6.0*M_PI/4.0, 7.0*M_PI/4.0};
//   IPChamberHolesRmin = {BPipeR-BPipeThickness, 76.505 *mm, 76.505 *mm, BPipeR-BPipeThickness, 76.505 *mm, 76.505 *mm};
//   IPChamberHolesPhi = {0.0, M_PI/3.0, 2.0*M_PI/3.0, 3.0*M_PI/3.0, 4.0*M_PI/3.0, 5.0*M_PI/3.0};
  IPChamberPipeLength = 47.15 *mm;
  IPChamberPipeFlangedR = 20.0 *mm;
  IPChamberPipeFlangeThickness = 10.0 *mm;
  IPChamberBottomThickness = 22.23 *mm;
  IPChamberTopThickness = 12.7 *mm;
  IPChamberTopR = 514.35 *mm;
  // This redefinition is to adapt beam pipe length to the new size of the interaction chamber
  IPContainerZ = 2.0 * (IPRmax + IPChamberPipeLength);
  IPChamberSupportH = -0.5*IPHight - IPChamberBottomThickness - FloorSurfaceYpos;
  IPChamberSupportR = 40.0 *mm;
  IPChamberSupportRpos = 0.7*IPRmin;
  
//   TypMBMagnetX = 900.0 *mm;
  TypMBMagnetX = 1196.0 *mm;
  TypMBMagnetY = 672.0 *mm;
  TypMBMagnetZ = 1520.0 *mm;
  TypMBMagnetCoreZ = 1396.0 *mm;
//   TypMBMagnetHoleX = 326.0 *mm;
  TypMBMagnetHoleX = 330.0 *mm;
  TypMBMagnetEndCapHoleX = 371.0 *mm;
  TypMBMagnetEndCapHoleOutX = 428.74 *mm;
  TypMBMagnetHoleY = 108.0 *mm;
  TypMBMagnetEndCapZ = 50.0 *mm;
  TypMBMagnetEndCapW = 804.0 *mm;
  TypMBFieldLength = 1029 *mm;

  TypMBWireMaterial = "Copper";
  TypMBWireGap = 20.0 *mm;
  TypMBWireXW = 100.0 *mm;
  TypMBWireXT = 50.0 *mm;
  
  TypMBCoreL = 1080.0 *mm;
  TypMBWireL = TypMBCoreL + 2.0*TypMBWireGap;
  TypMBWireH = 330.0 *mm;
  TypMBWireYShift = 100.0*mm;
  TypMBCoreW = TypMBWireH + 2.0*(TypMBWireXW + TypMBWireGap);
  TypMBCoreT = 0.5*(TypMBMagnetY - TypMBMagnetHoleY);
  TypMBWedgeLFraction = 0.3;
  
  TypMBSupportH = -FloorSurfaceYpos - TypMBMagnetY/2.0;
  TypMBSupportV = -FloorSurfaceYpos - TypMBMagnetX/2.0;
  TypMBSupportR = 40.0 *mm;
  TypMBSupportPosDz = 800.0 *mm;
  TypMBSupportPosDx = 0.7*TypMBMagnetX;

  FlashMMiddleCutY = 310.32 *mm;
  FlashMMiddleCutX = 570.0 *mm;
  FlashMFieldGapY = 59.92 *mm;
  FlashMagnetX = 980.0 *mm;
  FlashMagnetY = 700.0 *mm;
  FlashMagnetZ = 1200.0 *mm;
  FlashMagnetCutD = 50.0 *mm;
  FlashMagnetCoilZ = 1440.0 *mm;
  FlashMagnetCoilD = 100.0 *mm;
  FlashMagnetCoilH = 110.0 *mm;
  FlashMagnetCoilGapH = 10.0 *mm;
  FlashMagnetCoreX = FlashMMiddleCutX - 2.0*FlashMagnetCoilH - 4.0*FlashMagnetCoilGapH;
  FlashMagnetCoreZ = 1200.0 *mm;
  FlashMFieldX = FlashMMiddleCutX;
  FlashMagneteffY = FlashMFieldGapY;
  FlashMFieldLength = FlashMagnetCoreZ;

  IPMAgnetBeamPipeXGap = 1.0 *mm;
  VacChamberMaterial = BeamPipeMaterial;
  VacChamberSideMaterial = BeamPipeMaterial;
  VacChamberWindowMaterial = "Aluminium"; //"G4_KAPTON";
  VacChamberXZWallThickness = 20.0 *mm;
  VacChamberDetXZWidth = 580.0*mm; //according to 3D CAD //700.0*mm; //OPPPTrackerX + OPPPDetXPos;
  VacChamberMagXZWidth = FlashMFieldX - 2.0*IPMAgnetBeamPipeXGap; //TypMBMagnetEndCapHoleOutX;
  VacChambertoOPPPDetZGap = 20.0 *mm;
  VacChamberHight = TypMBMagnetHoleY; //FlashMFieldGapY;
  VacChamberWindowThickness = 0.30 *mm; //2.0 *mm;
  VacChamberWindowX = 500.0 *mm;
  VacChamberWindowY = 20.0 *mm;
  VacChamberSideWallThickness = 20.0 *mm;

  VacChamberFrontThickness = 20.0 *mm;
  VacChamberFrontCutX = 500.0 *mm;
  VacChamberFrontCutY = 20.0 *mm;
  VacChamberFrontCutXpos = 45.0 *mm; //distance from the cut rectangular edge to the beam axis
  VacChamberWindowPanelX = 556.0 *mm;
  VacChamberWindowPanelY = 76.0 *mm;
  VacChamberWindowPanelZ = 10.0 *mm;
  VacChamberWindowPanelXpos = 0.5 * 16.5 *mm;
  VacChamberWindowPanelCutY = 20.0 *mm;
  VacChamberWindowPanelCutX = 500.0 *mm;
  VacChamberWindowPanelCutXpos = 28 *mm; //distance from the cut rectangular edge to the edge of the panel
  VacChamberReinforceBarX = 0.5*65.0 *mm;
  VacChamberReinforceBarZ = 30.0 *mm;
  VacChamberReinforceFilletXZ = 15.0 *mm;
  VacChamberFrameLX = 30.0 *mm;
  VacChamberFrameLY = 60.0 *mm;
  VacChamberFrameHX = 305.0 *mm;
  VacChamberFrameHY = 60.0 *mm;
  VacChamberFrameHZ = 30.0 *mm;

  TAUIChamberMaterial = "StainlessSteel";
  TAUIChamberBottomX = 1620.0 *mm;
  TAUIChamberBottomY = 45.0 *mm;
  TAUIChamberBottomZ = 2120.0 *mm;
  TAUIChamberX = 1500.0 *mm;
  TAUIChamberY = 550.0 *mm;
  TAUIChamberZ = 2000.0 *mm;
  TAUIChamberTopX = TAUIChamberBottomX;
  TAUIChamberTopZ = TAUIChamberBottomZ;
  TAUIChamberTopY = TAUIChamberBottomY;
  TAUIChamberSideThickness = 30.0 *mm;
  TAUIChamberFrontThickness = 30.0 *mm;
  TAUIChamberXpos = 0.5*TAUIChamberX - TAUIChamberSideThickness - 920.0 *mm;
  TAUIChamberYpos = -0.5*TAUIChamberY - TAUIChamberBottomY + 270.0 *mm;
  TAUIChamberZpos = 0.5*TAUIChamberZ - TAUIChamberSideThickness - 807.52 *mm;
  TAUIChamberBeamBottomY = TAUIChamberY/2.0 + TAUIChamberBottomY - TAUIChamberYpos;
  TAUIChamberSupportR = 0.5* 168.30 *mm;
  TAUIChamberSupportH = -TAUIChamberBeamBottomY - FloorSurfaceYpos;
  TAUIChamberPipeLength = (100.0 + 20.0)*mm;
  TAUIChamberPipeFlangedH = 65.0 *mm;
  TAUIChamberPipeFlangedL = 20.0 *mm;
  IPContainerZ = TAUIChamberZ + 2.0*TAUIChamberPipeLength;
  TAUIChamberElPipeRIn = 100.0 *mm;
  TAUIChamberElPipeROut = 115.0 *mm;

  TAUIChamberBBoardX = TAUIChamberX - 2.0*TAUIChamberSideThickness - 100.0 *mm;
  TAUIChamberBBoardY = 12.6 *mm;
  TAUIChamberBBoardZ = TAUIChamberZ - 2.0*TAUIChamberSideThickness - 100.0 *mm;
  TAUIChamberBBoardGapY = 10.0 *mm;
  TAUIChamberBBoardMaterial = "StainlessSteel";
  TAUIChamberMirrorMaterial = "G4_SILICON_DIOXIDE";

  TAUIChamberMirrorR = {110.0*mm, 110.0*mm, 110.0*mm, 110.0*mm, 82.5*mm, 82.5*mm};
  TAUIChamberMirrorD = {28.0*mm, 28.0*mm, 28.0*mm, 28.0*mm, 21.0*mm, 21.0*mm};
  TAUIChamberMirrorTheta = {5.0*M_PI/4.0, 5.0*M_PI/4.0, -5.0*M_PI/4.0, M_PI/4.0, M_PI, 0.0};
  TAUIChamberMirrorXpos = {-720.0*mm, -370.0*mm, -145.0*mm, 165.0*mm, -160.0*mm, 160.0*mm};
  TAUIChamberMirrorZpos = {-525.0*mm, -225.0*mm, -525.0*mm, -225.0*mm, 505.0*mm, -505.0*mm};
  TAUIChamberMirrorName = {"_0", "_1", "_2", "_3", "_4", "_5"};
  TAUIChamberMirrorHolderMaterial = "Aluminium";
  TAUIChamberMirrorRingH = 10.0 *mm;
  TAUIChamberMirrorRingMaterial = "Aluminium";
  ICBBoardTableGapY = 155.0 *mm;
  ICBBoardTableX = 1220.0 *mm;
  ICBBoardTableY = 555.0 *mm;
  ICBBoardTableZ = 1600.0 *mm;
  ICBBoardTableLegR = 145.0 *mm;

  IPVolumeZ = 0.025 *mm;
  IPVolumeR = 50.0 *mm;
  
  OPPPDetOffsetZ = 5.0 *mm;  // offset of the position of the first stave with respect of the supporting plate. 
  OPPPDetSupportMaterial = "StainlessSteel";
  OPPPDetBottomSupportX = 580.0 *mm;
  OPPPDetBottomSupportY = 20.0 *mm;
  OPPPDetBottomSupportZ = 560.0 *mm;
  OPPPDetSupportBallH = 28.28 *mm;
  OPPPDetTopPlateX = 580.0 *mm;
  OPPPDetTopPlateY = 15.0 *mm;
  OPPPDetTopPlateZ = 360.0 *mm;
  OPPPAngleHolderX = 8.0 *mm;
  OPPPAngleHolderY = 20.0 *mm;
  OPPPAngleHolderZ = 25.0 *mm;
  OPPPAngleHolderH = 5.0 *mm;
  OPPPAngleHolderT = 3.0 *mm;
  OPPPAngleHolderBaseX = 13.0 *mm;
  OPPPAngleHolderBaseY = 10.0 *mm;
  OPPPAngleHolderBaseT = 5.0 *mm;
  OPPPAngleHolderLCutX = 4 *mm;
  OPPPAngleHolderLCutY = 8 *mm;
  OPPPAngleHolderRCutR = 2 *mm;
  OPPPAngleHolderTopCutXY = 1.5 *mm;

  OPPPColdPlateMaterial = "ColdPlateMaterial";
  OPPPFrameMaterial = "CarbonFiber";
  OPPPStaveSideMaterial = "StaveEndCapMaterial";
  OPPPStaveHolderMaterial = "Aluminium";
  
  OPPPColdPlateCarbonFiberThikness = 0.165 *mm;
  OPPPColdPlateEpoxyThikness = 0.100 *mm;
  
  OPPPCoolingPipeR = (0.5*1.024 + 0.025) *mm;
  OPPPCoolingPipeWall = 0.025 *mm;
  OPPPCoolingPipeYpos = 3.0 *mm;
  OPPPCoolingPipeMaterial = "G4_KAPTON";
  OPPPCoolingMaterial = "G4_WATER";
  OPPPColdPlateX = 290.0 *mm;
  OPPPColdPlateY = 15.0 *mm;
  OPPPTrussX = 15.0 *mm;
  OPPPTrussZ = 4.9 *mm;
  OPPPTrussTopXY = 0.7 *mm;
  OPPPTrussTopThickness = 0.3 *mm;
  OPPPNTruss = 19;
  OPPPStaveSideZ = 2.0 *mm;
  OPPPStaveSideZpos = 0.2 *mm;
  OPPPStaveLiftY = 10.0 *mm;
  OPPPStaveX = 318.0 *mm;
  OPPPStaveInOutZ = 12.0 *mm;
  OPPPStaveInOutX = 40.0 *mm;
  OPPPTrackerInterLayerZ = 100.0 *mm;
  OPPPSesorFPCMaterial1 = "G4_KAPTON";
  OPPPSesorFPCMaterial2 = "Aluminium";
  OPPPSesorFPCLayer1Z = 0.1 *mm;
  OPPPSesorFPCLayer2Z = 0.05 *mm;
  OPPPNStaveSensors = 9;
  OPPPStaveSensorsGapX = 0.1 *mm;
  OPPPSensorX = 30.0 *mm;
  OPPPSensorY = 15.0 *mm;
  OPPPSensorZ = 0.05 *mm;
  OPPPSensorNCellX = 1024;
  OPPPSensorNCellY = 512;
  OPPPSensorPixelX = 0.02924 *mm;
  OPPPSensorPixelY = 0.02688 *mm;
  OPPPSensorPixelZ = 0.025 *mm;
  OPPPTrackerActiveX =  OPPPNStaveSensors * (OPPPSensorX + OPPPStaveSensorsGapX) - OPPPStaveSensorsGapX; //270.8 *mm;

  OPPPHexapodMaterial = "StainlessSteel";
  OPPPHexapodUpR = 250.0/2.0 *mm;
  OPPPHexapodUpH = 21.5 *mm;
  OPPPHexapodDownR = 348.0/2.0 *mm;
  OPPPHexapodDownH = 30.0 *mm;
  OPPPHexapodY = 328.0 *mm;
  OPPPHexapodThickness = 10.0 *mm;

  OPPPBasePlateMaterial = "StainlessSteel";
  OPPPBasePlateX = 2.0 * (OPPPDetBottomSupportX + OPPPDetXPos);
  OPPPBasePlateY = 100.0 *mm;
  OPPPBasePlateZ = OPPPDetBottomSupportZ;
  OPPPBasePlateR = 50.0 *mm;

  OPPPServiceSupportX = 72.0 *mm;
  OPPPServiceSupportY = OPPPDetTopPlateY;
  OPPPServiceSupportZ = 410.0 *mm;
  OPPPServiceSupportTubeX = 181.0 *mm;
  OPPPServiceSupportTubeRin = 4.5 *mm;
  OPPPServiceSupportTubeThickness = 1.5 *mm;
  OPPPServiceSupportTubeMaterial = "StainlessSteel";

  OPPPServiceCoolPipeRout = 1.0 *mm;
  OPPPServiceCoolPipeThick = 0.35 *mm;
  OPPPServiceCoolPipeYdist = 5.0 *mm;
  OPPPServiceCoolPipeCnctrL = 4.0 *mm;
  OPPPServiceCoolPipeCnctrRin = 0.5 *mm;
  OPPPServiceCoolPipeCnctrRout = 1.0 *mm;
  OPPPServiceCoolPipeCnctrSL = 6.0 *mm;
  OPPPServiceCoolPipeCnctrSRin = 0.5 *mm;
  OPPPServiceCoolPipeCnctrSRout = OPPPServiceCoolPipeRout - OPPPServiceCoolPipeThick;
  OPPPServiceCoolPipeCnctrMaterial = "StaveEndCapMaterial";
  OPPPServiceCoolPipeMaterial = "G4_POLYVINYL_CHLORIDE";
  OPPPInnerServiceCoolPipeDZ = 30.0 *mm;
  OPPPInnerServiceCoolPipeBentR = 300.0 *mm;

  OPPPServiceCablMaterial = "G4_KAPTON";
  OPPPInnerServiceCablSec0X = 50.0 *mm;
  OPPPInnerServiceCablSec0Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec0Z = 2.0 *mm;
  OPPPInnerServiceCablSec0Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  OPPPInnerServiceCablSec1X = 50.0 *mm;
  OPPPInnerServiceCablSec1Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec1Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  OPPPInnerServiceCablSec20X = 80.0 *mm;
  OPPPInnerServiceCablSec20Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec20Thick = OPPPSensorZ;
  OPPPInnerServiceCablSec21X = 80.0 *mm;
  OPPPInnerServiceCablSec21Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec21Z = 17.0 *mm;
  OPPPInnerServiceCablSec21Thick = OPPPSesorFPCLayer1Z;
  OPPPInnerServiceCablSec22X = 80.0 *mm;
  OPPPInnerServiceCablSec22Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec22Z = 20.0 *mm;
  OPPPInnerServiceCablSec22Thick = OPPPSesorFPCLayer2Z;
  OPPPInnerServiceCablSec31Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec31Thick = OPPPSesorFPCLayer1Z;
  OPPPInnerServiceCablSec32Y = OPPPColdPlateY;
  OPPPInnerServiceCablSec32Thick = OPPPSesorFPCLayer2Z;

  OPPPOuterServiceCablSec0X = 40.0 *mm;
  OPPPOuterServiceCablSec0Y = OPPPColdPlateY;
  OPPPOuterServiceCablSec0Z = 3.0 *mm;
  OPPPOuterServiceCablSec0Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  OPPPOuterServiceCablSec1Y = OPPPColdPlateY;
  OPPPOuterServiceCablSec1Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;

  OPPPSideOutServiceCablSec0X = 50.0 *mm;
  OPPPSideOutServiceCablSec0Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec0Z = 3.0 *mm;
  OPPPSideOutServiceCablSec0Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  OPPPSideOutServiceCablSec1X = 30.0 *mm;
  OPPPSideOutServiceCablSec1Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec1Thick = OPPPSensorZ + OPPPSesorFPCLayer1Z + OPPPSesorFPCLayer2Z;
  OPPPSideOutServiceCablSec20X = 145.0 *mm;
  OPPPSideOutServiceCablSec20Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec20Thick = OPPPSensorZ;
  OPPPSideOutServiceCablSec21X = 50.0 *mm;
  OPPPSideOutServiceCablSec21Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec21Z = 3.0 *mm;
  OPPPSideOutServiceCablSec21Thick = OPPPSesorFPCLayer1Z;
  OPPPSideOutServiceCablSec22X = 50.0 *mm;
  OPPPSideOutServiceCablSec22Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec22Z = 15.0 *mm;
  OPPPSideOutServiceCablSec22Thick = OPPPSesorFPCLayer2Z;
  OPPPSideOutServiceCablSec31X = 95.0 *mm;
  OPPPSideOutServiceCablSec31Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec31Thick = OPPPSesorFPCLayer1Z;
  OPPPSideOutServiceCablSec32X = 80.0 *mm;
  OPPPSideOutServiceCablSec32Y = OPPPColdPlateY;
  OPPPSideOutServiceCablSec32Thick = OPPPSesorFPCLayer2Z;

  OPPPSideInServiceCablSec01X = 185.0 *mm;
  OPPPSideInServiceCablSec01Y = OPPPColdPlateY;
  OPPPSideInServiceCablSec01Thick = OPPPSesorFPCLayer1Z;
  OPPPSideInServiceCablSec02X = 185.0 *mm;
  OPPPSideInServiceCablSec02Y = OPPPColdPlateY;
  OPPPSideInServiceCablSec02Thick = OPPPSesorFPCLayer2Z;

  OPPPServiceCablTermBase1X = 25.0 *mm;
  OPPPServiceCablTermBase1Z = 12.0 *mm;
  OPPPServiceCablTermMid1X = 25.0 *mm;
  OPPPServiceCablTermMid1Y = OPPPColdPlateY;
  OPPPServiceCablTermMid1Z = 5.0 *mm;
  OPPPServiceCablTermBase2X = 25.0 *mm;
  OPPPServiceCablTermBase2Z = 25.0 *mm;
  OPPPServiceCablTermMid2X = 25.0 *mm;
  OPPPServiceCablTermMid2Y = OPPPColdPlateY;
  OPPPServiceCablTermMid2Z = 5.0 *mm;
  OPPPServiceCablTermMaterial = "Aluminium";

  TrackerElectronicsRackMaterial = "Lead";
  TrackerElectronicsRackX = 600.0 *mm;
  TrackerElectronicsRackY = 900.0 *mm;
  TrackerElectronicsRackZ = 600.0 *mm;
  TrackerElectronicsRackThick = 50.0 *mm;
  TrackerElectronicsPCBThick = 3.0 *mm;
  TrackerElectronicsRackXpos = 100.0 *mm;  //distance in x to the tracker supporting plate

  ECalX = 550.0 *mm;
  ECalY = 55.0 *mm;
  ECalZ = 200.0 *mm;
  ECalNLayers = 21;
  ECalSensorPixelX = 5.0 *mm;
  ECalSensorPixelY = 5.0 *mm;
  ECalSensorNCellX = static_cast<G4int>(ECalX/ECalSensorPixelX + 0.5);
  ECalSensorNCellY = static_cast<G4int>(ECalY/ECalSensorPixelY + 0.5);
  OPPPTrackerECalZ = 30.0 *mm;
  ECalTopPlateZ = OPPPDetBottomSupportZ - OPPPDetTopPlateZ - OPPPTrackerECalZ;
  ECalXpos = BPipeR + 15.0*mm;   //+ 5.0*mm; // actually the distance from the beam axis to the detector
  ECalCasingTopY = 10.0 *mm;
  ECalCasingSideX = 10.0 *mm;
  ECalCasingBackZ = 10.0 *mm;
  ECalCasingGapZ = 20.0 *mm;
  ECalCasingMaterial = "StainlessSteel";
  ECalPipeShieldX = 5.0 *mm;
  ECalPipeShieldY = 6.0 * BPipeR;
  ECalPipeShieldZ = 1000.0 *mm;
  ECalPipeShieldShiftZ = 35.0 *mm < VacChambertoOPPPDetZGap ? 35.0 *mm : VacChambertoOPPPDetZGap/2.0 ;
  ECalPipeShieldMaterial = "Tungsten"; //"Lead";

  OPPPECalPCBX = ECalX;
  OPPPECalPCBY = 30.0 *mm;
  OPPPECalPCBZ = 3.0 *mm;
  OPPPECalPCBMaterial = "ECAL_FR4";

  ECalDumpShieldMaterialInner = "Iron";
  ECalDumpShieldMaterialOuter = "G4_POLYETHYLENE";
  ECalDumpShieldX = OPPPDetTopPlateX;
  ECalDumpShieldY = 300.0 *mm;
  ECalDumpShieldInnerZ = 400.0 *mm;
  ECalDumpShieldOuterZ = 150.0 *mm;
  ECalDumpShieldXpos = OPPPDetXPos;
  ECalDumpShieldYpos = 0.0;
  ECalDumpShieldZpos = 200.0 *mm;  // actuall the distance form the ECal to the

  HICSDumpFrontZPos = 7.0 *m;
  HICSDumpFrontXPos = (9.0 + 1.2) *cm;  // it is absolut value, actually x<0 but it is accounted in construction.
  HICSDumpR = 15.0 *cm;
  HICSDumpZ = 75.0 *cm;
  HICSDumpAlInsertR = 4.0 *cm;
  HICSDumpAlInsertZ = 30.0 *cm;
  HICSDumpAirInsertR = 2.0 *cm;
  HICSDumpAirInsertZ = 15.0 *cm;
  HICSDumpMaterial = "Lead"; //"Copper";
  HICSDumpAlInsertMaterial = "Aluminium";
  HICSDumpProductionCut = 1.0 *cm;
  HICSDump2TPosX = 9.5 *cm;  // This si displacement of the center of the beam dump in case 2T field with respectr to 1T

  HICSShieldingSideMaterial = "Aluminium";
  HICSShieldingMiddleMaterial = "Lead"; //"Copper";
  HICSShieldingTopMaterial = "Iron";
  HICSShieldingSupportMaterial = "ShieldingConcrete";
  HICSShieldingXPos = -60.0 *cm;
  HICSShieldingX = 2.0 *m;
//   HICSShieldingY = 1.0 *m;
  HICSShieldingY = 2.4 * HICSDumpR;
  HICSShieldingSideZ = 20.0 *cm;
  HICSShieldingMiddleZ = 30.0 *cm;
  HICSShieldingGapX = 1.0 *m;
  HICSShieldingGapY = 2.0 *cm;
  HICSShieldingSupportZ = 1.5 *m;
  HICSShieldingSupportY = 3.0 *m;

  HICSNeutronAbsorberMaterial = "BoratedPolyethyleneConcrete";
  HICSNeutronAbsorberX = HICSShieldingX;
  HICSNeutronAbsorberY = HICSShieldingY;
  HICSNeutronAbsorberZ = 170.0 *mm;

  //HICSScintilatorXPos = 47 *mm;   // Displacement at 1.49 m from the 1T magnet (m): 0.0374821
  //HICSScintilatorXPos = 42 *mm;   // Displacement at 1.3 m from the magnet (m): 0.0339292
  //HICSScintilatorXPos = 80 *mm;     // Displacement required for 2T magnet to avoid 16.5 GeV beam 

  ComptonElBackshift =  OPPPDetZtoMagnet + 1.5*m; // defined from edge of B-field to start of Scintillator 

  ScintCameraYpos = 45.*cm;
//   ScintCameraZpos = IPMagnetZpos + TypMBFieldLength/2.0 + ComptonElBackshift
//                     + (ScintBaseZ + ScintPhosphorZ + ScintFinishZ - OPPPDetTopPlateZ)/2.0 + 1.*m;
  ScintCameraZpos = IPMagnetZpos + FlashMFieldLength/2.0 + ComptonElBackshift
                    + (ScintBaseZ + ScintPhosphorZ + ScintFinishZ - OPPPDetTopPlateZ)/2.0 + 1.*m;

  HICSScintilatorXPosMag1T = (36+30+10) *mm;   // Displacement required for 1T magnet to accept ~15.5 GeV
  // HICSScintilatorXPosMag2T = 70 *mm;     // Displacement required for 2T magnet to avoid 16.5 GeV beam - no longer true after 1.5m backshift, no longer a useful variable
  
  HICSCherenkovXPosMag1T = (20.25 + 30 + 10)*mm;     // Displacement at 1.69 + 1.5 m from the 1T magnet (m)
  // HICSCherenkovXPosMag2T = 68.25 *mm;     // Displacement at 1.69 m for 2T magnet

  ComptonElectronBeamtoStageY = 15 * cm;

  HICSScintSupportBallH = 50 *mm; // 28.28 *mm;
  
  HICSDetSupportMaterial = "StainlessSteel";

  HICSDetBottomSupportX = 2.0 * 580.0 *mm;
  HICSDetBottomSupportY = 20.0 *mm;
  HICSDetBottomSupportZ = 625. * mm; // 560.0 *mm;
   
  HICSBasePlateX = 2.0 * (OPPPDetBottomSupportX + OPPPDetXPos);
  HICSBasePlateY = 100.0 *mm;
  HICSBasePlateZ = OPPPDetBottomSupportZ;
  HICSBasePlateR = 50.0 *mm;

  BremBasePlateMaterial = "StainlessSteel";
  BremBasePlateX = 2.0 * (OPPPDetBottomSupportX + OPPPDetXPos);
  BremBasePlateY = 100.0 *mm;
  BremBasePlateZ = 2 * OPPPDetBottomSupportZ;
  BremBasePlateR = 50.0 *mm;

  HICSBasePlateMaterial = "StainlessSteel";

  BeamProfilerContanerMaterial = "Air20";
  BeamProfilerContanerX = 200.0 *mm;
  BeamProfilerContanerY = 150.0 *mm;
  BeamProfilerContanerZ = 150.0 *mm;
  BeamProfilerContanerD = 3.0 *mm;
  BeamProfilerContanerWallMaterial = "Aluminium";
  BeamProfilerWindowX = 20.0 *mm;
  BeamProfilerWindowY = 20.0 *mm;
  BeamProfilerWindowZ = 0.005 *mm;
  BeamProfilerWindowMaterial = "G4_KAPTON";
  
  BeamProfilerSensorX = 25.0 *mm;
  BeamProfilerSensorY = 25.0 *mm;
  BeamProfilerSensorZ = 0.1 *mm;
  BeamProfilerMetalizationZ = 0.005 *mm;
  G4double bp2zpos = 11800*mm; //GammaCaloZpos - 0.5*GammaCaloZ - 100*mm;
  BeamProfilerZpos =  {bp2zpos, bp2zpos + 15.0*mm}; //{6700.0*mm, 6710.0*mm, bp2zpos, bp2zpos + 10.0*mm};
  BeamProfilerSensorMaterial = "G4_ALUMINUM_OXIDE";  //Sapphire
  BeamProfilerNCellX = 5000;  //5um
  BeamProfilerNCellY = 5000;  //5um
  BeamPipe2ProfilerZ = bp2zpos - 100.0*mm - ComptonLysoZpos;
  BeamProfilerSensorMetalization = "Aluminium";

  BeamProfilerV408Material = "StainlessSteel";
  BeamProfilerV408X = 105.0 *mm;
  BeamProfilerV408Y = 25.0 *mm;
  BeamProfilerV408Z = 80.0 *mm;
  BeamProfilerMotorGap = 40 *mm;
  BeamProfilerQ545Material = "StainlessSteel";
  BeamProfilerQ545X = 45.0 *mm;
  BeamProfilerQ545Y = 63.0 *mm;
  BeamProfilerQ545Z = 15.0 *mm;
  BeamProfilerPCBX = 60.0 *mm;
  BeamProfilerPCBY = 60.0 *mm;
  BeamProfilerPCBZ = 3.0 *mm;
  BeamProfilerPCBCutX = 22.0 *mm;
  BeamProfilerPCBCutY = 22.0 *mm;
  BeamProfilerPCBHolderX = 50.0 *mm;
  BeamProfilerPCBHolderY = 50.0 *mm;
  BeamProfilerPCBHolderZ = 7.0 *mm;
  BeamProfilerPCBHolderCutR = 15.0 *mm;
  BeamProfilerPCBHolderMaterial = "Aluminium";

  BeamProfilerAngleSupportX = 30.0 *mm;
  BeamProfilerAngleSupportD = 5.0 *mm;
  BeamProfilerAngleSupportMaterial = "Aluminium";

  GamVacChamberGap = 20.0 *mm;
  BeamPipeLidMaterial = "G4_KAPTON";
  BeamPipeLidThickness = 0.2 *mm;

  BSMCaloMaterial = "G4_PbWO4";
  BSMCaloX = 2005.0 *mm;
  BSMCaloY = 2005.0 *mm;
  BSMCaloLayerZ = 0.5 * 8.903 *mm; // X0 = 8.903 mm
  BSMCaloNLayers = 40;
  BSMCaloZPos = 16.0 *m;
  BSMCaloNCellX = 401;
  BSMCaloNCellY = 401;

////// Telescope //////////
  WoodSupportZpos = 165.0 *mm;
  WoodSupportZ = 18.0 *mm;
  TelescopeFloorZpos = WoodSupportZpos + WoodSupportZ + 45.0 *mm;
  FloorZ = 100.0 *mm;
  FloorGapZ = 450.0 *mm;

  Collimator1X = 238.0*mm;
  Collimator1Y = 100.0*mm;
  Collimator1Z = 10.0*mm;
  Collimator1HoleR = 1.0*mm;

  SourceContainerRin = 20.0*mm;
  SourceContainerRout = 30.0*mm;
  SourceContainerZ = 100.0*mm;
  SourceContainerTopZ = 20.0*mm;

  Sr90SourceR = 20.4/2.0 *mm;
  Sr90SourceZ = 1.0 *mm;
  Sr90GapZ = 13.0 *mm;
  SrAlPlugR = Sr90SourceR;
  SrAlPlugZ = 2.0 *mm;
  SrSupportRingRin = Sr90SourceR;
  SrSupportRingRout = Sr90SourceR + 1.0 *mm;
  SrSupportRingZ = 3.18 *mm;

  CollimatorLeadZ = 10.0*mm;
  CollimatorLeadHoleR = 6.3/2.0 *mm;

  Collimator2X = 90.0 *mm;
  Collimator2Y = 90.0 *mm;

  CarrierPCBX = 80.0 *mm;
  CarrierPCBY = 70.0 *mm;
  CarrierPCBZ = 1.5 *mm;
  CarrierPCBCutX = OPPPSensorX - 1.0*mm;
  CarrierPCBCutY = OPPPSensorY - 1.0*mm;

  ROPCBX = 100.0 *mm;
  ROPCBY = 100.0 *mm;
  ROPCBZ = 1.5 *mm;
  CarrierROGap = 17.0 *mm;

  TeleFrameX = 150.0 *mm;
  TeleFrameY = 220.0 *mm;
  TeleFrameZ = 8.0 *mm;
  TeleFrameCutAlpideX = 50.0 *mm;
  TeleFrameCutAlpideY = 50.0 *mm;
  TeleFrameCutAlpideGapY = 15.5 *mm;
  TeleFramePCBGap = 27.0 *mm;
  TeleFramePCBShift = 6.0 *mm;
  TeleFrameHolderY = 36.0 *mm;


  TelescopeSensorZpos = {29.0*mm, (29.0+25.8)*mm, (29.0 + 2.0*25.8)*mm};

  TelescopeShieldingMaterial = "G4_Pb";
  TeleShieldingX = 360.0 *mm;
  TeleShieldingY = 465.0 *mm;
  TeleShieldingZ = 180.0 *mm;
  TeleShieldingThickess = 50.0 *mm;
  ShieldingGapY = 70.0 *mm;

  FrameHolderX = 52.0 *mm;
  FrameHolderY = 20.0 *mm;
  FrameHolderZ = 330.0 *mm;
  SupportPoleX = 60.0 *mm;
  SupportPoleY = 60.0 *mm;
  SupportPoleZ = 1100.0 *mm;
  SupportPoleThickness = 2.0 *mm;

  ShieldingBoxBottomX = 350.0 *mm;
  ShieldingBoxBottomY = 6.0 *mm;
  ShieldingBoxBottomZ = 450.0 *mm;

}



LXSetUp *LXSetUp::Instance() 
{
  if (!flxSetup) {
    flxSetup = new LXSetUp();
  }
  return flxSetup;
}




