
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

//  CeilingSurfaceYpos = 1816.0 *mm; // approx. for Scint screen camera (IP) ceiling mount. downstream (lower) part of chamber

  BPipeR = 24.13 *mm;          // DN40;

  BeamCrossingAngle = 0.3 *rad;
  BeamPipeVacuumMaterial = "XFELVacuum"; //"Galactic";

  BeamDumpProductionCut = 5.0 *cm;

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

  BTargetX = 10.0 *cm;
  BTargetY = 10.0 *cm;
  BTargetZ = 15.0*mm; //0.035 *mm;

  OPPPSensorX = 30.0 *mm;
  OPPPSensorY = 15.0 *mm;
  OPPPSensorZ = 0.05 *mm;
  OPPPSensorNCellX = 1024;
  OPPPSensorNCellY = 512;
  OPPPSensorPixelX = 0.02924 *mm;
  OPPPSensorPixelY = 0.02688 *mm;
  OPPPSensorPixelZ = 0.025 *mm;

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




