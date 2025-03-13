
#ifndef LXSETUP_h
#define LXSETUP_h 1

#include "globals.hh"
#include "G4Cache.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


class LXSetUp
{
public:
  static LXSetUp *Instance();
  ~LXSetUp();

protected:
  LXSetUp();


  static LXSetUp *flxSetup;

public:
  enum tTargetType : G4int {tfoil, twire};

public:
  G4bool OverlapTest;

  G4String EnvironmentMaterial;
  G4double FloorSurfaceYpos;
  G4double FloorY;

  G4double CeilingSurfaceYpos;

  G4double BPipeR;
  G4String BeamPipeVacuumMaterial;

  G4double BeamDumpProductionCut;
  G4double BeamCrossingAngle;


//   G4double OPPPDetZpos;
  G4double OPPPDetZtoMagnet;
  G4double OPPPDetXPos;
  G4double OPPPDetYLowPos;
  G4double OPPPDetX;
  G4double OPPPDetY;
  G4double OPPPDetZ;
  G4String OPPPDetMaterial;
  
  G4double OPPPTrackerX;
  G4double OPPPTrackerY;
  G4double OPPPTrackerZ;
  G4double OPPPTrackerDZ;
  G4int OPPPTrackerNLayers;
  G4String OPPPTrackerMaterial;
  G4int OPPPTrackerNCellX;
  G4int OPPPTrackerNCellY;

  G4double BTargetZpos;
  G4double BTargetX;
  G4double BTargetY;
  G4double BTargetZ;

  G4double OPPPSensorX;
  G4double OPPPSensorY;
  G4double OPPPSensorZ;
  G4double OPPPSensorNCellX;
  G4double OPPPSensorNCellY;
  G4double OPPPSensorPixelX;
  G4double OPPPSensorPixelY;
  G4double OPPPSensorPixelZ;

////// Telescope //////////
  G4double WoodSupportZpos;
  G4double WoodSupportZ;
  G4double TelescopeFloorZpos;
  G4double FloorZ;
  G4double FloorGapZ;

  G4double Collimator1X;
  G4double Collimator1Y;
  G4double Collimator1Z;
  G4double Collimator1HoleR;
  G4double SourceContainerTopZ;

  G4double SourceContainerRin;
  G4double SourceContainerRout;
  G4double SourceContainerZ;

  G4double Sr90SourceR;
  G4double Sr90SourceZ;
  G4double Sr90GapZ;
  G4double SrAlPlugR;
  G4double SrAlPlugZ;
  G4double SrSupportRingRin;
  G4double SrSupportRingRout;
  G4double SrSupportRingZ;

  G4double CollimatorLeadZ;
  G4double CollimatorLeadHoleR;

  G4double Collimator2X;
  G4double Collimator2Y;

  G4double CarrierPCBX;
  G4double CarrierPCBY;
  G4double CarrierPCBZ;
  G4double CarrierPCBCutX;
  G4double CarrierPCBCutY;
  G4double ROPCBX;
  G4double ROPCBY;
  G4double ROPCBZ;
  G4double CarrierROGap;

  G4double TeleFrameX;
  G4double TeleFrameY;
  G4double TeleFrameZ;
  G4double TeleFrameCutAlpideX;
  G4double TeleFrameCutAlpideY;
  G4double TeleFrameCutAlpideGapY;
  G4double TeleFramePCBGap;
  G4double TeleFramePCBShift;
  G4double TeleFrameHolderY;

  std::vector<G4double> TelescopeSensorZpos;

  G4String TelescopeShieldingMaterial;
  G4double TeleShieldingX;
  G4double TeleShieldingY;
  G4double TeleShieldingZ;
  G4double TeleShieldingThickess;
  G4double ShieldingGapY;

  G4double FrameHolderX;
  G4double FrameHolderY;
  G4double FrameHolderZ;
  G4double SupportPoleX;
  G4double SupportPoleY;
  G4double SupportPoleZ;
  G4double SupportPoleThickness;

  G4double ShieldingBoxBottomX;
  G4double ShieldingBoxBottomY;
  G4double ShieldingBoxBottomZ;

};

#endif
