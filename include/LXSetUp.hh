
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

  void SetGammaBeamDumpZ(const G4double size_z)
  {
    GammaBeamDumpZpos -= 0.5 * (GammaBeamDumpZ - size_z);    // keep front face in nomimal place
    GammaBeamDumpZ = size_z;
    GammaCaloZpos = GammaBeamDumpZpos - 0.5*(GammaBeamDumpZ + GammaCaloZ) - 100.0*mm;
    CAbsorberZpos = GammaBeamDumpZpos - 0.5*(GammaBeamDumpZ + CAbsorberThickness) - 20.0*mm;
  }

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

  G4double IPBoxX;
  G4double IPBoxY;
  G4double IPContainerZ;
  G4double IPBoxZ;
  G4double IPBoxThickness;

  G4double BPipeR;
  G4double BPipeThickness;
  G4String BeamPipeMaterial;
  G4String BeamPipeWindowMaterial;
  G4String BeamPipeVacuumMaterial;

  G4double LaserPipeR;
  G4double LaserPipeThickness;

  G4double BeamCrossingAngle;

  G4String MagnetMaterial;
  G4double DumpMagnetThickness;
  G4double DumpMagnetXpos;
  G4double DumpMagnetYpos;
  G4double DumpMagnetZpos;
  G4double DumpMagnetX;
  G4double DumpMagnetY;
  G4double DumpMagnetZ;
  G4double DMBPipeWindow;

  G4bool ScintCerenkovPhysics;
  G4double ScintPhysicsYield;
  G4double ScintAngle;
  G4double ScintXpos;
  G4double ScintZpos;
  G4double ScintFrameThickness;
  G4String ScintFrameMaterial;
  G4double ScintFrameBeamLoopX;
  G4double BremScintFrameWidthOffset;

  G4double CerenkovAngle;
  G4double CerenkovXpos;
  G4double CerenkovZpos;
  G4double ScintX;
  G4double ScintY;
  G4double ScintBaseZ;
  G4double ScintPhosphorZ;
  G4double ScintFinishZ;

  // G4double CerenkovWallWidth;
  G4double CerenkovBoxThickness;
  // G4double CerenkovWindowThickness;
  // G4double CerenkovMirrorThickness;
  // G4double CerenkovDeviceOuterThickness;
  G4double CerenkovBeamWindowThickness;
  // G4double CerenkovchannelX;
   G4double CerenkovchannelY;
  // G4double CerenkovchannelZ;
  G4int CerenkovChannels;
  G4double CerenkovStrawHeight;
  G4int CerenkovChannelLayers;
  G4double CerenkovStrawInnerRadius;
  G4double CerenkovStrawGraphiteLayerThickness;
  G4double CerenkovStrawAlLayerThickness;
  G4double CerenkovStrawKaptonLayerThickness;
  G4double CerenkovStrawPolyurethaneLayerThickness;
  G4String CerenkovStrawGraphiteLayerMaterial;
  G4String CerenkovStrawAlLayerMaterial;
  G4String CerenkovStrawKaptonLayerMaterial;
  G4String CerenkovStrawPolyurethaneLayerMaterial;
  G4double CerenkovShieldingPlateThickness;

  G4double CerenkovStrawXFrequency;
  G4double CerenkovStrawZFrequency;
  G4double CerenkovStrawLayerOffset;
  G4double CerenkovStrawPlateBuffer;
  G4double CerenkovElectronicsBoardThickness; 
  G4double CerenkovSupportThickness;
  G4double CerenkovSupportPosYtoStage;
  G4double CerenkovBeamWindowY;
  G4double CerenkovTotalBoxHeight;  G4double CerenkovTotalBoxLength;
  G4double CerenkovUpperSupporttoElectronicsBoard;
  
  G4double CerenkovLegHeight;
  G4String ScintBaseMaterial;
  G4String ScintPhosphorMaterial;
  G4String ScintFinishMaterial;
  
  G4String CerenkovMetal;
  G4String CerenkovMedium;
  G4String CerenkovElectronicsBoardMaterial;
  G4String CerenkovShieldingPlateMaterial;
  G4double ScintCameraX;
  G4double ScintCameraY;
  G4double ScintCameraZ;
  G4double ScintCameraAngle;
  G4double ScintCameraApertureZ;
  G4double ScintCameraApertureOuter;
  G4double ScintCameraApertureInner;
  G4double ScintCameraShieldThickness;
  G4double ScintCameraSupportThickness;
  G4double ScintCameraPlatformZ;
  G4bool StickScintScreentoBeamWindow;

  G4String ScintCameraMaterial;
  G4String ScintCameraShieldMaterial;
  G4String ScintCameraApertureMaterial;
  G4String ScintCameraSupportMaterial;
  G4double ScintCameraYpos;
  G4double ScintCameraZpos;

  G4double DumpMagFieldY;
  G4double IPMagFieldY;
  
  G4double BeamDumpXpos;
  G4double BeamDumpYpos;
  G4double BeamDumpZpos;
  G4double BeamDumpR;
  G4double BeamDumpZ;
  G4double BeamDumpAngle;
  G4double DumpBeamHoleRin;
  G4double DumpBeamHoleRout;
  G4String BeamDumpMaterial;
  G4double BeamDumpProductionCut;
  
  G4double BeamDumpFrontZpos;
  G4double BeamDumpFrontXpos;

  G4double BeamDumpInsertR;
  G4double BeamDumpInsertZ;
  G4String BeamDumpInsertMaterial;

  G4double BeamDumpPipeR;
  G4double BeamDumpPipeThickness;
  G4String BeamDumpPipeMaterial;
  
  G4String ShieldingMaterial;
  G4double ShieldingZpos;
  G4double ShieldingX;
  G4double ShieldingY;
  G4double ShieldingZ;
  G4double ShieldingDeepZ;
  G4double ShieldingDXBeamWall;
  G4double ShieldingDeepMargine;
  G4double ShieldingProductionCut;
  G4double ShieldingAbsorberX;
  G4double ShieldingAbsorberZ;
  G4String ShieldingAbsorberMaterial;
  G4double ShieldingAbsorberSlitX;
  G4double ShieldingAbsorberSlitZ;
  G4double ShieldingAbsorberTopZ;
  
  G4double BeamPipeSplitterZ;
  
  G4double IPMagnetXpos;
  G4double IPMagnetYpos;
  G4double IPMagnetZpos;

  G4double GammaDumpZpos;
  G4double GammaDumpZ;
  G4double GammaDumpR;
  G4String GammaDumpMaterial;
  G4double GammaDumpProductionCut;
  
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

  G4double OPPPCaloX;
  G4double OPPPCaloY;
  G4double OPPPCaloZ;
  G4String OPPPCaloMaterial;
  
  G4String BTargetMaterial;
  G4double BTargetZpos;
  G4double BTargetX;
  G4double BTargetY;
  G4double BTargetZ;
  G4double BTargetChamberMagnetGapZ;

  G4double TargetChamberX;
  G4double TargetChamberY;
  G4double TargetChamberZ;
  G4double TargetChamberXPos;
  G4double TargetChamberThickness;
  G4String TargetChamberMaterial;
  G4double ChamberTargetVolX;
  G4double ChamberTargetVolY;
  G4double ChamberTargetVolZ;
  G4String TargetChamberMotorMaterial;
  G4double TargetChamberMotorX;
  G4double TargetChamberMotorY;
  G4double TargetChamberMotorZ;
  G4double TargetFrameD;
  G4double TargetFrameZ;
  G4String TargetFrameMaterial;
  G4double TargetFrameHolderR;
  G4double TargetFrameHolderGroveY;

  G4double DumpCleanMagZpoz;
  G4double DumpCleanMagZ;
  G4double DumpCleanMagX;
  G4double DumpCleanMagY;
  G4double DumpCleanMagThickness;
  G4double DumpCleanMagField;

  G4double GMagnetXpos;
  G4double GMagnetYpos;
  G4double GMagnetZpos;
  G4double GMagnetX;
  G4double GMagnetY;
  G4double GMagnetZ;
  G4double GMagFieldY;
  G4double GMagCoilThickness;
  G4double GMagCoilH;
  G4double GMagnetCutX;
  G4double GMagnetCutY;
  G4double GMagneteffY;
  G4double GMagnetCoreX;
  G4double GMFieldX;
  G4double GMFieldLength;

  G4double GTargetX;
  G4double GTargetY;
  G4double GTargetZ;
  G4double GTargetZpos;
  G4String GTargetMaterial;
  tTargetType GTargetType;
  G4double GTargetContainerX;
  G4double GTargetContainerY;
  G4double GTargetContainerZ;
  G4double QBeamPipeContainerX;
  G4double GammaBPipeWindowThickness;
  G4double WideBeamPipeContainerX;
  G4double GChamberTBWallThickness;
  G4double GammaMagnetBeamPipeXGap;

  G4double ComptonLysoX;
  G4double ComptonLysoY;
  G4double ComptonLysoZ;
//  G4double ComptonLysoDZ;
  G4double ComptonLysoNCellX;
  G4double ComptonLysoNCellY;
  G4String ComptonLysoMaterial;
  G4double ComptonLysoZpos;
  G4double ComptonLysoXpos;
  G4double LYSODetTopPlateX;
  G4double LYSODetTopPlateZ;
  G4double LYSODetBottomSupportX;
  G4double LYSODetBottomSupportZ;
  G4double LYSODetSupportBallH;
  G4String LYSOPedestaMaterial;
  G4double LYSOBasePlateX;
  G4double LYSOBasePlateY;
  G4double LYSOBasePlateZ;

  G4double ComptonDetX;
  G4double ComptonDetY;
  G4double ComptonDetZ;
  G4double ComptonDetXpos;
  G4double ComptonDetYpos;
  G4double ComptonDetZpos;
  G4String ComptonDetMaterial;
  
  G4double ComptonTrackerX;
  G4double ComptonTrackerY;
  G4double ComptonTrackerZ;
  G4double ComptonTrackerDZ;
  G4int ComptonTrackerNLayers;
  G4String ComptonTrackerMaterial;
  G4int ComptonTrackerNCellX; 
  G4int ComptonTrackerNCellY; 

  G4String CollimatorMaterial;
  G4String CollimatorSupportMaterial;
  G4double CollimatorZ; 
  G4double CollimatorZpos;
  G4double CollimatorRout;
  G4double CollimatorRin;
  G4double ComptonCaloX;
  G4double ComptonCaloY;
  G4double ComptonCaloZ;
  G4String ComptonCaloMaterial;
  
  G4String GammaCaloFoilMaterial;
  G4double GammaCaloZ;
  G4double GammaCaloZpos;
  G4String GammaMonitorDetMaterial;
  G4String GammaCaloMaterial;
  G4double LeadGlassX;
  G4double LeadGlassY;
  G4double LeadGlassZ;
  G4double GammaCaloFoilThickness;
  G4double BPipeRLG;
  G4String CalSupportMaterial;
  G4double CalSupportR;
  G4double CalSupporthickness;
  G4String CAbsorberMaterial;
  G4double CAbsorberThickness;
  G4double CAbsorberZpos;
  
  G4double GammaBeamDumpZpos;
  G4double GammaBeamDumpZ;
  G4double GammaBeamDumpR;
  G4String GammaBeamDumpMaterial;
  G4double GammaBeamDumpInsertR ;
  G4double GammaBeamDumpInsertZ;
  G4String GammaBeamDumpInsertMaterial;
  G4String GammaBeamDumpSupportMaterial;
  G4String GammaMonitorSupportMaterial;
  G4String GammaMonitorHolderMaterial;
  
  G4String ComptShieldingMaterial;
  G4String ComptShieldingMaterial2;
  G4String ComptShieldingPlateMaterial;
  G4double ComptShieldingZpos;
  G4double ComptShieldingX;
  G4double ComptShieldingY;
  G4double ComptShieldingZ;
  G4double ComptShieldingDeepZ;
  G4double ComptShieldingDeepMargine;
  G4double ComptShieldingProductionCut;
  G4double ComptShieldingSuppY;
  
  G4String IPChamberMaterial;
  G4double IPRmin;
  G4double IPRmax;
  G4double IPHight;
  G4double IPRHoleCut;
  std::vector<G4double> IPChamberHolesRmin;
  std::vector<G4double> IPChamberHolesPhi; 
  G4double IPChamberPipeLength;
  G4double IPChamberPipeFlangedR;
  G4double IPChamberPipeFlangeThickness;
  G4double IPChamberBottomThickness;
  G4double IPChamberTopThickness;
  G4double IPChamberTopR;
  G4double IPChamberSupportH;
  G4double IPChamberSupportR;
  G4double IPChamberSupportRpos;
  
  G4double TypMBMagnetX;
  G4double TypMBMagnetY;
  G4double TypMBMagnetZ;
  G4double TypMBMagnetCoreZ;
  G4double TypMBMagnetEndCapZ;
  G4double TypMBMagnetEndCapW;
  G4double TypMBMagnetEndCapHoleX;
  G4double TypMBMagnetEndCapHoleOutX;
  G4double TypMBFieldLength;
  
  G4double TypMBMagnetHoleX;
  G4double TypMBMagnetHoleY;

  G4String TypMBWireMaterial;
  G4double TypMBWireGap;
  G4double TypMBWireXW; 
  G4double TypMBWireXT;
  G4double TypMBWireL;
  G4double TypMBWireH;
  G4double TypMBWireYShift;

  G4double TypMBCoreL;
  G4double TypMBCoreW;
  G4double TypMBCoreT;
  G4double TypMBWedgeLFraction;

  G4double TypMBSupportH;
  G4double TypMBSupportV;
  G4double TypMBSupportR;
  G4double TypMBSupportPosDz;
  G4double TypMBSupportPosDx;

  G4double FlashMMiddleCutY;
  G4double FlashMMiddleCutX;
  G4double FlashMFieldGapY;
  G4double FlashMagnetX;
  G4double FlashMagnetY;
  G4double FlashMagnetZ;
  G4double FlashMagnetCutD;
  G4double FlashMagnetCoilZ;
  G4double FlashMagnetCoilD;
  G4double FlashMagnetCoilH;
  G4double FlashMagnetCoilGapH;
  G4double FlashMagnetCoreX;
  G4double FlashMagnetCoreZ;
  G4double FlashMFieldX;
  G4double FlashMagneteffY;
  G4double FlashMFieldLength;

  G4double BeamPipe2ProfilerZ;

  G4String VacChamberMaterial;
  G4String VacChamberWindowMaterial;
  G4String VacChamberSideMaterial;
  G4double VacChamberXZWallThickness;
  G4double VacChamberDetXZWidth;
  G4double VacChamberMagXZWidth;
  G4double VacChambertoOPPPDetZGap;
  G4double VacChamberHight;
  G4double VacChamberWindowThickness;
  G4double VacChamberWindowX;
  G4double VacChamberWindowY;
  G4double VacChamberSideWallThickness;
  G4double IPMAgnetBeamPipeXGap;

  G4double VacChamberFrontThickness;
  G4double VacChamberFrontCutX;
  G4double VacChamberFrontCutY;
  G4double VacChamberFrontCutXpos;
  G4double VacChamberWindowPanelX;
  G4double VacChamberWindowPanelY;
  G4double VacChamberWindowPanelZ;
  G4double VacChamberWindowPanelXpos;
  G4double VacChamberWindowPanelCutY;
  G4double VacChamberWindowPanelCutX;
  G4double VacChamberWindowPanelCutXpos;
  G4double VacChamberReinforceBarX;
  G4double VacChamberReinforceBarZ;
  G4double VacChamberReinforceFilletXZ;
  G4double VacChamberFrameLX;
  G4double VacChamberFrameLY;
  G4double VacChamberFrameHX;
  G4double VacChamberFrameHY;
  G4double VacChamberFrameHZ;

  G4String TAUIChamberMaterial;
  G4double TAUIChamberBottomX;
  G4double TAUIChamberBottomZ;
  G4double TAUIChamberBottomY;
  G4double TAUIChamberX;
  G4double TAUIChamberZ;
  G4double TAUIChamberY;
  G4double TAUIChamberTopX;
  G4double TAUIChamberTopZ;
  G4double TAUIChamberTopY;
  G4double TAUIChamberSideThickness;
  G4double TAUIChamberFrontThickness;
  G4double TAUIChamberBeamBottomY;
  G4double TAUIChamberSupportR;
  G4double TAUIChamberSupportH;
  G4double TAUIChamberPipeLength;
  G4double TAUIChamberPipeFlangedH;
  G4double TAUIChamberPipeFlangedL;

  G4double TAUIChamberElPipeRIn;
  G4double TAUIChamberElPipeROut;
  G4double TAUIChamberXpos;
  G4double TAUIChamberYpos;
  G4double TAUIChamberZpos;
  G4double TAUIChamberBBoardX;
  G4double TAUIChamberBBoardY;
  G4double TAUIChamberBBoardZ;
  G4double TAUIChamberBBoardGapY;
  G4String TAUIChamberBBoardMaterial;
  G4String TAUIChamberMirrorMaterial;
  std::vector<G4double> TAUIChamberMirrorR;
  std::vector<G4double> TAUIChamberMirrorD;
  std::vector<G4double> TAUIChamberMirrorTheta;
  std::vector<G4double> TAUIChamberMirrorXpos;
  std::vector<G4double> TAUIChamberMirrorZpos;
  std::vector<G4String> TAUIChamberMirrorName;
  G4double TAUIChamberMirrorRingH;
  G4String TAUIChamberMirrorRingMaterial;
  G4String TAUIChamberMirrorHolderMaterial;

  G4double ICBBoardTableGapY;
  G4double ICBBoardTableX;
  G4double ICBBoardTableY;
  G4double ICBBoardTableZ;
  G4double ICBBoardTableLegR;

  G4double IPVolumeZ;
  G4double IPVolumeR;

  G4double OPPPDetOffsetZ;
  G4String OPPPDetSupportMaterial;
  G4double OPPPDetBottomSupportX;
  G4double OPPPDetBottomSupportY;
  G4double OPPPDetBottomSupportZ;
  G4double OPPPDetSupportBallH;
  G4double OPPPDetTopPlateX;
  G4double OPPPDetTopPlateY;
  G4double OPPPDetTopPlateZ;
  G4double OPPPAngleHolderX;
  G4double OPPPAngleHolderY;
  G4double OPPPAngleHolderZ;
  G4double OPPPAngleHolderBaseX;
  G4double OPPPAngleHolderBaseY;
  G4double OPPPAngleHolderBaseT;
  G4double OPPPAngleHolderH;
  G4double OPPPAngleHolderT;
  G4double OPPPAngleHolderLCutX;
  G4double OPPPAngleHolderLCutY;
  G4double OPPPAngleHolderRCutR;
  G4double OPPPAngleHolderTopCutXY;

  G4double OPPPCoolingPipeR;
  G4double OPPPCoolingPipeWall;
  G4double OPPPColdPlateCarbonFiberThikness;
  G4double OPPPColdPlateEpoxyThikness;
  G4String OPPPFrameMaterial;
  G4String OPPPColdPlateMaterial;
  G4String OPPPCoolingPipeMaterial;
  G4String OPPPCoolingMaterial;
  G4double OPPPCoolingPipeYpos;
  G4double OPPPColdPlateX;
  G4double OPPPColdPlateY;
  G4double OPPPTrussX;
  G4double OPPPTrussZ;
  G4double OPPPTrussTopXY;
  G4double OPPPTrussTopThickness;
  G4int    OPPPNTruss;
  G4double OPPPStaveX;
  G4double OPPPStaveSideZ;
  G4double OPPPStaveSideZpos;
  G4double OPPPStaveLiftY;
  G4double OPPPStaveInOutZ;
  G4double OPPPStaveInOutX;
  G4double OPPPTrackerActiveX;
  G4double OPPPTrackerInterLayerZ;
  G4String OPPPSesorFPCMaterial1;
  G4String OPPPSesorFPCMaterial2;
  G4String OPPPStaveSideMaterial;
  G4String OPPPStaveHolderMaterial;

  G4double OPPPSesorFPCLayer1Z;
  G4double OPPPSesorFPCLayer2Z;
  G4int OPPPNStaveSensors;
  G4double OPPPStaveSensorsGapX;
  G4double OPPPSensorX;
  G4double OPPPSensorY;
  G4double OPPPSensorZ;
  G4double OPPPSensorNCellX;
  G4double OPPPSensorNCellY;
  G4double OPPPSensorPixelX;
  G4double OPPPSensorPixelY;
  G4double OPPPSensorPixelZ;
  
  G4String OPPPHexapodMaterial;
  G4double OPPPHexapodUpR;
  G4double OPPPHexapodUpH;
  G4double OPPPHexapodDownR;
  G4double OPPPHexapodDownH;
  G4double OPPPHexapodY;
  G4double OPPPHexapodThickness;

  G4String OPPPBasePlateMaterial;
  G4double OPPPBasePlateX;
  G4double OPPPBasePlateY;
  G4double OPPPBasePlateZ;
  G4double OPPPBasePlateR;
  
  G4double OPPPServiceSupportX;
  G4double OPPPServiceSupportY;
  G4double OPPPServiceSupportZ;
  G4double OPPPServiceSupportTubeX;
  G4double OPPPServiceSupportTubeRin;
  G4double OPPPServiceSupportTubeThickness;
  G4String OPPPServiceSupportTubeMaterial;

  G4double OPPPServiceCoolPipeRout;
  G4double OPPPServiceCoolPipeThick;
  G4double OPPPServiceCoolPipeYdist;
  G4double OPPPServiceCoolPipeCnctrL;
  G4double OPPPServiceCoolPipeCnctrRin;
  G4double OPPPServiceCoolPipeCnctrRout;
  G4double OPPPServiceCoolPipeCnctrSL;
  G4double OPPPServiceCoolPipeCnctrSRin;
  G4double OPPPServiceCoolPipeCnctrSRout;
  G4String OPPPServiceCoolPipeCnctrMaterial;
  G4String OPPPServiceCoolPipeMaterial;
  G4double OPPPInnerServiceCoolPipeDZ;
  G4double OPPPInnerServiceCoolPipeBentR;

  G4String OPPPServiceCablMaterial;
  G4double OPPPInnerServiceCablSec0X;
  G4double OPPPInnerServiceCablSec0Y;
  G4double OPPPInnerServiceCablSec0Z;
  G4double OPPPInnerServiceCablSec0Thick;
  G4double OPPPInnerServiceCablSec1X;
  G4double OPPPInnerServiceCablSec1Y;
  G4double OPPPInnerServiceCablSec1Thick;
  G4double OPPPInnerServiceCablSec20X;
  G4double OPPPInnerServiceCablSec20Y;
  G4double OPPPInnerServiceCablSec20Thick;
  G4double OPPPInnerServiceCablSec21X;
  G4double OPPPInnerServiceCablSec21Y;
  G4double OPPPInnerServiceCablSec21Z;
  G4double OPPPInnerServiceCablSec21Thick;
  G4double OPPPInnerServiceCablSec22X;
  G4double OPPPInnerServiceCablSec22Y;
  G4double OPPPInnerServiceCablSec22Z;
  G4double OPPPInnerServiceCablSec22Thick;
  G4double OPPPInnerServiceCablSec31Y;
  G4double OPPPInnerServiceCablSec31Thick;
  G4double OPPPInnerServiceCablSec32Y;
  G4double OPPPInnerServiceCablSec32Thick;

  G4double OPPPOuterServiceCablSec0X;
  G4double OPPPOuterServiceCablSec0Y;
  G4double OPPPOuterServiceCablSec0Z;
  G4double OPPPOuterServiceCablSec0Thick;
  G4double OPPPOuterServiceCablSec1Y;
  G4double OPPPOuterServiceCablSec1Thick;

  G4double OPPPSideOutServiceCablSec0X;
  G4double OPPPSideOutServiceCablSec0Y;
  G4double OPPPSideOutServiceCablSec0Z;
  G4double OPPPSideOutServiceCablSec0Thick;
  G4double OPPPSideOutServiceCablSec1X;
  G4double OPPPSideOutServiceCablSec1Y;
  G4double OPPPSideOutServiceCablSec1Thick;
  G4double OPPPSideOutServiceCablSec20X;
  G4double OPPPSideOutServiceCablSec20Y;
  G4double OPPPSideOutServiceCablSec20Thick;
  G4double OPPPSideOutServiceCablSec21X;
  G4double OPPPSideOutServiceCablSec21Y;
  G4double OPPPSideOutServiceCablSec21Z;
  G4double OPPPSideOutServiceCablSec21Thick;
  G4double OPPPSideOutServiceCablSec22X;
  G4double OPPPSideOutServiceCablSec22Y;
  G4double OPPPSideOutServiceCablSec22Z;
  G4double OPPPSideOutServiceCablSec22Thick;
  G4double OPPPSideOutServiceCablSec31X;
  G4double OPPPSideOutServiceCablSec31Y;
  G4double OPPPSideOutServiceCablSec31Thick;
  G4double OPPPSideOutServiceCablSec32X;
  G4double OPPPSideOutServiceCablSec32Y;
  G4double OPPPSideOutServiceCablSec32Thick;

  G4double OPPPSideInServiceCablSec01X;
  G4double OPPPSideInServiceCablSec01Y;
  G4double OPPPSideInServiceCablSec01Thick;
  G4double OPPPSideInServiceCablSec02X;
  G4double OPPPSideInServiceCablSec02Y;
  G4double OPPPSideInServiceCablSec02Thick;

  G4double OPPPServiceCablTermBase1X;
  G4double OPPPServiceCablTermBase1Z;
  G4double OPPPServiceCablTermMid1X;
  G4double OPPPServiceCablTermMid1Y;
  G4double OPPPServiceCablTermMid1Z;
  G4double OPPPServiceCablTermBase2X;
  G4double OPPPServiceCablTermBase2Z;
  G4double OPPPServiceCablTermMid2X;
  G4double OPPPServiceCablTermMid2Y;
  G4double OPPPServiceCablTermMid2Z;
  G4String OPPPServiceCablTermMaterial;

  G4String TrackerElectronicsRackMaterial;
  G4double TrackerElectronicsRackX;
  G4double TrackerElectronicsRackY;
  G4double TrackerElectronicsRackZ;
  G4double TrackerElectronicsRackThick;
  G4double TrackerElectronicsPCBThick;
  G4double TrackerElectronicsRackXpos;

  G4double ECalX;
  G4double ECalY;
  G4double ECalZ;
  G4int    ECalNLayers;
  G4double ECalSensorPixelX;
  G4double ECalSensorPixelY;
  G4double ECalSensorNCellX;
  G4double ECalSensorNCellY;
  G4double OPPPTrackerECalZ;
  G4double ECalTopPlateZ;
  G4double ECalXpos;
  G4double ECalCasingTopY;
  G4double ECalCasingSideX;
  G4double ECalCasingBackZ;
  G4double ECalCasingGapZ;
  G4String ECalCasingMaterial;
  G4double ECalPipeShieldX;
  G4double ECalPipeShieldY;
  G4double ECalPipeShieldZ;
  G4double ECalPipeShieldShiftZ;
  G4String ECalPipeShieldMaterial;

  G4double OPPPECalPCBX;
  G4double OPPPECalPCBY;
  G4double OPPPECalPCBZ;
  G4String OPPPECalPCBMaterial;

  G4String ECalDumpShieldMaterialInner;
  G4String ECalDumpShieldMaterialOuter;
  G4double ECalDumpShieldX;
  G4double ECalDumpShieldY;
  G4double ECalDumpShieldInnerZ;
  G4double ECalDumpShieldOuterZ;
  G4double ECalDumpShieldXpos;
  G4double ECalDumpShieldYpos;
  G4double ECalDumpShieldZpos;

  G4double HICSDumpFrontZPos;
  G4double HICSDumpFrontXPos;
  G4double HICSDumpR;
  G4double HICSDumpZ;
  G4double HICSDumpAlInsertR;
  G4double HICSDumpAlInsertZ;
  G4double HICSDumpAirInsertR;
  G4double HICSDumpAirInsertZ;
  G4String HICSDumpMaterial;
  G4String HICSDumpAlInsertMaterial;
  G4double HICSDumpProductionCut;
  G4double HICSDump2TPosX;

  G4String HICSShieldingSideMaterial;
  G4String HICSShieldingMiddleMaterial;
  G4String HICSShieldingTopMaterial;
  G4String HICSShieldingSupportMaterial;
  G4double HICSShieldingXPos;
  G4double HICSShieldingX;
  G4double HICSShieldingY;
  G4double HICSShieldingSideZ;
  G4double HICSShieldingMiddleZ;
  G4double HICSShieldingGapX;
  G4double HICSShieldingGapY;
  G4double HICSShieldingSupportZ;
  G4double HICSShieldingSupportY;

  G4String HICSNeutronAbsorberMaterial;
  G4double HICSNeutronAbsorberX;
  G4double HICSNeutronAbsorberY;
  G4double HICSNeutronAbsorberZ;

  G4double ComptonElBackshift;

  G4double HICSScintilatorXPosMag1T;
  G4double HICSCherenkovXPosMag1T;
  G4double HICSScintilatorXPosMag2T;
  G4double HICSCherenkovXPosMag2T;
  G4double ComptonElectronBeamtoStageY;
  G4double HICSScintSupportBallH;

  G4String BremBasePlateMaterial;
  G4double BremBasePlateX;
  G4double BremBasePlateY;
  G4double BremBasePlateZ;
  G4double BremBasePlateR;

  G4String HICSDetSupportMaterial;
  G4double HICSBasePlateX;
  G4double HICSBasePlateY;
  G4double HICSBasePlateZ;
  G4double HICSBasePlateR;

  G4double HICSDetBottomSupportX;
  G4double HICSDetBottomSupportY;
  G4double HICSDetBottomSupportZ;
   
  G4String HICSBasePlateMaterial;

  G4String BeamProfilerContanerMaterial;
  G4double BeamProfilerContanerX;
  G4double BeamProfilerContanerY;
  G4double BeamProfilerContanerZ;
  G4double BeamProfilerContanerD;
  G4String BeamProfilerContanerWallMaterial;
  G4double BeamProfilerWindowX;
  G4double BeamProfilerWindowY;
  G4double BeamProfilerWindowZ;
  G4String BeamProfilerWindowMaterial;

  G4double BeamProfilerSensorX;
  G4double BeamProfilerSensorY;
  G4double BeamProfilerSensorZ;
  G4double BeamProfilerMetalizationZ;
  std::vector<G4double> BeamProfilerZpos;
  G4int BeamProfilerNCellX;
  G4int BeamProfilerNCellY;
  G4String BeamProfilerSensorMaterial;
  G4String BeamProfilerSensorMetalization;

  G4String BeamProfilerV408Material;
  G4double BeamProfilerV408X;
  G4double BeamProfilerV408Y;
  G4double BeamProfilerV408Z;
  G4double BeamProfilerMotorGap;
  G4String BeamProfilerQ545Material;
  G4double BeamProfilerQ545X;
  G4double BeamProfilerQ545Y;
  G4double BeamProfilerQ545Z;
  G4double BeamProfilerPCBX;
  G4double BeamProfilerPCBY;
  G4double BeamProfilerPCBZ;
  G4double BeamProfilerPCBCutX;
  G4double BeamProfilerPCBCutY;
  G4double BeamProfilerPCBHolderX;
  G4double BeamProfilerPCBHolderY;
  G4double BeamProfilerPCBHolderZ;
  G4double BeamProfilerPCBHolderCutR;
  G4String BeamProfilerPCBHolderMaterial;

  G4double BeamProfilerAngleSupportX;
  G4double BeamProfilerAngleSupportD;
  G4String BeamProfilerAngleSupportMaterial;

  G4String BeamPipeLidMaterial;
  G4double BeamPipeLidThickness;
  G4double GamVacChamberGap;
  G4double GamVacChamberSideWallThickness;
  G4double GamVacChamberTopBotThickness;

  G4String BSMCaloMaterial;
  G4double BSMCaloX;
  G4double BSMCaloY;
  G4double BSMCaloLayerZ;
  G4int    BSMCaloNLayers;
  G4double BSMCaloZPos;
  G4int BSMCaloNCellX;
  G4int BSMCaloNCellY;

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
