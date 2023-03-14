//
/// \brief Implementation of the DetectorMessenger class
//

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),
 fAbsMaterCmd(0),
 fAbsThickCmd(0),
 fAbsSizXYCmd(0),
 fAbsZposCmd(0),
 fWorldMaterCmd(0),
 fWorldZCmd(0),
 fWorldXYCmd(0),
 fAbsTypeCmd(0)
{
  fTestemDir = new G4UIdirectory("/lxphoton/");
  fTestemDir->SetGuidance("UI commands specific to this example.");

  fDetDir = new G4UIdirectory("/lxphoton/det/");
  fDetDir->SetGuidance("detector construction commands");

  fAbsMaterCmd = new G4UIcmdWithAString("/lxphoton/det/setAbsMat",this);
  fAbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  fAbsMaterCmd->SetParameterName("choice",false);
  fAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsMaterCmd->SetToBeBroadcasted(false);

  fAbsTypeCmd = new G4UIcmdWithAString("/lxphoton/det/setAbsType",this);
  fAbsTypeCmd->SetGuidance("Set tareget shape: (foil/wire).");
  fAbsTypeCmd->SetParameterName("TargetType",false);
  fAbsTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsTypeCmd->SetToBeBroadcasted(false);

  fWorldMaterCmd = new G4UIcmdWithAString("/lxphoton/det/setWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select Material of the World.");
  fWorldMaterCmd->SetParameterName("wchoice",false);
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldMaterCmd->SetToBeBroadcasted(false);

  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsThick",this);
  fAbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  fAbsThickCmd->SetParameterName("SizeZ",false);
  fAbsThickCmd->SetRange("SizeZ>0.");
  fAbsThickCmd->SetUnitCategory("Length");
  fAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsThickCmd->SetToBeBroadcasted(false);

  fAbsSizXYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsXY",this);
  fAbsSizXYCmd->SetGuidance("Set sizeXY of the Absorber");
  fAbsSizXYCmd->SetParameterName("SizeXY",false);
  fAbsSizXYCmd->SetRange("SizeXY>0.");
  fAbsSizXYCmd->SetUnitCategory("Length");
  fAbsSizXYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsSizXYCmd->SetToBeBroadcasted(false);

  fAbsZposCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setAbsZpos",this);
  fAbsZposCmd->SetGuidance("Set X pos. of the Absorber");
  fAbsZposCmd->SetParameterName("Xpos",false);
  fAbsZposCmd->SetUnitCategory("Length");
  fAbsZposCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fAbsZposCmd->SetToBeBroadcasted(false);

  fWorldZCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setWorldZ",this);
  fWorldZCmd->SetGuidance("Set X size of the World");
  fWorldZCmd->SetParameterName("WSizeX",false);
  fWorldZCmd->SetRange("WSizeX>0.");
  fWorldZCmd->SetUnitCategory("Length");
  fWorldZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldZCmd->SetToBeBroadcasted(false);

  fWorldXYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/setWorldXY",this);
  fWorldXYCmd->SetGuidance("Set sizeXY of the World");
  fWorldXYCmd->SetParameterName("WSizeXY",false);
  fWorldXYCmd->SetRange("WSizeXY>0.");
  fWorldXYCmd->SetUnitCategory("Length");
  fWorldXYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldXYCmd->SetToBeBroadcasted(false);

  fDetBFieldDir = new G4UIdirectory("/lxphoton/det/magnet_field/");
  fDetBFieldDir->SetGuidance("commands for setting magnetic fields");

  fSetBFieldValue = new G4UIcommand("/lxphoton/det/magnet_field/value", this);
  fSetBFieldValue->SetGuidance("Set magnetic field amplitude value");
  fSetBFieldValue->AvailableForStates(G4State_PreInit);
  fSetBFieldValue->SetToBeBroadcasted(false);
  fSetBFieldValue->SetParameter(new G4UIparameter('s'));

  fSetBFieldDistrib = new G4UIcommand("/lxphoton/det/magnet_field/distribution", this);
  fSetBFieldDistrib->SetGuidance("Set magnetic field distribution: MagnetID B_component distrib_coordinate model params...");
  fSetBFieldDistrib->AvailableForStates(G4State_PreInit);
  fSetBFieldDistrib->SetToBeBroadcasted(false);
  fSetBFieldDistrib->SetParameter(new G4UIparameter('s'));

  fGammaBeamDumpZCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/det/GammaBeamDumpZ",this);
  fGammaBeamDumpZCmd->SetGuidance("Set Z size of the GammaBeamDump");
  fGammaBeamDumpZCmd->SetParameterName("GammaBeamDumpZ",false);
  fGammaBeamDumpZCmd->SetRange("GammaBeamDumpZ>0.");
  fGammaBeamDumpZCmd->SetUnitCategory("Length");
  fGammaBeamDumpZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fGammaBeamDumpZCmd->SetToBeBroadcasted(false);

  fGammaBeamDumpMaterialCmd = new G4UIcmdWithAString("/lxphoton/det/GammaBeamDumpMaterial",this);
  fGammaBeamDumpMaterialCmd->SetGuidance("Set Material of the GammaBeamDump.");
  fGammaBeamDumpMaterialCmd->SetParameterName("GammaBeamDumpMaterial",false);
  fGammaBeamDumpMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fGammaBeamDumpMaterialCmd->SetToBeBroadcasted(false);
}



DetectorMessenger::~DetectorMessenger()
{
  delete fAbsMaterCmd; 
  delete fAbsThickCmd; 
  delete fAbsSizXYCmd;  
  delete fAbsZposCmd; 
  delete fWorldMaterCmd;
  delete fWorldZCmd;
  delete fWorldXYCmd;
  delete fDetDir;  
  delete fTestemDir;
  delete fAbsTypeCmd;
  delete fSetBFieldValue;
  delete fSetBFieldDistrib;
  delete fGammaBeamDumpZCmd;
  delete fGammaBeamDumpMaterialCmd;
  
}



void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if ( command == fAbsMaterCmd )
   {fDetector->SetAbsorberMaterial(newValue);}

  if ( command == fAbsTypeCmd )
   {fDetector->SetAbsorberType(newValue);}

  if ( command == fWorldMaterCmd )
   {fDetector->SetWorldMaterial(newValue);}

  if ( command == fAbsThickCmd )
  {fDetector->SetAbsorberThickness(fAbsThickCmd->GetNewDoubleValue(newValue));}

  if ( command == fAbsSizXYCmd )
   {fDetector->SetAbsorberSizeXY(fAbsSizXYCmd->GetNewDoubleValue(newValue));}

  if ( command == fAbsZposCmd )
   {fDetector->SetAbsorberZpos(fAbsZposCmd->GetNewDoubleValue(newValue));}

  if ( command == fWorldZCmd )
   {fDetector->SetWorldSizeZ(fWorldZCmd->GetNewDoubleValue(newValue));}

  if ( command == fWorldXYCmd )
   {fDetector->SetWorldSizeXY(fWorldXYCmd->GetNewDoubleValue(newValue));}

  if ( command == fSetBFieldValue ) {
    std::istringstream istr(newValue);
    G4String magid, fcomp;
    istr >> magid >> fcomp;
    G4String params(newValue.data() + istr.tellg());
    fDetector->AddBFieldModel(magid, std::make_tuple(fcomp, "bvalue", "bvalue", params));
  }

  if ( command == fSetBFieldDistrib ) {
    std::istringstream istr(newValue);
    G4String magid, fcomp, coord, fmodel;
    istr >> magid >> fcomp >> coord >> fmodel;
    G4String params(newValue.data() + istr.tellg());
//     params.strip(G4String::leading);
//     G4cout << "Parsed: !!" << magid << "!!  !!" << fcomp << "!!  !!" 
//            << coord << "!!  !!" << fmodel << "!!  !!" << params << "!!" << G4endl;
    fDetector->AddBFieldModel(magid, std::make_tuple(fcomp, coord, fmodel, params));
  }

  if ( command == fGammaBeamDumpZCmd )
   {fDetector->SetGammaBeamDumpLength(fGammaBeamDumpZCmd->GetNewDoubleValue(newValue));}

  if ( command == fGammaBeamDumpMaterialCmd )
   {fDetector->SetGammaBeamDumpMaterial(newValue);}

}

