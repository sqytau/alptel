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
 fWorldXYCmd(0)
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
}



void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if ( command == fAbsMaterCmd )
   {fDetector->SetAbsorberMaterial(newValue);}

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

}

