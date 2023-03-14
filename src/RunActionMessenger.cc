//
/// \brief Implementation of the RunMessenger class
//

#include <sstream>

#include "RunActionMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"


RunActionMessenger::RunActionMessenger(RunAction* ra)
:G4UImessenger(), fRunAction(ra),
 fRunSettingsDir(0),      
 fAddInterceptVolCmd(0)
{
  fRunSettingsDir = new G4UIdirectory("/luxe/run/");
  fRunSettingsDir->SetGuidance("Run settings");
 
  fAddInterceptVolCmd = new G4UIcmdWithAString("/luxe/run/add_intercept_volume",this);
  fAddInterceptVolCmd->SetGuidance("Tracks crossing the volume boundary are saved to the output tree");
  fAddInterceptVolCmd->SetParameterName("InterceptVolume",false);
  fAddInterceptVolCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fAddSensitiveVolCmd = new G4UIcmdWithAString("/luxe/run/add_sensitive_volume",this);
  fAddSensitiveVolCmd->SetGuidance("Sensitive volume, energy deposited in this volume will be saved in output tree.");
  fAddSensitiveVolCmd->SetParameterName("Sensitivevolume",false);
  fAddSensitiveVolCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDumpGeometryCmd = new G4UIcmdWithABool("/luxe/run/dump_geometry",this);
  fDumpGeometryCmd->SetGuidance("Save geometry in gdml file.");
  fDumpGeometryCmd->SetParameterName("DumpGeometry",false);
  fDumpGeometryCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSaveTrajectoryCmd = new G4UIcmdWithABool("/luxe/run/save_primary_trajectory",this);
  fSaveTrajectoryCmd->SetGuidance("Save trajectory of the primary particle.");
  fSaveTrajectoryCmd->SetParameterName("PrimaryTrajectory",false);
  fSaveTrajectoryCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fAddInterceptVolECutCmd = new G4UIcommand("/luxe/run/add_intercept_energy_cut", this);
  fAddInterceptVolECutCmd->SetGuidance("Energy cut for tracks crossing intercept volume.");
  fAddInterceptVolECutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4UIparameter*  sParamVol = new G4UIparameter('s');
  G4UIparameter*  dblParamECut = new G4UIparameter('d');
  G4UIparameter*  untParam = new G4UIparameter('s');
  untParam->SetParameterName("Unit");
  fAddInterceptVolECutCmd->SetParameter(sParamVol);
  fAddInterceptVolECutCmd->SetParameter(dblParamECut);
  fAddInterceptVolECutCmd->SetParameter(untParam);

  fExternaRunIDCmd = new G4UIcmdWithAnInteger("/luxe/run/run_id", this);
  fExternaRunIDCmd->SetGuidance("Set runID which will idetify the data in output TTrees.");
  fExternaRunIDCmd->SetParameterName("ExternalRunID", true);
  fExternaRunIDCmd->SetDefaultValue(-1);

}



RunActionMessenger::~RunActionMessenger()
{
  delete fRunSettingsDir;
  delete fAddInterceptVolCmd;
  delete fAddSensitiveVolCmd;
  delete fDumpGeometryCmd;
  delete fSaveTrajectoryCmd;
  delete fAddInterceptVolECutCmd;
  delete fExternaRunIDCmd;
}



void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fAddInterceptVolCmd)
    {fRunAction->AddInterceptVolume(newValue);}
    
  if (command == fAddSensitiveVolCmd)
    {fRunAction->AddSensitiveVolume(newValue);}
    
  if (command == fDumpGeometryCmd)
    {fRunAction->SetDumpGeometry(fDumpGeometryCmd->GetNewBoolValue(newValue));}
  
  if (command == fSaveTrajectoryCmd)
    {fRunAction->SetSaveTrajectory(fSaveTrajectoryCmd->GetNewBoolValue(newValue));}

  if (command == fAddInterceptVolECutCmd) {
     std::istringstream istr(newValue);
     G4String vname, units;
     G4double ecut;
     istr >> vname >> ecut >> units;

     G4double value_given = G4UIcommand::ValueOf(units);
     ecut *= value_given;
     fRunAction->AddInterceptVolumeECut(vname, ecut);
  }

  if (command == fExternaRunIDCmd)
    {fRunAction->SetExternalRunID(fExternaRunIDCmd->GetNewIntValue(newValue));}

}



