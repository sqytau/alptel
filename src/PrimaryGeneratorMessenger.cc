//
/// \brief Implementation of the PrimaryGeneratorMessenger class
//

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"



PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                             PrimaryGeneratorAction* Gun)
:G4UImessenger(),fAction(Gun),
 fGunDir(0),
 fDefaultCmd(0),
 fBeamTypeCmd(0),
 fBeamSigmaXCmd(0),
 fBeamSigmaYCmd(0),
 fBeamPosZCmd(0),
 fMCFileCmd(0),
 fMCFileListCmd(0),
 fMCWightScaleCmd(0),
 fMCSelectParticleCmd(0)

{
  fGunDir = new G4UIdirectory("/lxphoton/gun/");
  fGunDir->SetGuidance("gun control");

  fDefaultCmd = new G4UIcmdWithoutParameter("/lxphoton/gun/setDefault",this);
  fDefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
  fDefaultCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamTypeCmd = new G4UIcmdWithAString("/lxphoton/gun/beamType",this);
  fBeamTypeCmd->SetGuidance("Distribution of primary particles of the beam. Supported types: gaussian, mono, mc");
  fBeamTypeCmd->SetParameterName("beamType",false);
  fBeamTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSpectraFileCmd = new G4UIcmdWithAString("/lxphoton/gun/SpectraFile",this);
  fSpectraFileCmd->SetGuidance("File contaning tabulated PDF of primary particles of the beam.");
  fSpectraFileCmd->SetParameterName("SpectraFile",false);
  fSpectraFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamSigmaXCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/gun/setSigmaX",this);
  fBeamSigmaXCmd->SetGuidance("Set sigma for beam distribution in X direction at IP");
  fBeamSigmaXCmd->SetParameterName("SigmaX",false);
  fBeamSigmaXCmd->SetRange("SigmaX>0.");
  fBeamSigmaXCmd->SetUnitCategory("Length");
  fBeamSigmaXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamSigmaYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/gun/setSigmaY",this);
  fBeamSigmaYCmd->SetGuidance("Set sigma for beam distribution in Y direction at IP");
  fBeamSigmaYCmd->SetParameterName("SigmaY",false);
  fBeamSigmaYCmd->SetRange("SigmaY>0.");
  fBeamSigmaYCmd->SetUnitCategory("Length");
  fBeamSigmaYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamPosZCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/gun/setPosZ",this);
  fBeamPosZCmd->SetGuidance("Set z position of the beam");
  fBeamPosZCmd->SetParameterName("PosZ",false);
  fBeamPosZCmd->SetUnitCategory("Length");
  fBeamPosZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMCFileCmd = new G4UIcmdWithAString("/lxphoton/gun/MCParticlesFile",this);
  fMCFileCmd->SetGuidance("File contaning tabulated list of primary particles to simulate.");
  fMCFileCmd->SetParameterName("MCParticlesFile",false);
  fMCFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMCFileListCmd = new G4UIcmdWithAString("/lxphoton/gun/MCParticlesFileList",this);
  fMCFileListCmd->SetGuidance("File contaning list of files with primary particles to simulate.");
  fMCFileListCmd->SetParameterName("MCParticlesFileList",false);
  fMCFileListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMCWightScaleCmd = new G4UIcmdWithAnInteger("/lxphoton/gun/MCWeightScale", this);
  G4String sstr("If weight of MC event is greater than MCWeightScale, ");
  sstr += G4String("then MCWeightScale particles are generated and the weight is scaled as Weight/MCWeightScale");
  fMCWightScaleCmd->SetGuidance(sstr);
  fMCWightScaleCmd->SetParameterName("MCWeightScale",false);
  fMCWightScaleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMCSelectParticleCmd = new G4UIcmdWithAnInteger("/lxphoton/gun/MCSelectParticle", this);
  G4String sstrps("If something is specified, then only these slected particles from MC file will be processed, ");
  sstrps += G4String("otherwise all known will be processed without filters.");
  fMCSelectParticleCmd->SetGuidance(sstrps);
  fMCSelectParticleCmd->SetParameterName("MCSelectParticle",false);
  fMCSelectParticleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamPositionCmd = new G4UIcmdWith3VectorAndUnit("/lxphoton/gun/setPosition", this);
  fBeamPositionCmd->SetGuidance("Set primary particle position, relevant for mono beam type only.");
  fBeamPositionCmd->SetParameterName("ParticlePositionX", "ParticlePositionY", "ParticlePositionZ", false);
  fBeamPositionCmd->SetUnitCategory("Length");
  fBeamPositionCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMCSkipEventsCmd = new G4UIcmdWithAnInteger("/lxphoton/gun/skipMCEvents", this);
  fMCSkipEventsCmd->SetGuidance("Set the number of events which heed to be skipped in the beginning of MC file.");
  fMCSkipEventsCmd->SetParameterName("skipMCEvents", false);
  fMCSkipEventsCmd->SetDefaultValue(0);
  fMCSkipEventsCmd->SetRange("skipMCEvents >= 0");
  fMCSkipEventsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}



PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fDefaultCmd;
  delete fBeamTypeCmd;
  delete fSpectraFileCmd;
  delete fGunDir;
  delete fBeamSigmaXCmd;
  delete fBeamSigmaYCmd;
  delete fBeamPosZCmd;
  delete fMCFileCmd;
  delete fMCFileListCmd;
  delete fMCWightScaleCmd;
  delete fMCSelectParticleCmd;
  delete fBeamPositionCmd;
  delete fMCSkipEventsCmd;
}



void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  if (command == fDefaultCmd)
    {fAction->SetDefaultKinematic();}

  if (command == fBeamTypeCmd)
   { fAction->SetBeamType(newValue);}

  if (command == fSpectraFileCmd)
   { fAction->SetSpectraFile(newValue);}

  if (command == fBeamSigmaXCmd)
   { fAction->SetSigmaX(fBeamSigmaXCmd->GetNewDoubleValue(newValue));}

  if (command == fBeamSigmaYCmd)
   { fAction->SetSigmaY(fBeamSigmaYCmd->GetNewDoubleValue(newValue));}

  if (command == fBeamPosZCmd)
   { fAction->SetBeamPosZ(fBeamPosZCmd->GetNewDoubleValue(newValue));}

  if (command == fMCFileCmd)
   { fAction->SetMCParticleFile(newValue);}

   if (command == fMCFileListCmd)
   { fAction->SetMCParticleFile(newValue, true);}

   if (command == fMCWightScaleCmd)
   { fAction->SetMCWeightScale(fMCWightScaleCmd->GetNewIntValue(newValue));}

   if (command == fMCSelectParticleCmd)
   { fAction->AddMCSelectParticle(fMCSelectParticleCmd->GetNewIntValue(newValue));}

   if (command == fBeamPositionCmd)
   { fAction->SetBeamPosition(fBeamPositionCmd->GetNew3VectorValue(newValue));}

   if (command == fMCSkipEventsCmd)
   { fAction->SetSkipEvents(fMCSkipEventsCmd->GetNewIntValue(newValue));}

}



