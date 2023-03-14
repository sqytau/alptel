
/// \brief Implementation of the HistoMessenger class
//

#include "HistoMessenger.hh"

#include "HistoManager.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


HistoMessenger::HistoMessenger(HistoManager* histoM)
:G4UImessenger(), fHistManager(histoM), fTreeCutXCmd(0), fTreeCutYCmd(0)
{ 
  fTreeCutXCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/ttree_cut_x",this);
  fTreeCutXCmd->SetGuidance("Set cut for track abs(x) for output TTree");
  fTreeCutXCmd->SetParameterName("mTackCutX",false);
  fTreeCutXCmd->SetRange("mTackCutX>0.");
  fTreeCutXCmd->SetUnitCategory("Length");
  fTreeCutXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);   

  fTreeCutYCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/ttree_cut_y",this);
  fTreeCutYCmd->SetGuidance("Set cut for track abs(y) for output TTree");
  fTreeCutYCmd->SetParameterName("mTackCutY",false);
  fTreeCutYCmd->SetRange("mTackCutY>0.");
  fTreeCutYCmd->SetUnitCategory("Length");
  fTreeCutYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);   
   
  fTreeParticleTypeCmd = new G4UIcmdWithAString("/lxphoton/ttree_particle",this);
  fTreeParticleTypeCmd->SetGuidance("Particle type to save in output ttree");
  fTreeParticleTypeCmd->SetParameterName("ttreeParticle",false);
  fTreeParticleTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);   
  
}


HistoMessenger::~HistoMessenger()
{
  delete fTreeCutXCmd;
  delete fTreeCutYCmd;
  delete fTreeParticleTypeCmd;
}


void HistoMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fTreeCutXCmd)
    { fHistManager->SetTreeCutX(fTreeCutXCmd->GetNewDoubleValue(newValue));}
  if (command == fTreeCutYCmd)
    { fHistManager->SetTreeCutY(fTreeCutYCmd->GetNewDoubleValue(newValue));}
  if (command == fTreeParticleTypeCmd)
    { fHistManager->SetTreeParticle(newValue);}
}

