
/// \brief Implementation of the StepMaxMessenger class
//

#include "StepMaxMessenger.hh"

#include "StepMax.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


StepMaxMessenger::StepMaxMessenger(StepMax* stepM)
:G4UImessenger(),fStepMax(stepM),fStepMaxCmd(0)
{ 
  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/lxphoton/stepMax",this);
  fStepMaxCmd->SetGuidance("Set max allowed step length");
  fStepMaxCmd->SetParameterName("mxStep",false);
  fStepMaxCmd->SetRange("mxStep>0.");
  fStepMaxCmd->SetUnitCategory("Length");
}


StepMaxMessenger::~StepMaxMessenger()
{
  delete fStepMaxCmd;
}


void StepMaxMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fStepMaxCmd)
    { fStepMax->SetMaxStep(fStepMaxCmd->GetNewDoubleValue(newValue));}
}

