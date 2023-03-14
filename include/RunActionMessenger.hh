//
/// \brief Definition of the RunActionMessenger class
//

#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;


class RunActionMessenger: public G4UImessenger
{
  public:
    RunActionMessenger(RunAction* ra);
   ~RunActionMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    RunAction*    fRunAction;
    
    G4UIdirectory*             fRunSettingsDir;      
    G4UIcmdWithAString*        fAddInterceptVolCmd;
    G4UIcmdWithAString*        fAddSensitiveVolCmd;
    G4UIcmdWithABool*          fDumpGeometryCmd;
    G4UIcmdWithABool*          fSaveTrajectoryCmd;
    G4UIcommand*               fAddInterceptVolECutCmd;
    G4UIcmdWithAnInteger*      fExternaRunIDCmd;
};



#endif

