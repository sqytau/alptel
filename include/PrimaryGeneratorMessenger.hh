//
/// \brief Definition of the PrimaryGeneratorMessenger class
//

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;


class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
   ~PrimaryGeneratorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    PrimaryGeneratorAction*    fAction;
    
    G4UIdirectory*             fGunDir;      
    G4UIcmdWithoutParameter*   fDefaultCmd;
    G4UIcmdWithAString*        fBeamTypeCmd; 
    G4UIcmdWithAString*        fSpectraFileCmd;
    
    G4UIcmdWithADoubleAndUnit* fBeamSigmaXCmd;
    G4UIcmdWithADoubleAndUnit* fBeamSigmaYCmd;
    G4UIcmdWithADoubleAndUnit* fBeamPosZCmd;

    G4UIcmdWithAString*        fMCFileCmd;
    G4UIcmdWithAString*        fMCFileListCmd;
    G4UIcmdWithAnInteger*      fMCWightScaleCmd;
    G4UIcmdWithAnInteger*      fMCSelectParticleCmd;
    G4UIcmdWith3VectorAndUnit* fBeamPositionCmd;

    G4UIcmdWithAnInteger*      fMCSkipEventsCmd;
};



#endif

