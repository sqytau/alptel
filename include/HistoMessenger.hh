
/// \brief Definition of the HistoMessenger class
//
// $Id: HistoMessenger.hh
//

#ifndef HistoMessenger_h
#define HistoMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HistoManager;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;


class HistoMessenger: public G4UImessenger
{
  public:
    HistoMessenger(HistoManager* histoM);
   ~HistoMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    HistoManager* fHistManager;
    G4UIcmdWithADoubleAndUnit  *fTreeCutXCmd;
    G4UIcmdWithADoubleAndUnit  *fTreeCutYCmd;
    G4UIcmdWithAString         *fTreeParticleTypeCmd;
    
};

#endif
