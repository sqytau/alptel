//
/// \file electromagnetic/TestEm5/include/TrackingAction.hh
/// \brief Definition of the TrackingAction class
//

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;

class TrackingAction : public G4UserTrackingAction {

  public:  
    TrackingAction(DetectorConstruction*,EventAction*);
   ~TrackingAction() {};
   
    virtual void  PreUserTrackingAction(const G4Track*);   
    virtual void PostUserTrackingAction(const G4Track*);
    
    void SetVerbose(const G4int v = 0) { fVerbose = v; };
    
  private:
    DetectorConstruction* fDetector;
    EventAction*          fEventAction;
    
    G4double fZstartAbs, fZendAbs, fZendMagnet;
    G4double fPrimaryCharge;
    
    G4int  fVerbose;
};

#endif
