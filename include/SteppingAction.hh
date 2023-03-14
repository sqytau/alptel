//
/// \brief Definition of the SteppingAction class
//
// $Id: SteppingAction.hh 73026 2013-08-15 09:12:46Z gcosmo $
//

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class RunAction;
class EventAction;
class G4StepPoint;

class SteppingAction : public G4UserSteppingAction
{
  public:
   SteppingAction(DetectorConstruction*,EventAction*);
  ~SteppingAction();

   virtual void UserSteppingAction(const G4Step*);

   void ProcessInAbsorber(const G4Step* aStep);

  protected:
   G4int CheckPointHistory(const G4StepPoint *sp, const G4String vtname) const;

   void ProcessInDetector(const G4Step* aStep, const std::vector<G4int> &dettype, const G4int vdepth);
   void ProcessInTracker(const G4Step* aStep, const G4String &dname, const G4int dettype,
                                      const G4int detdepth, const G4int layerdepthd);
   void ProcessPrimaryTrack(const G4Step* aStep);
   void ProcessScintCerenkov(const G4Step* aStep);

  private:
    DetectorConstruction* fDetector;
    EventAction*          fEventAction;
};


#endif
