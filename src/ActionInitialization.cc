//
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
#include "SteppingVerbose.hh"



ActionInitialization::ActionInitialization(DetectorConstruction* det)
 : G4VUserActionInitialization(),fDetector(det)
{ }



ActionInitialization::~ActionInitialization()
{ }



void ActionInitialization::BuildForMaster() const
{
 SetUserAction(new RunAction(fDetector));
}




void ActionInitialization::Build() const
{
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(fDetector);
  
  SetUserAction(primary);
 
  RunAction* runaction = new RunAction(fDetector,primary);
  SetUserAction(runaction); 
  
  EventAction* eventaction = new EventAction();
  SetUserAction(eventaction);

  SetUserAction(new TrackingAction(fDetector,eventaction));

  SetUserAction(new SteppingAction(fDetector,eventaction));

  SetUserAction(new StackingAction(eventaction));
}  



G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
{
  return new SteppingVerbose();
}  


