//
//
/// \brief Definition of the PhysicsList class
//

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "OpticalPhysics.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;
class G4VDiscreteProcess;



class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
 ~PhysicsList();

  virtual void ConstructParticle();
        
  void AddPhysicsList(const G4String& name);
    
  virtual void ConstructProcess();    
  void AddDecay();
  void AddStepMax();

protected:
  void AlterPhysicsList(const G4String& name, G4VPhysicsConstructor* phys);
  
private:

  PhysicsListMessenger* fMessenger; 

  G4String fEmName;
  G4VPhysicsConstructor*  fEmPhysicsList;

  G4VDiscreteProcess*  fMaxStepLimit;
//   G4VPhysicsConstructor*  fMaxStepLimit;    // alternative solution
};



#endif

