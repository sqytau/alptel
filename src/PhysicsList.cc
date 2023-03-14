//
/// \brief Implementation of the PhysicsList class
//

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "PhysListEmStandardSSM.hh"
#include "PhysListEmStandardGS.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsSS.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4Decay.hh"
#include "StepMax.hh"

#include "G4LossTableManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4VDiscreteProcess.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4Electron.hh"
#include "OpticalPhysics.hh"
#include "LXSetUp.hh"




PhysicsList::PhysicsList() : G4VModularPhysicsList(),
 fMessenger(0),fEmPhysicsList(0)
{  
  fMessenger = new PhysicsListMessenger(this); 
  SetVerboseLevel(1);
     
  // EM physics
  fEmName = G4String("local");
  fEmPhysicsList = new PhysListEmStandard(fEmName);
  RegisterPhysics(fEmPhysicsList);
  
  G4LossTableManager::Instance();
  SetDefaultCutValue(1*um);  
  
  fMaxStepLimit = new StepMax();
// The following alternative must be used with fMaxStepLimit->ConstructProcess(); 
// in PhysicsList::ConstructProcess() and corresponding setting in DetectorConstruction
//   fMaxStepLimit = new G4StepLimiterPhysics();   
}



PhysicsList::~PhysicsList()
{
//   delete fEmPhysicsList;
  delete fMessenger;  
}



void PhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();

    LXSetUp *lxs = LXSetUp::Instance();

  if(lxs->ScintCerenkovPhysics){
     G4VPhysicsConstructor* opticalPhysics = new OpticalPhysics();
     opticalPhysics->ConstructParticle();}

}



void PhysicsList::ConstructProcess()
{
//   AddTransportation();
  G4VModularPhysicsList::ConstructProcess();
  AddDecay();  
  AddStepMax();
//   fMaxStepLimit->ConstructProcess();
      LXSetUp *lxs = LXSetUp::Instance();

  if(lxs->ScintCerenkovPhysics){
    G4VPhysicsConstructor* opticalPhysics = new OpticalPhysics();
    opticalPhysics->ConstructProcess();}

    
}



void PhysicsList::AddDecay()
{
  // Add Decay Process

  G4Decay* fDecayProcess = new G4Decay();

  G4ParticleTable::G4PTblDicIterator *theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) { 

      pmanager->AddProcess(fDecayProcess);

      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(fDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(fDecayProcess, idxAtRest);

    }
  }
}



void PhysicsList::AddStepMax()
{
  StepMax* stepMaxProcess = dynamic_cast<StepMax*>(fMaxStepLimit);

  G4ParticleTable::G4PTblDicIterator *theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
        {
          pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
  }
}



void PhysicsList::AlterPhysicsList(const G4String& name, G4VPhysicsConstructor* phys)
{
    fEmName = name;
      RemovePhysics(fEmPhysicsList);
    fEmPhysicsList = phys;
//     fEmPhysicsList->SetVerboseLevel(3);

//     G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
//     ph->SetVerboseLevel(4);
    RegisterPhysics(fEmPhysicsList);
//     ph->DumpOrdingParameterTable();
}


void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>-1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "local") {
    AlterPhysicsList(name, new PhysListEmStandard(name));

  } else if (name == "emstandard_opt0") {
    AlterPhysicsList(name, new G4EmStandardPhysics());

  } else if (name == "emstandard_opt1") {
    AlterPhysicsList(name, new G4EmStandardPhysics_option1());

  } else if (name == "emstandard_opt2") {
    AlterPhysicsList(name, new G4EmStandardPhysics_option2());
    
  } else if (name == "emstandard_opt3") {
    AlterPhysicsList(name, new G4EmStandardPhysics_option3());
    
  } else if (name == "emstandard_opt4") {
    AlterPhysicsList(name, new G4EmStandardPhysics_option4());
        
  } else if (name == "emstandardSS") {
    AlterPhysicsList(name, new G4EmStandardPhysicsSS());

  } else if (name == "standardSSM") {
    AlterPhysicsList(name, new PhysListEmStandardSSM());

  } else if (name == "emstandardWVI") {
    AlterPhysicsList(name, new G4EmStandardPhysicsWVI());

  } else if (name == "standardGS") {
    AlterPhysicsList(name, new PhysListEmStandardGS());

  } else if (name == "empenelope"){
    AlterPhysicsList(name, new G4EmPenelopePhysics());

  } else if (name == "emlowenergy"){
    AlterPhysicsList(name, new G4EmLowEPPhysics());

  } else if (name == "emlivermore"){
    AlterPhysicsList(name, new G4EmLivermorePhysics());
 
  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}




