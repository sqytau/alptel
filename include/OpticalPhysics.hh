#ifndef OpticalPhysics_h
#define OpticalPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4LossTableManager.hh"
//~ #include "G4EmConfigurator.hh"
//~ #include "G4BraggIonGasModel.hh"
//~ #include "G4BetheBlochIonGasModel.hh"
//~ #include "G4IonFluctuations.hh"
//~ #include "G4IonParametrisedLossModel.hh"
//~ #include "G4UniversalFluctuation.hh"
#include "globals.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpWLS.hh"

class OpticalPhysics : public G4VPhysicsConstructor {
	
	public:
		
		OpticalPhysics(const G4String& name = "OpticalPhysics");
		virtual ~OpticalPhysics();
		void ConstructParticle();
		void ConstructProcess();
};

#endif