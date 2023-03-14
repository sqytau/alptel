//
/// \brief Implementation of the TrackingAction class
//

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"


TrackingAction::TrackingAction(DetectorConstruction* DET, EventAction* EA)
:G4UserTrackingAction(),fDetector(DET), fEventAction(EA), fZstartAbs(0.0), fZendAbs(0.0), fZendMagnet(0.0),
 fPrimaryCharge(0.0), fVerbose(0) 
{ }
 

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  // few initialisations
  //
    
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
    
  if (aTrack->GetParentID() == 0) {
    fPrimaryCharge = aTrack->GetDefinition()->GetPDGCharge();
    
    G4ThreeVector position = aTrack->GetPosition();
    G4ThreeVector vertex   = aTrack->GetVertexPosition();
    G4ThreeVector momentum   = aTrack->GetMomentum();
    
    G4double evweght = 1.0;
    const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    if (evinf) evweght = evinf->GetWeight();
    
    analysisManager->FillH1(53, position.x(), evweght);
    analysisManager->FillH1(54, position.y(), evweght);
    analysisManager->FillH1(55, position.z(), evweght);
    
    analysisManager->FillH1(56, momentum.x(), evweght);
    analysisManager->FillH1(57, momentum.y(), evweght);
    analysisManager->FillH1(58, momentum.z(), evweght);
    
    G4double zfocus = 0.0;  // in IP
    G4double focusx = momentum.x()/momentum.z() * (zfocus - position.z()) + position.x();
    G4double focusy = momentum.y()/momentum.z() * (zfocus - position.z()) + position.y();
    analysisManager->FillH1(59, focusx, evweght);
    analysisManager->FillH1(60, focusy, evweght);

    analysisManager->FillH2(10, position.x(), momentum.x()/momentum.z(), evweght);
    analysisManager->FillH2(11, focusx, momentum.x()/momentum.z(), evweght);
//     analysisManager->FillH2(11, position.y(), momentum.y()/momentum.z());
  }
}



void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 

  G4double evweght = 1.0;
  const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
  if (evinf) evweght = evinf->GetWeight();

  G4ThreeVector position = aTrack->GetPosition();
  G4ThreeVector vertex   = aTrack->GetVertexPosition();  
  G4double charge        = aTrack->GetDefinition()->GetPDGCharge();

  G4bool transmit = ((position.z() >= fZendAbs) && (vertex.z() < fZendAbs));
  G4bool reflect  =  (position.z() <= fZstartAbs);
  G4bool notabsor = (transmit || reflect);

  //transmitted + reflected particles counter
  //
  G4int flag = 0;
  if (charge == fPrimaryCharge)  flag = 1;
  if (aTrack->GetTrackID() == 1) flag = 2;
  if (transmit) fEventAction->SetTransmitFlag(flag);
  if (reflect)  fEventAction->SetReflectFlag(flag);
      
  //
  //histograms
  //
  G4bool charged  = (charge != 0.);
  G4bool neutral = !charged;

  //energy spectrum at exit
  //
  G4int id = 0;  
       if (transmit && charged) id = 10;
  else if (transmit && neutral) id = 20;
  else if (reflect  && charged) id = 30;
  else if (reflect  && neutral) id = 40;

  if (id>0)   
    analysisManager->FillH1(id, aTrack->GetKineticEnergy(), evweght);
  if (id == 10) fEventAction->AddKinEnegyCharged(aTrack->GetKineticEnergy());
  if (id == 20) fEventAction->AddKinEnegyGamma(aTrack->GetKineticEnergy());
  //energy leakage
  //
  if (notabsor) {
    G4int trackID = aTrack->GetTrackID();
    G4int index = 0; if (trackID > 1) index = 1;    //primary=0, secondaries=1
    G4double eleak = aTrack->GetKineticEnergy();
    if ((aTrack->GetDefinition() == G4Positron::Positron()) && (trackID > 1))
      eleak += 2*electron_mass_c2;
    run->AddEnergyLeak(eleak,index);  
  }

  //space angle distribution at exit : dN/dOmega
  //
  G4ThreeVector direction = aTrack->GetMomentumDirection();    
  id = 0;   
       if (transmit && charged) id = 12;
  else if (transmit && neutral) id = 22;
  else if (reflect  && charged) id = 32;
  else if (reflect  && neutral) id = 42;

  if (id>0) {
    G4double theta  = std::acos(direction.z());
    if (theta > 0.0) {
      G4double dteta  = analysisManager->GetH1Width(id);
      G4double unit   = analysisManager->GetH1Unit(id);    
      G4double weight = (unit*unit)/(twopi*std::sin(theta)*dteta);
      /*
      G4cout << "theta, dteta, unit, weight: " 
             << theta << "  "   
             << dteta << "  "   
             << unit << "  "   
             << weight << G4endl;   
      */
      analysisManager->FillH1(id,theta, evweght*weight);
    }
  }
  
  //energy fluence at exit : dE(MeV)/dOmega
  //
  id = 0;  
  if (transmit && charged) id = 11;
  else if (transmit && neutral) id = 21;
  else if (reflect  && charged) id = 31;
  else if (reflect  && neutral) id = 41;

  if (id>0) {
    G4double theta  = std::acos(direction.z());
    if (theta > 0.0) {
      G4double dteta  = analysisManager->GetH1Width(id);
      G4double unit   = analysisManager->GetH1Unit(id);    
      G4double weight = (unit*unit)/(twopi*std::sin(theta)*dteta);
      weight *= (aTrack->GetKineticEnergy()/MeV); 
      analysisManager->FillH1(id,theta, evweght*weight);
    } 
  }
  
  //projected angles distribution at exit
  //
  id = 0;   
  if (transmit && charged) id = 13;
  else if (transmit && neutral) id = 23;
  else if (reflect  && charged) id = 33;
  else if (reflect  && neutral) id = 43;

  if (id>0) {
    if (direction.z() != 0.0) {
      G4double tet = std::atan(direction.x()/std::fabs(direction.z()));
      analysisManager->FillH1(id,tet, evweght);
      if (transmit && (flag == 2)) run->AddMscProjTheta(tet);

      tet = std::atan(direction.y()/std::fabs(direction.z()));
//       analysisManager->FillH1(id,tet, evweght);
      if (transmit && (flag == 2)) run->AddMscProjTheta(tet);
    }
  }

  //projected position and radius at exit
  //
  id = 0;   
  if (transmit && charged) id = 14;
  else if (transmit && neutral) id = 24;
//   if (reflect && charged) id = 14;
//   else if (reflect && neutral) id = 24;
  
  if (id>0) {
    G4double y = position.y(), x = position.x();
    G4double r = std::sqrt(y*y + x*x);
    analysisManager->FillH1(id,   x, evweght);
    analysisManager->FillH1(id+1, r, evweght);
    if (id == 14) analysisManager->FillH2(2, x, y, evweght);
    if (id == 24) analysisManager->FillH2(3, x, y, evweght);
  }
  
  //x-vertex of charged secondaries
  //
  if ((aTrack->GetParentID() == 1) && charged) {
    G4double xVertex = (aTrack->GetVertexPosition()).x();
    analysisManager->FillH1(6, xVertex, evweght);
    if (notabsor) analysisManager->FillH1(7, xVertex, evweght); 
  }
  
//energy spectrum vs position at exit
  //
  id = -1;  
  if (transmit && charged) id = 0;

  if (id>=0) {
    G4double x = position.x();
    analysisManager->FillH2(id, x, aTrack->GetKineticEnergy(), evweght);
//     G4double y = position.y();
//     G4double r = std::sqrt(y*y + x*x);
//     analysisManager->FillH1(id,   y, evweght);
//     analysisManager->FillH1(id+1, r, evweght);
  }

  //z-vertex of charged secondaries
  //
  id = -1;
  fZendMagnet = fDetector->GetMagnetZend();
  if (charged && vertex.z() > fZendMagnet) id = 51;
  if (neutral && vertex.z() > fZendMagnet) id = 52;
  if (id>=0) {
    if ((aTrack->GetParentID() == 1) ) {
      analysisManager->FillH1(id, vertex.z());
//       if (id == 51) analysisManager->FillH2(4, vertex.z(), position.x());
//       if (id == 52) analysisManager->FillH2(5, vertex.z(), position.x());
      if (id == 51) analysisManager->FillH2(4, vertex.z(), aTrack->GetKineticEnergy());
      if (id == 52) analysisManager->FillH2(5, vertex.z(), aTrack->GetKineticEnergy());
    }
  }

  G4double z_threshold = 0.95*0.5*fDetector->GetWorldSizeZ();
  if (position.z() <= z_threshold) return;
  
  G4int pdgid = aTrack->GetDefinition()->GetPDGEncoding(); 
  G4ThreeVector momentum = aTrack->GetMomentum();
// electron at exit
  if (pdgid == G4Electron::Electron()->GetPDGEncoding() && position.z() > z_threshold) {
    analysisManager->FillH1(61, aTrack->GetGlobalTime(), evweght);
    analysisManager->FillH2(6, aTrack->GetKineticEnergy(), momentum.getTheta(), evweght);
    analysisManager->FillH1(65, aTrack->GetKineticEnergy(), evweght);
    analysisManager->FillH1(68, position.x(), evweght);
//     analysisManager->FillH2(11, position.x(), momentum.x()/momentum.z(), evweght);
  }

// positron at exit
  if (pdgid == G4Positron::Positron()->GetPDGEncoding() && position.z() > z_threshold) {
    analysisManager->FillH1(67, aTrack->GetKineticEnergy(), evweght);
    analysisManager->FillH1(70, position.x(), evweght);
    analysisManager->FillH2(12, aTrack->GetKineticEnergy(), momentum.getTheta(), evweght);
  }
  
// photon at exit
  if (pdgid == G4Gamma::Gamma()->GetPDGEncoding() && position.z() > z_threshold) {
    analysisManager->FillH1(62, aTrack->GetGlobalTime(), evweght);
    analysisManager->FillH2(7, aTrack->GetKineticEnergy(), momentum.getTheta(), evweght);
    analysisManager->FillH2(8, aTrack->GetKineticEnergy(), momentum.getTheta(), evweght);
//     analysisManager->FillH2(9, position.perp(), momentum.getTheta(), evweght);
    analysisManager->FillH2(9, position.x(),position.y(), evweght);

    analysisManager->FillH1(66, aTrack->GetKineticEnergy(), evweght);
    analysisManager->FillH1(69, position.x(), evweght);
  }
  
  const PrimaryGeneratorAction *gen = dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4ThreeVector source_center = G4ThreeVector( gen->GetX0(), gen->GetY0(), fabs(0.5 * fDetector->GetWorldSizeZ() - gen->GetZ0()) );

  G4double gtime = aTrack->GetGlobalTime();
  G4double time0 = source_center.z() / c_light;
        
  G4ThreeVector pos_time = position;
  if (aTrack->GetTotalEnergy() > 0.0) {
    pos_time += momentum/aTrack->GetTotalEnergy() * c_light*(time0-gtime);
  }
  pos_time -= G4ThreeVector(0.0, 0.0, fDetector->GetWorldSizeZ()/2.0);

  if (fVerbose) {
    G4cout << "Particle : " << pdgid << G4endl;
    G4cout << "Position : " << position << G4endl;
    G4cout << "Energy, Momentum : " << aTrack->GetKineticEnergy() << "   " << momentum << G4endl;
    G4cout << "Speed : " << momentum/aTrack->GetKineticEnergy() << G4endl;
    G4cout << "gtime, time : " << gtime << "   " << time0 << G4endl;
    G4cout << "Position new : " << pos_time << G4endl;
  }

// Tree output  
   G4double txcut = run->GetTreeCutX();
   G4double tycut = run->GetTreeCutY();
   G4ParticleDefinition  *tparticle = run->GetTreeParticle();
   
  if (   ( (!tparticle) ||  (tparticle && pdgid == tparticle->GetPDGEncoding()) )
      && (position.z() > z_threshold) && (abs(pos_time.x()) < txcut && abs(pos_time.y()) < tycut) ) {
          
      analysisManager->FillNtupleDColumn(0, 0, aTrack->GetKineticEnergy()/GeV);
      analysisManager->FillNtupleDColumn(0, 1, pos_time.x());
      analysisManager->FillNtupleDColumn(0, 2, pos_time.y());
      analysisManager->FillNtupleDColumn(0, 3, pos_time.z());

//     analysisManager->FillNtupleDColumn(0, 1, position.x());
//     analysisManager->FillNtupleDColumn(0, 2, position.y());
//     analysisManager->FillNtupleDColumn(0, 3, position.z());

      analysisManager->FillNtupleDColumn(0, 4, momentum.x()/GeV);
      analysisManager->FillNtupleDColumn(0, 5, momentum.y()/GeV);
      analysisManager->FillNtupleDColumn(0, 6, momentum.z()/GeV);
    
      analysisManager->FillNtupleDColumn(0, 7, momentum.getTheta());
      analysisManager->FillNtupleDColumn(0, 8, momentum.getPhi());

      analysisManager->FillNtupleIColumn(0, 9, pdgid);
      
      G4int ptype = 0, psubtype = 0;
      const G4VProcess *pproc = aTrack->GetCreatorProcess();
      if (pproc) {
        ptype = pproc->GetProcessType();
        psubtype = pproc->GetProcessSubType();
      }
      
      analysisManager->FillNtupleIColumn(0, 10, ptype*1000 + psubtype);
      analysisManager->FillNtupleDColumn(0, 11, evweght);
      
      analysisManager->AddNtupleRow();
  }
}



