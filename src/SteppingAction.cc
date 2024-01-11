// 
/// \brief Implementation of the SteppingAction class
//

#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NavigationHistory.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4StepPoint.hh"
#include "G4Event.hh"
#include "G4Exception.hh"

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "LXSetUp.hh"

#include "G4Step.hh"

#include "G4AnalysisManager.hh"


SteppingAction::SteppingAction(DetectorConstruction* DET, EventAction* EA)
:G4UserSteppingAction(),fDetector(DET), fEventAction(EA)
{ }



SteppingAction::~SteppingAction()
{ }



void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  if (preStepPoint->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber()) ProcessInAbsorber(aStep);

  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
  if (!postStepPoint->GetTouchableHandle()->GetVolume()) {
    return;
  }

#ifdef LX_TRACK_TIMING_CUT
  if (postStepPoint->GetGlobalTime() > LX_TRACK_TIMING_CUT) {
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    return;
  }
#endif

  LXSetUp *lxs = LXSetUp::Instance();

  if (lxs->ScintCerenkovPhysics){  // && (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "opticalphoton")){
    ProcessScintCerenkov(aStep);
    if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "opticalphoton"){return;}
  }

//   const G4NavigationHistory *tchistps = postStepPoint->GetTouchable()->GetHistory();
//   const G4NavigationHistory *tchistpr = preStepPoint->GetTouchable()->GetHistory();
//   G4int vdepthps = tchistps->GetDepth();
//   G4int vdepthpr = tchistpr->GetDepth();
//   std::cout << "PreStepPoint : " << tchistpr->GetVolume(vdepthpr)->GetName() << "  depth: " << vdepthpr
//             << "  copy: " << tchistpr->GetVolume(vdepthpr)->GetCopyNo()
//             << "          PostStepPoint : " << tchistps->GetVolume(vdepthps)->GetName() << "  depth: " << vdepthps
//             << "  copy: " << tchistps->GetVolume(vdepthps)->GetCopyNo()
//             << "  TrackId: " << aStep->GetTrack()->GetTrackID() << std::endl;

  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  const auto &intercetvolmap = run->GetInterceptVolumes();

  for (const auto &volid : intercetvolmap) {
    //Check the particle enters OpppDetContainer
    G4int vdps = CheckPointHistory(postStepPoint, volid.first);
    if (vdps) {
      G4int vdpr = CheckPointHistory(preStepPoint, volid.first);
      G4int vdpscpn(0), vdprcpn(0);
      if (vdpr) {
        vdpscpn = postStepPoint->GetTouchableHandle()->GetHistory()->GetVolume(vdps)->GetCopyNo();
        vdprcpn = preStepPoint->GetTouchableHandle()->GetHistory()->GetVolume(vdpr)->GetCopyNo();
      }

      if (!vdpr || (vdpscpn != vdprcpn)) {
//         std::cout << "Container: " << preStepPoint->GetTouchableHandle()->GetVolume()->GetName()
//                   << "  " << postStepPoint->GetTouchableHandle()->GetVolume()->GetName()
//                   << "  " << aStep->GetTrack()->GetTrackID() << std::endl;

#ifdef LX_FAST_SIM
        if (volid.first=="BeamDumpAssembly" || volid.first=="Shielding" || volid.first=="GammaBeamDumpAssembly"
             || volid.first.contains("HICSDumpAssembly") || volid.first.contains("HICSDump2T")
             || volid.first.contains("ComptShieldingConcrete") || volid.first.contains("ComptShieldingPlate")
            ) {
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          break;
        }
#endif
        ProcessInDetector(aStep, volid.second, vdps);
      }
    }
  }

  const auto &sensitivevolmap = run->GetSensitiveVolumes();

  G4TouchableHandle aTouchable = preStepPoint->GetTouchableHandle();
//     G4TouchableHandle aTouchable = postStepPoint->GetTouchableHandle();
  G4String touchableVolName = aTouchable->GetVolume()->GetName();
  for (const auto &volid : sensitivevolmap) {
    if (touchableVolName == volid.first) {
//       std::cout << "Sensitive volume: " << volid.first << " Position: "
//                 << touchableVolName << std::endl;
      ProcessInTracker(aStep, volid.first, std::get<0>(volid.second), std::get<1>(volid.second), std::get<2>(volid.second));
    }
  }

  if (aStep->GetTrack()->GetTrackID() == 1) {
    if (run->GetSaveTrajectory()) ProcessPrimaryTrack(aStep);
  }

}



G4int SteppingAction::CheckPointHistory(const G4StepPoint *sp, const G4String vtname) const
{
  const G4NavigationHistory *tchist = sp->GetTouchableHandle()->GetHistory();
  G4int vdepth = tchist->GetDepth();
  G4String vname = tchist->GetVolume(vdepth)->GetName();
//   G4bool psopppin = vname.contains(vtname);
  G4bool psopppin = (vname == vtname);
  while (vdepth > 0 && !psopppin) {
    --vdepth;
    vname = tchist->GetVolume(vdepth)->GetName();
//     psopppin = vname.contains(vtname);
    psopppin = (vname == vtname);
  }
  return vdepth;
}



void SteppingAction::ProcessInAbsorber(const G4Step* aStep)
{
  fEventAction->AddEnergy (aStep->GetTotalEnergyDeposit());

  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  if (charge != 0.) {
    fEventAction->AddTrakLenCharg(aStep->GetStepLength());
    fEventAction->CountStepsCharg();
  } else {
    fEventAction->AddTrakLenNeutr(aStep->GetStepLength());
    fEventAction->CountStepsNeutr();
  }
}



void SteppingAction::ProcessInDetector(const G4Step* aStep, const std::vector<G4int> &dettype, const G4int vdepth)
{
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4StepPoint* stepp = aStep->GetPostStepPoint();
  G4Track* aTrack  = aStep->GetTrack();
  const G4TrackVector *sectrv = aStep->GetSecondary();
  G4int nsecondaries = sectrv->size();
  G4double esecondaries = 0.0;
//  std::for_each(sectrv->begin(), sectrv->end(), [&](const G4Track *t) {esecondaries += t->GetVertexKineticEnergy();});
  std::for_each(sectrv->begin(), sectrv->end(), [&](const G4Track *t) {esecondaries += t->GetKineticEnergy();});

//   G4ThreeVector position = aTrack->GetPosition();
//   G4ThreeVector momentum   = aTrack->GetMomentum();
//   G4double gtime = aTrack->GetGlobalTime();
  G4ThreeVector position = stepp->GetPosition();
  G4ThreeVector momentum   = stepp->GetMomentum();
  G4double gtime = stepp->GetGlobalTime();

  G4ThreeVector vertex   = aTrack->GetVertexPosition();
  G4int pdgid = aTrack->GetDefinition()->GetPDGEncoding();
  G4int trckid = aTrack->GetTrackID();
  G4int ptrackid = aTrack->GetParentID();
  if (trckid==1) {
    const EventInfo *evinf = dynamic_cast<EventInfo*>(
                             G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
    if (evinf) ptrackid = evinf->GetMCTrackId();
    else ptrackid = -1;
  }

  G4TouchableHandle theTouchable = aStep->GetPostStepPoint()->GetTouchableHandle();
  G4ThreeVector local = theTouchable->GetHistory()->GetTransform(vdepth).TransformPoint(position);

  G4int detid = dettype[0] + theTouchable->GetHistory()->GetVolume(vdepth)->GetCopyNo();
  G4int bdec = 1000;
  for (auto iv = dettype.cbegin()+1; iv != dettype.cend(); ++iv) {
    detid += bdec * theTouchable->GetHistory()->GetVolume(vdepth - (*iv))->GetCopyNo();
    bdec *= 10;
  }

//   std::cout << "ProcessInDetector: " << dettype[0] << "  depth " << vdepth << "  copy number: "
//             << theTouchable->GetHistory()->GetVolume(vdepth)->GetCopyNo() << "  Touchable copy number: "
//             << theTouchable->GetCopyNumber() << "  final detid: " << detid << std::endl;

  G4double ivecut = run->GetInterceptVolECut(dettype[0]);
  if ( aTrack->GetKineticEnergy() > ivecut) {
   // if ( true) {
    G4int ptype = 0, psubtype = 0;
    const G4VProcess *pproc = aTrack->GetCreatorProcess();
    if (pproc) {
      ptype = pproc->GetProcessType();
      psubtype = pproc->GetProcessSubType();
//      std::cout << "Process " << pproc->GetProcessName() << "  type/subtype: " <<  ptype << "   " << psubtype
//                << "  event: " << run->GetNumberOfEvent() << std::endl;
    }

    fEventAction->AddVolumeTrack(std::make_tuple(
                    aTrack->GetKineticEnergy()/GeV,
                    position.x(),
                    position.y(),
                    position.z(),
                    gtime,
                    vertex.x(),
                    vertex.y(),
                    vertex.z(),
                    momentum.x()/GeV,
                    momentum.y()/GeV,
                    momentum.z()/GeV,
                    momentum.getTheta(),
                    momentum.getPhi(),
                    local.x(),
                    local.y(),
                    local.z(),
                    aTrack->GetTrackID(),
                    detid,
                    pdgid,
                    ptype*1000 + psubtype,
                    ptrackid,
                    nsecondaries,
                    esecondaries/GeV));
  }
}



void SteppingAction::ProcessInTracker(const G4Step* aStep, const G4String &dname, const G4int dettype,
                                      const G4int detdepth, const G4int layerdepth)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  // make sure you aren't using geantinos
  if (edep <= 0.0 /*&& aStep->GetTrack()->GetDefinition()->GetParticleType() != "geantino"*/) {
//     std::cout << "No energy deposit!\n";
    return;
  }

  const G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  const G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
  G4Track* aTrack  = aStep->GetTrack();
  G4ThreeVector globalHitPos = 0.5 * (preStepPoint->GetPosition() + postStepPoint->GetPosition());
//   G4ThreeVector globalHitPos = postStepPoint->GetPosition();

  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector local = theTouchable->GetHistory()->GetTopTransform().TransformPoint(globalHitPos);

////////////////////////////////////////
//   G4AffineTransform transform = theTouchable->GetHistory()->GetTopTransform();
//   G4RotationMatrix mtrx = transform.NetRotation();
//   G4ThreeVector trns = transform.NetTranslation();
//   G4cout << "ProcessInTracker Translation " << trns << G4endl;
//   G4cout << "ProcessInTracker Rotation, Euler Angles  phi/theta/psi: " << mtrx.getPhi()
//          << "  " << mtrx.getTheta() << "  " << mtrx.getPsi() << G4endl;
//   G4cout << "Local position : " << local << "   Global: " << globalHitPos << G4endl;
//   G4cout << "Affine Transformation:  " << transform << G4endl;
////////////////////////////////////////

  G4int vdepth = theTouchable->GetHistory()->GetDepth();
  G4int layerid  = theTouchable->GetHistory()->GetVolume(vdepth - layerdepth)->GetCopyNo();
  G4int pdetid = dettype + theTouchable->GetHistory()->GetVolume(vdepth - detdepth)->GetCopyNo();

//   G4ThreeVector momentum   = aTrack->GetMomentum();
//   G4double gtime = aTrack->GetGlobalTime();
//   G4ThreeVector momentum   = stepp->GetMomentum();
//   G4double gtime = stepp->GetGlobalTime();
//   G4ThreeVector vertex   = aTrack->GetVertexPosition();
//   G4int pdgid = aTrack->GetDefinition()->GetPDGEncoding();

  G4int cellx(-1), celly(-1);

  const auto &sseg = fDetector->GetSensorSegmentation(dname);
  G4double xrange = std::get<0>(sseg);
  G4double yrange = std::get<1>(sseg);
  G4int ncellx = std::get<2>(sseg);
  G4int ncelly = std::get<3>(sseg);

  cellx = static_cast<int>((0.5 * xrange + local.x()) / xrange * ncellx);
  celly = static_cast<int>((0.5 * yrange + local.y()) / yrange * ncelly);

//   std::cout << "Added e_dep: " << pdetid << "  " << layerid << "  " << cellx << "  "
//             << celly << "  " << edep << std::endl;
  fEventAction->AddOpppTrackerCellDepposit(pdetid, layerid, cellx, celly, edep, aTrack);
}



void SteppingAction::ProcessPrimaryTrack(const G4Step* aStep)
{
//   const G4StepPoint *stepp = aStep->GetPostStepPoint();
  const G4StepPoint *stepp = aStep->GetPreStepPoint();
  G4Track* aTrack = aStep->GetTrack();
  G4ThreeVector position = stepp->GetPosition();
  G4double gtime = stepp->GetGlobalTime();
//   G4ThreeVector momentum = aTrack->GetMomentum();
//   G4double energyp = aTrack->GetKineticEnergy();
  G4ThreeVector momentum = stepp->GetMomentum();
  G4double energyp = stepp->GetKineticEnergy();
  G4int pdgid = aTrack->GetDefinition()->GetPDGEncoding();
  G4int ptrackid = aTrack->GetParentID();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4int nsecondaries = aStep->GetNumberOfSecondariesInCurrentStep();

  G4int ptype = 0, psubtype = 0;
  const G4VProcess *pproc = stepp->GetProcessDefinedStep();
  if (pproc) {
    ptype = pproc->GetProcessType();
    psubtype = pproc->GetProcessSubType();
  }

  G4int trcktype = 0;

  fEventAction->AddTrajectoryStep(aTrack->GetTrackID(),
                    std::make_tuple(position.x(), position.y(), position.z(), gtime,
                    energyp/GeV, momentum.x()/GeV, momentum.y()/GeV, momentum.z()/GeV,
                    edep, ptype*1000 + psubtype, nsecondaries),
                    std::make_tuple(trcktype, pdgid, ptrackid) );
}


void SteppingAction::ProcessScintCerenkov(const G4Step* aStep)
{ 
  if(!aStep->IsFirstStepInVolume()){return;}

  // G4cout << "First step in IP lanex or Cerenkov: " << (theTouchable->GetVolume()->GetName() << G4endl;

      LXSetUp *lxs = LXSetUp::Instance();

 G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();

 //G4cout << "First step in lanex: " << theTouchable->GetVolume()->GetName() << G4endl;
							
 G4int dettype =-1; //0 for cherenkov channel, 1 for scint screen, 2 for scint. camera
 if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "CerenkovStrawInnerPhysical")) // && (aStep->GetTrack()->GetCurrentStepNumber()==1))
   {dettype=0;}
 else if ((G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "scintPhosphorPhysical") ||  G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "LysoCal"))) //&& (aStep->GetTrack()->GetCurrentStepNumber()==1))
   {dettype=1;}
 else if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "ScintCameraApertureInnerPhysical"))
   {dettype=2;}
 else {return;}


  //  G4cout << "Step in IP lanex or Cerenkov" << G4endl;

 G4ThreeVector globalHitPos = aStep->GetPreStepPoint()->GetPosition();
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 G4int CopyNo, MotherCopyNo, DetectorCopyNo;
 G4ThreeVector local;;
 // For Cherenkov inner refractive medium volume -> Straw -> Detector
 if (dettype==0){
   CopyNo = theTouchable->GetCopyNumber();
   theTouchable->MoveUpHistory();
   MotherCopyNo = theTouchable->GetCopyNumber(); 

   theTouchable->MoveUpHistory();
   theTouchable->MoveUpHistory();
local = theTouchable->GetHistory()->GetTopTransform().TransformPoint(globalHitPos);
   DetectorCopyNo = theTouchable->GetCopyNumber();}
 
 else if  (dettype==1){

   local = theTouchable->GetHistory()->GetTopTransform().TransformPoint(globalHitPos);

   if(G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "LysoCal")){DetectorCopyNo = theTouchable->GetCopyNumber() + 2;}
   else {
   theTouchable->MoveUpHistory();
   if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "Brem")){DetectorCopyNo=0;}
   theTouchable->MoveUpHistory();
   if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "HICS")){DetectorCopyNo=1;}
   }
 }
 else if  (dettype==2){

   theTouchable->MoveUpHistory();
      MotherCopyNo = theTouchable->GetCopyNumber();  

   if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "Brem")){DetectorCopyNo=0;}
   if (G4StrUtil::contains(theTouchable->GetVolume()->GetName(), "HICS")){DetectorCopyNo=1; MotherCopyNo -= 2;}
 }


   // cherenkov det copy 0 = brem, 1 = hics ip, 2 & 3 = gamma spectrometer

 G4double evweight = 1.0;
 const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
 if (evinf){ evweight = evinf->GetWeight();
   //   normalising for event weight
 }

  if(aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
    if(aStep->GetTrack()->GetCurrentStepNumber()==1){
 if( dettype == 0){ // cherenkov
   G4int norm_mcn;
   if (MotherCopyNo<lxs->CerenkovChannels/lxs->CerenkovChannelLayers){norm_mcn = lxs->CerenkovChannelLayers*MotherCopyNo+3;} 
   else if (MotherCopyNo>(lxs->CerenkovChannels/lxs->CerenkovChannelLayers - 1) && MotherCopyNo<2*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)){norm_mcn =lxs->CerenkovChannelLayers *(MotherCopyNo-lxs->CerenkovChannels/lxs->CerenkovChannelLayers) + 2;} 
   else if (MotherCopyNo>2*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)-1 && MotherCopyNo<3*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)){norm_mcn = lxs->CerenkovChannelLayers*(MotherCopyNo-2*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)) + 1;} 
   else if (MotherCopyNo>3*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)-1 && MotherCopyNo<4*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers)){norm_mcn = lxs->CerenkovChannelLayers*(MotherCopyNo-3*(lxs->CerenkovChannels/lxs->CerenkovChannelLayers));} 

     analysisManager->FillH1(71+DetectorCopyNo, norm_mcn, evweight);
     //     G4cout << "cherenkov emission! " << local.x() << ", " << local.y() << G4endl;
     analysisManager->FillH2(17+DetectorCopyNo, local.x(), local.y(), evweight);

 }
 
 if( dettype == 1){ // scintillator
   if (abs(local.y())<5.0){  // selecting only central band =/- 5mm
      analysisManager->FillH1(75+DetectorCopyNo, local.x(), evweight/(lxs->ScintPhysicsYield));}
      analysisManager->FillH2(13+DetectorCopyNo, local.x(), local.y(), evweight/(lxs->ScintPhysicsYield));
 }}
 if( dettype == 2){ // camera
   analysisManager->FillH1(79+DetectorCopyNo, MotherCopyNo, evweight/(lxs->ScintPhysicsYield));}
  }

  
  else{      // for normal particles (non-optical photons)

    G4double E = aStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy();

    if( dettype == 0){ // cherenkov

    analysisManager->FillH2(25+DetectorCopyNo, local.x(), E, evweight);}

 if( dettype == 1){ // scintillator

   if (abs(local.y())<5.0){  // selecting only central band =/- 5mm
     analysisManager->FillH2(21+DetectorCopyNo, local.x(), E, evweight);}
 }

  }

}
 

