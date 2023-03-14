//
/// \brief Implementation of the EventAction class
//
#include <numeric>
#include <bitset>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"
#include "LXSetUp.hh"


EventAction::EventAction()
:G4UserEventAction(),
 fEnergyDeposit(0.),
 fTrakLenCharged(0.), fTrakLenNeutral(0.),
 fNbStepsCharged(0), fNbStepsNeutral(0),
 fTransmitFlag(0), fReflectFlag(0), fNHits(0)
{
  const int nDetIdFields = 4;
  fBSize.resize(nDetIdFields, 16);
  fBSizeSum.resize(nDetIdFields);
//   LXSetUp *lxs = LXSetUp::Instance();
//   fBSize[2] = static_cast<unsigned int>(ceil(log(lxs->OPPPTrackerNCellX)/log(2.0)));
//   fBSize[3] = static_cast<unsigned int>(ceil(log(lxs->OPPPTrackerNCellY)/log(2.0)));

  std::partial_sum (fBSize.begin(), fBSize.end(), fBSizeSum.begin());

  if (fBSizeSum[nDetIdFields-1] > DET_ENCODING_BIT_SIZE) {
    G4String msgstr("Number of bits required to encode sensitive detector is greater than the DET_ENCODING_BIT_SIZE!\n");
    G4Exception("EventAction::", "EventAction()", FatalException, msgstr.c_str());
  }   
}



EventAction::~EventAction()
{ }



int EventAction::TestSegmentationEncoding(const int maxx, const int maxy)
{
  if (fBSize.size() < 2) return -1;
  int yind = fBSize.size() - 1;
  int xind = yind - 1;
  if (maxx >= pow(2.0, fBSize[xind])) return 1;
  if (maxy >= pow(2.0, fBSize[yind])) return 1;
  return 0;
}



void EventAction::BeginOfEventAction(const G4Event* )
{
 // initialisation per event
  fEnergyDeposit  = 0.;
  fTrakLenCharged = fTrakLenNeutral = 0.;
  fNbStepsCharged = fNbStepsNeutral = 0;
  fTransmitFlag   = fReflectFlag    = 0;

  fKinetikEnergyGamma.clear();
  fKinetikEnergyCharged.clear();

//   ClearData();
}



void EventAction::EndOfEventAction(const G4Event *event)
{
   Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->AddEnergy(fEnergyDeposit);
  run->AddTrakLenCharg(fTrakLenCharged);
  run->AddTrakLenNeutr(fTrakLenNeutral);

  run->CountStepsCharg(fNbStepsCharged);
  run->CountStepsNeutr(fNbStepsNeutral);

  run->CountTransmit (fTransmitFlag);
  run->CountReflect  (fReflectFlag);

  G4AnalysisManager::Instance()->FillH1(1, fEnergyDeposit);

  G4double evweght = 1.0;
  const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
  if (evinf) evweght = evinf->GetWeight();
  G4int nvtxp = event->GetPrimaryVertex()->GetNumberOfParticle();
  if (nvtxp) G4AnalysisManager::Instance()->FillH1(0, 0.0, evweght*nvtxp);

  G4double npart = static_cast<G4double>(fKinetikEnergyCharged.size() * fKinetikEnergyGamma.size());
  for (std::vector<G4double>::const_iterator itre = fKinetikEnergyCharged.begin(); itre != fKinetikEnergyCharged.end(); ++itre) {
    for (std::vector<G4double>::const_iterator itrg = fKinetikEnergyGamma.begin(); itrg != fKinetikEnergyGamma.end(); ++itrg) {
      G4AnalysisManager::Instance()->FillH2(1, *itre, *itrg, 1.0/npart);
    }
  }

  if(fHitId.size() > 0) {
    UpdateOPPPTracker();
  }

  if (fVolumeTracks.size() > 0) {
    UpdateVolumeTracks();
  }

  if (fTrajTracks.size() > 0) {
    UpdateTrajectories(event);
  }

  ClearData();

}



void  EventAction::AddOpppTrackerCellDepposit(const G4int opppdet, const G4int layerid, const G4int cellx, const G4int celly,
                                              const G4double edep, const G4Track* aTrack)
{
  std::bitset<DET_ENCODING_BIT_SIZE> detcellid(0);
  detcellid |= ((1UL << fBSize[0]) - 1) & opppdet ;
  detcellid |= ( ((1UL << fBSize[1]) - 1) & layerid ) << fBSizeSum[0];
  detcellid |= ( ((1UL << fBSize[2]) - 1) & cellx ) << fBSizeSum[1];
  detcellid |= ( ((1UL << fBSize[3]) - 1) & celly ) << fBSizeSum[2];

//   std::cout << "AddOpppTrackerCellDepposit: " << opppdet << "  " << layerid << "  " << cellx << "  " << celly << "  " <<  edep << std::endl;
// 
//   int opppdet1 = ( detcellid & std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[0]) - 1) ).to_ulong();
//   int layerid1 = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[1]) - 1) & (detcellid >> fBSizeSum[0])).to_ulong();
//   int cellx1 = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[2]) - 1) & (detcellid >> fBSizeSum[1])).to_ulong();
//   int celly1 = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[3]) - 1) & (detcellid >> fBSizeSum[2])).to_ulong();
// 
//   std::cout << "Updating tracker tree: " << opppdet1 << "  " << layerid1 << "  " << cellx1 << "  " << celly1 << std::endl;

  // Update/insert hit
  G4int trackid = aTrack->GetTrackID();
  auto eitr = fHitId.find(detcellid.to_string());
  if (eitr != fHitId.end()) {
    G4int hitId  = eitr->second;
    fHitEDep[hitId] += edep;
    const auto htrkitr = std::find(fHitTrackId[hitId].begin(), fHitTrackId[hitId].end(), trackid);
    if (htrkitr == fHitTrackId[hitId].end()) {
      fHitTrackGPos[hitId].push_back(aTrack->GetPosition());
      fHitTrackId[hitId].push_back(trackid);
      fHitTrackGTime[hitId].push_back(aTrack->GetGlobalTime());
      fHitTrackEDep[hitId].push_back(edep);
    } else {
      fHitTrackEDep[hitId].at(htrkitr-fHitTrackId[hitId].begin()) += edep;
    }

  } else {
    fHitId[detcellid.to_string()] = fNHits; 
    fHitEDep[fNHits] = edep;
    fHitTrackGPos[fNHits].push_back(aTrack->GetPosition());
    fHitTrackId[fNHits].push_back(trackid);
    fHitTrackGTime[fNHits].push_back(aTrack->GetGlobalTime());
    fHitTrackEDep[fNHits].push_back(edep);
    ++fNHits;
  }

  // Update/insert track
  //if ( aTrack->GetKineticEnergy() > 0.1*GeV ) {
  if ( aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy() > -10.0*GeV ) {

    if (std::find(ftrackId.begin(), ftrackId.end(), trackid) == ftrackId.end()) {
      ftrackId.push_back(trackid);
      ftrackVtxPos.push_back(aTrack->GetVertexPosition());
      ftrackMomentum.push_back(aTrack->GetStep()->GetPreStepPoint()->GetMomentum()/GeV);
      ftrackE.push_back(aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy()/GeV);
//       ftrackMomentum.push_back(aTrack->GetVertexMomentumDirection()/GeV);
//       ftrackE.push_back(aTrack->GetVertexKineticEnergy()/GeV);
      ftrackPDG.push_back(aTrack->GetDefinition()->GetPDGEncoding());
      G4int ptrackid = aTrack->GetParentID();
      if (trackid==1) {
        const EventInfo *evinf = dynamic_cast<EventInfo*>(
                                 G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
        if (evinf) ptrackid = evinf->GetMCTrackId();
        else ptrackid = -1;
      }
      ftrackPTD.push_back(ptrackid);

      G4int ptype = 0, psubtype = 0;
      const G4VProcess *pproc = aTrack->GetCreatorProcess();
      if (pproc) {
        ptype = pproc->GetProcessType();
        psubtype = pproc->GetProcessSubType();
      }
      ftrackPProc.push_back(ptype*1000 + psubtype);
    }
  }
}



void  EventAction::UpdateOPPPTracker()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  HistoManager *histManager = run->GetHistoManager();

  G4double evweght = 1.0;
  const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
  if (evinf) evweght = evinf->GetWeight();

//   std::cout << "End of event. Updating tracker tree: \n";
  for (auto eitr = fHitId.begin(); eitr != fHitId.end(); ++eitr) {
    std::bitset<DET_ENCODING_BIT_SIZE> dd(eitr->first);
    G4int hitId = eitr->second;

    int opppdet = ( dd & std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[0]) - 1) ).to_ulong();
    int layerid = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[1]) - 1) & (dd >> fBSizeSum[0])).to_ulong();
    int cellx = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[2]) - 1) & (dd >> fBSizeSum[1])).to_ulong();
    int celly = ( std::bitset<DET_ENCODING_BIT_SIZE>((1UL << fBSize[3]) - 1) & (dd >> fBSizeSum[2])).to_ulong();

    if (!opppdet) opppdet = -1;  // This is just to write it to the Hits tree with the same semantic as in Tracks tree.
//    std::cout << "Updating tracker tree: " << opppdet << "  " << layerid << "  " << cellx << "  " << celly << "  " << eitr->second << std::endl;

    analysisManager->FillNtupleIColumn(2, 0, run->GetNumberOfEvent() + run->GetSkipEvents());
    analysisManager->FillNtupleIColumn(2, 1, opppdet);
    analysisManager->FillNtupleIColumn(2, 2, layerid);
    analysisManager->FillNtupleIColumn(2, 3, cellx);
    analysisManager->FillNtupleIColumn(2, 4, celly);
    analysisManager->FillNtupleDColumn(2, 5, fHitEDep[hitId]/GeV);
    analysisManager->FillNtupleIColumn(2, 6, hitId);
    histManager->SetHitTrackList(fHitTrackId[hitId]);
    histManager->SetHitTrackPosition(fHitTrackGPos[hitId]);
    histManager->SetHitTrackTime(fHitTrackGTime[hitId]);
    std::vector<G4double> tredep(fHitTrackEDep[hitId].size());
    std::transform(fHitTrackEDep[hitId].begin(), fHitTrackEDep[hitId].end(), tredep.begin(), 
                                                        std::bind2nd(std::divides<double>(), GeV));
    histManager->SetHitTrackEDep(tredep);
    analysisManager->FillNtupleDColumn(2, 13, evweght);
    analysisManager->FillNtupleIColumn(2, 14, run->GetRunID());
    analysisManager->AddNtupleRow(2);
  }

  if(fHitId.size() > 0){
    analysisManager->FillNtupleIColumn(3, 0, run->GetNumberOfEvent() + run->GetSkipEvents());
    histManager->SetTracks(ftrackId);
    histManager->SetTracksVtx(ftrackVtxPos);
    histManager->SetTracksMomentum(ftrackMomentum);
    histManager->SetTracksE(ftrackE);
    histManager->SetTracksPDG(ftrackPDG);
    histManager->SetTracksPhysProc(ftrackPProc);
    histManager->SetTracksPTId(ftrackPTD);
    analysisManager->FillNtupleDColumn(3, 12, evweght);
    analysisManager->FillNtupleIColumn(3, 13, run->GetRunID());
    analysisManager->AddNtupleRow(3);
  }
}



void EventAction::UpdateVolumeTracks()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  HistoManager *histManager = run->GetHistoManager();
  
  G4double evweght = 1.0;
  const EventInfo *evinf = dynamic_cast<EventInfo*>(G4RunManager::GetRunManager()->GetCurrentEvent()->GetUserInformation());
  if (evinf) evweght = evinf->GetWeight();
  
//   std::cout << "End of event. Updating volume tracks tree: \n";
  histManager->SetVolumeTracks(fVolumeTracks);
  analysisManager->FillNtupleIColumn(1, 0, run->GetNumberOfEvent() + run->GetSkipEvents());
  analysisManager->FillNtupleDColumn(1, 21, evweght);
  analysisManager->FillNtupleIColumn(1, 25, run->GetRunID());
  analysisManager->AddNtupleRow(1);
  fVolumeTracks.clear();
}



void EventAction::AddTrajectoryStep(const G4int trckid, const HistoManager::TrajectoryTuple &trjpoint,
                                         const std::tuple<G4int, G4int, G4int> &trckinfo)
{
  const auto &itr = fTrajectories.find(trckid);
  if (itr != fTrajectories.end()) {
    itr->second.push_back(trjpoint);
  } else {
    fTrajectories[trckid] = std::vector<HistoManager::TrajectoryTuple>(1, trjpoint);
    fTrajTracks[trckid] = trckinfo;
  }
}



void EventAction::UpdateTrajectories(const G4Event *event)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  HistoManager *histManager = run->GetHistoManager();

  G4double evweght = 1.0;
  const EventInfo *evinf = dynamic_cast<EventInfo*>(event->GetUserInformation());
  if (evinf) evweght = evinf->GetWeight();
  
  for (const auto &itr : fTrajTracks) {
    const G4int trckid = itr.first;

    analysisManager->FillNtupleIColumn(5, 0, run->GetNumberOfEvent() + run->GetSkipEvents());
    analysisManager->FillNtupleIColumn(5, 1, trckid);
    analysisManager->FillNtupleIColumn(5, 2, std::get<0>(itr.second));
    analysisManager->FillNtupleIColumn(5, 3, std::get<1>(itr.second));
    analysisManager->FillNtupleIColumn(5, 4, std::get<2>(itr.second));
    analysisManager->FillNtupleDColumn(5, 5, evweght);
    analysisManager->FillNtupleIColumn(5, 17, run->GetRunID());
    histManager->SetTrajectory(fTrajectories[trckid]);
    analysisManager->AddNtupleRow(5);
  }

  fTrajectories.clear();
  fTrajTracks.clear();
}



void  EventAction::ClearData()
{
  fHitId.clear();
  fHitEDep.clear();
  for (auto mitr = fHitTrackGPos.begin(); mitr != fHitTrackGPos.end(); ++mitr) {
    mitr->second.clear();
  }
  fHitTrackGPos.clear();
  
  for (auto mitr = fHitTrackId.begin(); mitr != fHitTrackId.end(); ++mitr) {
    mitr->second.clear();  
  }
  fHitTrackId.clear();
  
  for (auto mitr = fHitTrackGTime.begin(); mitr != fHitTrackGTime.end(); ++mitr) {
    mitr->second.clear();  
  }
  fHitTrackGTime.clear();
  
  for (auto mitr = fHitTrackEDep.begin(); mitr != fHitTrackEDep.end(); ++mitr) {
    mitr->second.clear();  
  }
  fHitTrackEDep.clear();

  ftrackId.clear();
  ftrackVtxPos.clear();
  ftrackMomentum.clear();
  ftrackE.clear();
  ftrackPDG.clear();
  ftrackPTD.clear();
  ftrackPProc.clear();

  fTrajectories.clear();
  fTrajTracks.clear();

  fNHits = 0;

//   fVolumeTracks.clear();
}

