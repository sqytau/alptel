//
/// \brief Definition of the EventAction class
//

#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "EventAction.hh"
#include "globals.hh"

#include "G4VUserEventInformation.hh"
#include "HistoManager.hh"

#include <vector>
#include <map>
#include <tuple>
#include <bitset>

#define DET_ENCODING_BIT_SIZE 64

class G4Track;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event *event);

    void AddEnergy      (G4double edep)   {fEnergyDeposit  += edep;}
    void AddTrakLenCharg(G4double length) {fTrakLenCharged += length;}
    void AddTrakLenNeutr(G4double length) {fTrakLenNeutral += length;}

    void CountStepsCharg ()               {fNbStepsCharged++ ;}
    void CountStepsNeutr ()               {fNbStepsNeutral++ ;}

    void SetTransmitFlag (G4int flag)
                           {if (flag > fTransmitFlag) fTransmitFlag = flag;}
    void SetReflectFlag  (G4int flag)
                           {if (flag > fReflectFlag)   fReflectFlag = flag;}

    void AddKinEnegyCharged (const G4double ek) { fKinetikEnergyCharged.push_back(ek); }
    void AddKinEnegyGamma (const G4double ek) { fKinetikEnergyGamma.push_back(ek); }

    void AddOpppTrackerCellDepposit(const G4int opppdet, const G4int layerid, const G4int cellx, const G4int celly,
                                    const G4double edep, const G4Track* aTrack);
    void AddVolumeTrack(const HistoManager::VolumeTrackTuple &ttrck) {fVolumeTracks.push_back(ttrck);}
    void AddTrajectoryStep(const G4int trckid, const HistoManager::TrajectoryTuple &trjpoint,
                                         const std::tuple<G4int, G4int, G4int> &trckinfo);
    void UpdateOPPPTracker();
    void UpdateVolumeTracks();
    void UpdateTrajectories(const G4Event *event);
    void ClearData();
    
    int TestSegmentationEncoding(const int maxx, const int maxy);

  private:
    G4double fEnergyDeposit;
    G4double fTrakLenCharged, fTrakLenNeutral;
    G4int    fNbStepsCharged, fNbStepsNeutral;
    G4int    fTransmitFlag,   fReflectFlag;

    std::vector<G4double> fKinetikEnergyGamma;
    std::vector<G4double> fKinetikEnergyCharged;

    int fNHits;

    std::map<std::string, G4int> fHitId;
    std::map<G4int, G4double> fHitEDep;
    std::map<G4int, std::vector<G4ThreeVector> > fHitTrackGPos;
    std::map<G4int, std::vector<G4int> > fHitTrackId;
    std::map<G4int, std::vector<G4double> > fHitTrackGTime;
    std::map<G4int, std::vector<G4double> > fHitTrackEDep;

    std::vector<G4int> ftrackId;
    std::vector<G4ThreeVector> ftrackVtxPos;
    std::vector<G4ThreeVector> ftrackMomentum;
    std::vector<G4double> ftrackE;
    std::vector<G4int> ftrackPDG;
    std::vector<G4int> ftrackPProc;
    std::vector<G4int> ftrackPTD;     // Parent track id
    std::vector<HistoManager::VolumeTrackTuple> fVolumeTracks;
    std::map<G4int, std::vector<HistoManager::TrajectoryTuple> > fTrajectories;
    std::map<G4int, std::tuple<G4int, G4int, G4int> > fTrajTracks;

    std::vector<unsigned int> fBSize;      // Number of bits to encode sensitive detector cells
    std::vector<unsigned int> fBSizeSum;   // std::partial_sum of fBSize, Number of bis to shift

};



class EventInfo : public G4VUserEventInformation
{
  public:
    EventInfo() : G4VUserEventInformation(), weight(0.0), mctrackid(-1) {};
    EventInfo(const G4double w, const G4int mctid) : G4VUserEventInformation(), weight(w), mctrackid(mctid) {};
    virtual void Print() const { G4cout << "EventInfo: event weight: " << weight << " track: " << mctrackid <<  G4endl; }

    void SetWeight(const G4double w) {weight = w;}
    void SetMCTrackId(const G4int mctid) {mctrackid = mctid;}

    G4double GetWeight() const {return weight; }
    G4int GetMCTrackId() const {return mctrackid; }

  protected:
    G4double weight;
    G4int    mctrackid;
};

#endif

 
