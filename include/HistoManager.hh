//
/// \brief Definition of the HistoManager class
//

#ifndef HistoManager_h
#define HistoManager_h 1

#include <algorithm>
#include <tuple>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>


#include "g4root.hh"
//#include "g4xml.hh"

class HistoMessenger;

class HistoManager
{
  public:
    HistoManager();
   ~HistoManager();
    
    void SetTreeCutX(const G4double cx) { fTreeCutX = cx; }
    void SetTreeCutY(const G4double cy) { fTreeCutY = cy; }
    void SetTreeParticle(const G4String pstr) {fTreeParticle = pstr; }

    G4double GetTreeCutX() const { return fTreeCutX; }
    G4double GetTreeCutY() const { return fTreeCutY; }
    G4String GetTreeParticle() const { return fTreeParticle;}
    
    void SetHitTrackList(const std::vector<G4int> &vht) { fvHitTrackList = vht; }
    void SetHitTrackPosition(const std::vector<G4ThreeVector> &vht);
    void SetHitTrackTime(const std::vector<G4double> &vht) { ftrackt = vht; }
    void SetHitTrackEDep(const std::vector<G4double> &vht) { ftrackedep = vht; }
    void SetTracks(const std::vector<G4int> &vht)       { fvTracks = vht; }
    void SetTracksVtx(const std::vector<G4ThreeVector> &vht);
    void SetTracksMomentum(const std::vector<G4ThreeVector> &vht);
    void SetTracksE(const std::vector<G4double> &vht) { fE = vht; }
    void SetTracksPDG(const std::vector<G4int> &vht) { fPDG = vht; }
    void SetTracksPhysProc(const std::vector<G4int> &vht) { fPhysProc = vht; }
    void SetTracksPTId(const std::vector<G4int> &vht) { fPTId = vht; }

    typedef std::tuple<G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, 
                       G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, 
                       G4int, G4int, G4int, G4int, G4int, G4int, G4double> VolumeTrackTuple;
    void SetVolumeTracks(const std::vector<VolumeTrackTuple> &tvtracks);
    
    typedef std::tuple<G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double, 
                        G4double, G4int, G4int> TrajectoryTuple;

    void SetTrajectory(const std::vector<TrajectoryTuple> &tvtracks);
 
  private:
    void Book();
    void AllocateVolumeTracksMap();
    
    G4String  fFileName;
    
    HistoMessenger *fHMessanger;
    G4double  fTreeCutX, fTreeCutY;
    G4String  fTreeParticle;
    
    std::vector<G4int> fvHitTrackList;
    std::vector<G4int> fvTracks;

    std::vector<G4double> fVtxx, fVtxy, fVtxz;
    std::vector<G4double> ftrackx, ftracky, ftrackz, ftrackt, ftrackedep;
    std::vector<G4double> fPx, fPy, fPz, fE;
    std::vector<G4int> fPDG;
    std::vector<G4int> fPhysProc;
    std::vector<G4int> fPTId;

    std::map <G4int, std::vector<G4double> > fvolTrackDMap;
    std::map <G4int, std::vector<G4int> > fvolTrackIMap;
    
    std::map <G4int, std::vector<G4double> > fTrajectoryD;
    std::map <G4int, std::vector<G4int> > fTrajectoryI;
    std::vector<G4String> fTrajectoryS;

};



#endif

