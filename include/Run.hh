//
/// \brief Definition of the Run class
//

#ifndef Run_h
#define Run_h 1

#include <vector>
#include <map>
#include <tuple>

#include "G4Run.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "globals.hh"

class DetectorConstruction;
class G4ParticleDefinition;
class HistoManager;
class G4VPhysicalVolume;
class RunMessenger;


class Run : public G4Run
{
 public:
   Run(DetectorConstruction*);
  ~Run();

 public:

   void SetPrimary(G4ParticleDefinition* particle, G4double energy);

   void AddEnergy (G4double edep)
              {fEnergyDeposit += edep; fEnergyDeposit2 += edep*edep;};

   void AddTrakLenCharg (G4double length)
              {fTrakLenCharged += length; fTrakLenCharged2 += length*length;};

   void AddTrakLenNeutr (G4double length)
              {fTrakLenNeutral += length; fTrakLenNeutral2 += length*length;};

   void AddMscProjTheta (G4double theta)
              {if (std::abs(theta) <= fMscThetaCentral) { fMscEntryCentral++;
                 fMscProjecTheta += theta;  fMscProjecTheta2 += theta*theta;}
              };

   void CountStepsCharg (G4int nSteps)
              {fNbStepsCharged += nSteps; fNbStepsCharged2 += nSteps*nSteps;};

   void CountStepsNeutr (G4int nSteps)
              {fNbStepsNeutral += nSteps; fNbStepsNeutral2 += nSteps*nSteps;};

   void CountParticles (G4ParticleDefinition* part)
              {     if (part == G4Gamma::Gamma())       fNbGamma++ ;
               else if (part == G4Electron::Electron()) fNbElect++ ;
               else if (part == G4Positron::Positron()) fNbPosit++ ; };

   void CountTransmit (G4int flag)
              {     if (flag == 1)  fTransmit[0]++;
               else if (flag == 2) {fTransmit[0]++; fTransmit[1]++; }};

   void CountReflect (G4int flag)
              {     if (flag == 1)  fReflect[0]++;
               else if (flag == 2) {fReflect[0]++; fReflect[1]++; }};

   void AddEnergyLeak (G4double eleak, G4int index)
            {fEnergyLeak[index] += eleak; fEnergyLeak2[index] += eleak*eleak;};

   G4double ComputeMscHighland();

   virtual void Merge(const G4Run*);

   void EndOfRun();

   void SetHistoManager (HistoManager *hmng);
   HistoManager *GetHistoManager () { return fHistoManager; }
   G4double GetTreeCutX(void) const { return fTreeCutX; }
   G4double GetTreeCutY(void) const { return fTreeCutY; }
   G4ParticleDefinition *GetTreeParticle(void) const { return fTreeParticle; }

   void ScanGeometry(const G4VPhysicalVolume *topvol, int verbose = 0);
   void SaveDetectorTransformation(const G4VPhysicalVolume *pv, std::vector<const G4VPhysicalVolume*> pvol, const int l);
   void SaveVolumeTransformation(const G4VPhysicalVolume *pv, std::vector<const G4VPhysicalVolume*> pvol, const int l);
   void SetInterceptVolumes(const std::map<G4String, std::vector<G4int> > &vmap) { fInterceptVol = vmap; }
   void SetSensitiveVolumes(const std::map<G4String, std::tuple<G4int, G4int, G4int> > &vmap) { fSensitiveVol = vmap; }
   void SetInterceptVolECut(const std::map<G4int, G4double> &mevc) { fInterceptVolECut = mevc; }
   void SetSaveTrajectory(const G4bool ff) { fSaveTrajectory = ff; }
   void SetSkipEvents(const size_t skipeve) { fSkipEvents = skipeve; }

   const std::map<G4String, std::vector<G4int> > &GetInterceptVolumes() const { return fInterceptVol; }
   const std::map<G4String, std::tuple<G4int, G4int, G4int> > &GetSensitiveVolumes() const { return fSensitiveVol; }
   G4int GetLayerDepth(const G4String &str) const { return std::get<1>(fSensitiveVol.at(str)); }
   G4int GetSubDetDepth(const G4String &str) const { return std::get<2>(fSensitiveVol.at(str)); }
   G4double GetInterceptVolECut(const G4int volid) const { const auto& itr = fInterceptVolECut.find(volid);
                 if (itr != fInterceptVolECut.end()) return itr->second; else return 0.0; }
   G4bool GetSaveTrajectory() { return fSaveTrajectory; }
   size_t GetSkipEvents() const { return fSkipEvents; }

 private:
    DetectorConstruction*  fDetector;
    G4ParticleDefinition*  fParticle;
    G4double  fEkin;

    G4double fEnergyDeposit,  fEnergyDeposit2;
    G4double fTrakLenCharged, fTrakLenCharged2;
    G4double fTrakLenNeutral, fTrakLenNeutral2;
    G4double fNbStepsCharged, fNbStepsCharged2;
    G4double fNbStepsNeutral, fNbStepsNeutral2;
    G4double fMscProjecTheta, fMscProjecTheta2;
    G4double fMscThetaCentral;
    
    G4int    fNbGamma, fNbElect, fNbPosit;
    G4int    fTransmit[2],   fReflect[2];
    G4int    fMscEntryCentral;

    G4double fEnergyLeak[2],  fEnergyLeak2[2];

    HistoManager *fHistoManager;
    G4double  fTreeCutX, fTreeCutY;
    G4ParticleDefinition  *fTreeParticle;

    std::map<G4String, std::vector<G4int> > fInterceptVol;     //VolName VolID, depth for copy number
    std::map<G4int, G4double> fInterceptVolECut;   //VolID  Tracks_ECut
    std::map<G4String, std::tuple<G4int, G4int, G4int> > fSensitiveVol;  //VolName, VolID, SubDetDepth, LayerDepth

    G4bool  fSaveTrajectory;
    size_t  fSkipEvents;
};


#endif

